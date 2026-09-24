"""
Bioinformatics File Analyzer  (FASTA / FASTQ / VCF / SAM)
=========================================================

Run:
    pip install streamlit pandas numpy plotly reportlab
    streamlit run app.py                                  # default 200 MB upload limit
    streamlit run app.py --server.maxUploadSize 1024      # allow larger files

Structure (each layer is independent so the core can move into a FastAPI backend):
    1. Core   - detection, parsers, analyzers, report builders (no UI code)
    2. Charts - Plotly figure builders (take DataFrames, return figures)
    3. UI     - Streamlit rendering only

`process_file(raw_bytes) -> dict` is the single entry point of the core: it is a
pure function of the file bytes, so a future FastAPI endpoint can call it directly.
"""
from __future__ import annotations

import gzip
import io
import os
import re
import sys
from collections import namedtuple
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Iterator

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.platypus import Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle
from xml.sax.saxutils import escape

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #
UNKNOWN = "Unknown"
RECORD_LABELS = {"FASTA": "sequences", "FASTQ": "reads", "VCF": "variants", "SAM": "alignments"}
MAX_WARNINGS = 25            # parser warnings kept per file (rest are only counted)
HEAD_CHARS = 262_144         # amount of text inspected for format detection
PDF_PREVIEW_ROWS = 25        # rows of the detail table shown in the PDF
UI_TABLE_ROWS = 1_000        # rows of the detail table shown in the browser
HIST_SAMPLE_LIMIT = 200_000  # max points sent to the browser for a raw histogram
MAX_VCF_BARS = 40            # max contigs drawn in the VCF bar chart
PLOT_TEMPLATE = "plotly_white"

# Lightweight record types: namedtuples use far less memory than dicts on big files
FastaRecord = namedtuple("FastaRecord", ["id", "length", "gc_count", "n_count"])
FastqRecord = namedtuple("FastqRecord", ["id", "length", "qual_sum"])  # qual_sum = sum of Phred scores
VcfRecord = namedtuple("VcfRecord", ["chrom", "pos", "id", "ref", "alt", "qual", "filter"])
SamRecord = namedtuple("SamRecord", ["qname", "flag", "rname", "pos", "mapq"])


@dataclass
class ParseResult:
    """Output of every parser: records plus non-fatal problems found on the way."""
    records: list = field(default_factory=list)
    warnings: list = field(default_factory=list)
    skipped: int = 0                      # malformed lines / records that were dropped
    meta: dict = field(default_factory=dict)

    def warn(self, message: str, skipped: bool = False) -> None:
        """Register a warning (capped) and optionally count a dropped record."""
        if skipped:
            self.skipped += 1
        if len(self.warnings) < MAX_WARNINGS:
            self.warnings.append(message)
        elif len(self.warnings) == MAX_WARNINGS:
            self.warnings.append("Further warnings were suppressed.")


# --------------------------------------------------------------------------- #
# Small helpers
# --------------------------------------------------------------------------- #
def iter_lines(text: str) -> Iterator[str]:
    """Yield lines lazily (no giant list in memory); handles both \\n and \\r\\n."""
    start, n = 0, len(text)
    while start < n:
        end = text.find("\n", start)
        if end == -1:
            end = n
        yield text[start:end].rstrip("\r")
        start = end + 1


def _short(text: str, limit: int = 40) -> str:
    """Trim long strings for use in messages."""
    return text if len(text) <= limit else text[: limit - 3] + "..."


def human_size(num_bytes: float) -> str:
    """Format a byte count, e.g. 1536 -> '1.5 KB'."""
    size = float(num_bytes)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if size < 1024 or unit == "TB":
            return f"{size:.0f} {unit}" if unit == "B" else f"{size:.1f} {unit}"
        size /= 1024
    return f"{size:.1f} TB"


def format_metric(value: Any) -> str:
    """Uniform display of metric values (UI, CSV header and PDF)."""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return "N/A"
    if isinstance(value, (int, np.integer)):
        return f"{int(value):,}"
    if isinstance(value, (float, np.floating)):
        return f"{float(value):,.2f}"
    return str(value)


def chrom_sort_key(name: str) -> tuple:
    """Natural chromosome order: 1..22, X, Y, M/MT, then other contigs alphabetically."""
    core = name[3:] if name.lower().startswith("chr") else name
    if core.isdigit():
        return (0, int(core), "")
    special = {"X": 1, "Y": 2, "M": 3, "MT": 3}
    if core.upper() in special:
        return (1, special[core.upper()], "")
    return (2, 0, core)


def n50(lengths: np.ndarray) -> int:
    """N50: length L such that sequences >= L contain at least half of all bases."""
    ordered = np.sort(lengths)[::-1]
    cumulative = np.cumsum(ordered)
    return int(ordered[np.searchsorted(cumulative, cumulative[-1] / 2)])


# --------------------------------------------------------------------------- #
# 1. File type detection (content based, extension is ignored)
# --------------------------------------------------------------------------- #
_SAM_HEADER_RE = re.compile(r"^@(?:(?:HD|SQ|RG|PG)\t[A-Za-z][A-Za-z0-9]:|CO\t)")
_SEQ_LINE_RE = re.compile(r"[A-Za-z*\-.]+")


def _looks_like_sam_record(line: str) -> bool:
    """True if a line has >= 11 tab fields with numeric FLAG, POS and MAPQ."""
    f = line.split("\t")
    return len(f) >= 11 and f[1].isdigit() and f[3].isdigit() and f[4].isdigit()


def detect_file_type(content: str) -> str:
    """Return 'FASTA', 'FASTQ', 'VCF', 'SAM' or 'Unknown' by inspecting the first lines."""
    lines = []
    for line in iter_lines(content[:HEAD_CHARS]):
        if line.strip():
            lines.append(line)
        if len(lines) >= 20:
            break
    if not lines:
        return UNKNOWN
    first = lines[0]

    # VCF: '##fileformat=VCF...' or a '#CHROM' column header line
    if first.startswith("##fileformat=VCF") or any(l.startswith("#CHROM") for l in lines if l[0] == "#"):
        return "VCF"
    # SAM: '@HD/@SQ/@RG/@PG/@CO' header lines
    if _SAM_HEADER_RE.match(first):
        return "SAM"
    # FASTQ: '@' header and a '+' separator on the third line (SAM headers handled above)
    if first.startswith("@"):
        if len(lines) >= 3 and lines[2].startswith("+"):
            return "FASTQ"
        if len(lines) == 2 and _SEQ_LINE_RE.fullmatch(lines[1].strip()):
            return "FASTQ"  # truncated file: reported later as incomplete
        return UNKNOWN
    # FASTA: '>' header followed only by sequence-like lines
    if first.startswith(">"):
        body = [l.strip() for l in lines[1:] if not l.startswith(">")]
        return "FASTA" if all(_SEQ_LINE_RE.fullmatch(l) for l in body) else UNKNOWN
    # SAM without header: first line already looks like an alignment record
    if _looks_like_sam_record(first):
        return "SAM"
    return UNKNOWN


# --------------------------------------------------------------------------- #
# 2. Parsers  (single pass, tolerant of malformed input)
# --------------------------------------------------------------------------- #
_FASTA_SEQ_RE = _SEQ_LINE_RE


def _finish_fasta_record(result: ParseResult, seq_id, length: int, gc: int, n: int) -> None:
    """Store a completed FASTA record; header-only records are dropped with a warning."""
    if seq_id is None:
        return
    if length == 0:
        result.warn(f"Sequence '{_short(seq_id)}' has a header but no sequence data; skipped.", skipped=True)
        return
    result.records.append(FastaRecord(seq_id, length, gc, n))


def parse_fasta(content: str) -> ParseResult:
    """Parse FASTA. Sequences are not kept in memory: only length, G+C and N counts."""
    result = ParseResult()
    seq_id, length, gc, n = None, 0, 0, 0
    for line_no, raw in enumerate(iter_lines(content), start=1):
        line = raw.strip()
        if not line:
            continue
        if line[0] == ">":  # new record: close the previous one
            _finish_fasta_record(result, seq_id, length, gc, n)
            header = line[1:].strip()
            seq_id = header.split(None, 1)[0] if header else f"unnamed_{len(result.records) + 1}"
            length = gc = n = 0
            continue
        if seq_id is None:
            result.warn(f"Line {line_no}: sequence data before the first '>' header ignored.", skipped=True)
            continue
        if not _FASTA_SEQ_RE.fullmatch(line):
            result.warn(f"Line {line_no}: invalid characters in '{_short(seq_id)}'; line skipped.", skipped=True)
            continue
        upper = line.upper()
        length += len(upper)
        gc += upper.count("G") + upper.count("C")
        n += upper.count("N")
    _finish_fasta_record(result, seq_id, length, gc, n)
    return result


def parse_fastq(content: str) -> ParseResult:
    """Parse FASTQ (strict 4-line records). Stores length and the Phred+33 score sum per read."""
    result = ParseResult()
    lines = iter_lines(content)
    for header in lines:
        if not header.strip():
            continue
        if header[0] != "@":
            result.warn(f"Line outside a FASTQ record skipped: '{_short(header)}'.", skipped=True)
            continue
        seq, plus, qual = next(lines, None), next(lines, None), next(lines, None)
        header_tokens = header[1:].split(None, 1)
        read_id = header_tokens[0] if header_tokens else f"read_{len(result.records) + 1}"
        if qual is None:  # file ended in the middle of a record
            result.warn(f"Read '{_short(read_id)}' is truncated (incomplete final record).", skipped=True)
            break
        seq, qual = seq.strip(), qual.strip()
        if not plus.startswith("+"):
            result.warn(f"Read '{_short(read_id)}': missing '+' separator line; skipped.", skipped=True)
            continue
        if len(seq) != len(qual):
            result.warn(f"Read '{_short(read_id)}': sequence and quality lengths differ; skipped.", skipped=True)
            continue
        q_bytes = qual.encode("ascii", "replace")
        if q_bytes and min(q_bytes) < 33:  # below '!' is not valid Phred+33
            result.warn(f"Read '{_short(read_id)}': invalid quality characters; skipped.", skipped=True)
            continue
        result.records.append(FastqRecord(read_id, len(seq), sum(q_bytes) - 33 * len(q_bytes)))
    return result


def parse_vcf(content: str) -> ParseResult:
    """Parse VCF data lines (first 8 mandatory columns) and read sample names from #CHROM."""
    result = ParseResult(meta={"samples": [], "has_header": False})
    for line_no, line in enumerate(iter_lines(content), start=1):
        if not line.strip() or line.startswith("##"):
            continue
        if line.startswith("#"):
            cols = line.lstrip("#").rstrip().split("\t")
            if cols[0] == "CHROM":
                result.meta["has_header"] = True
                result.meta["samples"] = cols[9:]
            continue
        f = line.split("\t", 8)  # CHROM POS ID REF ALT QUAL FILTER INFO [rest]
        if len(f) < 8 or not f[0]:
            result.warn(f"Line {line_no}: fewer than 8 tab-separated VCF columns; skipped.", skipped=True)
            continue
        try:
            pos = int(f[1])
        except ValueError:
            result.warn(f"Line {line_no}: POS '{_short(f[1], 15)}' is not an integer; skipped.", skipped=True)
            continue
        try:
            qual = None if f[5] == "." else float(f[5])
        except ValueError:
            qual = None
        result.records.append(VcfRecord(sys.intern(f[0]), pos, f[2], f[3], f[4], qual, f[6]))
    if result.records and not result.meta["has_header"]:
        result.warn("No '#CHROM' header line found; sample names are unavailable.")
    return result


def parse_sam(content: str) -> ParseResult:
    """Parse SAM alignment lines (11 mandatory fields) and count @SQ reference lines."""
    result = ParseResult(meta={"references": 0})
    for line_no, line in enumerate(iter_lines(content), start=1):
        if not line.strip():
            continue
        if line[0] == "@":
            if line.startswith("@SQ"):
                result.meta["references"] += 1
            continue
        f = line.split("\t", 11)  # only the first 11 fields are needed
        if len(f) < 11:
            result.warn(f"Line {line_no}: fewer than 11 tab-separated SAM fields; skipped.", skipped=True)
            continue
        try:
            flag, pos, mapq = int(f[1]), int(f[3]), int(f[4])
        except ValueError:
            result.warn(f"Line {line_no}: FLAG/POS/MAPQ are not integers; skipped.", skipped=True)
            continue
        if flag < 0 or not 0 <= mapq <= 255:
            result.warn(f"Line {line_no}: FLAG or MAPQ out of range; skipped.", skipped=True)
            continue
        result.records.append(SamRecord(f[0], flag, sys.intern(f[2]), pos, mapq))
    return result


# --------------------------------------------------------------------------- #
# 3. Analysis  (each returns (metrics dict, detail DataFrame))
# --------------------------------------------------------------------------- #
def analyze_fasta(parsed: ParseResult) -> tuple[dict, pd.DataFrame]:
    """Sequence count, length statistics and GC content (G+C over non-N bases)."""
    df = pd.DataFrame(parsed.records, columns=FastaRecord._fields).rename(columns={"id": "sequence_id"})
    valid_bases = (df["length"] - df["n_count"]).replace(0, np.nan)
    df["gc_percent"] = (df["gc_count"] / valid_bases * 100).round(2)
    lengths = df["length"].to_numpy()
    total_bases = int(lengths.sum())
    total_valid = total_bases - int(df["n_count"].sum())
    metrics = {
        "Sequences": len(df),
        "Total bases": total_bases,
        "Shortest": int(lengths.min()),
        "Longest": int(lengths.max()),
        "Mean length": float(lengths.mean()),
        "Median length": float(np.median(lengths)),
        "N50": n50(lengths),
        "Overall GC (%)": float(df["gc_count"].sum() / total_valid * 100) if total_valid else float("nan"),
        "Mean GC / seq (%)": float(df["gc_percent"].mean()),
    }
    return metrics, df[["sequence_id", "length", "gc_count", "n_count", "gc_percent"]]


def analyze_fastq(parsed: ParseResult) -> tuple[dict, pd.DataFrame]:
    """Read count, length distribution and average Phred+33 quality."""
    df = pd.DataFrame(parsed.records, columns=FastqRecord._fields)
    total_bases = int(df["length"].sum())
    read_mean_q = df["qual_sum"] / df["length"].replace(0, np.nan)  # mean quality per read
    metrics = {
        "Reads": len(df),
        "Total bases": total_bases,
        "Min read length": int(df["length"].min()),
        "Mean read length": float(df["length"].mean()),
        "Median read length": float(df["length"].median()),
        "Max read length": int(df["length"].max()),
        "Avg quality (Phred+33)": float(df["qual_sum"].sum() / total_bases) if total_bases else float("nan"),
        "Reads with mean Q>=30 (%)": float((read_mean_q >= 30).mean() * 100),
    }
    # Compact length distribution (scales to millions of reads)
    dist = df.groupby("length", as_index=False).agg(read_count=("id", "size"), qual_sum=("qual_sum", "sum"))
    dist["mean_quality"] = (dist["qual_sum"] / (dist["length"] * dist["read_count"]).replace(0, np.nan)).round(2)
    table = dist.rename(columns={"length": "read_length"})[["read_length", "read_count", "mean_quality"]]
    return metrics, table


def classify_variant(ref: str, alt: str) -> str:
    """SNP (1 bp -> 1 bp), Indel (length change) or Other (MNP, symbolic, breakend, '*')."""
    alts = alt.split(",")
    if any(a in ("*", ".") or a.startswith("<") or "[" in a or "]" in a for a in alts):
        return "Other"
    if len(ref) == 1 and all(len(a) == 1 for a in alts):
        return "SNP"
    if all(len(a) == len(ref) for a in alts):
        return "Other"  # multi-nucleotide substitution
    return "Indel"


def analyze_vcf(parsed: ParseResult) -> tuple[dict, pd.DataFrame]:
    """Total variants and per-chromosome distribution with SNP/Indel/Other breakdown."""
    df = pd.DataFrame(parsed.records, columns=VcfRecord._fields)
    df["type"] = [classify_variant(r, a) for r, a in zip(df["ref"], df["alt"])]
    counts = df.groupby(["chrom", "type"]).size().unstack(fill_value=0)
    counts = counts.reindex(columns=["SNP", "Indel", "Other"], fill_value=0)
    counts = counts.loc[sorted(counts.index, key=chrom_sort_key)]
    counts["variants"] = counts.sum(axis=1)
    counts.columns.name = None
    table = counts.reset_index().rename(columns={"chrom": "chromosome", "SNP": "snps", "Indel": "indels", "other": "other", "Other": "other"})
    table["percent"] = (table["variants"] / len(df) * 100).round(2)
    table = table[["chromosome", "variants", "percent", "snps", "indels", "other"]]
    metrics = {
        "Total variants": len(df),
        "Chromosomes/contigs": int(df["chrom"].nunique()),
        "SNPs": int((df["type"] == "SNP").sum()),
        "Indels": int((df["type"] == "Indel").sum()),
        "Other": int((df["type"] == "Other").sum()),
        "PASS variants": int((df["filter"] == "PASS").sum()),
        "Mean QUAL": float(pd.to_numeric(df["qual"], errors="coerce").mean()),
        "Samples": len(parsed.meta.get("samples", [])),
    }
    return metrics, table


def analyze_sam(parsed: ParseResult) -> tuple[dict, pd.DataFrame]:
    """Alignment totals and MAPQ statistics (mapped alignments, MAPQ 255 = 'unavailable' excluded)."""
    df = pd.DataFrame(parsed.records, columns=SamRecord._fields)
    mapped = (df["flag"] & 4) == 0                      # FLAG bit 0x4 = segment unmapped
    usable = mapped & (df["mapq"] != 255)
    mapq = df.loc[usable, "mapq"]
    total, mapped_n = len(df), int(mapped.sum())
    metrics = {
        "Total alignments": total,
        "Mapped": mapped_n,
        "Unmapped": total - mapped_n,
        "Mapped (%)": mapped_n / total * 100,
        "Mean MAPQ": float(mapq.mean()),
        "Median MAPQ": float(mapq.median()),
        "Std dev MAPQ": float(mapq.std(ddof=0)),
        "Min MAPQ": float(mapq.min()) if len(mapq) else float("nan"),
        "Max MAPQ": float(mapq.max()) if len(mapq) else float("nan"),
        "MAPQ>=30 (%)": float((mapq >= 30).mean() * 100) if len(mapq) else float("nan"),
        "References (@SQ)": int(parsed.meta.get("references", 0)),
    }
    unavailable = int((mapped & (df["mapq"] == 255)).sum())
    if unavailable:
        metrics["MAPQ 255 (n/a)"] = unavailable
    table = mapq.value_counts().sort_index().rename_axis("mapq").reset_index(name="alignments")
    table["percent"] = (table["alignments"] / max(int(table["alignments"].sum()), 1) * 100).round(2)
    return metrics, table


# Dispatch table: format -> (parser, analyzer). Add new formats here.
PIPELINES = {
    "FASTA": (parse_fasta, analyze_fasta),
    "FASTQ": (parse_fastq, analyze_fastq),
    "VCF": (parse_vcf, analyze_vcf),
    "SAM": (parse_sam, analyze_sam),
}


# --------------------------------------------------------------------------- #
# 4. Reports (CSV via pandas, PDF via reportlab)
# --------------------------------------------------------------------------- #
def build_csv(file_type: str, metrics: dict, table: pd.DataFrame) -> bytes:
    """CSV = '# metric: value' summary lines followed by the detail table.
    Read back with pandas.read_csv(path, comment='#')."""
    buf = io.StringIO()
    buf.write(f"# file_type: {file_type}\n")
    for key, value in metrics.items():
        buf.write(f"# {key}: {format_metric(value)}\n")
    table.to_csv(buf, index=False)
    return buf.getvalue().encode("utf-8")


def build_preview(table: pd.DataFrame, n: int = PDF_PREVIEW_ROWS) -> list:
    """First n rows of the detail table as strings (header row first) for the PDF."""
    def cell(v: Any) -> str:
        if v is None or (isinstance(v, float) and np.isnan(v)):
            return "-"
        return f"{v:.2f}" if isinstance(v, (float, np.floating)) else str(v)
    head = table.head(n)
    return [[str(c) for c in head.columns]] + [[cell(v) for v in row] for row in head.itertuples(index=False)]


def _pdf_clip(text: str, max_chars: int) -> str:
    """Make text safe for the built-in PDF fonts (Latin-1) and short enough for its column."""
    text = text.encode("latin-1", "replace").decode("latin-1")
    return text if len(text) <= max_chars else text[: max(max_chars - 3, 1)] + "..."


def _pdf_table(data: list, width: float, header: bool, ratios=None, font_size: int = 8) -> Table:
    """Styled reportlab table with column-aware text clipping."""
    ncols = len(data[0])
    widths = [width * r for r in (ratios or [1 / ncols] * ncols)]
    clipped = [[_pdf_clip(str(c), int(w / (font_size * 0.52))) for c, w in zip(row, widths)] for row in data]
    table = Table(clipped, colWidths=widths, repeatRows=1 if header else 0)
    style = [
        ("FONTSIZE", (0, 0), (-1, -1), font_size),
        ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#c8ced6")),
        ("ROWBACKGROUNDS", (0, 1 if header else 0), (-1, -1), [colors.white, colors.HexColor("#f4f6f8")]),
        ("TOPPADDING", (0, 0), (-1, -1), 3), ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
    ]
    if header:
        style += [("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#1f4e79")),
                  ("TEXTCOLOR", (0, 0), (-1, 0), colors.white), ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold")]
    else:
        style.append(("FONTNAME", (0, 0), (0, -1), "Helvetica-Bold"))
    table.setStyle(TableStyle(style))
    return table


@st.cache_data(show_spinner=False, max_entries=32)
def build_pdf_report(filename: str, size_bytes: int, file_type: str, n_records: int,
                     metrics: dict, preview: list, warnings: list) -> bytes:
    """Formatted one-file PDF summary: file info, statistics, detail preview, warnings."""
    buf = io.BytesIO()
    doc = SimpleDocTemplate(buf, pagesize=A4, leftMargin=18 * mm, rightMargin=18 * mm,
                            topMargin=18 * mm, bottomMargin=18 * mm, title="Bioinformatics analysis report")
    styles = getSampleStyleSheet()
    info = [
        ["File name", filename], ["File size", human_size(size_bytes)], ["Detected format", file_type],
        [RECORD_LABELS[file_type].capitalize(), f"{n_records:,}"],
        ["Generated", datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")],
    ]
    story = [
        Paragraph("Bioinformatics File Analysis Report", styles["Title"]), Spacer(1, 4 * mm),
        Paragraph("File information", styles["Heading2"]),
        _pdf_table(info, doc.width, header=False, ratios=[0.3, 0.7], font_size=9), Spacer(1, 5 * mm),
        Paragraph("Summary statistics", styles["Heading2"]),
        _pdf_table([[k, format_metric(v)] for k, v in metrics.items()], doc.width, header=False,
                   ratios=[0.5, 0.5], font_size=9),
        Spacer(1, 5 * mm),
    ]
    if len(preview) > 1:
        story += [Paragraph(f"Detailed results (first {len(preview) - 1} rows)", styles["Heading2"]),
                  _pdf_table(preview, doc.width, header=True), Spacer(1, 5 * mm)]
    if warnings:
        story.append(Paragraph("Parsing warnings", styles["Heading2"]))
        for w in warnings[:10]:
            safe = escape(w.encode("latin-1", "replace").decode("latin-1"))
            story.append(Paragraph(f"&bull; {safe}", styles["BodyText"]))
    doc.build(story)
    return buf.getvalue()


# --------------------------------------------------------------------------- #
# 5. Core entry point: bytes -> analysis result (cached, UI independent)
# --------------------------------------------------------------------------- #
@st.cache_data(show_spinner=False, max_entries=16)
def process_file(raw: bytes) -> dict:
    """Detect, parse and analyze one file. Never raises: problems land in result['error']."""
    result = {"file_type": UNKNOWN, "error": None, "warnings": [], "skipped": 0, "n_records": 0,
              "metrics": {}, "table": None, "csv": None, "preview": []}
    try:
        if not raw or raw.isspace():
            result["error"] = "The file is empty."
            return result
        if raw[:2] == b"\x1f\x8b":  # transparently support .gz uploads
            try:
                raw = gzip.decompress(raw)
            except (OSError, EOFError):
                result["error"] = "The gzip archive is corrupted or incomplete."
                return result
            if not raw or raw.isspace():
                result["error"] = "The compressed file is empty."
                return result
        if b"\x00" in raw[:4096]:
            result["error"] = "The file looks binary (e.g. BAM/CRAM). Convert it to text (SAM/VCF/FASTA/FASTQ) first."
            return result

        content = raw.decode("utf-8-sig", errors="replace")   # decode once, reuse everywhere
        file_type = detect_file_type(content)
        result["file_type"] = file_type
        if file_type == UNKNOWN:
            result["error"] = "Unsupported or unrecognized format. Supported: FASTA, FASTQ, VCF and SAM."
            return result

        parse_fn, analyze_fn = PIPELINES[file_type]
        parsed = parse_fn(content)
        result.update(warnings=parsed.warnings, skipped=parsed.skipped, n_records=len(parsed.records))
        if not parsed.records:
            result["error"] = f"Detected {file_type}, but no valid records could be parsed (file may be empty, truncated or malformed)."
            return result

        metrics, table = analyze_fn(parsed)
        result.update(metrics=metrics, table=table, csv=build_csv(file_type, metrics, table),
                      preview=build_preview(table))
    except Exception as exc:  # last line of defense: the UI must never crash
      result["error"] = f"Unexpected error while processing the file: {exc}"
    return result


# --------------------------------------------------------------------------- #
# 6. Charts (Plotly only)
# --------------------------------------------------------------------------- #
def _finish(fig: go.Figure, title: str, xlabel: str, ylabel: str, bargap: float = 0.05) -> go.Figure:
    """Shared, responsive figure styling."""
    fig.update_layout(title=title, xaxis_title=xlabel, yaxis_title=ylabel, template=PLOT_TEMPLATE,
                      bargap=bargap, autosize=True, margin=dict(l=10, r=10, t=50, b=10))
    return fig


def _raw_histogram(values: pd.Series, name: str, title: str, xlabel: str, color: str) -> go.Figure:
    """Histogram of raw values, down-sampled if huge so the browser stays responsive."""
    note = ""
    if len(values) > HIST_SAMPLE_LIMIT:
        values, note = values.sample(HIST_SAMPLE_LIMIT, random_state=0), f" (sample of {HIST_SAMPLE_LIMIT:,})"
    fig = px.histogram(values.to_frame(name=name), x=name, nbins=40, color_discrete_sequence=[color])
    fig.update_traces(marker_line_width=1, marker_line_color="white")
    return _finish(fig, title + note, xlabel, "Number of sequences")


def charts_fasta(table: pd.DataFrame) -> list:
    return [
        ("length", _raw_histogram(table["length"], "length", "Sequence length distribution", "Length (bp)", "#1f77b4")),
        ("gc", _raw_histogram(table["gc_percent"].dropna(), "gc_percent", "GC content per sequence", "GC (%)", "#2ca02c")),
    ]


def charts_fastq(table: pd.DataFrame) -> list:
    # Weighted histogram built from the compact (length, count) table -> exact and fast
    fig = px.histogram(table, x="read_length", y="read_count", histfunc="sum",
                       nbins=int(min(60, max(len(table), 1))), color_discrete_sequence=["#ff7f0e"])
    fig.update_traces(marker_line_width=1, marker_line_color="white")
    return [("length", _finish(fig, "Read length distribution", "Read length (bp)", "Number of reads"))]


def charts_sam(table: pd.DataFrame) -> list:
    if table.empty:
        return []
    fig = px.histogram(table, x="mapq", y="alignments", histfunc="sum", color_discrete_sequence=["#9467bd"])
    top = float(table["mapq"].max())
    fig.update_traces(xbins=dict(start=-0.5, end=top + 0.5, size=1), marker_line_width=0.5, marker_line_color="white")
    fig.update_xaxes(range=[-1, top + 1])
    return [("mapq", _finish(fig, "Mapping quality (MAPQ) distribution", "MAPQ", "Number of alignments"))]


def charts_vcf(table: pd.DataFrame) -> list:
    data, note = table, ""
    if len(table) > MAX_VCF_BARS:  # keep the chart readable for assemblies with thousands of contigs
        data = table.nlargest(MAX_VCF_BARS, "variants")
        data = data.loc[sorted(data.index, key=lambda i: chrom_sort_key(data.at[i, "chromosome"]))]
        note = f" (top {MAX_VCF_BARS} contigs)"
    fig = px.bar(data, x="chromosome", y="variants", hover_data=["percent", "snps", "indels", "other"],
                 color_discrete_sequence=["#d62728"])
    fig.update_xaxes(type="category", categoryorder="array", categoryarray=list(data["chromosome"]))
    return [("dist", _finish(fig, "Variants per chromosome" + note, "Chromosome / contig", "Number of variants", 0.2))]


CHARTS = {"FASTA": charts_fasta, "FASTQ": charts_fastq, "SAM": charts_sam, "VCF": charts_vcf}
CHART_NOTES = {
    "FASTA": "GC content = (G+C) / non-N bases. N50 uses sequence lengths.",
    "FASTQ": "Quality scores decoded as Phred+33. Histogram is built from the exact length distribution.",
    "SAM": "MAPQ statistics and histogram cover mapped alignments only (MAPQ 255 = 'not available' is excluded).",
    "VCF": "Variant types: SNP = 1 bp substitution, Indel = length change, Other = MNP / symbolic / breakend.",
}


# --------------------------------------------------------------------------- #
# 7. Streamlit UI
# --------------------------------------------------------------------------- #
EXAMPLES = {
    "FASTA": ("example.fasta", ">seq1 example sequence\nATGCGTACGTTAGCTAGCTA\nGGCTTAAGCC\n>seq2\nGGCTTAAGCCTAGGATCC"),
    "FASTQ": ("example.fastq", "@read1\nGATTACAGATTACAGATTAC\n+\nIIIIHHHGGGFFFEEEDDDC\n@read2\nACGTACGTAC\n+\nIIIIIIIIII"),
    "VCF": ("example.vcf", "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                           "chr1\t12345\trs001\tA\tG\t50\tPASS\tDP=20\nchr2\t67890\t.\tAT\tA\t30\tPASS\tDP=15"),
    "SAM": ("example.sam", "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:248956422\n"
                           "read1\t0\tchr1\t100\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n"
                           "read2\t16\tchr1\t250\t42\t10M\t*\t0\t0\tTTGGCCAATT\tIIIIIIIIII"),
}


def render_sidebar() -> None:
    """Instructions, supported formats and example snippets."""
    with st.sidebar:
        st.header("🧬 Bio File Analyzer")
        st.subheader("Instructions")
        st.markdown("1. Upload one or more files (drag & drop).\n"
                    "2. The format is detected from the **content**, not the extension.\n"
                    "3. Review metrics and interactive charts for each file.\n"
                    "4. Download a CSV or PDF report per file.\n\n"
                    "Gzip-compressed files (`.gz`) are supported.")
        st.subheader("Supported formats")
        st.markdown("- **FASTA**: sequences, lengths, GC content\n- **FASTQ**: reads, length, Phred+33 quality\n"
                    "- **VCF**: variants per chromosome\n- **SAM**: alignments, MAPQ statistics")
        st.subheader("Example formats")
        for name, (filename, snippet) in EXAMPLES.items():
            with st.expander(name):
                st.code(snippet, language=None)
                st.download_button("Download example", snippet + "\n", file_name=filename,
                                   mime="text/plain", key=f"example_{name}")


def render_metrics(metrics: dict, per_row: int = 4) -> None:
    """Metric cards laid out in rows of `per_row` columns."""
    items = list(metrics.items())
    for start in range(0, len(items), per_row):
        cols = st.columns(per_row)
        for col, (label, value) in zip(cols, items[start:start + per_row]):
            col.metric(label, format_metric(value))


def render_overview(files: list, results: list) -> None:
    """Cross-file summary cards shown above the per-file sections."""
    ok = [r for r in results if not r["error"]]
    st.subheader("Overview")
    cols = st.columns(5)
    cols[0].metric("Files uploaded", len(files))
    cols[1].metric("Files analyzed", len(ok))
    cols[2].metric("Files with issues", len(results) - len(ok))
    cols[3].metric("Records parsed", f"{sum(r['n_records'] for r in ok):,}")
    cols[4].metric("Total size", human_size(sum(f.size for f in files)))


def render_file_section(idx: int, name: str, size: int, res: dict) -> None:
    """One bordered card per uploaded file."""
    with st.container(border=True):
        st.subheader(f"📄 {name}")
        c1, c2, c3 = st.columns([3, 1, 1])
        c1.markdown(f"**File name**  \n`{name}`")
        c2.markdown(f"**Size**  \n{human_size(size)}")
        c3.markdown(f"**Detected type**  \n{res['file_type']}")

        if res["error"]:
            st.error(res["error"])
            return

        if res["skipped"]:
            st.warning(f"{res['skipped']:,} malformed line(s)/record(s) were skipped; results cover the valid data only.")
        if res["warnings"]:
            with st.expander(f"Parsing notes ({len(res['warnings'])})"):
                for w in res["warnings"]:
                    st.write(f"- {w}")

        st.markdown("#### Summary")
        render_metrics(res["metrics"])

        st.markdown("#### Visualizations")
        charts = CHARTS[res["file_type"]](res["table"])
        if charts:
            for col, (key, fig) in zip(st.columns(len(charts)), charts):
                with col:
                    st.plotly_chart(fig, key=f"chart_{idx}_{key}")
        else:
            st.info("No chart available: there is no data to plot.")
        st.caption(CHART_NOTES[res["file_type"]])

        with st.expander("Detailed results"):
            st.dataframe(res["table"].head(UI_TABLE_ROWS), hide_index=True)
            if len(res["table"]) > UI_TABLE_ROWS:
                st.caption(f"Showing the first {UI_TABLE_ROWS:,} of {len(res['table']):,} rows. The CSV contains all rows.")

        st.markdown("#### Download report")
        stem = os.path.splitext(name)[0] or "report"
        d1, d2 = st.columns(2)
        d1.download_button("⬇️ CSV (analysis results)", res["csv"], file_name=f"{stem}_analysis.csv",
                           mime="text/csv", key=f"csv_{idx}")
        try:
            pdf = build_pdf_report(name, size, res["file_type"], res["n_records"],
                                   res["metrics"], res["preview"], res["warnings"])
            d2.download_button("⬇️ PDF (formatted summary)", pdf, file_name=f"{stem}_report.pdf",
                               mime="application/pdf", key=f"pdf_{idx}")
        except Exception as exc:
            d2.error(f"PDF could not be generated: {exc}")
        st.caption("CSV: summary lines start with '#'; read with pandas.read_csv(path, comment='#').")


def run_pipeline(files: list) -> list:
    """Process every uploaded file independently, with a progress bar."""
    results = []
    progress = st.progress(0.0, text="Starting analysis...")
    for i, f in enumerate(files, start=1):
        progress.progress((i - 1) / len(files), text=f"Analyzing {f.name} ({i}/{len(files)})...")
        results.append(process_file(f.getvalue()))   # file bytes read once; result cached by content
    progress.progress(1.0, text="Analysis complete")
    progress.empty()
    return results


def main() -> None:
    st.set_page_config(page_title="Bioinformatics File Analyzer", page_icon="🧬", layout="wide")
    render_sidebar()
    st.title("🧬 Bioinformatics File Analyzer")
    st.caption("Upload FASTA, FASTQ, VCF or SAM files for instant statistics, interactive charts and downloadable reports.")

    files = st.file_uploader("Upload one or more files", accept_multiple_files=True,
                             help="Format is detected from file content. Plain text or .gz.")
    if not files:
        st.info("👆 Upload one or more files to get started. See the sidebar for supported formats and examples.")
        return

    results = run_pipeline(files)
    render_overview(files, results)
    st.divider()
    st.header("File analyses")
    for idx, (f, res) in enumerate(zip(files, results)):
        render_file_section(idx, f.name, f.size, res)


if __name__ == "__main__":
    main()
