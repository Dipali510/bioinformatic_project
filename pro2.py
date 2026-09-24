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
    Read back with pandas.read_csv with the correct header parameter.
    """
    summary = f"# File type: {file_type}\n"
    for k, v in metrics.items():
        summary += f"# {k}: {v}\n"
    # Make sure index name is in header for CSV
    table.index.name = table.index.name or "index"
    return (summary + table.to_csv()).encode("utf-8", errors="replace")

def build_pdf(file_type: str, metrics: dict, table: pd.DataFrame, filename: str) -> io.BytesIO:
    """Build PDF report via reportlab, returning an in-memory binary buffer."""
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.styles import getSampleStyleSheet
    from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle
    from reportlab.lib import colors

    buffer = io.BytesIO()
    doc = SimpleDocTemplate(buffer, pagesize=A4)
    styles = getSampleStyleSheet()
    story = []

    # Title
    title = Paragraph(f"Bioinformatics File Studio – {file_type} Report", styles["Title"])
    story.append(title)
    story.append(Spacer(1, 12))

    # Hero banner-style metadata box
    md_text = f"""
    <div style="background-color:#0f1f38; border:1px solid #38bdf8; border-radius:12px; padding:16px;">
        <h3 style="color:#38bdf8; margin:0 0 12px 0;">File Analysis Summary</h3>
        <p style="color:#94a3b8; margin:0 0 6px 0;"><strong>File type:</strong> {file_type}</p>
        <p style="color:#94a3b8; margin:0 0 6px 0;"><strong>Total entries:</strong> {len(table) if isinstance(table.index, pd.RangeIndex) else len(table)}</p>
        <p style="color:#94a3b8; margin:0;"><strong>Generated on:</strong> {datetime.now().isoformat()}</p>
    </div>
    """
    story.append(Paragraph(md_text, styles["Normal"]))
    story.append(Spacer(1, 18))

    # Metrics section
    story.append(Paragraph("<b>Metrics</b>", styles["Heading3"]))
    story.append(Spacer(1, 10))
    metric_rows = [["Metric", "Value"]]
    for k, v in metrics.items():
        metric_rows.append([k, str(v)])
    m_table = Table(metric_rows, colWidths=[150, None])
    m_table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.rgba(56, 189, 248, 0.15)),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("ALIGN", (0, 0), (-1, -1), "LEFT"),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTNAME", (0, 1), (-1, -1), "Helvetica"),
        ("BOTTOMPADDING", (0, 0), (-1, 0), 12),
        ("BACKGROUND", (0, 1), (-1, -1), colors.rgba(15, 23, 42, 0.5)),
        ("GRID", (0, 0), (-1, -1), 1, colors.rgba(56, 189, 248, 0.3)),
        ("LEFTPADDING", (0, 0), (-1, -1), 8),
        ("RIGHTPADDING", (0, 0), (-1, -1), 8),
    ]))
    story.append(m_table)
    story.append(Spacer(1, 24))

    # Data table
    story.append(Paragraph("<b>Detailed Data</b>", styles["Heading3"]))
    story.append(Spacer(1, 10))
    d_rows = [table.columns.tolist()] + table.values.tolist()
    d_table = Table(d_rows, colWidths=[80] * len(table.columns))
    d_table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.rgba(56, 189, 248, 0.25)),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("ALIGN", (0, 0), (-1, -1), "LEFT"),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTNAME", (0, 1), (-1, -1), "Helvetica"),
        ("BOTTOMPADDING", (0, 0), (-1, 0), 8),
        ("BACKGROUND", (0, 1), (-1, -1), colors.rgba(15, 23, 42, 0.35)),
        ("GRID", (0, 0), (-1, -1), 1, colors.rgba(56, 189, 248, 0.25)),
    ]))
    story.append(d_table)

    doc.build(story)
    buffer.seek(0)
    return buffer

# --------------------------------------------------------------------------- #
# 5. Streamlit User Interface
# --------------------------------------------------------------------------- #
st.markdown('<div class="molecular-hero"><h1>🧬 Biomolecule Multi-Format Studio</h1><p>Analyze and visualize genetic sequences and alignment data (FASTA, FASTQ, VCF, SAM) with instant summaries and downloadable reports.</p></div>', unsafe_allow_html=True)

# Sidebar
with st.sidebar:
    st.markdown("### 📂 Upload File")
    uploaded_file = st.file_uploader(
        "Choose a file (FASTA, FASTQ, VCF, SAM)",
        type=["fasta", "fa", "fastq", "fq", "vcf", "sam"],
        help="Supported formats: FASTA, FASTQ, VCF, and SAM"
    )

    if uploaded_file:
        st.sidebar.markdown("**File uploaded successfully!**", unsafe_allow_html=True)

        # Peek at raw bytes for binary detection (first 8 KB)
        raw_bytes = uploaded_file.getvalue()[:8192]
        is_binary, bin_msg = is_binary_file(raw_bytes)

        if is_binary:
            st.error(f"❌ Binary format detected: {bin_msg}")
            st.info("Please use a standard bioinformatics tool (e.g., samtools, bgzip, etc.) to convert your file to a plain-text format before uploading.")
        else:
            # Read content as string (safe, as it's not binary)
            content = uploaded_file.getvalue().decode("utf-8", errors="replace")
            file_ext = uploaded_file.name.split(".")[-1].upper()

            # Detect format from content
            detected_format, reason = detect_file_format(content)

            if detected_format == "Empty File":
                st.error("⚠️ The uploaded file is empty or contains no valid records.")
            elif detected_format == "Unknown/Unsupported":
                st.error(f"❌ Unsupported format: {reason}. Please upload FASTA, FASTQ, VCF, or SAM.")
            else:
                # Use detected format, but allow user to override for ambiguous cases (e.g. .fq vs .fastq)
                fmt_label = detected_format
                if file_ext == "FQ" and detected_format == "FASTQ":
                    fmt_label = "FASTQ (.fq)"
                st.success(f"✅ Format detected: <b>{fmt_label}</b> ({reason})", unsafe_allow_html=True)
                st.markdown(f"**File name:** {