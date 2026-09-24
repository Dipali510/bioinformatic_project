import sys
import streamlit as st
import pandas as pd

# -----------------------------------------------------------------------------
# 1. Page Configuration & Molecular Theme Styling
# -----------------------------------------------------------------------------
st.set_page_config(
    page_title="Biomolecule Multi-Format Studio",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

def inject_molecular_theme():
    """Injects custom CSS for a modern molecular / biotech aesthetic."""
    st.markdown("""
    <style>
        /* Main background and base font */
        @import url('https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@400;600&family=Plus+Jakarta+Sans:wght@400;600;700;800&display=swap');

        html, body, [class*="css"] {
            font-family: 'Plus Jakarta Sans', -apple-system, BlinkMacSystemFont, sans-serif;
        }

        code, pre, [data-testid="stCode"] {
            font-family: 'JetBrains Mono', monospace !important;
        }

        /* Hero Banner Container */
        .molecular-hero {
            background: radial-gradient(circle at 10% 20%, rgba(0, 210, 255, 0.12) 0%, rgba(13, 27, 42, 0.95) 90%),
                        linear-gradient(135deg, #09111e 0%, #0f1f38 50%, #0b1526 100%);
            border: 1px solid rgba(0, 210, 255, 0.25);
            border-radius: 18px;
            padding: 26px 30px;
            margin-bottom: 25px;
            box-shadow: 0 10px 30px -5px rgba(0, 195, 255, 0.15);
        }

        .molecular-hero h1 {
            margin: 0;
            font-size: 2.2rem;
            font-weight: 800;
            letter-spacing: -0.5px;
        }

        .molecular-hero .gradient-text {
            background: linear-gradient(90deg, #00d2ff 0%, #3a7bd5 40%, #00f2fe 80%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
        }

        .molecular-hero p {
            color: #94a3b8;
            margin: 8px 0 0 0;
            font-size: 1.02rem;
        }

        /* Molecular Cards */
        .molecular-card {
            background: rgba(15, 23, 42, 0.65);
            border: 1px solid rgba(56, 189, 248, 0.2);
            border-radius: 14px;
            padding: 18px 22px;
            margin-bottom: 18px;
            backdrop-filter: blur(12px);
            transition: transform 0.2s ease, border-color 0.2s ease;
        }

        .molecular-card:hover {
            border-color: rgba(56, 189, 248, 0.45);
            transform: translateY(-2px);
        }

        /* Metric Cards Styling */
        div[data-testid="stMetric"] {
            background: rgba(15, 23, 42, 0.6);
            border: 1px solid rgba(56, 189, 248, 0.18);
            border-radius: 12px;
            padding: 12px 18px;
            box-shadow: 0 4px 14px rgba(0, 0, 0, 0.15);
        }

        div[data-testid="stMetricValue"] {
            color: #38bdf8 !important;
            font-weight: 700;
            font-size: 1.6rem !important;
        }

        div[data-testid="stMetricLabel"] {
            color: #94a3b8 !important;
            font-size: 0.88rem;
            font-weight: 600;
        }

        /* Pill Badge */
        .format-badge {
            display: inline-block;
            padding: 6px 14px;
            border-radius: 9999px;
            font-size: 0.95rem;
            font-weight: 700;
            letter-spacing: 0.5px;
            text-transform: uppercase;
        }

        .badge-fasta { background: rgba(16, 185, 129, 0.18); color: #34d399; border: 1px solid #10b981; }
        .badge-fastq { background: rgba(6, 182, 212, 0.18); color: #22d3ee; border: 1px solid #06b6d4; }
        .badge-vcf   { background: rgba(168, 85, 247, 0.18); color: #c084fc; border: 1px solid #a855f7; }
        .badge-sam   { background: rgba(245, 158, 11, 0.18); color: #fbbf24; border: 1px solid #f59e0b; }
        .badge-unsupported { background: rgba(239, 68, 68, 0.18); color: #f87171; border: 1px solid #ef4444; }
    </style>
    """, unsafe_allow_html=True)


# -----------------------------------------------------------------------------
# 2. Validation & Binary Inspection Engine
# -----------------------------------------------------------------------------
def is_binary_file(raw_bytes: bytes) -> tuple[bool, str]:
    """
    Checks if the uploaded file is binary or compressed (e.g., BAM, BCF, GZIP, ZIP).
    Returns (is_binary, description).
    """
    if not raw_bytes:
        return False, ""

    if raw_bytes.startswith(b"\x1f\x8b"):
        return True, "GZIP Compressed Archive (.gz). Please decompress the file to plain text before analysis."
    if raw_bytes.startswith(b"BAM\x01"):
        return True, "Binary Alignment Map (BAM) file detected. Please convert to plain text SAM format."
    if raw_bytes.startswith(b"BCF"):
        return True, "Binary Call Format (BCF) file detected. Please convert to plain text VCF format."
    if raw_bytes.startswith(b"PK\x03\x04"):
        return True, "ZIP Archive detected. Please extract the uncompressed data file."

    # General null byte check for arbitrary binary files
    sample = raw_bytes[:4096]
    if b"\x00" in sample:
        return True, "Binary / Non-text data detected. This studio accepts standard plain-text bioinformatics formats."

    return False, ""


# -----------------------------------------------------------------------------
# 3. Format Auto-Detection Engine
# -----------------------------------------------------------------------------
def detect_file_format(content: str) -> tuple[str, str]:
    """
    Detects whether the uploaded file content corresponds to FASTA, FASTQ, VCF, or SAM
    based on internal structural syntax.
    Returns: (format_name, reason)
    """
    if not content or not content.strip():
        return "Empty File", "The uploaded file contains 0 bytes or has no readable text lines."

    lines = [line.strip() for line in content.splitlines() if line.strip()]

    if not lines:
        return "Empty File", "The file contains only blank lines or whitespace."

    first_line = lines[0]

    # 1. VCF (Variant Call Format)
    if first_line.startswith("##fileformat=VCF") or (
        first_line.startswith("##") and any(l.startswith("#CHROM") for l in lines[:50])
    ):
        return (
            "VCF",
            "Identified VCF header starting with '##fileformat=VCF' / metadata tags followed by '#CHROM' variant columns."
        )

    # 2. SAM (Sequence Alignment/Map) with header
    sam_header_tags = ("@HD", "@SQ", "@RG", "@PG", "@CO")
    if any(first_line.startswith(tag) for tag in sam_header_tags):
        return (
            "SAM",
            f"Found standard SAM header record starting with '{first_line[:3]}'."
        )

    # 3. FASTQ (Fast Quality Format)
    if len(lines) >= 4:
        if lines[0].startswith("@") and lines[2].startswith("+"):
            seq_len = len(lines[1])
            qual_len = len(lines[3])
            if abs(seq_len - qual_len) <= 1:
                return (
                    "FASTQ",
                    "Detected 4-line repeating FASTQ read block: '@' header, bases, '+' separator, and matching quality scores."
                )

    # 4. FASTA
    if first_line.startswith(">"):
        return (
            "FASTA",
            "Header starts with '>' definition line followed by biological sequence data."
        )

    # 5. SAM without header (Headerless alignment records)
    fields = first_line.split("\t")
    if len(fields) >= 11:
        if fields[1].isdigit() and fields[3].isdigit() and fields[4].isdigit():
            return (
                "SAM",
                "Found tab-separated alignment records matching 11 standard SAM fields without headers."
            )

    return (
        "Unsupported Format",
        "The file content structure does not match standard FASTA, FASTQ, VCF, or SAM signatures."
    )


# -----------------------------------------------------------------------------
# 4. Bioinformatics Utility Functions
# -----------------------------------------------------------------------------
def get_sequence_length_distribution(lengths: list[int], num_bins: int = 10) -> pd.DataFrame:
    """Computes a binned sequence length distribution frequency table."""
    if not lengths:
        return pd.DataFrame(columns=["Length Interval", "Count"]).set_index("Length Interval")

    if min(lengths) == max(lengths):
        return pd.DataFrame([{"Length Interval": f"{lengths[0]} bp", "Count": len(lengths)}]).set_index("Length Interval")

    unique_vals = sorted(list(set(lengths)))
    if len(unique_vals) <= 8:
        counts = pd.Series(lengths).value_counts().sort_index()
        return pd.DataFrame({
            "Length Interval": [f"{val} bp" for val in counts.index],
            "Count": counts.values
        }).set_index("Length Interval")

    s = pd.Series(lengths)
    actual_bins = min(num_bins, len(unique_vals))
    binned = pd.cut(s, bins=actual_bins, right=True)
    counts = binned.value_counts().sort_index()

    formatted_data = []
    for interval, count in counts.items():
        formatted_data.append({
            "Length Interval": f"{int(interval.left)}-{int(interval.right)} bp",
            "Count": count
        })
    return pd.DataFrame(formatted_data).set_index("Length Interval")


def calculate_base_composition(sequences: list[str]) -> dict:
    """Counts nucleotide occurrences (A, T, C, G, Other)."""
    comp = {"A": 0, "T": 0, "C": 0, "G": 0, "Other": 0}
    for seq in sequences:
        for char in seq.upper():
            if char in comp:
                comp[char] += 1
            else:
                comp["Other"] += 1
    return comp


def natural_sort_chromosomes(chrom_counts: dict) -> pd.DataFrame:
    """Sorts chromosomes naturally (chr1, chr2, ... chr22, chrX, chrY, chrM)."""
    def chrom_key(chrom_name):
        c = str(chrom_name).lower().replace("chr", "")
        if c.isdigit():
            return (0, int(c))
        return (1, c)

    sorted_items = sorted(chrom_counts.items(), key=lambda x: chrom_key(x[0]))
    return pd.DataFrame(sorted_items, columns=["Chromosome", "Variant Count"]).set_index("Chromosome")


# -----------------------------------------------------------------------------
# 5. Modular Parsers
# -----------------------------------------------------------------------------

# --- FASTA Parser ---
def calculate_gc_content(sequence: str) -> float:
    """Calculates GC content percentage."""
    if not sequence:
        return 0.0
    seq_upper = sequence.upper()
    gc_count = seq_upper.count("G") + seq_upper.count("C")
    return round((gc_count / len(seq_upper)) * 100, 2)


def parse_fasta(content: str) -> tuple[list[dict], list[str]]:
    """Parses FASTA content into sequence records with validation."""
    records = []
    warnings = []
    current_id = ""
    current_desc = ""
    current_seq_parts = []
    line_num = 0

    for line in content.splitlines():
        line_num += 1
        line = line.strip()
        if not line:
            continue

        if line.startswith(">"):
            if current_id or current_seq_parts:
                full_seq = "".join(current_seq_parts).upper()
                if not full_seq:
                    warnings.append(f"Header '{current_id}' on line {line_num} had an empty sequence.")
                else:
                    records.append({
                        "id": current_id,
                        "description": current_desc,
                        "sequence": full_seq,
                        "length": len(full_seq),
                        "gc_content": calculate_gc_content(full_seq)
                    })

            header = line[1:].strip()
            current_desc = header
            current_id = header.split()[0] if header else "Unnamed_Seq"
            current_seq_parts = []
        else:
            if not current_id:
                warnings.append(f"Line {line_num} contains sequence data before any '>' header line.")
            current_seq_parts.append(line)

    if current_id or current_seq_parts:
        full_seq = "".join(current_seq_parts).upper()
        if not full_seq:
            warnings.append(f"Final header '{current_id}' had an empty sequence.")
        else:
            records.append({
                "id": current_id,
                "description": current_desc,
                "sequence": full_seq,
                "length": len(full_seq),
                "gc_content": calculate_gc_content(full_seq)
            })

    return records, warnings


# --- FASTQ Parser ---
def parse_fastq(content: str, phred_offset: int = 33) -> tuple[dict, list[str]]:
    """Parses FASTQ content and computes length and quality statistics."""
    lines = [line.strip() for line in content.splitlines() if line.strip()]
    reads = []
    warnings = []
    all_quality_scores = []
    lengths = []
    sequences = []

    if len(lines) % 4 != 0:
        warnings.append(f"Total line count ({len(lines)}) is not a multiple of 4. Some records may be incomplete.")

    for i in range(0, len(lines) - 3, 4):
        header = lines[i]
        seq = lines[i + 1]
        plus = lines[i + 2]
        qual_str = lines[i + 3]
        record_idx = (i // 4) + 1

        if not header.startswith("@"):
            warnings.append(f"Read #{record_idx}: Header does not begin with '@'.")
            continue

        if not plus.startswith("+"):
            warnings.append(f"Read #{record_idx}: Separator line does not begin with '+'.")
            continue

        if len(seq) != len(qual_str):
            warnings.append(
                f"Read #{record_idx}: Sequence length ({len(seq)}) does not match quality score length ({len(qual_str)})."
            )

        read_id = header[1:].split()[0] if len(header) > 1 else f"Read_{record_idx}"
        q_scores = [ord(char) - phred_offset for char in qual_str]

        if any(q < 0 or q > 93 for q in q_scores):
            warnings.append(f"Read #{record_idx}: Quality scores out of standard Phred range (Q0 - Q93).")

        all_quality_scores.extend(q_scores)
        seq_len = len(seq)
        lengths.append(seq_len)
        sequences.append(seq)
        avg_read_q = round(sum(q_scores) / len(q_scores), 2) if q_scores else 0.0

        reads.append({
            "id": read_id,
            "header": header[1:],
            "sequence": seq,
            "quality_string": qual_str,
            "length": seq_len,
            "avg_quality": avg_read_q
        })

    num_reads = len(reads)
    if num_reads == 0:
        return {
            "num_reads": 0,
            "reads": [],
            "lengths": [],
            "min_length": 0,
            "max_length": 0,
            "avg_length": 0.0,
            "total_bases": 0,
            "all_quality_scores": [],
            "base_composition": {"A": 0, "T": 0, "C": 0, "G": 0, "Other": 0},
            "quality_stats": {
                "mean_quality": 0.0,
                "min_quality": 0,
                "max_quality": 0,
                "q20_percentage": 0.0,
                "q30_percentage": 0.0,
            }
        }, warnings

    total_bases = sum(lengths)
    avg_len = round(total_bases / num_reads, 2)
    min_len = min(lengths)
    max_len = max(lengths)
    base_comp = calculate_base_composition(sequences)

    if all_quality_scores:
        mean_q = round(sum(all_quality_scores) / len(all_quality_scores), 2)
        min_q = min(all_quality_scores)
        max_q = max(all_quality_scores)
        q20_count = sum(1 for q in all_quality_scores if q >= 20)
        q30_count = sum(1 for q in all_quality_scores if q >= 30)
        q20_pct = round((q20_count / len(all_quality_scores)) * 100, 2)
        q30_pct = round((q30_count / len(all_quality_scores)) * 100, 2)
    else:
        mean_q, min_q, max_q, q20_pct, q30_pct = 0.0, 0, 0, 0.0, 0.0

    return {
        "num_reads": num_reads,
        "reads": reads,
        "lengths": lengths,
        "min_length": min_len,
        "max_length": max_len,
        "avg_length": avg_len,
        "total_bases": total_bases,
        "all_quality_scores": all_quality_scores,
        "base_composition": base_comp,
        "quality_stats": {
            "mean_quality": mean_q,
            "min_quality": min_q,
            "max_quality": max_q,
            "q20_percentage": q20_pct,
            "q30_percentage": q30_pct,
        }
    }, warnings


# --- VCF Parser ---
def parse_vcf(content: str) -> tuple[dict, list[str]]:
    """Parses VCF content and calculates variant statistics."""
    meta_lines = []
    header_cols = []
    raw_records = []
    warnings = []
    line_num = 0

    for line in content.splitlines():
        line_num += 1
        line = line.strip()
        if not line:
            continue
        if line.startswith("##"):
            meta_lines.append(line)
        elif line.startswith("#CHROM"):
            header_cols = line.split("\t")
        elif header_cols and not line.startswith("#"):
            fields = line.split("\t")
            if len(fields) < 8:
                warnings.append(f"Line {line_num}: Malformed variant record with only {len(fields)} columns (min 8 required).")
                continue
            raw_records.append(fields)

    if not header_cols and raw_records:
        warnings.append("Missing standard '#CHROM' column header line in VCF file.")

    total_variants = len(raw_records)
    samples = header_cols[9:] if len(header_cols) > 9 else []

    snps = 0
    insertions = 0
    deletions = 0
    others = 0
    pass_count = 0
    filter_counts = {}
    chrom_counts = {}
    parsed_rows = []

    for fields in raw_records:
        chrom = fields[0]
        pos = fields[1]
        var_id = fields[2]
        ref = fields[3].upper()
        alt = fields[4].upper()
        qual = fields[5]
        filter_status = fields[6]

        chrom_counts[chrom] = chrom_counts.get(chrom, 0) + 1
        filter_counts[filter_status] = filter_counts.get(filter_status, 0) + 1
        if filter_status.upper() == "PASS":
            pass_count += 1

        if len(ref) == 1 and len(alt) == 1 and ref in "ACGT" and alt in "ACGT":
            snps += 1
            vtype = "SNP"
        elif len(ref) < len(alt):
            insertions += 1
            vtype = "Insertion"
        elif len(ref) > len(alt):
            deletions += 1
            vtype = "Deletion"
        else:
            others += 1
            vtype = "Complex/Other"

        parsed_rows.append({
            "Chrom": chrom,
            "Pos": pos,
            "ID": var_id,
            "Ref": ref,
            "Alt": alt,
            "Type": vtype,
            "Qual": qual,
            "Filter": filter_status
        })

    pass_pct = round((pass_count / total_variants) * 100, 2) if total_variants > 0 else 0.0

    return {
        "total_variants": total_variants,
        "meta_lines_count": len(meta_lines),
        "samples": samples,
        "snps": snps,
        "insertions": insertions,
        "deletions": deletions,
        "others": others,
        "pass_count": pass_count,
        "pass_pct": pass_pct,
        "filter_counts": filter_counts,
        "chrom_counts": chrom_counts,
        "records": parsed_rows
    }, warnings


# --- SAM Parser ---
def parse_sam(content: str) -> tuple[dict, list[str]]:
    """Parses SAM content and calculates alignment statistics."""
    header_lines = []
    alignments = []
    mapq_scores = []
    warnings = []
    mapped_count = 0
    unmapped_count = 0
    forward_strand = 0
    reverse_strand = 0
    ref_counts = {}
    line_num = 0

    for line in content.splitlines():
        line_num += 1
        line = line.strip()
        if not line:
            continue
        if line.startswith("@"):
            header_lines.append(line)
        else:
            fields = line.split("\t")
            if len(fields) < 11:
                warnings.append(f"Line {line_num}: Incomplete alignment record with {len(fields)} fields (expected ≥11).")
                continue

            qname = fields[0]
            try:
                flag = int(fields[1])
            except ValueError:
                flag = 0
                warnings.append(f"Line {line_num}: Non-integer FLAG '{fields[1]}'. Defaulted to 0.")

            rname = fields[2]
            pos = fields[3]
            try:
                mapq = int(fields[4])
                mapq_scores.append(mapq)
            except ValueError:
                mapq = 0
                warnings.append(f"Line {line_num}: Non-integer MAPQ '{fields[4]}'. Defaulted to 0.")

            cigar = fields[5]

            is_unmapped = bool(flag & 0x4)
            if is_unmapped or rname == "*":
                unmapped_count += 1
            else:
                mapped_count += 1
                ref_counts[rname] = ref_counts.get(rname, 0) + 1

            if flag & 0x10:
                reverse_strand += 1
            else:
                forward_strand += 1

            alignments.append({
                "QNAME": qname,
                "FLAG": flag,
                "RNAME": rname,
                "POS": pos,
                "MAPQ": mapq,
                "CIGAR": cigar,
                "Status": "Unmapped" if is_unmapped else "Mapped"
            })

    total_alignments = len(alignments)
    avg_mapq = round(sum(mapq_scores) / len(mapq_scores), 2) if mapq_scores else 0.0
    mapping_rate = round((mapped_count / total_alignments) * 100, 2) if total_alignments > 0 else 0.0

    return {
        "total_alignments": total_alignments,
        "header_lines_count": len(header_lines),
        "mapped_count": mapped_count,
        "unmapped_count": unmapped_count,
        "mapping_rate": mapping_rate,
        "avg_mapq": avg_mapq,
        "forward_strand": forward_strand,
        "reverse_strand": reverse_strand,
        "ref_counts": ref_counts,
        "mapq_scores": mapq_scores,
        "alignments": alignments
    }, warnings


# -----------------------------------------------------------------------------
# 6. Molecular Dashboard Renderers
# -----------------------------------------------------------------------------
def display_validation_warnings(warnings: list[str]):
    """Renders validation alerts if data anomalies occurred."""
    if warnings:
        with st.expander(f"⚠️ Structural Warnings ({len(warnings)})", expanded=False):
            st.warning("The following structural inconsistencies were identified in the file:")
            for w in warnings[:15]:
                st.markdown(f"- `{w}`")
            if len(warnings) > 15:
                st.caption(f"... and {len(warnings) - 15} additional warnings.")


# --- Render FASTA ---
def render_fasta_dashboard(content: str):
    records, warnings = parse_fasta(content)
    display_validation_warnings(warnings)

    num_records = len(records)
    if num_records == 0:
        st.error("❌ No valid FASTA sequences could be parsed.")
        return

    lengths = [r["length"] for r in records]
    total_bases = sum(lengths)
    avg_length = round(total_bases / num_records, 2)
    avg_gc = round(sum(r["gc_content"] for r in records) / num_records, 2)

    tab1, tab2, tab3 = st.tabs(["📊 Genomic Metrics", "📈 Molecular Distributions", "📋 Sequence Inspector"])

    with tab1:
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Total Sequences", f"{num_records:,}")
        c2.metric("Total Base Pairs", f"{total_bases:,} bp")
        c3.metric("Average Length", f"{avg_length:,} bp")
        c4.metric("Mean GC Content", f"{avg_gc}%")

        c5, c6 = st.columns(2)
        c5.metric("Min Sequence Length", f"{min(lengths):,} bp")
        c6.metric("Max Sequence Length", f"{max(lengths):,} bp")

    with tab2:
        g1, g2 = st.columns(2)
        with g1:
            st.markdown("##### 📏 Sequence Length Distribution")
            df_len = get_sequence_length_distribution(lengths)
            st.bar_chart(df_len)

        with g2:
            st.markdown("##### 🧬 Nucleotide Base Frequencies (A, T, C, G)")
            base_comp = calculate_base_composition([r["sequence"] for r in records])
            df_bases = pd.DataFrame(list(base_comp.items()), columns=["Base", "Count"]).set_index("Base")
            st.bar_chart(df_bases)

    with tab3:
        st.markdown("##### 📋 Parsed Sequences Overview")
        df_records = pd.DataFrame([
            {
                "Sequence ID": r["id"],
                "Length (bp)": r["length"],
                "GC (%)": f"{r['gc_content']}%",
                "Description": r["description"]
            }
            for r in records
        ])
        st.dataframe(df_records, use_container_width=True)

        st.markdown("##### 🔬 Detailed Sequence Inspection")
        for i, rec in enumerate(records[:25], start=1):
            with st.expander(f"🧬 #{i} • {rec['id']} ({rec['length']:,} bp | GC: {rec['gc_content']}%)"):
                st.markdown(f"**Header:** `{rec['description']}`")
                st.text_area(f"Sequence {i}", value=rec["sequence"], height=100, key=f"fa_seq_{i}")

        if num_records > 25:
            st.caption(f"Showing 25 of {num_records:,} sequences.")


# --- Render FASTQ ---
def render_fastq_dashboard(content: str):
    fastq_data, warnings = parse_fastq(content)
    display_validation_warnings(warnings)

    num_reads = fastq_data["num_reads"]
    if num_reads == 0:
        st.error("❌ No valid FASTQ reads could be parsed.")
        return

    q_stats = fastq_data["quality_stats"]
    lengths = fastq_data["lengths"]

    tab1, tab2, tab3 = st.tabs(["📊 Read Quality & Stats", "📈 Quality & Length Distributions", "📋 Read Inspector"])

    with tab1:
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Total Reads", f"{num_reads:,}")
        c2.metric("Total Bases", f"{fastq_data['total_bases']:,} bp")
        c3.metric("Avg Read Length", f"{fastq_data['avg_length']} bp")
        c4.metric("Length Range", f"{fastq_data['min_length']}-{fastq_data['max_length']} bp")

        q1, q2, q3, q4 = st.columns(4)
        q1.metric("Mean Quality (Q)", f"{q_stats['mean_quality']}")
        q2.metric("Quality Range", f"Q{q_stats['min_quality']} - Q{q_stats['max_quality']}")
        q3.metric("Bases ≥ Q20 (99% Acc)", f"{q_stats['q20_percentage']}%")
        q4.metric("Bases ≥ Q30 (99.9% Acc)", f"{q_stats['q30_percentage']}%")

    with tab2:
        g1, g2 = st.columns(2)
        with g1:
            st.markdown("##### 📏 Read Length Distribution")
            df_len = get_sequence_length_distribution(lengths)
            st.bar_chart(df_len)

        with g2:
            st.markdown("##### 🎯 Quality Score Tier Distribution")
            all_q = fastq_data["all_quality_scores"]
            q_tiers = {
                "High (Q ≥ 30)": sum(1 for q in all_q if q >= 30),
                "Moderate (Q 20-29)": sum(1 for q in all_q if 20 <= q < 30),
                "Low (Q < 20)": sum(1 for q in all_q if q < 20)
            }
            df_q = pd.DataFrame(list(q_tiers.items()), columns=["Quality Tier", "Base Count"]).set_index("Quality Tier")
            st.bar_chart(df_q)

    with tab3:
        st.markdown("##### 📋 Reads Overview Table")
        df_reads = pd.DataFrame([
            {
                "Read ID": r["id"],
                "Length (bp)": r["length"],
                "Avg Q": r["avg_quality"],
                "Quality Preview": r["quality_string"][:25] + ("..." if len(r["quality_string"]) > 25 else "")
            }
            for r in fastq_data["reads"]
        ])
        st.dataframe(df_reads, use_container_width=True)

        st.markdown("##### 🔬 Individual Read Inspector")
        for i, r in enumerate(fastq_data["reads"][:25], start=1):
            with st.expander(f"🔬 Read #{i} • {r['id']} ({r['length']} bp | Avg Q: {r['avg_quality']})"):
                st.text_input(f"Sequence {i}", value=r["sequence"], key=f"fq_s_{i}")
                st.text_input(f"Quality String {i}", value=r["quality_string"], key=f"fq_q_{i}")

        if num_reads > 25:
            st.caption(f"Showing 25 of {num_reads:,} reads.")


# --- Render VCF ---
def render_vcf_dashboard(content: str):
    vcf_data, warnings = parse_vcf(content)
    display_validation_warnings(warnings)

    total_vars = vcf_data["total_variants"]
    if total_vars == 0:
        st.error("❌ No variant records found in this VCF file.")
        return

    tab1, tab2, tab3 = st.tabs(["📊 Variant Metrics", "📈 Variant Count Analytics", "📋 Variant Call Table"])

    with tab1:
        v1, v2, v3, v4 = st.columns(4)
        v1.metric("Total Variants", f"{total_vars:,}")
        v2.metric("SNPs", f"{vcf_data['snps']:,}")
        v3.metric("Insertions", f"{vcf_data['insertions']:,}")
        v4.metric("Deletions", f"{vcf_data['deletions']:,}")

        v5, v6, v7, v8 = st.columns(4)
        v5.metric("Complex / Others", f"{vcf_data['others']:,}")
        v6.metric("PASS Rate", f"{vcf_data['pass_pct']}%")
        v7.metric("Sample Count", f"{len(vcf_data['samples']):,}")
        v8.metric("Meta Lines", f"{vcf_data['meta_lines_count']:,}")

    with tab2:
        g1, g2 = st.columns(2)
        with g1:
            st.markdown("##### 🧬 Variant Counts by Mutation Type")
            vtypes = {
                "SNPs": vcf_data["snps"],
                "Insertions": vcf_data["insertions"],
                "Deletions": vcf_data["deletions"],
                "Complex/Other": vcf_data["others"]
            }
            df_vtypes = pd.DataFrame(list(vtypes.items()), columns=["Type", "Count"]).set_index("Type")
            st.bar_chart(df_vtypes)

        with g2:
            st.markdown("##### 🗺️ Genomic Distribution (Variants per Chromosome)")
            if vcf_data["chrom_counts"]:
                df_chrom = natural_sort_chromosomes(vcf_data["chrom_counts"])
                st.bar_chart(df_chrom)
            else:
                st.info("No chromosome data available.")

    with tab3:
        st.markdown("##### 📋 Variant Calls Table")
        df_vars = pd.DataFrame(vcf_data["records"])
        st.dataframe(df_vars, use_container_width=True)

        if vcf_data["samples"]:
            st.info(f"**Identified Genotyped Samples:** {', '.join(vcf_data['samples'])}")


# --- Render SAM ---
def render_sam_dashboard(content: str):
    sam_data, warnings = parse_sam(content)
    display_validation_warnings(warnings)

    total_alns = sam_data["total_alignments"]
    if total_alns == 0:
        st.error("❌ No alignment records found in this SAM file.")
        return

    tab1, tab2, tab3 = st.tabs(["📊 Alignment Metrics", "📈 Alignment Distribution", "📋 Alignment Table"])

    with tab1:
        s1, s2, s3, s4 = st.columns(4)
        s1.metric("Total Alignments", f"{total_alns:,}")
        s2.metric("Mapped Reads", f"{sam_data['mapped_count']:,}")
        s3.metric("Unmapped Reads", f"{sam_data['unmapped_count']:,}")
        s4.metric("Mapping Rate", f"{sam_data['mapping_rate']}%")

        s5, s6, s7, s8 = st.columns(4)
        s5.metric("Mean MAPQ Score", f"{sam_data['avg_mapq']}")
        s6.metric("Forward Strand", f"{sam_data['forward_strand']:,}")
        s7.metric("Reverse Strand", f"{sam_data['reverse_strand']:,}")
        s8.metric("SAM Header Lines", f"{sam_data['header_lines_count']:,}")

    with tab2:
        g1, g2 = st.columns(2)
        with g1:
            st.markdown("##### 🎯 Mapping & Strand Breakdown")
            mapping_breakdown = {
                "Mapped Reads": sam_data["mapped_count"],
                "Unmapped Reads": sam_data["unmapped_count"],
                "Forward Strand": sam_data["forward_strand"],
                "Reverse Strand": sam_data["reverse_strand"]
            }
            df_map = pd.DataFrame(list(mapping_breakdown.items()), columns=["Category", "Count"]).set_index("Category")
            st.bar_chart(df_map)

        with g2:
            st.markdown("##### 📌 Top Reference Contigs Alignment Counts")
            if sam_data["ref_counts"]:
                df_ref = pd.DataFrame(
                    list(sam_data["ref_counts"].items()),
                    columns=["Contig", "Count"]
                ).sort_values(by="Count", ascending=False).head(10).set_index("Contig")
                st.bar_chart(df_ref)
            else:
                st.info("No mapped reference contig records available.")

    with tab3:
        st.markdown("##### 📋 Alignment Records Table")
        df_alns = pd.DataFrame(sam_data["alignments"])
        st.dataframe(df_alns, use_container_width=True)


# -----------------------------------------------------------------------------
# 7. Main Application
# -----------------------------------------------------------------------------
def main():
    inject_molecular_theme()

    # Molecular Hero Header
    st.markdown("""
    <div class="molecular-hero">
        <h1 style="color: #ffffff;"><span class="gradient-text">🧬 Biomolecule Multi-Format Studio</span></h1>
        <p>Unified computational workbench for FASTA, FASTQ, VCF, and SAM high-throughput genomic data.</p>
    </div>
    """, unsafe_allow_html=True)

    # Sidebar: File Upload & Guidance
    with st.sidebar:
        st.header("🔬 Input Workspace")
        st.markdown("Upload any standard genomics file to automatically classify its format and extract structural intelligence.")

        uploaded_file = st.file_uploader(
            "Choose a genomics file",
            type=None,
            help="Supported formats: FASTA (.fa, .fasta), FASTQ (.fq, .fastq), VCF (.vcf), SAM (.sam)"
        )

        st.markdown("---")
        st.markdown("### 🧬 Supported Formats")
        st.markdown("""
        - **FASTA**: Biological nucleotide/protein sequences
        - **FASTQ**: Sequencing reads with Phred quality scores
        - **VCF**: Genomic variant calls and annotations
        - **SAM**: Sequence Alignment/Map coordinates
        """)

    if uploaded_file is not None:
        raw_bytes = uploaded_file.read()

        # 1. Empty file validation
        if not raw_bytes or len(raw_bytes) == 0:
            st.error("⚠️ **Empty File Detected**: The uploaded file contains 0 bytes. Please upload a valid genomics file.")
            return

        # 2. Binary / Compressed file detection
        binary_detected, binary_reason = is_binary_file(raw_bytes)
        if binary_detected:
            st.error(f"⚠️ **Incompatible Binary Format**: {binary_reason}")
            st.info("💡 **Note**: Please supply plain-text FASTA, FASTQ, VCF, or SAM files. (e.g. convert BAM to SAM).")
            return

        # Decode content safely
        try:
            content = raw_bytes.decode("utf-8")
        except UnicodeDecodeError:
            content = raw_bytes.decode("latin-1", errors="replace")

        # 3. Empty text validation
        if not content.strip():
            st.warning("⚠️ **Blank File**: The uploaded file contains only whitespace or blank lines.")
            return

        # Format detection
        detected_format, detection_reason = detect_file_format(content)

        # File general metrics
        file_name = uploaded_file.name
        file_size_bytes = uploaded_file.size
        file_extension = file_name.split(".")[-1].upper() if "." in file_name else "NO EXT"
        line_count = len(content.splitlines())
        char_count = len(content)

        if file_size_bytes < 1024:
            formatted_size = f"{file_size_bytes} Bytes"
        elif file_size_bytes < 1024 * 1024:
            formatted_size = f"{file_size_bytes / 1024:.2f} KB"
        else:
            formatted_size = f"{file_size_bytes / (1024 * 1024):.2f} MB"

        badge_classes = {
            "FASTA": "badge-fasta",
            "FASTQ": "badge-fastq",
            "VCF": "badge-vcf",
            "SAM": "badge-sam",
            "Unsupported Format": "badge-unsupported"
        }
        badge_class = badge_classes.get(detected_format, "badge-unsupported")

        # Top Metadata Banner
        st.markdown(f"""
        <div class="molecular-card">
            <div style="display: flex; justify-content: space-between; align-items: center; flex-wrap: wrap; gap: 10px;">
                <div>
                    <span style="color: #94a3b8; font-size: 0.85rem; font-weight: 600;">CLASSIFIED FORMAT</span><br>
                    <span class="format-badge {badge_class}">{detected_format}</span>
                </div>
                <div style="text-align: right;">
                    <span style="color: #94a3b8; font-size: 0.85rem; font-weight: 600;">ACTIVE FILE</span><br>
                    <span style="color: #f8fafc; font-weight: 700; font-size: 1.1rem;">{file_name}</span>
                </div>
            </div>
            <div style="margin-top: 12px; color: #cbd5e1; font-size: 0.93rem;">
                <strong>Detection Logic:</strong> {detection_reason}
            </div>
        </div>
        """, unsafe_allow_html=True)

        # General file statistics cards
        m1, m2, m3, m4 = st.columns(4)
        m1.metric("File Name", file_name)
        m2.metric("Extension", f".{file_extension}" if file_extension != "NO EXT" else "None")
        m3.metric("Size", formatted_size)
        m4.metric("Total Lines", f"{line_count:,}")

        # Dynamic Format Dashboard Dispatch
        if detected_format == "FASTA":
            render_fasta_dashboard(content)
        elif detected_format == "FASTQ":
            render_fastq_dashboard(content)
        elif detected_format == "VCF":
            render_vcf_dashboard(content)
        elif detected_format == "SAM":
            render_sam_dashboard(content)
        else:
            st.error("⚠️ **Unsupported Genomics Format**")
            st.warning("""
            The uploaded file structure did not match any of the supported formats:
            - **FASTA**: Must start with `>` followed by sequence letters.
            - **FASTQ**: Must have 4-line read blocks starting with `@` and matching quality strings.
            - **VCF**: Must have `##fileformat=VCF` metadata and `#CHROM` variant table columns.
            - **SAM**: Must have `@HD`/`@SQ` headers or 11 standard tab-delimited alignment columns.
            """)

        # Raw Stream Expander
        st.divider()
        with st.expander("📄 Raw File Data Stream", expanded=(detected_format == "Unsupported Format")):
            st.text_area("Raw Content", value=content, height=280, help="Complete raw text representation.")

    else:
        st.info("👆 Please upload a bioinformatics file from the sidebar or click browse above to begin.")


if __name__ == "__main__":
    if not (hasattr(st, "runtime") and st.runtime.exists()):
        from streamlit.web import cli as stcli
        sys.argv = ["streamlit", "run", __file__]
        sys.exit(stcli.main())
    else:
        main()
