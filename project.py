import streamlit as st

st.set_page_config(page_title="Bioinformatics File Format Analyzer")

st.title("🧬 Bioinformatics File Format Analyzer")

st.markdown(
    "<h3 style='color:red;'><u>Upload your bioinformatics file below</u></h3>",
    unsafe_allow_html=True
)

uploaded_file = st.file_uploader(
    "Choose a bioinformatics file",
    type=None
)

# ----------------------------------------------------
# Detect File Type
# ----------------------------------------------------
def detect_file_type(content):

    lines = content.splitlines()

    if not lines:
        return "Empty File"

    first_line = lines[0].strip()

    if first_line.startswith(">"):
        return "FASTA"

    if len(lines) >= 4:
        if lines[0].startswith("@") and lines[2].startswith("+"):
            return "FASTQ"

    if first_line.startswith("##") or first_line.startswith("#CHROM"):
        return "VCF"

    if first_line.startswith("@HD") or first_line.startswith("@SQ") \
       or first_line.startswith("@RG") or first_line.startswith("@PG") \
       or first_line.startswith("@CO"):
        return "SAM"

    return "Unknown"


# ----------------------------------------------------
# Parse FASTA File
# ----------------------------------------------------
def parse_fasta(content):

    sequences = []

    seq_id = ""
    sequence = ""

    for line in content.splitlines():

        line = line.strip()

        if not line:
            continue

        if line.startswith(">"):

            if seq_id != "":
                sequences.append((seq_id, sequence))

            seq_id = line[1:]
            sequence = ""

        else:
            sequence += line

    if seq_id != "":
        sequences.append((seq_id, sequence))

    return sequences


# ----------------------------------------------------
# GC Content
# ----------------------------------------------------
def gc_content(sequence):

    sequence = sequence.upper()

    if len(sequence) == 0:
        return 0

    gc = sequence.count("G") + sequence.count("C")

    return round((gc / len(sequence)) * 100, 2)


# ----------------------------------------------------
# Main Program
# ----------------------------------------------------
if uploaded_file is not None:

    content = uploaded_file.read().decode("utf-8", errors="ignore")

    file_type = detect_file_type(content)

    st.success("✅ File Uploaded Successfully")

    st.write("### File Information")

    st.write("**File Name:**", uploaded_file.name)
    st.write("**File Size:**", uploaded_file.size, "bytes")
    st.write("**Detected Type:**", file_type)

    st.write("**Lines:**", len(content.splitlines()))
    st.write("**Characters:**", len(content))

    # FASTA Analysis
    if file_type == "FASTA":

        st.header("🧬 FASTA Analysis")

        sequences = parse_fasta(content)

        for i, (seq_id, sequence) in enumerate(sequences, start=1):

            st.subheader(f"Sequence {i}")

            st.write("**Sequence ID:**", seq_id)
            st.write("**Sequence Length:**", len(sequence))
            st.write("**GC Content:**", f"{gc_content(sequence)} %")

            st.text_area(
                f"Sequence {i}",
                sequence,
                height=120
            )

    st.header("Raw File Content")

    st.text_area(
        "",
        content,
        height=300
    )

else:

    st.info("Please upload a bioinformatics file.")