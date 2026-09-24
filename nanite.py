import streamlit as st

# Set page configuration
st.set_page_config(
    page_title="Bioinformatics Multi-Format File Analyzer",
    page_icon="🧬",
    layout="wide"
)

# App Title & Description
st.title("🧬 Bioinformatics Multi-Format File Analyzer")
st.markdown("Upload your bioinformatics data files to inspect basic file details and preview content.")

st.divider()

# File Uploader
uploaded_file = st.file_uploader(
    "Choose a file",
    type=None,
    help="Upload any bioinformatics or plain text data file (e.g., .fasta, .fastq, .vcf, .sam, .bed, .txt, etc.)"
)

if uploaded_file is not None:
    # Read file content
    raw_bytes = uploaded_file.read()
    
    # Try decoding content as UTF-8
    try:
        content = raw_bytes.decode("utf-8")
    except UnicodeDecodeError:
        content = raw_bytes.decode("latin-1", errors="replace")

    # File basic statistics
    file_name = uploaded_file.name
    file_size_bytes = uploaded_file.size
    file_extension = file_name.split(".")[-1].upper() if "." in file_name else "Unknown"
    line_count = len(content.splitlines())
    char_count = len(content)

    # Human-readable file size format
    if file_size_bytes < 1024:
        formatted_size = f"{file_size_bytes} Bytes"
    elif file_size_bytes < 1024 * 1024:
        formatted_size = f"{file_size_bytes / 1024:.2f} KB"
    else:
        formatted_size = f"{file_size_bytes / (1024 * 1024):.2f} MB"

    st.success("File uploaded and read successfully!")

    # Display Basic Information
    st.subheader("📋 Basic File Information")
    
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric(label="File Name", value=file_name)
    with col2:
        st.metric(label="File Extension", value=f".{file_extension}" if file_extension != "Unknown" else "None")
    with col3:
        st.metric(label="File Size", value=formatted_size)
    with col4:
        st.metric(label="Total Lines", value=f"{line_count:,}")

    st.markdown(f"**Total Character Count:** {char_count:,} characters")

    # Content Display
    st.subheader("📄 File Content Preview")
    st.text_area(
        label="Raw Content",
        value=content,
        height=350,
        help="Scroll to view the complete content of the uploaded file."
    )

else:
    st.info("👆 Please upload a file to begin.")
