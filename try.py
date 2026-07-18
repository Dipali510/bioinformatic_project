#!/usr/bin/env python3
"""
fasta_analyzer.py
==================
A production-quality CLI tool for analyzing FASTA files.

Features:
    - Accepts single files, multiple files, or directories (recursive scan)
    - Streams large files instead of loading them fully into memory
    - Detects DNA vs Protein sequences automatically (or via override)
    - Computes GC content (DNA) or amino acid composition (Protein)
    - Validates sequence characters and reports anomalies
    - Prints a clean formatted table and optionally exports to CSV

Usage:
    python fasta_analyzer.py file1.fasta
    python fasta_analyzer.py file1.fa file2.faa
    python fasta_analyzer.py ./data_dir --type DNA
    python fasta_analyzer.py ./data_dir --csv results.csv
"""

import argparse
import csv
import os
import sys
from collections import Counter
from dataclasses import dataclass, field
from typing import Iterator, List, Optional


# --------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------

VALID_FASTA_EXTENSIONS = {".fasta", ".fa", ".fna", ".faa"}

DNA_VALID_CHARS = set("ATGCN")
PROTEIN_VALID_CHARS = set("ACDEFGHIKLMNPQRSTVWYX")  # 20 amino acids + X (unknown)


# --------------------------------------------------------------------------
# Data model
# --------------------------------------------------------------------------

@dataclass
class SequenceRecord:
    """Holds parsed information and computed stats for a single FASTA entry."""
    file_name: str
    seq_id: str
    description: str
    sequence: str
    seq_type: str = ""
    length: int = 0
    gc_content: Optional[float] = None
    aa_composition: dict = field(default_factory=dict)
    invalid_chars: set = field(default_factory=set)


# --------------------------------------------------------------------------
# File collection
# --------------------------------------------------------------------------

def collect_fasta_files(paths: List[str]) -> List[str]:
    """
    Given a list of file/directory paths, return a de-duplicated list of
    valid FASTA file paths. Directories are scanned recursively.
    """
    collected = []

    for path in paths:
        if not os.path.exists(path):
            print(f"[WARNING] Path does not exist, skipping: {path}", file=sys.stderr)
            continue

        if os.path.isfile(path):
            ext = os.path.splitext(path)[1].lower()
            if ext in VALID_FASTA_EXTENSIONS:
                collected.append(path)
            else:
                print(f"[WARNING] Not a recognized FASTA extension, skipping: {path}",
                      file=sys.stderr)

        elif os.path.isdir(path):
            for root, _dirs, files in os.walk(path):
                for fname in files:
                    ext = os.path.splitext(fname)[1].lower()
                    if ext in VALID_FASTA_EXTENSIONS:
                        collected.append(os.path.join(root, fname))

    # Remove duplicates while preserving order
    seen = set()
    unique_files = []
    for f in collected:
        real = os.path.realpath(f)
        if real not in seen:
            seen.add(real)
            unique_files.append(f)

    return unique_files


# --------------------------------------------------------------------------
# FASTA parsing (streaming, memory-conscious)
# --------------------------------------------------------------------------

def read_fasta(file_path: str) -> Iterator[SequenceRecord]:
    """
    Stream-parse a FASTA file, yielding one SequenceRecord per entry.
    Raises ValueError on malformed FASTA content.
    """
    seq_id = None
    description = ""
    seq_chunks: List[str] = []
    found_any_header = False
    file_name = os.path.basename(file_path)

    def _flush():
        """Build a SequenceRecord from accumulated chunks."""
        sequence = "".join(seq_chunks)
        return SequenceRecord(
            file_name=file_name,
            seq_id=seq_id,
            description=description,
            sequence=sequence,
        )

    with open(file_path, "r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()

            if not line:
                continue  # ignore empty lines

            if line.startswith(">"):
                found_any_header = True
                # Yield the previous record before starting a new one
                if seq_id is not None:
                    if not seq_chunks:
                        raise ValueError(
                            f"Malformed FASTA in {file_name}: header '{seq_id}' "
                            f"has no following sequence."
                        )
                    yield _flush()

                header_line = line[1:].strip()
                if not header_line:
                    raise ValueError(f"Malformed FASTA in {file_name}: empty header line.")

                tokens = header_line.split(maxsplit=1)
                seq_id = tokens[0]
                description = tokens[1] if len(tokens) > 1 else ""
                seq_chunks = []

            else:
                if seq_id is None:
                    raise ValueError(
                        f"Malformed FASTA in {file_name}: sequence data found "
                        f"before any header ('>')."
                    )
                seq_chunks.append(line)

    if not found_any_header:
        raise ValueError(f"No valid FASTA headers found in {file_name}. Empty or invalid file.")

    # Yield the final record
    if seq_id is not None:
        if not seq_chunks:
            raise ValueError(
                f"Malformed FASTA in {file_name}: header '{seq_id}' has no sequence."
            )
        yield _flush()


# --------------------------------------------------------------------------
# Sequence processing helpers
# --------------------------------------------------------------------------

def clean_sequence(sequence: str) -> str:
    """Uppercase the sequence and strip any stray whitespace."""
    return "".join(sequence.split()).upper()


def detect_sequence_type(sequence: str) -> str:
    """
    Auto-detect whether a sequence is DNA or PROTEIN based on its
    character composition. DNA sequences consist overwhelmingly of
    A, T, G, C, N; anything with a significantly wider alphabet is
    treated as protein.
    """
    if not sequence:
        return "UNKNOWN"

    seq_chars = set(sequence)
    dna_like_chars = seq_chars & DNA_VALID_CHARS
    non_dna_chars = seq_chars - DNA_VALID_CHARS

    # If (almost) every character belongs to the DNA alphabet, call it DNA
    dna_fraction = sum(1 for c in sequence if c in DNA_VALID_CHARS) / len(sequence)

    if dna_fraction >= 0.9 and len(non_dna_chars) <= 1:
        return "DNA"
    return "PROTEIN"


def find_invalid_chars(sequence: str, seq_type: str) -> set:
    """Return the set of characters not valid for the given sequence type."""
    valid_set = DNA_VALID_CHARS if seq_type == "DNA" else PROTEIN_VALID_CHARS
    return set(sequence) - valid_set


def calculate_gc(sequence: str) -> float:
    """Calculate GC content percentage for a DNA sequence."""
    if not sequence:
        return 0.0
    gc_count = sequence.count("G") + sequence.count("C")
    return round((gc_count / len(sequence)) * 100, 2)


def amino_acid_composition(sequence: str) -> dict:
    """Return amino acid composition as percentages, rounded to 2 decimals."""
    if not sequence:
        return {}
    counts = Counter(sequence)
    total = len(sequence)
    return {aa: round((count / total) * 100, 2) for aa, count in sorted(counts.items())}


# --------------------------------------------------------------------------
# Main analysis pipeline
# --------------------------------------------------------------------------

def generate_stats(record: SequenceRecord, forced_type: Optional[str]) -> SequenceRecord:
    """
    Populate a SequenceRecord with cleaned sequence, detected/forced type,
    length, GC content or amino acid composition, and invalid character info.
    """
    cleaned = clean_sequence(record.sequence)
    record.sequence = cleaned
    record.length = len(cleaned)

    record.seq_type = forced_type if forced_type else detect_sequence_type(cleaned)
    record.invalid_chars = find_invalid_chars(cleaned, record.seq_type)

    if record.seq_type == "DNA":
        record.gc_content = calculate_gc(cleaned)
    else:
        record.aa_composition = amino_acid_composition(cleaned)

    return record


def process_files(file_paths: List[str], forced_type: Optional[str]) -> List[SequenceRecord]:
    """Parse and analyze all given FASTA files, collecting results and reporting errors."""
    all_records: List[SequenceRecord] = []

    for path in file_paths:
        try:
            if os.path.getsize(path) == 0:
                print(f"[WARNING] Skipping empty file: {path}", file=sys.stderr)
                continue

            for raw_record in read_fasta(path):
                try:
                    stats_record = generate_stats(raw_record, forced_type)
                    all_records.append(stats_record)

                    if stats_record.invalid_chars:
                        print(
                            f"[WARNING] {stats_record.file_name} :: {stats_record.seq_id} "
                            f"contains invalid characters for {stats_record.seq_type}: "
                            f"{sorted(stats_record.invalid_chars)}",
                            file=sys.stderr,
                        )
                except Exception as inner_err:
                    print(f"[ERROR] Failed processing a record in {path}: {inner_err}",
                          file=sys.stderr)

        except ValueError as fasta_err:
            print(f"[ERROR] {fasta_err}", file=sys.stderr)
        except FileNotFoundError:
            print(f"[ERROR] File not found: {path}", file=sys.stderr)
        except Exception as unexpected_err:
            print(f"[ERROR] Unexpected error reading {path}: {unexpected_err}", file=sys.stderr)

    return all_records


# --------------------------------------------------------------------------
# Output formatting
# --------------------------------------------------------------------------

def print_results_table(records: List[SequenceRecord]) -> None:
    """Print a clean, aligned table of results to the console."""
    if not records:
        print("No sequences to display.")
        return

    header = f"{'File':<20} {'Seq ID':<20} {'Type':<10} {'Length':<10} {'GC %':<8}"
    separator = "-" * len(header)
    print(separator)
    print(header)
    print(separator)

    for rec in records:
        gc_display = f"{rec.gc_content:.2f}" if rec.gc_content is not None else "N/A"
        print(f"{rec.file_name:<20} {rec.seq_id:<20} {rec.seq_type:<10} "
              f"{rec.length:<10} {gc_display:<8}")

    print(separator)


def print_summary(records: List[SequenceRecord]) -> None:
    """Print aggregate statistics: total sequences, DNA vs protein counts."""
    total = len(records)
    dna_count = sum(1 for r in records if r.seq_type == "DNA")
    protein_count = sum(1 for r in records if r.seq_type == "PROTEIN")

    print("\nSummary")
    print("-------")
    print(f"Total sequences   : {total}")
    print(f"DNA sequences     : {dna_count}")
    print(f"Protein sequences : {protein_count}")


def export_csv(records: List[SequenceRecord], output_path: str) -> None:
    """Export all processed records into a single CSV file."""
    fieldnames = ["file_name", "seq_id", "description", "type", "length",
                  "gc_content", "amino_acid_composition", "invalid_chars"]

    try:
        with open(output_path, "w", newline="", encoding="utf-8") as csv_file:
            writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
            writer.writeheader()

            for rec in records:
                writer.writerow({
                    "file_name": rec.file_name,
                    "seq_id": rec.seq_id,
                    "description": rec.description,
                    "type": rec.seq_type,
                    "length": rec.length,
                    "gc_content": rec.gc_content if rec.gc_content is not None else "",
                    "amino_acid_composition": rec.aa_composition if rec.aa_composition else "",
                    "invalid_chars": ",".join(sorted(rec.invalid_chars)) if rec.invalid_chars else "",
                })

        print(f"\nResults exported to: {output_path}")

    except OSError as err:
        print(f"[ERROR] Could not write CSV file '{output_path}': {err}", file=sys.stderr)


# --------------------------------------------------------------------------
# CLI entry point
# --------------------------------------------------------------------------

def build_arg_parser() -> argparse.ArgumentParser:
    """Construct the argparse CLI parser."""
    parser = argparse.ArgumentParser(
        prog="fasta_analyzer.py",
        description="Analyze FASTA files: detect sequence type, compute GC content "
                    "or amino acid composition, and export results.",
    )

    parser.add_argument(
        "paths",
        nargs="+",
        help="One or more FASTA files and/or directories to analyze.",
    )
    parser.add_argument(
        "--type",
        choices=["DNA", "PROTEIN"],
        default=None,
        help="Force sequence type instead of auto-detecting.",
    )
    parser.add_argument(
        "--csv",
        metavar="OUTPUT.csv",
        default=None,
        help="Export results to the given CSV file path.",
    )

    return parser


def main() -> None:
    parser = build_arg_parser()
    args = parser.parse_args()

    # Step 1: Collect valid FASTA files from given paths
    fasta_files = collect_fasta_files(args.paths)

    if not fasta_files:
        print("[ERROR] No valid FASTA files found for the given input(s).", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(fasta_files)} FASTA file(s) to process.\n")

    # Step 2: Parse and analyze all sequences
    records = process_files(fasta_files, forced_type=args.type)

    if not records:
        print("[ERROR] No sequences were successfully processed.", file=sys.stderr)
        sys.exit(1)

    # Step 3: Display results
    print_results_table(records)
    print_summary(records)

    # Step 4: Optional CSV export
    if args.csv:
        export_csv(records, args.csv)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n[INFO] Interrupted by user.", file=sys.stderr)
        sys.exit(130)
    except Exception as fatal_err:
        print(f"[FATAL] Unexpected error: {fatal_err}", file=sys.stderr)
        sys.exit(1)