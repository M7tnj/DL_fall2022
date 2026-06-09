import argparse
import os
import sys
import traceback
from pathlib import Path
from datetime import datetime

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqIO import QualityIO
except ImportError:
    print("ERROR: biopython is required. Install with: pip install biopython")
    sys.exit(1)


# ─── Constants ──────────────────────────────────────────────────────────────

SANGER_EXTENSIONS = {".ab1", ".scf", ".phd.1", ".phd", ".seq"}
FORMATS_WITH_QUALITY = {".ab1", ".scf", ".phd.1", ".phd"}
FORMATS_NO_QUALITY = {".seq"}
DEFAULT_QUALITY_SCORE = 40  # Phred score assigned when no quality data exists


def detect_format(filepath: str) -> str:
    """
    Detect the sequencing file format from the file extension.

    Args:
        filepath: Path to the sequencing file.

    Returns:
        Format string recognized by Bio.SeqIO, or empty string if unknown.
    """
    path = Path(filepath)
    name = path.name.lower()

    # Check compound extensions first (e.g., .phd.1)
    if name.endswith(".phd.1"):
        return "phd"
    if name.endswith(".phd"):
        return "phd"

    ext = path.suffix.lower()
    format_map = {
        ".ab1": "abi",
        ".scf": "scf",
        ".seq": "fasta",  # .seq is usually just a FASTA sequence
    }
    return format_map.get(ext, "")


def has_quality_scores(filepath: str) -> bool:
    """Check if the file format contains quality scores."""
    name = Path(filepath).name.lower()
    if name.endswith(".phd.1") or name.endswith(".phd"):
        return True
    ext = Path(filepath).suffix.lower()
    return ext in FORMATS_WITH_QUALITY


def convert_with_quality(filepath: str, fmt: str, default_quality: int = DEFAULT_QUALITY_SCORE) -> SeqRecord:
    """
    Convert a file that has quality scores to a SeqRecord with letter_annotations.

    Args:
        filepath: Path to the input file.
        fmt: BioPython format string (abi, scf, phd).
        default_quality: Fallback quality if parsing fails.

    Returns:
        SeqRecord with quality scores populated.
    """
    try:
        record = SeqIO.read(filepath, fmt)
        # Verify quality scores are present
        if "phred_quality" in record.letter_annotations:
            quals = record.letter_annotations["phred_quality"]
            if len(quals) == len(record.seq):
                return record
            else:
                print(f"  ⚠ Quality length mismatch in {filepath}: "
                      f"seq={len(record.seq)}, qual={len(quals)}. "
                      f"Using default quality.")
        # Fall through to default quality
    except Exception as e:
        print(f"  ⚠ Could not parse quality from {filepath}: {e}")
        print(f"    Assigning default Phred quality of {default_quality}")

    # Read sequence without quality, then assign default
    return convert_without_quality(filepath, fmt, default_quality)


def convert_without_quality(filepath: str, fmt: str, default_quality: int = DEFAULT_QUALITY_SCORE) -> SeqRecord:
    """
    Convert a file without quality scores, assigning a uniform default quality.

    Args:
        filepath: Path to the input file.
        fmt: BioPython format string.
        default_quality: Default Phred quality score to assign.

    Returns:
        SeqRecord with default quality scores.
    """
    try:
        record = SeqIO.read(filepath, fmt)
    except Exception as e:
        # For .seq files, try plain text reading
        if fmt == "fasta":
            record = _read_plain_seq(filepath)
        else:
            raise ValueError(f"Cannot read {filepath}: {e}")

    # Assign default quality scores
    quals = [default_quality] * len(record.seq)
    record.letter_annotations["phred_quality"] = quals
    return record


def _read_plain_seq(filepath: str) -> SeqRecord:
    """
    Read a plain .seq file (just nucleotides, possibly with a header line).
    Many .seq files are simple FASTA or raw sequence.
    """
    try:
        # Try FASTA format first
        record = SeqIO.read(filepath, "fasta")
        return record
    except Exception:
        pass

    # Try reading as raw sequence (one line or multi-line, no header)
    try:
        with open(filepath, "r") as f:
            lines = f.readlines()

        # Filter out empty lines and comment lines
        seq_lines = []
        header = ""
        for line in lines:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith(">"):
                header = line[1:].strip()
                continue
            # Only keep valid nucleotide characters
            cleaned = "".join(c for c in line.upper() if c in "ACGTN")
            if cleaned:
                seq_lines.append(cleaned)

        if not seq_lines:
            raise ValueError(f"No sequence data found in {filepath}")

        sequence = "".join(seq_lines)
        seq_id = header or Path(filepath).stem

        return SeqRecord(
            Seq(sequence),
            id=seq_id,
            name=seq_id,
            description=f"Converted from {Path(filepath).name}"
        )
    except Exception as e:
        raise ValueError(f"Cannot parse .seq file {filepath}: {e}")


def convert_file(filepath: str, output_folder: str, default_quality: int = DEFAULT_QUALITY_SCORE) -> str:
    """
    Convert a single Sanger sequencing file to FASTQ format.

    Args:
        filepath: Path to the input file.
        output_folder: Directory to write the output .fastq file.
        default_quality: Default Phred quality for files without quality data.

    Returns:
        Path to the generated FASTQ file, or empty string on failure.
    """
    filepath = str(filepath)
    fmt = detect_format(filepath)

    if not fmt:
        print(f"  ✗ Skipping {filepath}: unrecognized file format")
        return ""

    # Read and convert
    if has_quality_scores(filepath):
        record = convert_with_quality(filepath, fmt, default_quality)
    else:
        record = convert_without_quality(filepath, fmt, default_quality)

    # Preserve original sequence name/ID, clean it up for FASTQ compatibility
    if not record.id or record.id == "<unknown id>":
        record.id = Path(filepath).stem
    if not record.name or record.name == "<unknown name>":
        record.name = Path(filepath).stem

    # Ensure the description includes source info
    source_note = f"Converted from {Path(filepath).name}"
    if has_quality_scores(filepath):
        source_note += " (real Phred quality scores)"
    else:
        source_note += f" (default Phred quality={default_quality}, no original quality data)"

    if not record.description or record.description == "<unknown description>":
        record.description = source_note
    else:
        # Append conversion note if not already there
        if "Converted from" not in record.description:
            record.description += f" {source_note}"

    # Write FASTQ
    out_name = Path(filepath).stem + ".fastq"
    out_path = os.path.join(output_folder, out_name)

    # Handle duplicate filenames
    counter = 1
    while os.path.exists(out_path):
        out_name = f"{Path(filepath).stem}_{counter}.fastq"
        out_path = os.path.join(output_folder, out_name)
        counter += 1

    with open(out_path, "w") as f:
        SeqIO.write(record, f, "fastq")

    return out_path


def convert_folder(input_folder: str, output_folder: str = "", default_quality: int = DEFAULT_QUALITY_SCORE) -> dict:
    """
    Convert all Sanger sequencing files in a folder to FASTQ.

    Args:
        input_folder: Directory containing sequencing files.
        output_folder: Directory for output FASTQ files (default: input_folder/fastq_output).
        default_quality: Default Phred quality for files without quality data.

    Returns:
        Summary dictionary with conversion results.
    """
    input_folder = str(input_folder)

    if not os.path.isdir(input_folder):
        raise FileNotFoundError(f"Input folder not found: {input_folder}")

    # Set default output folder
    if not output_folder:
        output_folder = os.path.join(input_folder, "fastq_output")

    os.makedirs(output_folder, exist_ok=True)

    # Find all convertible files
    files_to_convert = []
    for fname in sorted(os.listdir(input_folder)):
        fpath = os.path.join(input_folder, fname)
        if not os.path.isfile(fpath):
            continue
        name_lower = fname.lower()
        # Check compound extensions
        if name_lower.endswith(".phd.1") or name_lower.endswith(".phd"):
            files_to_convert.append(fpath)
            continue
        ext = Path(fname).suffix.lower()
        if ext in SANGER_EXTENSIONS:
            files_to_convert.append(fpath)

    if not files_to_convert:
        print(f"No convertible files found in {input_folder}")
        print(f"Supported formats: {', '.join(sorted(SANGER_EXTENSIONS))}")
        return {"total": 0, "success": 0, "failed": 0, "output_folder": output_folder, "results": []}

    print(f"Found {len(files_to_convert)} file(s) to convert in {input_folder}")
    print(f"Output folder: {output_folder}")
    print("-" * 60)

    results = []
    for fpath in files_to_convert:
        fname = os.path.basename(fpath)
        try:
            out_path = convert_file(fpath, output_folder, default_quality)
            if out_path:
                has_qual = has_quality_scores(fpath)
                qual_note = "real Phred scores" if has_qual else f"default Q={default_quality}"
                print(f"  ✓ {fname} → {os.path.basename(out_path)} ({qual_note})")
                results.append({
                    "input": fpath,
                    "output": out_path,
                    "has_quality": has_qual,
                    "status": "success"
                })
            else:
                print(f"  ✗ {fname}: conversion failed (unknown format)")
                results.append({"input": fpath, "output": "", "status": "failed", "error": "unknown format"})
        except Exception as e:
            print(f"  ✗ {fname}: {e}")
            results.append({"input": fpath, "output": "", "status": "failed", "error": str(e)})

    success = sum(1 for r in results if r["status"] == "success")
    failed = sum(1 for r in results if r["status"] == "failed")

    print("-" * 60)
    print(f"Conversion complete: {success} succeeded, {failed} failed out of {len(results)} files")

    return {
        "total": len(results),
        "success": success,
        "failed": failed,
        "output_folder": output_folder,
        "results": results
    }


def merge_fastq_files(fastq_folder: str, output_path: str = "") -> str:
    """
    Merge all FASTQ files in a folder into a single combined FASTQ file.
    This is useful when you have multiple Sanger reads covering different
    regions of the same virus genome.

    Args:
        fastq_folder: Directory containing .fastq files.
        output_path: Path for the merged file (default: fastq_folder/combined.fastq).

    Returns:
        Path to the merged FASTQ file.
    """
    if not output_path:
        output_path = os.path.join(fastq_folder, "combined.fastq")

    fastq_files = sorted([
        os.path.join(fastq_folder, f)
        for f in os.listdir(fastq_folder)
        if f.endswith(".fastq") and f != "combined.fastq"
    ])

    if not fastq_files:
        print("No FASTQ files found to merge.")
        return ""

    total_reads = 0
    with open(output_path, "w") as out_f:
        for fq_path in fastq_files:
            try:
                count = 0
                for record in SeqIO.parse(fq_path, "fastq"):
                    # Ensure unique read IDs by prepending source filename
                    source = Path(fq_path).stem
                    if not record.id.startswith(source):
                        record.id = f"{source}_{record.id}"
                        record.name = record.id
                    SeqIO.write(record, out_f, "fastq")
                    count += 1
                total_reads += count
                print(f"  Merged {os.path.basename(fq_path)} ({count} read(s))")
            except Exception as e:
                print(f"  ⚠ Error reading {fq_path}: {e}")

    print(f"\n  Combined {total_reads} reads into {output_path}")
    return output_path


def main():
    """Command-line entry point."""
    parser = argparse.ArgumentParser(
        description="Convert Sanger sequencing files (.ab1, .scf, .phd.1, .phd, .seq) to FASTQ format.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Convert all files in a folder
  python convert_to_fastq.py --input_folder ./my_sequences

  # Convert a single file
  python convert_to_fastq.py --input_file ./sample.ab1

  # Specify output folder and custom default quality
  python convert_to_fastq.py --input_folder ./sanger_data --output_folder ./fastq_output --default_quality 30

  # Convert and merge into one FASTQ (for multi-read Sanger runs)
  python convert_to_fastq.py --input_folder ./sanger_data --merge

Supported formats:
  .ab1   — Applied Biosystems chromatogram (with quality scores)
  .scf   — Standard Chromatogram Format (with quality scores)
  .phd.1 — Phred output (with quality scores)
  .phd   — Phred output (with quality scores)
  .seq   — Plain sequence file (NO quality scores → default Phred 40)

Important notes:
  - Sanger FASTQ uses Phred+33 encoding (same as Illumina 1.8+).
  - .ab1 and .scf files contain real quality scores from the instrument.
  - .seq files have NO quality data; a uniform default quality is assigned.
    Use these cautiously in variant-calling pipelines.
  - Sanger reads are typically 500-1200 bp (much longer than NGS reads).
"""
    )

    # Input options (mutually exclusive)
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--input_folder", "-i",
        type=str,
        help="Folder containing Sanger sequencing files to convert"
    )
    input_group.add_argument(
        "--input_file", "-f",
        type=str,
        help="Single Sanger sequencing file to convert"
    )

    # Output options
    parser.add_argument(
        "--output_folder", "-o",
        type=str,
        default="",
        help="Output folder for FASTQ files (default: <input_folder>/fastq_output)"
    )
    parser.add_argument(
        "--default_quality", "-q",
        type=int,
        default=DEFAULT_QUALITY_SCORE,
        help=f"Default Phred quality score for files without quality data (default: {DEFAULT_QUALITY_SCORE})"
    )

    # Merge option
    parser.add_argument(
        "--merge", "-m",
        action="store_true",
        help="Merge all converted FASTQ files into a single combined.fastq file"
    )

    args = parser.parse_args()

    print("=" * 60)
    print("Sanger Sequencing → FASTQ Converter")
    print("=" * 60)
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()

    # ─── Single file mode ──────────────────────────────────────────────
    if args.input_file:
        if not os.path.isfile(args.input_file):
            print(f"ERROR: File not found: {args.input_file}")
            sys.exit(1)

        output_folder = args.output_folder or os.path.join(
            os.path.dirname(args.input_file), "fastq_output"
        )
        os.makedirs(output_folder, exist_ok=True)

        print(f"Converting: {args.input_file}")
        print(f"Output to:  {output_folder}")
        print("-" * 60)

        try:
            out_path = convert_file(args.input_file, output_folder, args.default_quality)
            if out_path:
                has_qual = has_quality_scores(args.input_file)
                qual_note = "real Phred scores" if has_qual else f"default Q={args.default_quality}"
                print(f"  ✓ Converted: {os.path.basename(args.input_file)} → {os.path.basename(out_path)} ({qual_note})")
                print(f"\nFASTQ file: {out_path}")
            else:
                print(f"  ✗ Failed to convert {args.input_file}")
                sys.exit(1)
        except Exception as e:
            print(f"  ✗ Error: {e}")
            traceback.print_exc()
            sys.exit(1)

    # ─── Folder mode ───────────────────────────────────────────────────
    else:
        summary = convert_folder(args.input_folder, args.output_folder, args.default_quality)

        if summary["success"] > 0 and args.merge:
            print(f"\nMerging FASTQ files...")
            merged_path = merge_fastq_files(summary["output_folder"])
            if merged_path:
                print(f"\nMerged FASTQ: {merged_path}")
                summary["merged_file"] = merged_path

    print("\n" + "=" * 60)
    print("Done!")
    print("=" * 60)


if __name__ == "__main__":
    main()
