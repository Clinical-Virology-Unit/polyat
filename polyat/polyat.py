from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path
from typing import Iterable, Tuple


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="polyat",
        description=(
            "Quantify poly-A/T stretches (>=10/15/20 nt) across FASTQ reads "
            "and summarize counts per sample."
        ),
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Directory containing .fastq/.fastq.gz/.fq/.fq.gz files.",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Directory where the summary table will be written.",
    )
    return parser.parse_args(argv)


def resolve_directory(path_str: str, role: str) -> Path:
    path = Path(path_str).expanduser().resolve()
    if role == "input":
        if not path.exists():
            sys.exit(f"[error] Input path does not exist: {path}")
        if not path.is_dir():
            sys.exit(f"[error] Input path is not a directory: {path}")
    else:
        path.mkdir(parents=True, exist_ok=True)
    return path


def find_fastq_files(input_dir: Path) -> list[Path]:
    fastq_files: list[Path] = []
    for entry in sorted(input_dir.iterdir()):
        if entry.is_file() and has_fastq_suffix(entry.name):
            fastq_files.append(entry)
    return fastq_files


def has_fastq_suffix(filename: str) -> bool:
    valid_suffixes = (".fastq", ".fastq.gz", ".fq", ".fq.gz")
    return any(filename.endswith(suffix) for suffix in valid_suffixes)


def open_fastq(path: Path) -> Iterable[str]:
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, "r", encoding="utf-8")


def sanitize_sample_name(file_path: Path) -> str:
    name = file_path.name
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def count_poly_runs(file_path: Path) -> Tuple[int, int, int, int]:
    total = 0
    poly10 = 0
    poly15 = 0
    poly20 = 0
    with open_fastq(file_path) as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline().strip()
            handle.readline()  # +
            handle.readline()  # quality
            if not seq:
                continue
            total += 1
            longest = longest_poly_run(seq)
            if longest >= 10:
                poly10 += 1
            if longest >= 15:
                poly15 += 1
            if longest >= 20:
                poly20 += 1
    return total, poly10, poly15, poly20


def longest_poly_run(sequence: str) -> int:
    longest = 0
    current = 0
    prev_char = ""
    for base in sequence.upper():
        if base not in ("A", "T"):
            current = 0
            prev_char = ""
            continue
        if base == prev_char:
            current += 1
        else:
            current = 1
            prev_char = base
        if current > longest:
            longest = current
    return longest


def format_percent(count: int, total: int) -> str:
    if total == 0:
        return "0.00"
    return f"{(count * 100) / total:.2f}"


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    input_dir = resolve_directory(args.input, "input")
    output_dir = resolve_directory(args.output, "output")
    output_file = output_dir / "polyA_counts.txt"

    fastq_files = find_fastq_files(input_dir)
    if not fastq_files:
        sys.exit(f"[error] No FASTQ/FASTQ.GZ files found in {input_dir}")

    header = (
        "Sample\tTotal_Reads\tPolyA/T_10+\tPolyA/T_15+\tPolyA/T_20+\t"
        "Percent_10+\tPercent_15+\tPercent_20+"
    )

    with open(output_file, "w", encoding="utf-8") as out_handle:
        out_handle.write(header + "\n")
        for file_path in fastq_files:
            sample = sanitize_sample_name(file_path)
            total, poly10, poly15, poly20 = count_poly_runs(file_path)
            pct10 = format_percent(poly10, total)
            pct15 = format_percent(poly15, total)
            pct20 = format_percent(poly20, total)
            line = (
                f"{sample}\t{total}\t{poly10}\t{poly15}\t{poly20}\t"
                f"{pct10}\t{pct15}\t{pct20}"
            )
            out_handle.write(line + "\n")

    print(f"[polyat] Summary written to {output_file}")


if __name__ == "__main__":
    main()

