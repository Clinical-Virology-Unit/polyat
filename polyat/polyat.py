from __future__ import annotations

import argparse
import gzip
import sys
from html import escape
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


def write_html_summary(
    output_dir: Path, headers: list[tuple[str, str]], rows: list[list[str]]
) -> None:
    output_file = output_dir / "polyA_counts.html"
    css = (
        "body{font-family:Arial,sans-serif;margin:20px;background:#fefefe;}"
        "h1{margin-bottom:0.5em;}"
        "table{border-collapse:collapse;width:100%;font-family:Arial,sans-serif;}"
        "th,td{border:1px solid #ccc;padding:8px;text-align:center;}"
        "th{background-color:#f4f4f4;}"
        ".filters th{background-color:#fafafa;}"
        ".filters input{width:100%;box-sizing:border-box;padding:4px;}"
        "tr:nth-child(even){background:#fafafa;}"
    )
    script = (
        "const table=document.getElementById('polyat-table');"
        "const columnFilters=table.querySelectorAll('thead input[data-col]');"
        "function applyFilters(){"
        "table.querySelectorAll('tbody tr').forEach(row=>{"
        "let visible=true;"
        "columnFilters.forEach(input=>{"
        "if(!visible)return;"
        "const value=input.value.trim();"
        "if(!value)return;"
        "const col=parseInt(input.dataset.col,10);"
        "const type=input.dataset.type;"
        "const cell=row.children[col];"
        "if(!cell)return;"
        "const cellText=cell.innerText.trim();"
        "if(type==='number'){"
        "const cellValue=parseFloat(cellText);"
        "const filterValue=parseFloat(value);"
        "if(isNaN(filterValue)||isNaN(cellValue))return;"
        "if(cellValue<filterValue){visible=false;}"
        "}else{"
        "if(!cellText.toLowerCase().includes(value.toLowerCase())){"
        "visible=false;"
        "}"
        "}"
        "});"
        "row.style.display=visible?'':'none';"
        "});"
        "}"
        "columnFilters.forEach(input=>input.addEventListener('input',applyFilters));"
        "applyFilters();"
    )
    parts = [
        "<!DOCTYPE html>",
        "<html lang='en'>",
        "<head>",
        "<meta charset='utf-8' />",
        "<title>polyA Counts</title>",
        f"<style>{css}</style>",
        "</head>",
        "<body>",
        "<h1>polyA/T Summary</h1>",
        "<table id='polyat-table'>",
        "<thead>",
        "<tr>",
    ]
    for label, _type in headers:
        parts.append(f"<th>{escape(label)}</th>")
    parts.append("</tr>")
    parts.append("<tr class='filters'>")
    for idx, (_label, col_type) in enumerate(headers):
        if col_type == "number":
            placeholder = "min value"
            input_type = "number"
        else:
            placeholder = "text"
            input_type = "text"
        parts.append(
            "<th>"
            f"<input data-col='{idx}' data-type='{col_type}' "
            f"type='{input_type}' placeholder='{placeholder}' />"
            "</th>"
        )
    parts.append("</tr></thead><tbody>")
    for row in rows:
        parts.append("<tr>")
        parts.extend(f"<td>{escape(value)}</td>" for value in row)
        parts.append("</tr>")
    parts.extend(
        [
            "</tbody></table>",
            f"<script>{script}</script>",
            "</body>",
            "</html>",
        ]
    )
    output_file.write_text("\n".join(parts), encoding="utf-8")


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    input_dir = resolve_directory(args.input, "input")
    output_dir = resolve_directory(args.output, "output")
    output_file = output_dir / "polyA_counts.txt"

    fastq_files = find_fastq_files(input_dir)
    if not fastq_files:
        sys.exit(f"[error] No FASTQ/FASTQ.GZ files found in {input_dir}")

    headers = [
        ("Sample", "text"),
        ("Total_Reads", "number"),
        ("PolyA/T_10+", "number"),
        ("PolyA/T_15+", "number"),
        ("PolyA/T_20+", "number"),
        ("Percent_10+", "number"),
        ("Percent_15+", "number"),
        ("Percent_20+", "number"),
    ]
    header_line = "\t".join(label for label, _type in headers)
    html_rows: list[list[str]] = []

    with open(output_file, "w", encoding="utf-8") as out_handle:
        out_handle.write(header_line + "\n")
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
            html_rows.append(
                [
                    sample,
                    str(total),
                    str(poly10),
                    str(poly15),
                    str(poly20),
                    pct10,
                    pct15,
                    pct20,
                ]
            )

    write_html_summary(output_dir, headers, html_rows)

    print(f"[polyat] Summary written to {output_file}")
    print(f"[polyat] HTML summary written to {output_dir / 'polyA_counts.html'}")


if __name__ == "__main__":
    main()

