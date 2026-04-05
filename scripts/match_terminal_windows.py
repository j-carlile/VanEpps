from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import pandas as pd

# Easy-to-adjust default paths for this project.
DEFAULT_INPUT_TSV = Path("data/references/terminal_windows.tsv")
DEFAULT_FASTA_DIR = Path("data/raw_fastas/bacteria")
DEFAULT_OUTPUT_TSV = Path("data/references/terminal_windows_with_indices.tsv")

# This script reports 1-based coordinates to make the output easier to read in a
# biological context. The end index is inclusive.
INDEXING_NOTE = "1-based inclusive"


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Match terminal window sequences to genome FASTA files and record "
            "their positions."
        )
    )
    parser.add_argument(
        "--input-tsv",
        type=Path,
        default=DEFAULT_INPUT_TSV,
        help=f"Path to the input TSV file (default: {DEFAULT_INPUT_TSV})",
    )
    parser.add_argument(
        "--fasta-dir",
        type=Path,
        default=DEFAULT_FASTA_DIR,
        help=f"Directory containing FASTA files (default: {DEFAULT_FASTA_DIR})",
    )
    parser.add_argument(
        "--output-tsv",
        type=Path,
        default=DEFAULT_OUTPUT_TSV,
        help=f"Path for the output TSV file (default: {DEFAULT_OUTPUT_TSV})",
    )
    parser.add_argument(
        "--filters",
        nargs="*",
        default=[],
        help=(
            "Optional row filters in the form gene:end, for example "
            "--filters tufA:3 tufB:5"
        ),
    )
    return parser.parse_args()


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    translation_table = str.maketrans(
        {
            "A": "T",
            "T": "A",
            "C": "G",
            "G": "C",
            "N": "N",
            "a": "t",
            "t": "a",
            "c": "g",
            "g": "c",
            "n": "n",
        }
    )
    return sequence.translate(translation_table)[::-1].upper()


def normalize_sequence(value: object) -> str:
    """Clean sequence text so exact matching is consistent."""
    if pd.isna(value):
        return ""
    return "".join(str(value).split()).upper()


def load_fasta_sequences(fasta_path: Path) -> Dict[str, str]:
    """
    Read a FASTA file into a dictionary of {header: sequence}.

    The parser is intentionally simple and avoids extra dependencies.
    """
    contigs: Dict[str, List[str]] = {}
    current_header: str | None = None

    with fasta_path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()

            if not line:
                continue

            if line.startswith(">"):
                current_header = line[1:].strip()
                if not current_header:
                    raise ValueError(
                        f"Found an empty FASTA header in {fasta_path} at line {line_number}."
                    )
                contigs[current_header] = []
                continue

            if current_header is None:
                raise ValueError(
                    f"Found sequence data before the first FASTA header in {fasta_path}."
                )

            contigs[current_header].append(line.upper())

    return {header: "".join(sequence_parts) for header, sequence_parts in contigs.items()}


def find_all_exact_matches(sequence: str, query: str) -> List[Tuple[int, int]]:
    """
    Find all exact matches of query inside sequence.

    Returns a list of (start_index, end_index) tuples using 1-based inclusive
    coordinates.
    """
    matches: List[Tuple[int, int]] = []
    if not query:
        return matches

    search_start = 0
    while True:
        position = sequence.find(query, search_start)
        if position == -1:
            break

        start_index = position + 1
        end_index = start_index + len(query) - 1
        matches.append((start_index, end_index))

        # Advance by one base so repeated or overlapping exact matches are found.
        search_start = position + 1

    return matches


def resolve_column_name(columns: Iterable[str], expected_name: str, aliases: List[str]) -> str:
    """Match a required column name using a small set of aliases."""
    normalized_map = {column.strip().lower(): column for column in columns}

    for candidate in [expected_name, *aliases]:
        match = normalized_map.get(candidate.strip().lower())
        if match:
            return match

    raise KeyError(
        f"Could not find a column for '{expected_name}'. Checked aliases: {aliases}"
    )


def prepare_query(sequence: str, orientation: str) -> Tuple[str, bool]:
    """Return the sequence to search and whether reverse-complementing was used."""
    normalized_orientation = str(orientation).strip()

    if normalized_orientation == "-":
        return reverse_complement(sequence), True
    if normalized_orientation == "+":
        return sequence, False

    raise ValueError(
        f"Unexpected orientation value '{orientation}'. Expected '+' or '-'."
    )


def filter_terminal_windows(
    terminal_windows: pd.DataFrame,
    gene_col: str,
    end_col: str,
    filters: List[str],
) -> pd.DataFrame:
    """Keep only rows that match one or more gene:end filters."""
    if not filters:
        return terminal_windows

    mask = pd.Series(False, index=terminal_windows.index)

    for raw_filter in filters:
        if ":" not in raw_filter:
            raise ValueError(
                f"Invalid filter '{raw_filter}'. Use the format gene:end, for example tufA:3."
            )

        gene_value, end_value = raw_filter.split(":", 1)
        gene_value = gene_value.strip()
        end_value = end_value.strip()

        if not gene_value or not end_value:
            raise ValueError(
                f"Invalid filter '{raw_filter}'. Both gene and end must be provided."
            )

        current_mask = (
            terminal_windows[gene_col].astype(str).str.strip().eq(gene_value)
            & terminal_windows[end_col].astype(str).str.strip().eq(end_value)
        )
        mask = mask | current_mask

    filtered = terminal_windows.loc[mask].copy()
    if filtered.empty:
        raise ValueError(
            "The provided filters did not match any rows in the input TSV."
        )

    return filtered


def build_result_row(
    row: pd.Series,
    row_number: int,
    gene_col: str,
    end_col: str,
    original_sequence: str,
    searched_sequence: str,
    used_reverse_complement: bool,
    match_count: int,
    match_status: str,
    fasta_filename: object = pd.NA,
    contig_name: object = pd.NA,
    start_index: object = pd.NA,
    end_index: object = pd.NA,
) -> Dict[str, object]:
    """Create one output row while preserving the original TSV columns."""
    result = row.to_dict()
    result.update(
        {
            "row_id": row_number + 1,
            "gene_name": row[gene_col],
            "gene_end": row[end_col],
            "original_sequence": original_sequence,
            "searched_sequence": searched_sequence,
            "reverse_complement_used": used_reverse_complement,
            "fasta_filename": fasta_filename,
            "contig_name": contig_name,
            "start_index": start_index,
            "end_index": end_index,
            "match_count": match_count,
            "match_status": match_status,
            "indexing": INDEXING_NOTE,
        }
    )
    return result


def match_terminal_windows(
    terminal_windows: pd.DataFrame,
    fasta_dir: Path,
    gene_col: str,
    end_col: str,
    orientation_col: str,
    sequence_col: str,
) -> pd.DataFrame:
    """Match every terminal window sequence against every FASTA contig."""
    fasta_paths = (
        sorted(fasta_dir.glob("*.fa"))
        + sorted(fasta_dir.glob("*.fasta"))
        + sorted(fasta_dir.glob("*.fna"))
    )
    if not fasta_paths:
        raise FileNotFoundError(f"No FASTA files found in {fasta_dir}")

    prepared_rows: List[Dict[str, object]] = []
    matches_by_row: Dict[int, List[Dict[str, object]]] = {}

    for row_number, row in terminal_windows.iterrows():
        original_sequence = normalize_sequence(row[sequence_col])
        orientation = str(row[orientation_col]).strip()
        searched_sequence, used_reverse_complement = prepare_query(
            original_sequence,
            orientation,
        )

        prepared_rows.append(
            {
                "row_number": row_number,
                "row": row,
                "original_sequence": original_sequence,
                "searched_sequence": searched_sequence,
                "used_reverse_complement": used_reverse_complement,
            }
        )
        matches_by_row[row_number] = []

    for fasta_index, fasta_path in enumerate(fasta_paths, start=1):
        print(f"[{fasta_index}/{len(fasta_paths)}] Searching {fasta_path.name} ...")
        contigs = load_fasta_sequences(fasta_path)

        for contig_name, contig_sequence in contigs.items():
            for prepared in prepared_rows:
                searched_sequence = prepared["searched_sequence"]
                if not searched_sequence:
                    continue

                exact_matches = find_all_exact_matches(contig_sequence, searched_sequence)
                for start_index, end_index in exact_matches:
                    matches_by_row[prepared["row_number"]].append(
                        {
                            "fasta_filename": fasta_path.name,
                            "contig_name": contig_name,
                            "start_index": start_index,
                            "end_index": end_index,
                        }
                    )

    results: List[Dict[str, object]] = []

    for prepared in prepared_rows:
        row_number = prepared["row_number"]
        row_matches = matches_by_row[row_number]
        match_count = len(row_matches)

        if match_count == 0:
            match_status = "no_match"
        elif match_count == 1:
            match_status = "unique_match"
        else:
            match_status = "multiple_matches"

        if row_matches:
            for match in row_matches:
                results.append(
                    build_result_row(
                        row=prepared["row"],
                        row_number=row_number,
                        gene_col=gene_col,
                        end_col=end_col,
                        original_sequence=prepared["original_sequence"],
                        searched_sequence=prepared["searched_sequence"],
                        used_reverse_complement=prepared["used_reverse_complement"],
                        match_count=match_count,
                        match_status=match_status,
                        fasta_filename=match["fasta_filename"],
                        contig_name=match["contig_name"],
                        start_index=match["start_index"],
                        end_index=match["end_index"],
                    )
                )
        else:
            results.append(
                build_result_row(
                    row=prepared["row"],
                    row_number=row_number,
                    gene_col=gene_col,
                    end_col=end_col,
                    original_sequence=prepared["original_sequence"],
                    searched_sequence=prepared["searched_sequence"],
                    used_reverse_complement=prepared["used_reverse_complement"],
                    match_count=0,
                    match_status=match_status,
                )
            )

    return pd.DataFrame(results)


def main() -> None:
    args = parse_args()

    if not args.input_tsv.exists():
        raise FileNotFoundError(f"Input TSV not found: {args.input_tsv}")

    if not args.fasta_dir.exists():
        raise FileNotFoundError(f"FASTA directory not found: {args.fasta_dir}")

    terminal_windows = pd.read_csv(args.input_tsv, sep="\t")

    gene_col = resolve_column_name(terminal_windows.columns, "gene", ["gene_name"])
    end_col = resolve_column_name(
        terminal_windows.columns,
        "end",
        ["gene_end", "end_position"],
    )
    orientation_col = resolve_column_name(
        terminal_windows.columns,
        "orientation",
        ["strand"],
    )
    sequence_col = resolve_column_name(
        terminal_windows.columns,
        "sequence",
        ["seq", "window_sequence"],
    )
    terminal_windows = filter_terminal_windows(
        terminal_windows=terminal_windows,
        gene_col=gene_col,
        end_col=end_col,
        filters=args.filters,
    )

    matched_df = match_terminal_windows(
        terminal_windows=terminal_windows,
        fasta_dir=args.fasta_dir,
        gene_col=gene_col,
        end_col=end_col,
        orientation_col=orientation_col,
        sequence_col=sequence_col,
    )

    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    matched_df.to_csv(args.output_tsv, sep="\t", index=False)

    print(f"Input TSV: {args.input_tsv}")
    print(f"FASTA directory: {args.fasta_dir}")
    if args.filters:
        print(f"Filters applied: {', '.join(args.filters)}")
    print(f"Rows processed: {len(terminal_windows)}")
    print(f"Output written to: {args.output_tsv}")
    print(f"Coordinate system: {INDEXING_NOTE}")


if __name__ == "__main__":
    main()
