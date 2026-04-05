from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import pandas as pd

DEFAULT_INPUT_TSV = Path("data/references/terminal_windows.tsv")
DEFAULT_FASTA_PATH = Path("data/raw_fastas/bacteria/GCA_019815155.1.fasta")
DEFAULT_OUTPUT_TSV = Path("data/references/debug_GCA_019815155_matches.tsv")
INDEXING_NOTE = "1-based inclusive"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Debug exact sequence matching against a single FASTA file."
    )
    parser.add_argument(
        "--input-tsv",
        type=Path,
        default=DEFAULT_INPUT_TSV,
        help=f"Path to terminal_windows TSV (default: {DEFAULT_INPUT_TSV})",
    )
    parser.add_argument(
        "--fasta-path",
        type=Path,
        default=DEFAULT_FASTA_PATH,
        help=f"Path to one FASTA file to search (default: {DEFAULT_FASTA_PATH})",
    )
    parser.add_argument(
        "--filters",
        nargs="*",
        default=["tufA:3", "tufB:5"],
        help="Row filters in the form gene:end, for example tufA:3 tufB:5",
    )
    parser.add_argument(
        "--output-tsv",
        type=Path,
        default=DEFAULT_OUTPUT_TSV,
        help=f"Path for debug output (default: {DEFAULT_OUTPUT_TSV})",
    )
    return parser.parse_args()


def reverse_complement(sequence: str) -> str:
    translation_table = str.maketrans("ATCGNatcgn", "TAGCNtagcn")
    return sequence.translate(translation_table)[::-1].upper()


def normalize_sequence(value: object) -> str:
    if pd.isna(value):
        return ""
    return "".join(str(value).split()).upper()


def load_fasta_sequences(fasta_path: Path) -> Dict[str, str]:
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

    return {header: "".join(parts) for header, parts in contigs.items()}


def find_all_exact_matches(sequence: str, query: str) -> List[Tuple[int, int]]:
    matches: List[Tuple[int, int]] = []
    if not query:
        return matches

    start = 0
    while True:
        position = sequence.find(query, start)
        if position == -1:
            break

        one_based_start = position + 1
        one_based_end = one_based_start + len(query) - 1
        matches.append((one_based_start, one_based_end))
        start = position + 1

    return matches


def resolve_column_name(columns: Iterable[str], expected_name: str, aliases: List[str]) -> str:
    normalized_map = {column.strip().lower(): column for column in columns}
    for candidate in [expected_name, *aliases]:
        match = normalized_map.get(candidate.strip().lower())
        if match:
            return match
    raise KeyError(
        f"Could not find a column for '{expected_name}'. Checked aliases: {aliases}"
    )


def filter_terminal_windows(
    terminal_windows: pd.DataFrame,
    gene_col: str,
    end_col: str,
    filters: List[str],
) -> pd.DataFrame:
    mask = pd.Series(False, index=terminal_windows.index)

    for raw_filter in filters:
        if ":" not in raw_filter:
            raise ValueError(
                f"Invalid filter '{raw_filter}'. Use the format gene:end, for example tufA:3."
            )

        gene_value, end_value = raw_filter.split(":", 1)
        current_mask = (
            terminal_windows[gene_col].astype(str).str.strip().eq(gene_value.strip())
            & terminal_windows[end_col].astype(str).str.strip().eq(end_value.strip())
        )
        mask = mask | current_mask

    filtered = terminal_windows.loc[mask].copy()
    if filtered.empty:
        raise ValueError("The provided filters did not match any rows in the input TSV.")

    return filtered


def main() -> None:
    args = parse_args()

    terminal_windows = pd.read_csv(args.input_tsv, sep="\t")
    gene_col = resolve_column_name(terminal_windows.columns, "gene", ["gene_name"])
    end_col = resolve_column_name(terminal_windows.columns, "end", ["gene_end"])
    orientation_col = resolve_column_name(terminal_windows.columns, "orientation", ["strand"])
    sequence_col = resolve_column_name(terminal_windows.columns, "sequence", ["seq"])

    filtered = filter_terminal_windows(
        terminal_windows=terminal_windows,
        gene_col=gene_col,
        end_col=end_col,
        filters=args.filters,
    )

    contigs = load_fasta_sequences(args.fasta_path)
    results: List[Dict[str, object]] = []

    for row_number, row in filtered.iterrows():
        original_sequence = normalize_sequence(row[sequence_col])
        orientation = str(row[orientation_col]).strip()
        searched_sequence = (
            reverse_complement(original_sequence) if orientation == "-" else original_sequence
        )

        row_matches: List[Tuple[str, int, int]] = []
        for contig_name, contig_sequence in contigs.items():
            for start_index, end_index in find_all_exact_matches(contig_sequence, searched_sequence):
                row_matches.append((contig_name, start_index, end_index))

        if row_matches:
            for contig_name, start_index, end_index in row_matches:
                results.append(
                    {
                        "row_id": row_number + 1,
                        "gene": row[gene_col],
                        "end": row[end_col],
                        "orientation": orientation,
                        "original_sequence": original_sequence,
                        "searched_sequence": searched_sequence,
                        "reverse_complement_used": orientation == "-",
                        "fasta_filename": args.fasta_path.name,
                        "contig_name": contig_name,
                        "start_index": start_index,
                        "end_index": end_index,
                        "match_count": len(row_matches),
                        "match_status": "unique_match" if len(row_matches) == 1 else "multiple_matches",
                        "indexing": INDEXING_NOTE,
                    }
                )
        else:
            results.append(
                {
                    "row_id": row_number + 1,
                    "gene": row[gene_col],
                    "end": row[end_col],
                    "orientation": orientation,
                    "original_sequence": original_sequence,
                    "searched_sequence": searched_sequence,
                    "reverse_complement_used": orientation == "-",
                    "fasta_filename": args.fasta_path.name,
                    "contig_name": pd.NA,
                    "start_index": pd.NA,
                    "end_index": pd.NA,
                    "match_count": 0,
                    "match_status": "no_match",
                    "indexing": INDEXING_NOTE,
                }
            )

    results_df = pd.DataFrame(results)
    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    results_df.to_csv(args.output_tsv, sep="\t", index=False)

    print(f"Input TSV: {args.input_tsv}")
    print(f"FASTA searched: {args.fasta_path}")
    print(f"Filters applied: {', '.join(args.filters)}")
    print(f"Rows processed: {len(filtered)}")
    print(f"Output written to: {args.output_tsv}")
    print(results_df[['row_id', 'gene', 'end', 'match_count', 'match_status']].to_string(index=False))


if __name__ == "__main__":
    main()
