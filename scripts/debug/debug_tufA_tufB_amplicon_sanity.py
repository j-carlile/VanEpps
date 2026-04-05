from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import pandas as pd

INPUT_TSV = Path("data/references/terminal_windows.tsv")
FASTA_PATH = Path("data/raw_fastas/bacteria/GCA_019815155.1.fasta")
NOTEBOOK_PATH = Path("notebooks/VanEpps_workflow_testing.ipynb")
EXPECTED_XLSX_PATH = Path("data/ver_1_data/gene_pair_amplicon_summary.xlsx")
DEFAULT_TARGETS = [("tufA", "3"), ("tufB", "5")]
WINDOW_SIZES = [14, 15, 16, 18, 20]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Estimate candidate amplicon spans from primer-scale exact matches."
    )
    parser.add_argument("--gene1", default="tufA")
    parser.add_argument("--end1", default="3")
    parser.add_argument("--gene2", default="tufB")
    parser.add_argument("--end2", default="5")
    parser.add_argument("--expected-length", type=int, default=None)
    return parser.parse_args()


def print_header(title: str) -> None:
    print()
    print("=" * 90)
    print(title)
    print("=" * 90)


def reverse_complement(sequence: str) -> str:
    translation_table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return sequence.translate(translation_table)[::-1].upper()


def normalize_sequence(value: object) -> str:
    if pd.isna(value):
        return ""
    return "".join(str(value).split()).upper()


def resolve_column_name(columns: Iterable[str], expected_name: str, aliases: List[str]) -> str:
    normalized_map = {column.strip().lower(): column for column in columns}
    for candidate in [expected_name, *aliases]:
        match = normalized_map.get(candidate.strip().lower())
        if match:
            return match
    raise KeyError(
        f"Could not find a column for '{expected_name}'. Checked aliases: {aliases}"
    )


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
                    raise ValueError(f"Empty FASTA header in {fasta_path} line {line_number}")
                contigs[current_header] = []
                continue
            if current_header is None:
                raise ValueError(f"Sequence before first FASTA header in {fasta_path}")
            contigs[current_header].append(line.upper())

    return {header: "".join(parts) for header, parts in contigs.items()}


def find_all_exact_matches(sequence: str, query: str) -> List[Tuple[int, int]]:
    matches: List[Tuple[int, int]] = []
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


def relevant_range(sequence_length: int, end_value: str, size: int = 100) -> Tuple[int, int]:
    if str(end_value).strip() == "3":
        return max(1, sequence_length - size + 1), sequence_length
    return 1, min(size, sequence_length)


def find_candidate_hits(
    row: pd.Series,
    row_index: int,
    gene_col: str,
    end_col: str,
    orientation_col: str,
    sequence_col: str,
    accession_col: str,
    contigs: Dict[str, str],
) -> List[Dict[str, object]]:
    original_sequence = normalize_sequence(row[sequence_col])
    orientation = str(row[orientation_col]).strip()
    searched_sequence = reverse_complement(original_sequence) if orientation == "-" else original_sequence
    region_start, region_end = relevant_range(len(searched_sequence), str(row[end_col]))

    hits: List[Dict[str, object]] = []
    for window_size in WINDOW_SIZES:
        last_start = region_end - window_size + 1
        if last_start < region_start:
            continue

        for subwindow_start in range(region_start, last_start + 1):
            subwindow_end = subwindow_start + window_size - 1
            subwindow_sequence = searched_sequence[subwindow_start - 1 : subwindow_end]

            current_hits: List[Tuple[str, int, int]] = []
            for contig_name, contig_sequence in contigs.items():
                matches = find_all_exact_matches(contig_sequence, subwindow_sequence)
                for genome_start, genome_end in matches:
                    current_hits.append((contig_name, genome_start, genome_end))

            for contig_name, genome_start, genome_end in current_hits:
                hits.append(
                    {
                        "gene": row[gene_col],
                        "end": row[end_col],
                        "row_id": row_index + 1,
                        "orientation": orientation,
                        "anchor_accession": row[accession_col],
                        "searched_sequence_length": len(searched_sequence),
                        "relevant_range_start": region_start,
                        "relevant_range_end": region_end,
                        "window_size": window_size,
                        "subwindow_start": subwindow_start,
                        "subwindow_end": subwindow_end,
                        "subwindow_sequence": subwindow_sequence,
                        "contig_name": contig_name,
                        "genome_start": genome_start,
                        "genome_end": genome_end,
                        "match_count_in_genome": len(current_hits),
                    }
                )

    return hits


def choose_best_hits_per_gene(hits: List[Dict[str, object]]) -> List[Dict[str, object]]:
    if not hits:
        return []
    return sorted(
        hits,
        key=lambda hit: (
            int(hit["match_count_in_genome"]),
            -int(hit["window_size"]),
            abs(int(hit["subwindow_start"]) - int(hit["relevant_range_start"])),
            int(hit["genome_start"]),
        ),
    )


def estimate_pairings(
    tufa_hits: List[Dict[str, object]],
    tufb_hits: List[Dict[str, object]],
    expected_length: int | None,
) -> List[Dict[str, object]]:
    pairings: List[Dict[str, object]] = []

    for left_hit in tufa_hits:
        for right_hit in tufb_hits:
            if left_hit["contig_name"] != right_hit["contig_name"]:
                continue

            left_start = int(left_hit["genome_start"])
            left_end = int(left_hit["genome_end"])
            right_start = int(right_hit["genome_start"])
            right_end = int(right_hit["genome_end"])

            upstream_hit = left_hit
            downstream_hit = right_hit
            upstream_start = left_start
            upstream_end = left_end
            downstream_start = right_start
            downstream_end = right_end

            if right_start < left_start:
                upstream_hit = right_hit
                downstream_hit = left_hit
                upstream_start = right_start
                upstream_end = right_end
                downstream_start = left_start
                downstream_end = left_end

            gap_between_hits = downstream_start - upstream_end - 1
            total_span = downstream_end - upstream_start + 1

            pairings.append(
                {
                    "contig_name": left_hit["contig_name"],
                    "upstream_gene": upstream_hit["gene"],
                    "downstream_gene": downstream_hit["gene"],
                    "upstream_window_size": upstream_hit["window_size"],
                    "downstream_window_size": downstream_hit["window_size"],
                    "upstream_seq_range": f"{upstream_hit['subwindow_start']}-{upstream_hit['subwindow_end']}",
                    "downstream_seq_range": f"{downstream_hit['subwindow_start']}-{downstream_hit['subwindow_end']}",
                    "upstream_genome_range": f"{upstream_start}-{upstream_end}",
                    "downstream_genome_range": f"{downstream_start}-{downstream_end}",
                    "gap_between_hits": gap_between_hits,
                    "total_span": total_span,
                    "formula": f"{downstream_end} - {upstream_start} + 1 = {total_span}",
                    "combined_uniqueness_score": int(upstream_hit["match_count_in_genome"]) + int(downstream_hit["match_count_in_genome"]),
                }
            )

    return sorted(
        pairings,
        key=lambda pairing: (
            pairing["combined_uniqueness_score"],
            -pairing["upstream_window_size"] - pairing["downstream_window_size"],
            abs(pairing["total_span"] - expected_length) if expected_length is not None else pairing["total_span"],
        ),
    )


def find_expected_amplicon_length(pair_id: str) -> Tuple[int | None, str]:

    try:
        excel = pd.ExcelFile(EXPECTED_XLSX_PATH)
        for sheet_name in excel.sheet_names:
            df = pd.read_excel(EXPECTED_XLSX_PATH, sheet_name=sheet_name)
            if "pair_id" in df.columns and pair_id in set(df["pair_id"].astype(str)):
                row = df[df["pair_id"].astype(str) == pair_id].iloc[0]
                if "coverage" in df.columns and pd.notna(row["coverage"]):
                    return int(row["coverage"]), f"{EXPECTED_XLSX_PATH} (sheet {sheet_name})"
    except Exception as exc:
        excel_error = str(exc)
    else:
        excel_error = "not found in workbook"

    try:
        notebook = json.loads(NOTEBOOK_PATH.read_text(encoding="utf-8"))
        pattern = re.compile(
            re.escape(pair_id) + r"\s+\S+\s+\d+\s+\S+\s+\d+\s+(\d+)"
        )
        for cell in notebook.get("cells", []):
            for output in cell.get("outputs", []):
                text_chunks = output.get("text", [])
                joined = "".join(text_chunks)
                match = pattern.search(joined)
                if match:
                    return int(match.group(1)), f"{NOTEBOOK_PATH} (fallback from notebook parsed output)"

        raw_text = NOTEBOOK_PATH.read_text(encoding="utf-8")
        raw_match = re.search(
            re.escape(pair_id) + r"\s+\S+\s+\d+\s+\S+\s+\d+\s+(\d+)",
            raw_text,
        )
        if raw_match:
            return int(raw_match.group(1)), f"{NOTEBOOK_PATH} (fallback from raw notebook text)"
    except Exception as exc:
        return None, f"Could not read expected length from Excel or notebook. Excel error: {excel_error}; notebook error: {exc}"

    return None, f"Expected length not found. Excel error: {excel_error}"


def main() -> None:
    args = parse_args()
    targets = [(args.gene1, str(args.end1)), (args.gene2, str(args.end2))]
    pair_id = f"{args.gene1}_{args.end1}-{args.gene2}_{args.end2}"

    print_header("1. Inspect Workspace Paths And Expected Amplicon Source")
    print(f"terminal_windows.tsv: {INPUT_TSV.resolve()}")
    print(f"Genome FASTA: {FASTA_PATH.resolve()}")
    expected_length = args.expected_length
    if expected_length is not None:
        expected_source = "command-line override"
        print(f"Expected amplicon length provided manually: {expected_length}")
        print(f"Expected length source: {expected_source}")
    else:
        expected_length, expected_source = find_expected_amplicon_length(pair_id)
        if expected_length is None:
            print(f"Expected amplicon length source not found cleanly: {expected_source}")
        else:
            print(f"Expected amplicon length found: {expected_length}")
            print(f"Expected length source: {expected_source}")

    print_header("2. Load terminal_windows.tsv And Resolve Actual Column Names")
    terminal_windows = pd.read_csv(INPUT_TSV, sep="\t")
    print(f"Detected columns: {list(terminal_windows.columns)}")

    gene_col = resolve_column_name(terminal_windows.columns, "gene", ["gene_name"])
    end_col = resolve_column_name(terminal_windows.columns, "end", ["gene_end", "end_position"])
    orientation_col = resolve_column_name(terminal_windows.columns, "orientation", ["strand"])
    sequence_col = resolve_column_name(terminal_windows.columns, "sequence", ["seq", "window_sequence"])
    accession_col = resolve_column_name(terminal_windows.columns, "anchor_accession", ["accession"])

    print(f"Using gene column: {gene_col}")
    print(f"Using end column: {end_col}")
    print(f"Using orientation column: {orientation_col}")
    print(f"Using sequence column: {sequence_col}")
    print(f"Using accession column: {accession_col}")

    selected_rows = terminal_windows[
        terminal_windows.apply(
            lambda row: (
                (str(row[gene_col]).strip(), str(row[end_col]).strip()) in targets
            ),
            axis=1,
        )
    ].copy()

    print_header("3. Selected Rows")
    if selected_rows.empty:
        raise SystemExit("Did not find the target tufA/tufB rows in terminal_windows.tsv")

    for row_index, row in selected_rows.iterrows():
        cleaned_sequence = normalize_sequence(row[sequence_col])
        print(f"row index (0-based): {row_index}")
        print(f"row_id (1-based): {row_index + 1}")
        print(f"gene: {row[gene_col]}")
        print(f"end: {row[end_col]}")
        print(f"orientation: {row[orientation_col]}")
        print(f"anchor_accession: {row[accession_col]}")
        print(f"sequence length: {len(cleaned_sequence)}")
        print(f"sequence: {cleaned_sequence}")
        print("-" * 60)

    print_header("4. Parse Genome And Search Primer-Scale Exact Matches")
    contigs = load_fasta_sequences(FASTA_PATH)
    print(f"Loaded contigs from genome: {len(contigs)}")

    hits_by_target: Dict[Tuple[str, str], List[Dict[str, object]]] = {}
    for row_index, row in selected_rows.iterrows():
        target_key = (str(row[gene_col]).strip(), str(row[end_col]).strip())
        cleaned_sequence = normalize_sequence(row[sequence_col])
        searched_sequence = (
            reverse_complement(cleaned_sequence)
            if str(row[orientation_col]).strip() == "-"
            else cleaned_sequence
        )
        region_start, region_end = relevant_range(len(searched_sequence), str(row[end_col]))
        end_label = "end" if str(row[end_col]).strip() == "3" else "beginning"
        print(
            f"{target_key[0]} / end={target_key[1]}: searching near the biologically relevant "
            f"{end_label} of the searched sequence, positions {region_start}-{region_end}"
        )
        hits = find_candidate_hits(
            row=row,
            row_index=row_index,
            gene_col=gene_col,
            end_col=end_col,
            orientation_col=orientation_col,
            sequence_col=sequence_col,
            accession_col=accession_col,
            contigs=contigs,
        )
        hits_by_target[target_key] = hits
        print(f"Candidate exact hits found: {len(hits)}")
        for hit in choose_best_hits_per_gene(hits)[:10]:
            print(
                f"- window={hit['window_size']} bp | seq {hit['subwindow_start']}-{hit['subwindow_end']} | "
                f"contig={hit['contig_name']} | genome={hit['genome_start']}-{hit['genome_end']} | "
                f"match_count={hit['match_count_in_genome']}"
            )

    print_header("5. Pair tufA And tufB Hits To Estimate Amplicon Length")
    tufa_hits = choose_best_hits_per_gene(hits_by_target.get(targets[0], []))
    tufb_hits = choose_best_hits_per_gene(hits_by_target.get(targets[1], []))

    pairings = estimate_pairings(tufa_hits[:25], tufb_hits[:25], expected_length)
    if not pairings:
        print("No same-contig hit pairs were found from exact 14-20 bp matches.")
    else:
        print(f"Same-contig candidate pairings found: {len(pairings)}")
        for pairing in pairings[:15]:
            print(
                f"- contig={pairing['contig_name']} | upstream={pairing['upstream_gene']} {pairing['upstream_genome_range']} "
                f"| downstream={pairing['downstream_gene']} {pairing['downstream_genome_range']} | "
                f"gap_between_hits={pairing['gap_between_hits']} | total_span={pairing['total_span']} | "
                f"formula: {pairing['formula']}"
            )

    print_header("6. Compare Candidate Lengths Against Prior Expectation")
    if expected_length is None:
        print("No expected amplicon length was available for comparison.")
    elif not pairings:
        print(f"Expected amplicon length: {expected_length}")
        print("Observed candidate lengths: none")
    else:
        print(f"Expected amplicon length: {expected_length}")
        for pairing in pairings[:10]:
            difference = pairing["total_span"] - expected_length
            print(
                f"- contig={pairing['contig_name']} | observed total_span={pairing['total_span']} | "
                f"difference from expected={difference}"
            )

    print_header("7. Interpretation")
    if not pairings:
        print(
            "No confident estimate can be made from exact 14-20 bp matches alone for this genome. "
            "The hit sets do not yield a same-contig tufA/tufB pairing."
        )
    else:
        best_pairing = pairings[0]
        if expected_length is not None and abs(best_pairing["total_span"] - expected_length) <= 25:
            print(
                "Likely correct pairing found and the estimated amplicon length looks broadly consistent "
                "with prior expectations."
            )
        elif expected_length is not None:
            print(
                "Multiple or weak candidate pairings exist, and the best exact-hit span is not very close "
                "to the prior expected amplicon length."
            )
        else:
            print(
                "Candidate pairings exist, but without a readable expected-length file this remains a rough "
                "sanity check only."
            )


if __name__ == "__main__":
    main()
