from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import pandas as pd

INPUT_TSV = Path("data/references/terminal_windows.tsv")
FASTA_DIR = Path("data/raw_fastas/bacteria")
WINDOW_SIZES = [14, 15, 16, 18, 20]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Search full 400 bp windows for exact anchor hits across multiple genomes."
    )
    parser.add_argument("--gene1", default="rpoB")
    parser.add_argument("--end1", default="3")
    parser.add_argument("--gene2", default="rpoC")
    parser.add_argument("--end2", default="5")
    parser.add_argument("--expected-length", type=int, default=174)
    parser.add_argument(
        "--genomes",
        nargs="+",
        default=["GCA_019815155.1", "GCA_050293385.1", "GCA_050360595.1", "GCA_910592365.1"],
    )
    return parser.parse_args()


def print_header(title: str) -> None:
    print()
    print("=" * 100)
    print(title)
    print("=" * 100)


def normalize_sequence(value: object) -> str:
    if pd.isna(value):
        return ""
    return "".join(str(value).split()).upper()


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1].upper()


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
        matches.append((position + 1, position + len(query)))
        start = position + 1
    return matches


def terminal_offset(sequence_length: int, subwindow_start: int, subwindow_end: int, end_value: str) -> int:
    if str(end_value).strip() == "3":
        return sequence_length - subwindow_end
    return subwindow_start - 1


def direct_identity(seq1: str, seq2: str) -> Tuple[float, int]:
    matches = sum(base1 == base2 for base1, base2 in zip(seq1, seq2))
    return (matches / len(seq1) * 100.0), matches


def global_alignment_identity(seq1: str, seq2: str) -> Dict[str, object]:
    rows = len(seq1) + 1
    cols = len(seq2) + 1
    score = [[0] * cols for _ in range(rows)]
    pointer = [[""] * cols for _ in range(rows)]
    gap_penalty = -1
    match_score = 1
    mismatch_penalty = -1

    for i in range(1, rows):
        score[i][0] = i * gap_penalty
        pointer[i][0] = "up"
    for j in range(1, cols):
        score[0][j] = j * gap_penalty
        pointer[0][j] = "left"

    for i in range(1, rows):
        for j in range(1, cols):
            diagonal = score[i - 1][j - 1] + (
                match_score if seq1[i - 1] == seq2[j - 1] else mismatch_penalty
            )
            up = score[i - 1][j] + gap_penalty
            left = score[i][j - 1] + gap_penalty
            best = max(diagonal, up, left)
            score[i][j] = best
            pointer[i][j] = "diag" if best == diagonal else "up" if best == up else "left"

    aligned_seq1: List[str] = []
    aligned_seq2: List[str] = []
    i = len(seq1)
    j = len(seq2)
    while i > 0 or j > 0:
        direction = pointer[i][j]
        if direction == "diag":
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif direction == "up":
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append("-")
            i -= 1
        else:
            aligned_seq1.append("-")
            aligned_seq2.append(seq2[j - 1])
            j -= 1

    aligned_seq1.reverse()
    aligned_seq2.reverse()
    alignment_1 = "".join(aligned_seq1)
    alignment_2 = "".join(aligned_seq2)

    compared_positions = sum(
        base1 != "-" and base2 != "-" for base1, base2 in zip(alignment_1, alignment_2)
    )
    exact_matches = sum(
        base1 == base2 and base1 != "-" for base1, base2 in zip(alignment_1, alignment_2)
    )
    identity = (exact_matches / compared_positions * 100.0) if compared_positions else 0.0
    gaps = sum(base1 == "-" or base2 == "-" for base1, base2 in zip(alignment_1, alignment_2))
    return {
        "percent_identity": identity,
        "exact_matches": exact_matches,
        "compared_positions": compared_positions,
        "gaps": gaps,
    }


def collect_hits_for_row(
    row: pd.Series,
    row_index: int,
    gene_col: str,
    end_col: str,
    orientation_col: str,
    sequence_col: str,
    contigs: Dict[str, str],
) -> List[Dict[str, object]]:
    original_sequence = normalize_sequence(row[sequence_col])
    searched_sequence = (
        reverse_complement(original_sequence)
        if str(row[orientation_col]).strip() == "-"
        else original_sequence
    )
    hits: List[Dict[str, object]] = []

    for window_size in WINDOW_SIZES:
        last_start = len(searched_sequence) - window_size + 1
        for subwindow_start in range(1, last_start + 1):
            subwindow_end = subwindow_start + window_size - 1
            subwindow_sequence = searched_sequence[subwindow_start - 1 : subwindow_end]

            current_hits: List[Tuple[str, int, int]] = []
            for contig_name, contig_sequence in contigs.items():
                for genome_start, genome_end in find_all_exact_matches(contig_sequence, subwindow_sequence):
                    current_hits.append((contig_name, genome_start, genome_end))

            for contig_name, genome_start, genome_end in current_hits:
                hits.append(
                    {
                        "row_id": row_index + 1,
                        "gene": row[gene_col],
                        "end": str(row[end_col]).strip(),
                        "orientation": str(row[orientation_col]).strip(),
                        "searched_sequence": searched_sequence,
                        "window_size": window_size,
                        "subwindow_start": subwindow_start,
                        "subwindow_end": subwindow_end,
                        "subwindow_sequence": subwindow_sequence,
                        "contig_name": contig_name,
                        "genome_start": genome_start,
                        "genome_end": genome_end,
                        "match_count_in_genome": len(current_hits),
                        "terminal_offset": terminal_offset(
                            len(searched_sequence), subwindow_start, subwindow_end, str(row[end_col]).strip()
                        ),
                    }
                )

    preferred_center = len(searched_sequence) if str(row[end_col]).strip() == "3" else 1
    return sorted(
        hits,
        key=lambda hit: (
            int(hit["match_count_in_genome"]),
            -int(hit["window_size"]),
            abs(((int(hit["subwindow_start"]) + int(hit["subwindow_end"])) / 2) - preferred_center),
            int(hit["genome_start"]),
        ),
    )


def project_full_window(hit: Dict[str, object], contig_sequence: str) -> Dict[str, object]:
    searched_sequence = str(hit["searched_sequence"])
    projected_start = int(hit["genome_start"]) - (int(hit["subwindow_start"]) - 1)
    projected_end = int(hit["genome_end"]) + (len(searched_sequence) - int(hit["subwindow_end"]))
    if projected_start < 1 or projected_end > len(contig_sequence):
        raise ValueError("Projected full window falls outside contig bounds.")
    genome_segment = contig_sequence[projected_start - 1 : projected_end]
    return {
        "projected_start": projected_start,
        "projected_end": projected_end,
        "genome_segment": genome_segment,
    }


def pair_hits(hits1: List[Dict[str, object]], hits2: List[Dict[str, object]], expected_length: int) -> List[Dict[str, object]]:
    pairings: List[Dict[str, object]] = []
    for hit1 in hits1:
        for hit2 in hits2:
            if hit1["contig_name"] != hit2["contig_name"]:
                continue

            upstream = hit1 if int(hit1["genome_start"]) <= int(hit2["genome_start"]) else hit2
            downstream = hit2 if upstream is hit1 else hit1

            raw_gap = int(downstream["genome_start"]) - int(upstream["genome_end"]) - 1
            raw_span = int(downstream["genome_end"]) - int(upstream["genome_start"]) + 1
            adjusted_gap = raw_gap - int(upstream["terminal_offset"]) - int(downstream["terminal_offset"])
            adjusted_span = raw_span - int(upstream["terminal_offset"]) - int(downstream["terminal_offset"])

            pairings.append(
                {
                    "contig_name": upstream["contig_name"],
                    "upstream": upstream,
                    "downstream": downstream,
                    "raw_gap": raw_gap,
                    "raw_span": raw_span,
                    "adjusted_gap": adjusted_gap,
                    "adjusted_span": adjusted_span,
                    "distance_from_expected": abs(adjusted_gap - expected_length),
                    "combined_uniqueness": int(upstream["match_count_in_genome"]) + int(downstream["match_count_in_genome"]),
                }
            )

    return sorted(
        pairings,
        key=lambda pairing: (
            pairing["combined_uniqueness"],
            pairing["distance_from_expected"],
            -(int(pairing["upstream"]["window_size"]) + int(pairing["downstream"]["window_size"])),
        ),
    )


def main() -> None:
    args = parse_args()
    target_pairs = [(args.gene1, str(args.end1)), (args.gene2, str(args.end2))]

    print_header("Inputs")
    print(f"terminal_windows.tsv: {INPUT_TSV.resolve()}")
    print(f"FASTA directory: {FASTA_DIR.resolve()}")
    print(f"Pair: {args.gene1}_{args.end1}-{args.gene2}_{args.end2}")
    print(f"Expected length: {args.expected_length}")
    print(f"Genomes: {', '.join(args.genomes)}")

    terminal_windows = pd.read_csv(INPUT_TSV, sep="\t")
    gene_col = resolve_column_name(terminal_windows.columns, "gene", ["gene_name"])
    end_col = resolve_column_name(terminal_windows.columns, "end", ["gene_end"])
    orientation_col = resolve_column_name(terminal_windows.columns, "orientation", ["strand"])
    sequence_col = resolve_column_name(terminal_windows.columns, "sequence", ["seq"])

    selected_rows = terminal_windows[
        terminal_windows.apply(
            lambda row: ((str(row[gene_col]).strip(), str(row[end_col]).strip()) in target_pairs),
            axis=1,
        )
    ].copy()
    if len(selected_rows) != 2:
        raise SystemExit("Expected exactly two matching rows in terminal_windows.tsv.")

    for genome_id in args.genomes:
        fasta_path = FASTA_DIR / f"{genome_id}.fasta"
        print_header(f"Genome: {genome_id}")
        contigs = load_fasta_sequences(fasta_path)
        print(f"Contigs loaded: {len(contigs)}")

        hits_by_target: Dict[Tuple[str, str], List[Dict[str, object]]] = {}
        for row_index, row in selected_rows.iterrows():
            key = (str(row[gene_col]).strip(), str(row[end_col]).strip())
            hits = collect_hits_for_row(
                row=row,
                row_index=row_index,
                gene_col=gene_col,
                end_col=end_col,
                orientation_col=orientation_col,
                sequence_col=sequence_col,
                contigs=contigs,
            )
            hits_by_target[key] = hits
            print(f"{key[0]}_{key[1]} total full-window hits found: {len(hits)}")
            for hit in hits[:10]:
                print(
                    f"- window={hit['window_size']} | seq={hit['subwindow_start']}-{hit['subwindow_end']} "
                    f"| contig={hit['contig_name']} | genome={hit['genome_start']}-{hit['genome_end']} "
                    f"| offset_to_terminal={hit['terminal_offset']} | match_count={hit['match_count_in_genome']}"
                )

        pairings = pair_hits(hits_by_target[target_pairs[0]], hits_by_target[target_pairs[1]], args.expected_length)
        if not pairings:
            print("No same-contig pairings found even when searching the full 400 bp sequences.")
            continue

        print()
        print(f"Same-contig pairings found: {len(pairings)}")
        for pairing in pairings[:10]:
            upstream = pairing["upstream"]
            downstream = pairing["downstream"]
            print(
                f"- contig={pairing['contig_name']} | {upstream['gene']} {upstream['genome_start']}-{upstream['genome_end']} "
                f"({upstream['subwindow_sequence']}) -> {downstream['gene']} {downstream['genome_start']}-{downstream['genome_end']} "
                f"({downstream['subwindow_sequence']})"
            )
            print(
                f"  raw_gap={pairing['raw_gap']} | adjusted_gap={pairing['adjusted_gap']} | "
                f"raw_span={pairing['raw_span']} | adjusted_span={pairing['adjusted_span']} | "
                f"distance_from_expected={pairing['distance_from_expected']}"
            )

        best_pairing = pairings[0]
        print_header(f"Best-Pair Preservation Check: {genome_id}")
        for label, hit in [("upstream", best_pairing["upstream"]), ("downstream", best_pairing["downstream"])]:
            projected = project_full_window(hit, contigs[str(hit["contig_name"])])
            direct_pct, direct_matches = direct_identity(str(hit["searched_sequence"]), projected["genome_segment"])
            aligned = global_alignment_identity(str(hit["searched_sequence"]), projected["genome_segment"])
            print(f"{label}: {hit['gene']}_{hit['end']}")
            print(f"- contig: {hit['contig_name']}")
            print(f"- anchor positions in sequence: {hit['subwindow_start']}-{hit['subwindow_end']}")
            print(f"- anchor genome positions: {hit['genome_start']}-{hit['genome_end']}")
            print(f"- projected 400 bp region: {projected['projected_start']}-{projected['projected_end']}")
            print(f"- direct identity: {direct_matches}/400 = {direct_pct:.2f}%")
            print(
                f"- gap-allowing identity: {aligned['exact_matches']}/{aligned['compared_positions']} "
                f"= {aligned['percent_identity']:.2f}% with {aligned['gaps']} gaps"
            )


if __name__ == "__main__":
    main()
