from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import pandas as pd

INPUT_TSV = Path("data/references/terminal_windows.tsv")
FASTA_PATH = Path("data/raw_fastas/bacteria/GCA_019815155.1.fasta")
WINDOW_SIZES = [14, 15, 16, 18, 20]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Check how well terminal-window regions are preserved in one genome."
    )
    parser.add_argument("--gene1", default="rpoB")
    parser.add_argument("--end1", default="3")
    parser.add_argument("--gene2", default="rpoC")
    parser.add_argument("--end2", default="5")
    return parser.parse_args()


def print_header(title: str) -> None:
    print()
    print("=" * 90)
    print(title)
    print("=" * 90)


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
                    raise ValueError(f"Empty FASTA header at line {line_number}")
                contigs[current_header] = []
                continue
            if current_header is None:
                raise ValueError("Sequence encountered before first FASTA header")
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


def relevant_range(sequence_length: int, end_value: str, size: int = 100) -> Tuple[int, int]:
    if str(end_value).strip() == "3":
        return max(1, sequence_length - size + 1), sequence_length
    return 1, min(size, sequence_length)


def choose_anchor_hit(
    row: pd.Series,
    row_index: int,
    gene_col: str,
    end_col: str,
    orientation_col: str,
    sequence_col: str,
    contigs: Dict[str, str],
) -> Dict[str, object] | None:
    original_sequence = normalize_sequence(row[sequence_col])
    searched_sequence = (
        reverse_complement(original_sequence)
        if str(row[orientation_col]).strip() == "-"
        else original_sequence
    )
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
                for genome_start, genome_end in find_all_exact_matches(contig_sequence, subwindow_sequence):
                    current_hits.append((contig_name, genome_start, genome_end))

            for contig_name, genome_start, genome_end in current_hits:
                hits.append(
                    {
                        "row_id": row_index + 1,
                        "gene": row[gene_col],
                        "end": row[end_col],
                        "orientation": row[orientation_col],
                        "original_sequence": original_sequence,
                        "searched_sequence": searched_sequence,
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

    if not hits:
        return None

    return sorted(
        hits,
        key=lambda hit: (
            int(hit["match_count_in_genome"]),
            -int(hit["window_size"]),
            abs(int(hit["subwindow_start"]) - region_start),
            int(hit["genome_start"]),
        ),
    )[0]


def project_full_window(anchor_hit: Dict[str, object], contig_sequence: str) -> Dict[str, object]:
    searched_sequence = str(anchor_hit["searched_sequence"])
    query_start = int(anchor_hit["subwindow_start"])
    query_end = int(anchor_hit["subwindow_end"])
    genome_start = int(anchor_hit["genome_start"])
    genome_end = int(anchor_hit["genome_end"])

    projected_start = genome_start - (query_start - 1)
    projected_end = genome_end + (len(searched_sequence) - query_end)

    if projected_start < 1 or projected_end > len(contig_sequence):
        raise ValueError(
            f"Projected region {projected_start}-{projected_end} falls outside the contig bounds."
        )

    genome_segment = contig_sequence[projected_start - 1 : projected_end]
    return {
        "projected_start": projected_start,
        "projected_end": projected_end,
        "genome_segment": genome_segment,
    }


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
            pointer[i][j] = (
                "diag" if best == diagonal else "up" if best == up else "left"
            )

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
        base1 == base2 and base1 != "-"
        for base1, base2 in zip(alignment_1, alignment_2)
    )
    identity = (exact_matches / compared_positions * 100.0) if compared_positions else 0.0
    gaps = sum(base1 == "-" or base2 == "-" for base1, base2 in zip(alignment_1, alignment_2))

    return {
        "percent_identity": identity,
        "exact_matches": exact_matches,
        "compared_positions": compared_positions,
        "gaps": gaps,
        "alignment_1": alignment_1,
        "alignment_2": alignment_2,
    }


def mismatch_examples(seq1: str, seq2: str, limit: int = 10) -> List[str]:
    examples: List[str] = []
    for index, (base1, base2) in enumerate(zip(seq1, seq2), start=1):
        if base1 != base2:
            examples.append(f"pos {index}: terminal_windows={base1}, genome={base2}")
        if len(examples) >= limit:
            break
    return examples


def main() -> None:
    args = parse_args()
    target_pairs = [(args.gene1, str(args.end1)), (args.gene2, str(args.end2))]

    print_header("1. Inspect Workspace Inputs")
    print(f"terminal_windows.tsv: {INPUT_TSV.resolve()}")
    print(f"Genome FASTA: {FASTA_PATH.resolve()}")

    terminal_windows = pd.read_csv(INPUT_TSV, sep="\t")
    print(f"Detected columns: {list(terminal_windows.columns)}")

    gene_col = resolve_column_name(terminal_windows.columns, "gene", ["gene_name"])
    end_col = resolve_column_name(terminal_windows.columns, "end", ["gene_end"])
    orientation_col = resolve_column_name(terminal_windows.columns, "orientation", ["strand"])
    sequence_col = resolve_column_name(terminal_windows.columns, "sequence", ["seq"])

    selected = terminal_windows[
        terminal_windows.apply(
            lambda row: (
                (str(row[gene_col]).strip(), str(row[end_col]).strip()) in target_pairs
            ),
            axis=1,
        )
    ].copy()

    if selected.empty:
        raise SystemExit("Target rows were not found in terminal_windows.tsv.")

    contigs = load_fasta_sequences(FASTA_PATH)

    print_header("2. Choose Best Terminal Anchors")
    for row_index, row in selected.iterrows():
        anchor_hit = choose_anchor_hit(
            row=row,
            row_index=row_index,
            gene_col=gene_col,
            end_col=end_col,
            orientation_col=orientation_col,
            sequence_col=sequence_col,
            contigs=contigs,
        )

        gene_label = f"{row[gene_col]}_{row[end_col]}"
        if anchor_hit is None:
            print(f"{gene_label}: no primer-scale exact anchor hit found.")
            continue

        print(f"{gene_label}:")
        print(f"- row_id: {anchor_hit['row_id']}")
        print(f"- orientation: {anchor_hit['orientation']}")
        print(f"- chosen anchor size: {anchor_hit['window_size']} bp")
        print(f"- anchor sequence positions in window: {anchor_hit['subwindow_start']}-{anchor_hit['subwindow_end']}")
        print(f"- anchor sequence: {anchor_hit['subwindow_sequence']}")
        print(f"- contig: {anchor_hit['contig_name']}")
        print(f"- genome anchor positions: {anchor_hit['genome_start']}-{anchor_hit['genome_end']}")

        projected = project_full_window(
            anchor_hit=anchor_hit,
            contig_sequence=contigs[str(anchor_hit["contig_name"])],
        )

        direct_pct, direct_matches = direct_identity(
            str(anchor_hit["searched_sequence"]),
            projected["genome_segment"],
        )
        aligned = global_alignment_identity(
            str(anchor_hit["searched_sequence"]),
            projected["genome_segment"],
        )

        print_header(f"3. Preservation Check For {gene_label}")
        print(
            f"Projected genome region from anchor: {projected['projected_start']}-{projected['projected_end']}"
        )
        print(f"Projected genome segment length: {len(projected['genome_segment'])}")
        print(
            "This projection assumes the anchor sits in the same relative place within "
            "the approximately conserved 400 bp region."
        )
        print()
        print("Direct position-by-position comparison (no gaps):")
        print(
            f"- exact matches: {direct_matches}/{len(str(anchor_hit['searched_sequence']))}"
        )
        print(f"- percent identity: {direct_pct:.2f}%")
        print()
        print("Simple global alignment comparison (allows gaps):")
        print(f"- exact aligned matches: {aligned['exact_matches']}/{aligned['compared_positions']}")
        print(f"- percent identity: {aligned['percent_identity']:.2f}%")
        print(f"- alignment gaps introduced: {aligned['gaps']}")
        print()
        print("First mismatch examples from the direct projection:")
        examples = mismatch_examples(
            str(anchor_hit["searched_sequence"]),
            projected["genome_segment"],
        )
        if examples:
            for example in examples:
                print(f"- {example}")
        else:
            print("- none; projected region is an exact match")


if __name__ == "__main__":
    main()
