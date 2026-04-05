from __future__ import annotations

from collections import Counter
from difflib import SequenceMatcher
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import pandas as pd

INPUT_TSV = Path("data/references/terminal_windows.tsv")
FASTA_PATH = Path("data/raw_fastas/bacteria/GCA_019815155.1.fasta")
TARGET_GENE = "tufA"
TARGET_END = "3"
WINDOW_SIZES = [15, 20, 25, 30]
PRIMER_WINDOW_SIZES = [14, 15]


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


def sliding_window_matches(
    query_sequence: str,
    contigs: Dict[str, str],
    window_size: int,
    restrict_to_range: Tuple[int, int] | None = None,
) -> List[Dict[str, object]]:
    records: List[Dict[str, object]] = []

    if restrict_to_range is None:
        first_start = 1
        last_start = len(query_sequence) - window_size + 1
    else:
        region_start, region_end = restrict_to_range
        first_start = region_start
        last_start = region_end - window_size + 1

    if last_start < first_start:
        return records

    for start_1_based in range(first_start, last_start + 1):
        end_1_based = start_1_based + window_size - 1
        subwindow = query_sequence[start_1_based - 1 : end_1_based]

        hit_records: List[Tuple[str, int, int]] = []
        for contig_name, contig_sequence in contigs.items():
            matches = find_all_exact_matches(contig_sequence, subwindow)
            for genome_start, genome_end in matches:
                hit_records.append((contig_name, genome_start, genome_end))

        if hit_records:
            contig_summary = "; ".join(
                f"{contig_name} @ {genome_start}-{genome_end}"
                for contig_name, genome_start, genome_end in hit_records[:10]
            )
            if len(hit_records) > 10:
                contig_summary += f"; ... plus {len(hit_records) - 10} more hits"
        else:
            contig_summary = ""

        records.append(
            {
                "window_size": window_size,
                "subwindow_start": start_1_based,
                "subwindow_end": end_1_based,
                "subwindow_sequence": subwindow,
                "matched": bool(hit_records),
                "match_count": len(hit_records),
                "match_locations": contig_summary,
            }
        )

    return records


def region_label_from_position(start: int, end: int, sequence_length: int) -> str:
    midpoint = (start + end) / 2
    third = sequence_length / 3
    if midpoint <= third:
        return "beginning"
    if midpoint <= 2 * third:
        return "middle"
    return "end"


def summarize_window_results(records: List[Dict[str, object]], sequence_length: int) -> None:
    if not records:
        print("No windows were generated for this setting.")
        return

    total_windows = len(records)
    matched_records = [record for record in records if record["matched"]]
    matched_count = len(matched_records)
    print(f"Total windows checked: {total_windows}")
    print(f"Windows with at least one exact hit: {matched_count}")

    if not matched_records:
        print("No exact subwindow matches were found at this window size.")
        return

    region_counter = Counter(
        region_label_from_position(
            int(record["subwindow_start"]),
            int(record["subwindow_end"]),
            sequence_length,
        )
        for record in matched_records
    )
    print(f"Matched-window distribution: {dict(region_counter)}")

    match_count_values = [int(record["match_count"]) for record in matched_records]
    print(
        "Genome multiplicity among matched windows: "
        f"min={min(match_count_values)}, median={sorted(match_count_values)[len(match_count_values)//2]}, max={max(match_count_values)}"
    )

    longest_run = 0
    current_run = 0
    for record in records:
        if record["matched"]:
            current_run += 1
            longest_run = max(longest_run, current_run)
        else:
            current_run = 0

    print(f"Longest continuous run of matching windows: {longest_run}")

    print("Example matched windows:")
    for record in matched_records[:5]:
        print(
            f"- seq {record['subwindow_start']}-{record['subwindow_end']} | "
            f"match_count={record['match_count']} | {record['match_locations']}"
        )


def best_kmer_seed(records: List[Dict[str, object]]) -> Dict[str, object] | None:
    matched_records = [record for record in records if record["matched"]]
    if not matched_records:
        return None
    return min(matched_records, key=lambda record: (int(record["match_count"]), int(record["subwindow_start"])))


def extract_candidate_regions(
    contigs: Dict[str, str],
    seed_record: Dict[str, object],
    target_length: int,
    flank: int = 200,
) -> List[Dict[str, object]]:
    candidates: List[Dict[str, object]] = []
    seed_sequence = str(seed_record["subwindow_sequence"])

    for contig_name, contig_sequence in contigs.items():
        for seed_start, seed_end in find_all_exact_matches(contig_sequence, seed_sequence):
            region_start = max(1, seed_start - flank)
            region_end = min(len(contig_sequence), seed_end + flank)
            segment = contig_sequence[region_start - 1 : region_end]

            if len(segment) < target_length // 2:
                continue

            candidates.append(
                {
                    "contig_name": contig_name,
                    "seed_start": seed_start,
                    "seed_end": seed_end,
                    "region_start": region_start,
                    "region_end": region_end,
                    "segment": segment,
                }
            )

    return candidates


def local_alignment_report(query_sequence: str, candidates: List[Dict[str, object]]) -> None:
    if not candidates:
        print("No candidate regions were available for local alignment.")
        return

    best_candidate = None
    best_match = None

    for candidate in candidates:
        matcher = SequenceMatcher(None, query_sequence, candidate["segment"], autojunk=False)
        match = matcher.find_longest_match(
            0,
            len(query_sequence),
            0,
            len(candidate["segment"]),
        )
        if best_match is None or match.size > best_match.size:
            best_match = match
            best_candidate = candidate

    if best_match is None or best_candidate is None:
        print("Approximate local comparison did not return any candidate match.")
        return

    print(
        f"Best candidate contig: {best_candidate['contig_name']} | "
        f"candidate region {best_candidate['region_start']}-{best_candidate['region_end']}"
    )
    print(
        "Approximate local comparison result based on the longest exact shared block "
        "between the 400 bp query and the candidate genome region:"
    )
    print(f"- longest shared exact block length: {best_match.size}")
    print(
        f"- query coordinates (1-based): {best_match.a + 1}-{best_match.a + best_match.size}"
    )
    print(
        f"- candidate-region coordinates (1-based within extracted segment): "
        f"{best_match.b + 1}-{best_match.b + best_match.size}"
    )
    candidate_genome_start = best_candidate["region_start"] + best_match.b
    candidate_genome_end = candidate_genome_start + best_match.size - 1
    print(
        f"- approximate genome coordinates: {candidate_genome_start}-{candidate_genome_end}"
    )
    shared_block = query_sequence[best_match.a : best_match.a + best_match.size]
    print(f"- shared block sequence: {shared_block}")


def main() -> None:
    print_header("1. Inspect terminal_windows.tsv And Resolve Actual Column Names")
    print(f"Input TSV: {INPUT_TSV.resolve()}")
    terminal_windows = pd.read_csv(INPUT_TSV, sep="\t")
    print(f"Detected columns: {list(terminal_windows.columns)}")

    gene_col = resolve_column_name(terminal_windows.columns, "gene", ["gene_name"])
    end_col = resolve_column_name(terminal_windows.columns, "end", ["gene_end", "end_position"])
    orientation_col = resolve_column_name(terminal_windows.columns, "orientation", ["strand"])
    sequence_col = resolve_column_name(terminal_windows.columns, "sequence", ["seq", "window_sequence"])

    print("Resolved columns used for this debug run:")
    print(f"- gene column: {gene_col}")
    print(f"- end column: {end_col}")
    print(f"- orientation column: {orientation_col}")
    print(f"- sequence column: {sequence_col}")

    target_rows = terminal_windows[
        terminal_windows[gene_col].astype(str).str.strip().eq(TARGET_GENE)
        & terminal_windows[end_col].astype(str).str.strip().eq(TARGET_END)
    ].copy()

    print(f"Rows matching gene={TARGET_GENE!r} and end={TARGET_END!r}: {len(target_rows)}")
    if target_rows.empty:
        raise SystemExit("No matching row found.")

    print(target_rows.to_string())
    row_index = int(target_rows.index[0])
    row = target_rows.iloc[0]

    original_sequence = normalize_sequence(row[sequence_col])
    orientation = str(row[orientation_col]).strip()
    searched_sequence = reverse_complement(original_sequence) if orientation == "-" else original_sequence

    print_header("2. Print Target Row And Search Sequence")
    print(f"Selected row index (0-based): {row_index}")
    print(f"Selected row_id (1-based): {row_index + 1}")
    print(f"gene: {row[gene_col]}")
    print(f"end: {row[end_col]}")
    print(f"orientation: {orientation}")
    print(f"original sequence length: {len(original_sequence)}")
    print(f"searched sequence length: {len(searched_sequence)}")
    print(f"original sequence: {original_sequence}")
    print(f"searched sequence: {searched_sequence}")

    print_header("3. Confirm Full-Length Exact Match Against Single FASTA")
    print(f"FASTA path: {FASTA_PATH.resolve()}")
    contigs = load_fasta_sequences(FASTA_PATH)
    print(f"Contig count: {len(contigs)}")

    full_hits: List[Tuple[str, int, int]] = []
    for contig_name, contig_sequence in contigs.items():
        matches = find_all_exact_matches(contig_sequence, searched_sequence)
        for start_index, end_index in matches:
            full_hits.append((contig_name, start_index, end_index))

    if full_hits:
        print("Full searched sequence has an exact match.")
        for contig_name, start_index, end_index in full_hits:
            print(f"- {contig_name} | start={start_index} | end={end_index}")
    else:
        print("No exact full-length match was found. Continuing to sliding-window checks.")

    print_header("4. Sliding-Window Exact-Match Analysis Across The Full Sequence")
    all_window_results: Dict[int, List[Dict[str, object]]] = {}
    for window_size in WINDOW_SIZES:
        print()
        print(f"Window size: {window_size} bp")
        records = sliding_window_matches(searched_sequence, contigs, window_size)
        all_window_results[window_size] = records
        summarize_window_results(records, len(searched_sequence))

    print_header("5. Focused Check Near The Relevant Terminal End")
    if TARGET_END == "3":
        relevant_end_label = "end of the searched sequence"
        relevant_range = (max(1, len(searched_sequence) - 99), len(searched_sequence))
    else:
        relevant_end_label = "beginning of the searched sequence"
        relevant_range = (1, min(100, len(searched_sequence)))

    print(
        f"Using gene-end metadata, the likely primer-building region is near the {relevant_end_label} "
        f"for this {TARGET_END}' terminal window."
    )
    print(f"Focused sequence range checked: {relevant_range[0]}-{relevant_range[1]}")
    print(f"Focused sequence segment: {searched_sequence[relevant_range[0]-1:relevant_range[1]]}")

    for window_size in PRIMER_WINDOW_SIZES:
        print()
        print(f"Primer-scale window size: {window_size} bp")
        primer_records = sliding_window_matches(
            searched_sequence,
            contigs,
            window_size,
            restrict_to_range=relevant_range,
        )
        summarize_window_results(primer_records, len(searched_sequence))
        if any(record["matched"] for record in primer_records):
            print(
                "Warning: short exact matches can happen by chance in a bacterial genome, "
                "so these hits alone do not confirm the correct biological site."
            )

    print_header("6. Local Alignment Check Triggered By Best Exact Seed")
    alignment_seed = None
    for seed_window_size in [30, 25, 20, 15]:
        alignment_seed = best_kmer_seed(all_window_results.get(seed_window_size, []))
        if alignment_seed is not None:
            print(
                f"Using a {seed_window_size} bp exact seed from sequence positions "
                f"{alignment_seed['subwindow_start']}-{alignment_seed['subwindow_end']} "
                f"to anchor candidate regions."
            )
            break

    if alignment_seed is None:
        print("No exact seed window was found anywhere in the searched sequence.")
        print("Without an exact seed, this script cannot anchor a local alignment candidate region in the genome.")
    else:
        candidates = extract_candidate_regions(contigs, alignment_seed, len(searched_sequence))
        print(f"Candidate regions extracted from the genome: {len(candidates)}")
        local_alignment_report(searched_sequence, candidates)

    print_header("7. Interpretation")
    if full_hits:
        print("Interpretation: the 400 bp searched sequence is present as an exact substring in this genome.")
    else:
        matched_30 = sum(1 for record in all_window_results[30] if record["matched"])
        matched_25 = sum(1 for record in all_window_results[25] if record["matched"])
        matched_20 = sum(1 for record in all_window_results[20] if record["matched"])
        matched_15 = sum(1 for record in all_window_results[15] if record["matched"])

        if matched_30 == 0 and matched_25 == 0 and matched_20 == 0 and matched_15 == 0:
            print("Interpretation: there is no evidence here for even short exact conserved subsequences from this searched window in this genome.")
        elif matched_30 == 0 and matched_25 == 0:
            print("Interpretation: only short exact conserved subsequences appear to exist, not a long exact conserved block.")
        else:
            print("Interpretation: parts of a larger approximately matching region may be present, but the full 400 bp window is not an exact substring.")

        print(
            "Most likely conclusion for this case: treat this 400 bp window as an approximate conserved region "
            "or a source region from which smaller primer-relevant exact matches may need to be identified, "
            "rather than assuming the whole 400 bp sequence should match exactly in every genome."
        )


if __name__ == "__main__":
    main()
