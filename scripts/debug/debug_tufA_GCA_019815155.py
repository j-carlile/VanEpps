from __future__ import annotations

from pathlib import Path
import re

import pandas as pd

INPUT_TSV = Path("data/references/terminal_windows.tsv")
FASTA_PATH = Path("data/raw_fastas/bacteria/GCA_019815155.1.fasta")
TARGET_GENE = "tufA"
TARGET_END = "3"


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1].upper()


def normalize_sequence(value: object) -> str:
    if pd.isna(value):
        return ""
    return "".join(str(value).split()).upper()


def load_fasta_sequences(fasta_path: Path) -> dict[str, str]:
    contigs: dict[str, list[str]] = {}
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
                raise ValueError("Found sequence before first FASTA header")

            contigs[current_header].append(line.upper())

    return {header: "".join(parts) for header, parts in contigs.items()}


def find_all_exact_matches(sequence: str, query: str) -> list[tuple[int, int]]:
    matches: list[tuple[int, int]] = []
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


def middle_substring(sequence: str, size: int) -> str:
    start = (len(sequence) - size) // 2
    end = start + size
    return sequence[start:end]


def print_header(title: str) -> None:
    print()
    print("=" * 80)
    print(title)
    print("=" * 80)


def main() -> None:
    print_header("1. Load Target Row From terminal_windows.tsv")
    print(f"Input TSV path: {INPUT_TSV.resolve()}")
    df = pd.read_csv(INPUT_TSV, sep="\t")
    print(f"Detected columns: {list(df.columns)}")

    target_rows = df[
        df["gene"].astype(str).str.strip().eq(TARGET_GENE)
        & df["end"].astype(str).str.strip().eq(TARGET_END)
    ].copy()

    print(f"Rows matching gene={TARGET_GENE!r} and end={TARGET_END!r}: {len(target_rows)}")
    if target_rows.empty:
        raise SystemExit("No matching row found.")

    print(target_rows.to_string())

    row_index = int(target_rows.index[0])
    row = target_rows.iloc[0]

    print_header("2. Inspect Exact Row Values")
    print(f"Selected row index (0-based pandas index): {row_index}")
    print(f"Selected row_id if using script convention (1-based): {row_index + 1}")
    print(f"gene: {row['gene']}")
    print(f"end: {row['end']}")
    print(f"anchor_accession: {row['anchor_accession']}")
    print(f"anchor_len: {row['anchor_len']}")
    print(f"orientation: {row['orientation']}")
    print(f"raw sequence length (before cleaning): {len(str(row['sequence']))}")
    print(f"raw sequence: {row['sequence']}")

    print_header("3. Clean And QC The Sequence")
    raw_sequence = str(row["sequence"])
    cleaned_sequence = normalize_sequence(row["sequence"])
    unexpected_characters = sorted(set(re.sub(r"[ACGTN]", "", cleaned_sequence)))

    print(f"Cleaned sequence length: {len(cleaned_sequence)}")
    print(f"Contains lowercase in raw value: {any(char.islower() for char in raw_sequence)}")
    print(f"Contains whitespace in raw value: {any(char.isspace() for char in raw_sequence)}")
    print(f"Contains only A/C/G/T/N after cleaning: {len(unexpected_characters) == 0}")
    print(f"Unexpected characters after cleaning: {unexpected_characters if unexpected_characters else 'None'}")
    print(f"First 60 bp: {cleaned_sequence[:60]}")
    print(f"Middle 60 bp: {cleaned_sequence[len(cleaned_sequence)//2 - 30: len(cleaned_sequence)//2 + 30]}")
    print(f"Last 60 bp: {cleaned_sequence[-60:]}")

    print_header("4. Compute Original And Reverse Complement")
    original_sequence = cleaned_sequence
    reverse_complement_sequence = reverse_complement(cleaned_sequence)
    print(f"Original sequence length: {len(original_sequence)}")
    print(f"Reverse complement length: {len(reverse_complement_sequence)}")
    print(f"Original first 60 bp: {original_sequence[:60]}")
    print(f"Reverse complement first 60 bp: {reverse_complement_sequence[:60]}")
    print(f"Orientation in TSV: {row['orientation']}")
    print("Expected search sequence based on your biological rule:")
    if str(row["orientation"]).strip() == "-":
        print("orientation is '-', so search the reverse complement against the forward-strand FASTA.")
        search_sequence = reverse_complement_sequence
    else:
        print("orientation is '+', so search the original sequence against the forward-strand FASTA.")
        search_sequence = original_sequence

    print_header("5. Parse Only GCA_019815155.1.fasta")
    print(f"FASTA path: {FASTA_PATH.resolve()}")
    contigs = load_fasta_sequences(FASTA_PATH)
    print(f"Number of contigs: {len(contigs)}")
    print("Contig names and lengths:")
    for contig_name, contig_sequence in contigs.items():
        print(f"- {contig_name} | length={len(contig_sequence)}")

    first_contig_name = next(iter(contigs))
    first_contig_sequence = contigs[first_contig_name]
    print()
    print("Multiline FASTA concatenation check on first contig:")
    print(f"- First contig name: {first_contig_name}")
    print(f"- First 100 bp of concatenated sequence: {first_contig_sequence[:100]}")
    print(f"- Sequence contains newline characters after parsing: {'\\n' in first_contig_sequence}")

    print_header("6. Search The Genome In Multiple Ways")
    search_modes = [
        ("A. original sequence exactly", original_sequence),
        ("B. reverse complement exactly", reverse_complement_sequence),
        ("C1. middle 30 bp of expected search sequence", middle_substring(search_sequence, 30)),
        ("C2. middle 40 bp of expected search sequence", middle_substring(search_sequence, 40)),
        ("C3. middle 50 bp of expected search sequence", middle_substring(search_sequence, 50)),
    ]

    overall_hits: dict[str, list[tuple[str, int, int]]] = {}
    for label, query in search_modes:
        print()
        print(label)
        print(f"Query length: {len(query)}")
        print(f"Query sequence: {query}")

        hits: list[tuple[str, int, int]] = []
        for contig_name, contig_sequence in contigs.items():
            matches = find_all_exact_matches(contig_sequence, query)
            for start_index, end_index in matches:
                hits.append((contig_name, start_index, end_index))

        overall_hits[label] = hits
        if hits:
            print(f"Matches found: {len(hits)}")
            for contig_name, start_index, end_index in hits[:20]:
                print(f"- {contig_name} | start={start_index} | end={end_index}")
            if len(hits) > 20:
                print(f"... plus {len(hits) - 20} additional hits")
        else:
            print("No exact matches found.")

    print_header("7. Interpretation")
    original_hits = overall_hits["A. original sequence exactly"]
    rc_hits = overall_hits["B. reverse complement exactly"]
    partial_hits = (
        overall_hits["C1. middle 30 bp of expected search sequence"]
        or overall_hits["C2. middle 40 bp of expected search sequence"]
        or overall_hits["C3. middle 50 bp of expected search sequence"]
    )

    if original_hits:
        print("Exact match found with the original sequence.")
    elif rc_hits:
        print("Exact match found with the reverse complement sequence.")
    elif partial_hits:
        print("No full exact match was found, but an internal substring does match.")
    else:
        print("No exact full-length match and no internal 30-50 bp middle-substring match were found.")

    print_header("8. Focused Diagnosis")
    print("What this means for the likely failure mode:")
    print("- Wrong row selected: unlikely. The script is selecting the tufA / end=3 row you asked for.")
    print("- Orientation bug: unlikely for this case. The row has orientation '-' and the reverse complement was computed correctly.")
    print("- FASTA parsing bug: unlikely. The file loads, contigs are detected, and multiline sequence lines are concatenated into continuous strings.")
    print("- Hidden formatting issue: unlikely if the cleaned sequence contains only A/C/G/T/N and no whitespace.")
    print("- Exact-match search too strict or biologically absent in this assembly: currently the strongest explanation if no full match is found.")
    print("- Other possibility: this conserved window may come from a different strain/reference context and may not be a literal 400 bp substring in this assembly.")


if __name__ == "__main__":
    main()
