Project

Name: GenomicGapID

Goal

Analyze bacterial genome FASTA files to locate conserved source regions and primer-scale subsequences, identify expected amplicons, and support discrimination scoring for primer-set selection.

Current biological understanding
FASTA files are treated as forward-strand genomic sequence
Orientation matters for sequence comparison
Positive orientation sequences can be compared directly
Negative orientation sequences must be reverse complemented before searching/comparison
terminal_windows.tsv stores larger conserved genomic regions, not necessarily the exact final primers
Actual primer-like regions are much shorter, roughly 14 to 15 bp
Reverse primer handling must account for complement/reverse-complement logic correctly
Same-length amplicons can be indistinguishable in downstream sequencing/scoring
Current workflow priority
Define discrimination scoring logic
Define data structures
Build biological data extraction pipeline
Validate on small examples
Scale to full dataset
Key files
terminal_windows.tsv: conserved genomic window data, orientation, gene/end metadata
gene_pair_amplicon_summary or similar summary file: expected amplicon lengths per genome/pair
Genome FASTA files: bacterial genome sequences
Human genome FASTA: used for biological feasibility comparisons if needed
Current coding tasks
Build FASTA parsing utilities
Build sequence search function
Handle orientation transforms correctly
Validate matches against known amplicon lengths
Extract amplicon sequence between matched primer positions
Prepare scoring-ready output structure
Rules for Codex
Do not assume biological rules unless explicitly stated here
Keep implementation modular
Always propose a smallest test case first
When debugging, explain likely failure points in order
Distinguish clearly between “confirmed” and “assumed” logic

You are helping me implement and debug my GenomicGapID project in VS Code.

Project goal:
Analyze bacterial FASTA genomes, map conserved genomic regions from terminal_windows.tsv, identify primer-scale subsequences and their genomic locations, extract the expected amplicons, and support downstream discrimination scoring.

What I need from you:
- Help me write clean, modular Python code
- Help me debug errors step by step
- Suggest small test cases before scaling up
- Keep biological assumptions explicit and separate from coding logic
- Do not make biological decisions for me unless I state the rule explicitly
- When uncertain, flag assumptions clearly

Core biological rules currently in use:
- FASTA genomes are treated as forward-strand genomic sequence
- Sequences with negative orientation must be reverse complemented before comparison
- terminal_windows.tsv contains larger genomic conserved regions, not the final primer sequences
- Final primers are expected to be short subsequences, around 14 to 15 bp
- Reverse primer matching requires careful handling of reverse complement / inverse complement logic
- Amplicon lengths should be checked against gene_pair_amplicon_summary when validating matches

Coding preferences:
- Use Python
- Prioritize modular functions over one large script
- Separate parsing, sequence matching, validation, and scoring into different functions/files
- Add clear comments and docstrings
- Suggest simple print/debug checks when troubleshooting
- Prefer test-small-then-scale workflow

Current repository structure:
- data/
- scripts/
- notebooks/
- results/
- docs/
- .venv environment

Important implementation philosophy:
- Define scoring and data structures first
- Then build extraction pipeline
- Then optimize and scale

At the start of this session:
1. Read the session notes below
2. Summarize your understanding in 5 bullets
3. Identify the immediate coding task
4. Propose the smallest next implementation step


