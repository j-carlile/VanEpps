PROJECT CONTEXT - GenomicGapID (Current State)

Goal
Analyze bacterial genomes (FASTA), map conserved regions from terminal_windows.tsv, identify primer-binding locations, extract amplicons, and compute discrimination scores.

Confirmed Current State
1. terminal_windows.tsv sequences (~400 bp) are not reliable exact substrings of genomes.
   They should be treated as approximate conserved source regions rather than literal match targets.
2. Short exact matches in the ~14-20 bp range do exist inside some of these regions.
   These are anchors, not true biological boundaries.
3. Anchor-only matching is not sufficient to recover robust primer positions or amplicon lengths across genomes.
4. Full-window exact matching is too strict, while short exact matching alone is too permissive.
5. Raw anchor-to-anchor distances often disagree with expected amplicon lengths.
6. In some isolated cases, offset-adjusted anchor distances can move closer to expected lengths, but this does not generalize reliably.
7. Multi-genome checks showed:
   - some genomes have no valid same-contig pairings
   - some genomes have pairings but with implausible distances
   - rpoC anchors are often more stable than rpoB anchors
   - allowing full-window anchor search helped one genome but did not solve the problem generally

Default Assumptions For This Session
- terminal_windows.tsv entries are approximate conserved regions
- short exact matches are anchors, not final primer boundaries
- primer locations cannot currently be inferred by simple exact matching alone
- mapping anchors to true primer positions is still an open problem
- orientation still matters:
  - '+' means search the sequence directly
  - '-' means search the reverse complement against the forward-strand FASTA

Biological Interpretation
- The biologically relevant primer-binding sites are likely smaller conserved subsequences within the larger ~400 bp source regions.
- The rest of the conserved region may still exist in a genome but be mutated enough that exact full-window matching fails.
- Therefore, a short anchor can indicate approximate regional location without defining the true primer boundary.

Core Unresolved Problem
There is currently no reliable method to map:
approximate conserved regions -> precise primer-binding locations
in a way that is:
- consistent across genomes
- biologically sensible
- able to reproduce expected amplicon lengths

What Not To Assume Going Forward
- do not assume terminal_windows.tsv sequences are exact genome substrings
- do not assume primer locations can be recovered by simple exact substring matching
- do not assume anchor hits define the actual primer boundaries

Current Direction
Move toward:
- stronger pairing logic
- explicit biological constraints
- alignment-based or structure-aware regional matching
- methods that use anchors as seeds, then refine the surrounding region

Avoid relying on:
- raw full-window exact matching
- raw short-substring matching without downstream validation

Session Understanding Added By Codex
This session established that the project should no longer treat `terminal_windows.tsv`
rows as exact genome substrings. The correct working model is that these rows define
approximate conserved source regions, and short exact matches inside them are only anchor
signals. Those anchor signals can sometimes localize a region and occasionally produce an
adjusted distance that looks biologically plausible, but anchor-only logic does not
generalize across genomes and is not reliable enough to define primer positions or final
amplicons by itself.

The strongest conclusion from the debugging work is that the core problem is no longer
"find exact matches", but instead "map approximate conserved regions to biologically valid
primer-binding positions using stronger constraints." Exact 400 bp matching is too strict,
while isolated 14-20 bp matches are too permissive. The pipeline therefore needs an
intermediate method that uses short exact matches only as seeds and then reasons about the
surrounding region.

Best Next Step For The Next Session
Implement a seed-and-extend regional alignment workflow for one gene pair at a time.

Why this is the right next step
- It keeps the useful part of the current work: short exact anchors can still localize
  approximate regions.
- It fixes the current weakness: anchors alone do not define primer boundaries.
- It is more biologically realistic because it allows the conserved region to be present
  with substitutions and indels.
- It creates a path toward reproducible pairing logic that can later be scaled across
  genomes and fed into discrimination scoring.

Recommended Implementation For Next Session
1. Build a new focused debug/prototype script in `scripts/debug/` for one pair, likely
   `rpoB_3-rpoC_5`, using one or a few genomes first.
2. For each terminal window:
   - find short exact seed matches anywhere in the oriented 400 bp searched sequence
   - keep candidate seed hits with uniqueness and biological plausibility metadata
3. Around each candidate seed hit:
   - extract a projected genome region large enough to contain the full approximate
     conserved source region
   - align that projected region against the terminal-window searched sequence
   - compute alignment quality metrics such as percent identity, mismatch count, and gap
     count
4. For gene-pairing:
   - only compare same-contig candidate regions
   - rank candidate pairs using a combination of:
     - alignment quality for both source regions
     - expected gene order / direction
     - expected amplicon length proximity
     - hit uniqueness
5. Print a ranked table of candidate pairings instead of a single raw match result.

Minimum Useful Deliverable Next Session
A script that, for one pair and one genome:
- finds seed hits
- extends to candidate regional alignments
- scores candidate pairings
- reports the top-ranked candidate amplicon with enough detail to explain why it was
  chosen

Important Default Assumption For Next Session
Do not restart from exact-substring assumptions. Continue assuming:
- terminal windows are approximate conserved regions
- short exact matches are anchors, not boundaries
- mapping anchors to true primer positions is still an open algorithmic problem
