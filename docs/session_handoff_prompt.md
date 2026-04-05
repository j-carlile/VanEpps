# Session Handoff Prompt

Use this prompt at the end of a work session to create a durable handoff note and preserve the current project state for the next session.

```text
I want you to create a session handoff note that preserves today’s progress and becomes the default context for the next session.

Before writing anything:
1. Read all relevant context files in `docs/codex/` first.
2. Read the most recent relevant file(s) in `docs/daily_context/`.
3. Review the work we completed in this session.
4. Read any important local files we updated or relied on.

Then do all of the following:

1. Write a clean, structured session summary.
2. Save it to:
   docs/daily_context/<SESSION_NOTES_FILENAME>.md
3. After saving the file, paste the full summary back to me in the chat.

The note should include these sections:

1. Project Goal
- Briefly restate the current project goal in plain language.

2. What We Worked On
- Summarize the specific coding, debugging, and analysis tasks we performed this session.

3. Confirmed Findings
- List what we now know from actual testing or code inspection.
- Separate confirmed behavior from guesses.

4. Biological Assumptions
- List the biological rules or assumptions currently in use.
- Clearly distinguish assumptions from experimentally supported conclusions.

5. Coding / Pipeline Assumptions
- List the implementation assumptions that future work should use by default.

6. What Failed or Did Not Generalize
- Record unsuccessful approaches, dead ends, or methods that gave misleading results.
- Explicitly say what should NOT be assumed next time.

7. Unresolved Problems
- List the open technical and biological problems still blocking progress.

8. Current Best Interpretation
- Explain the current working model in plain language.
- This should become the default context for the next session.

9. Next Best Step
- Recommend the single best next technical step.
- Explain why it is better than the approaches already tried.

10. Concrete Implementation Plan For Next Session
- Give a short, practical plan for what to code next.
- Prefer a smallest-test-first approach.

11. Important Files
- List the main files created, modified, or relied on this session.

Important instructions:
- Do not revert to outdated assumptions from earlier in the project.
- Treat today’s conclusions as the default working context unless explicitly contradicted later.
- Be explicit about:
  - what is confirmed
  - what is inferred
  - what is still unresolved
- Keep the note readable and useful for restarting work later.
- Save the note, then show me the final contents in chat.

Here is the filename to use:
<SESSION_NOTES_FILENAME>
```

Daily context filename convention:

```text
Use this naming convention for daily context notes:
<M>_<D>_<YY>_notes.md

Rules:
- month is numeric
- day is numeric
- year is the last 2 digits
- no leading zeros are required
- always end with `_notes.md`

Examples:
4_4_26_notes.md
11_17_26_notes.md
```

Suggested filename example:

```text
4_4_26_notes.md
```
