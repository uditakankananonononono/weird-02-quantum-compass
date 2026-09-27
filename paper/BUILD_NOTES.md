# Manuscript build note

- Build: `cd paper && pdflatex -interaction=nonstopmode -halt-on-error main.tex` twice. Compiled and visually inspected on 2026-09-27.
- Source base: committed main tip `6a7ef5f15a4e9ac77eaed65cca97b938503a93d3`; working branch `paper-build` is not merged into main.
- Current PDF: 23 pages, A4, 12-point text. This draft is incomplete and does not reach the 50-page submission floor.
- Font substitution: `mathptmx` supplies Times-compatible text and math. Genuine licensed Times New Roman is not installed here, so this PDF must not be described as Times New Roman. Replace or rebuild on a system with the licensed font if that exact face is mandatory.
- The top rule is a design aid, not a science result. The data appendix is traceable to committed result paths named in the manuscript. Negative results and [PENDING] gates have not been softened.

- The exploratory post-hoc perturbation and geometric negative-control appendix transcribes committed attempt-5b files; it does not revise the locked gate.

- Rebuilt from science main 6a7ef5f on September 27; the main-branch manuscript is authoritative and the appendix does not reverse later negative gates.
- The 7d sampling robustness passed, but the later locked 7e Cry2 specificity control failed; the 7c mechanism-specific novelty claim remains downgraded.
- Attempt-4 gate ledger transcribes committed `results/h1_attempt4_gates.json`, preserving ALL=false and each sub-gate.
- Attempt-4 intermediate lifetime-to-score grid and three-species mapping transcribe committed `results/h1_attempt4_stage1.json`, without claiming measured kinetics.
- 16-point lifetime transfer grid and three-species input ledger are printed compactly on page 18; an earlier float-only build had a sparse extra page and was replaced.

- Attempt 7f and 7g falsifications are now folded into the science-main headline pivot and Section 17 of the appendix. Specificity stays closed at the tested descriptor level; the family-level canalization result and physical constraint are retained.

- The existing Udita Phookan byline is retained. PDF Author metadata is absent.
