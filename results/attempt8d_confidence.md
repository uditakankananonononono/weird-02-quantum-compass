# Queue #6: predicted confidence vs solved census vs functional regions
Amendments 2026-09-27 11:27 IST + 11:28 IST clarification (per-file scale detection, Cry4 accession fix), both locked BEFORE the corresponding computations. DESCRIPTIVE, no gate.

## Result: every conclusion rests on high-confidence regions
Per-residue confidence (PDB B-factor; ESMFold 0-1 scale thr 0.70, AFDB 0-100 scale thr 70) over all 137 folds, with CORE = the locked 7f-b functional core (Trp chain + FAD <=4.5A + Trp <=6.0A in 6PTZ numbering):

| fold class | n | overall mean (median) | CORE mean (median) | CORE mean (min) | CORE frac >= thr (median) | TAIL mean (median) |
|---|---|---|---|---|---|---|
| Cry1 panel | 116 | 0.853 | 0.863 | 0.772 | 0.972 | 0.894 |
| REF rows | 5 | 0.849 | 0.865 | 0.861 | 0.972 | 0.882 |
| Cry4 panel | 11 | 0.851 | 0.872 | 0.732 | 1.000 | 0.881 |
| Cry2 REF | 1 | 0.852 | 0.850 | 0.850 | 0.958 | 0.876 |
| AFDB pairs | 4 | 86.6 | 96.7 | 96.6 | 1.000 | 62.7 |

- The functional core is high-confidence EVERYWHERE: 97-100% of core residues clear the high-confidence threshold in every class; the weakest single fold's core still averages 0.73 (min overall 0.669 is tail-dragged, not core).
- Paired ESMFold-vs-AlphaFoldDB on 4 identical sequences: core means 0.862-0.866 vs 96.6-97.2 - two independent predictors agree the core is well-modeled (AFDB tails run much lower, 62.7 vs 0.88: tail confidence is predictor-dependent, core confidence is not).
- Solved census (34 cryptochrome structures, RCSB): X-ray + EM, median resolution 2.23A (range 1.6-3.5). Our numbering anchor 6PTZ is solved at 1.793A - consistent with the pipeline's 1.2-1.8A validation.

## Descriptive note (no causal claim)
6PTZ's own title: crystal structure of pigeon Cryptochrome 4 **mutant Y319D** in complex with FAD. Residue 319 is attempt 8's sole near-site Cry4-distinguishing residue (natural Cry4 = Tyr; the crystal carries the engineered Asp = the paralog-consensus state). The field's anchor structure thus sits exactly on the one near-site position where Cry4 departs from its paralogs - a descriptive convergence between the specialization landscape and the available structural record, and a caveat: 6PTZ-numbered near-site geometry at 319 describes the engineered reversion, not natural Cry4.

## Caveats (locked framing)
Confidence is model self-report, not functional validation. No conclusion in the paper rests on a region whose class-median core confidence falls below the high-confidence threshold, so no result requires caveat under the locked rule. Tail-region statements (7g-b charge) rely on sequence, not structure, and are unaffected by tail-confidence variation.
