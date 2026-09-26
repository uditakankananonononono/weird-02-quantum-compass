# W02 PREREG ATTEMPT-5A - locked 09:08 IST 2026-09-26, BEFORE any 5a metric evaluation
Source: ChatGPT redirection consult round 3, rank-1 (judge/consult_round3_attempt5_redirection.txt); parent approved 5a+5b as separate locks 09:07. Attempts 1-4 stand untouched as honest FAILs.

## Hypothesis (constraint result, framed positively)
A frozen, physically explicit structure-only lifetime model (M1.2: back-ET from terminal Trp, beta=1.4/A, Xu-anchored) CANNOT recover the experimentally reported avian Cry4 magnetic-sensitivity ordering. The positive contribution: a physics-constrained computational test rules out terminal-Trp geometry as a sufficient explanation and narrows the mechanism space to non-geometric sources (recombination kinetics, hyperfine environment, dynamics).

## Locked gates
- A1 (primary): Spearman rho between model sensitivity score (S at each species' locked M1.2 modeled tau, from attempt-4 stage-1 outputs, frozen) and Xu 2021 experimental ordering (erCry4 most sensitive; ClCry4/GgCry4 both less) is <= 0, AND bootstrap 95% CI (10k resamples over species pairs) excludes a positive correlation.
- A2: robin does NOT rank above both controls under the frozen model.
- A3: permutation test of rank agreement vs Xu ordering - HONEST RESOLUTION NOTE: with 3 species there are only 6 orderings, so the minimum achievable exact p is 1/6 = 0.167; the consult's p<0.05 is unreachable at n=3. Locked adaptation: full exact enumeration of all 6 permutations reported; gate = observed agreement is WORSE than random (agreement score in the lower half of the permutation distribution, i.e. robin strictly below both controls).
- A4: magnitude reported: S ratio controls/robin under the frozen model (no gate - descriptive).
Success = A1 AND A2 AND A3. Honest negative recorded if the model accidentally recovers the ordering (which would instead support geometry-sufficiency).
