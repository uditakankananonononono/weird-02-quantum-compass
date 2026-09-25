# W02 PRE-REGISTRATION AMENDMENT - H1 ATTEMPT 2 (locked 2026-09-25 23:52 IST, BEFORE any attempt-2 scoring)

Context: attempt-1 (structure-derived amplitude) FAILED the locked gate (AUC 0.3414 vs >=0.75). Diagnosis: chain geometry conserved across birds; amplitude cannot separate. This amendment locks the attempt-2 test BEFORE it is run. No attempt-2 scoring has been performed as of this lock.

## H1-A2 (locked)
Sequence-level discriminators from the LOCKED feature set v1 (116 features, frozen earlier today) separate migratory from sedentary birds' Cry sequences.

### Locked scoring protocol
1. Features (all from frozen feature set v1, no new features, no feature selection after seeing labels):
   a. Trp-tetrad completeness: count of chain-Trp columns (6PTZ W395/W372/W318/W369 MSA columns) occupied by Trp in the row (0-4).
   b. Terminal-Trp presence: binary, Trp at the W369-equivalent MSA column.
   c. FAD-pocket conservation: fraction of the locked pocket-column set identical to the 6PTZ residue.
   d. Trp-spacing profile: alignment-column distances between successive chain-Trp columns (3 numbers).
2. Classifier: L2-regularized logistic regression (C=1.0, fixed), standardized features, leave-one-out cross-validation over the frozen panel (mig vs sed ONLY; partial class excluded from attempt-2 and its absence diagnosed separately).
3. GATE: LOO-CV AUC >= 0.75 (same threshold family as attempt 1). Also report: label-shuffle null (1000 permutations) - locked expectation: null median ~0.5, gate requires observed AUC above the 95th percentile of the null.
4. No re-tuning after seeing results. Any further change = new amendment.
5. If H1-A2 FAILS: attempt-3 candidates (to be locked only if needed): (i) integrate inclination-range/migration-distance traits (AVONET) as labels instead of binary mig/sed; (ii) behavioral-data-anchored scoring (Hiscock 2016 per-species magnetic thresholds); (iii) escalate to parent if sequence-level angles exhaust.

## Panel hygiene locked with this amendment
- Partial-aptitude class: absent from attempt-1 scored rows (group-prefix parse). Attempt-2 will re-derive group labels from the frozen panel CSV (not MSA prefixes) and REPORT partial-class rows separately; they remain excluded from the AUC gate.
- 400-aa fold cap: irrelevant to attempt-2 (no folds used); any future structure attempt must handle C-terminal truncation explicitly.
