# Queue #8: ML benchmark vs physics descriptors (attempt8g) - REPORT ONLY

Amendment 2026-09-27 12:31 IST. Delegated-executor computation at ref 4845c5e, uncommitted.

## Headline
Neither generic classifier reaches the detectable-effect bound on the committed physics descriptor matrix:
- Logistic regression (C=1.0, LOO): AUC **0.4021** - reproduces the committed attempt-3a reference (0.402) exactly.
- Random forest (100 trees, LOO): AUC **0.6234**.
Both are BELOW the attempt8c MDE floor 0.6424, i.e. below the panel's detectable-effect bound.

## Setup (locked)
Features = frozen v1 physics descriptors, the attempt-3a set the 8c classifier_auc MDE is calibrated on: tail_after_ref, charge_pH7, pI, aromatic_frac, kd_mean, length, ident_erCry1, trp_total (data/cry1_features_v1.json). Panel = locked 116-species Cry1 panel (36 migratory / 80 sedentary, same as 7f-b/8c). Seed 260927. Leave-one-out CV, pooled held-out probabilities, single AUC. PLM-embedding arm: UNAVAILABLE (no local embeddings tool; multi-GB downloads barred) - disclosed, not run.

## Locked framing
This benchmark bounds descriptor information content only. No migratory-vs-sedentary prediction is claimed as a project result (locked anti-goal). The random forest's margin over logistic (0.6234 vs 0.4021) stays under the floor and earns no phenotype claim.
