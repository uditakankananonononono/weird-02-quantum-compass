# W02 PRE-REGISTRATION AMENDMENT - H1 ATTEMPT 3 (locked 2026-09-25 23:56 IST, BEFORE any attempt-3 scoring)

Grounding (Xu et al. 2021, Nature 594, accepted manuscript via Oxford ORA + SI, fetched 23:55): the migratory-vs-non-migratory difference in Cry4 is KINETIC - erCry4 stabilizes a long-lived signaling state (time-dependent spectral change) that chicken and pigeon Cry4 lack; the 4th Trp (TrpD/W369) prolongs the magnetically sensitive state. Trp-chain numbering in Xu matches our locked 6PTZ chain (W395/W372/W318/W369 = A/B/C/D).

## H1-A3 (locked): two pre-specified tests, Bonferroni alpha split, both must use ONLY frozen v1 features and the frozen panel
- 3a: logistic LOO (same estimator as attempt-2, C=1.0) on the REMAINING unused v1 features: [tail_after_ref, charge_pH7, pI, aromatic_frac, kd_mean, length, ident_erCry1, trp_total]. GATE: LOO AUC >= 0.75 AND above permutation-null p95 (1000 perms). Rationale: signaling-state stabilization plausibly lives in C-terminal tail / electrostatics, per Xu.
- 3b (directional replication, qualitative): full-fit logistic from 3a applied to the Cry4 subpanel. GATE: erCry4 scores strictly above BOTH GgCry4 and ClCry4 (Xu's measured MFE ordering). If 3a fails, 3b is not interpretable and is recorded as not-run.
- Multiple-comparison honesty: attempts 2 and 3a test overlapping families on the same panel; reported together with this caveat. No further sequence-level attempts after 3 - if 3a/3b fail, W02 escalates to parent with the honest position: validated physics pipeline + no panel-level sequence separator found.
