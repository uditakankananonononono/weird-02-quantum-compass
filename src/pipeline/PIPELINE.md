# W02 hypothesis-testing pipeline: sequences+labels in -> evidence scores out
Reproducibility map over committed artifacts (queue #17; amendment 2026-09-27 13:29 IST). Each stage: inputs -> script -> outputs -> locking amendment.

1. PANEL ASSEMBLY: NCBI/Ensembl Cry1+Cry4 sequences -> data/cry_aves_proteins.fasta, cry4_subpanel.fasta, cry4_subpanel_members.json (provenance ledger) -> original pre-registration.
2. FEATURE EXTRACTION: sequences -> data/cry1_features_v1.json (per-species: identity to erCry1/erCry4, pI, charge_pH7, kd_mean, aromatics, Trp geometry) -> locked v1.
3. BINARY SCREEN (attempts 1-3): features + committed labels -> results/h1_scoring.json (AUC 0.3414, 1000-perm nulls) -> pre-registered gates.
4. STRUCTURE ARMS (7a-7g): AlphaFold/6PTZ RMSD census, Cry2 control -> results/fold_rmsd_vs_6ptz.json, attempt7* -> 7-series amendments.
5. ML BENCHMARK (8g): src/attempt8g_*.py -> results/attempt8g_ml_benchmark.json (logistic 0.4021 reproduces 3a; RF 0.6234 < MDE 0.6424) -> 12:31 amendment.
6. STRUCTURE AGREEMENT (8f): -> results/attempt8f_structure_agreement.json (r=0.093) -> 12:31 amendment.
7. CONTINUOUS TRAITS (3'): AVONET join -> results/h_continuous_traits.json (18/18 ns) -> 13:22 amendment.
8. ELECTROSTATICS (12): -> results/h_electrostatic_map.json (all ns; tail localized) -> 13:14+13:16 amendments.
9. POWER (15): -> results/h_power_analysis.json (MDA 0.6453) -> 13:17 amendment.
10. SELECTION SCREEN (11): src/attempt8h_dnds.py + data/cds (sha256 ledger) -> results/attempt8h_dnds.json (p=0.3792) -> 12:59+13:15+13:28 amendments.
Every stage's amendment sits in PRE-REGISTRATION.md before its compute; every negative is verbatim in results/*.md.
