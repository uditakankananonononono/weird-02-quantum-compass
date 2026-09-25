# W02 PRE-REGISTRATION - H1 ATTEMPT 4: PHYSICS-FORWARD MFE REPLICATION (locked 2026-09-26 00:12 IST, BEFORE any attempt-4 computation)

Directive: user (12:08 AM, parent-relayed): redirect to positive angle; use ChatGPT for redirection/validation. Design sharpened by ChatGPT consult round 1 (judge/consult_round1_attempt4_design.txt). Attempts 1-3 remain honest record.

## Framing (locked)
Can an independently parameterized radical-pair model, using only protein structure + physically motivated kinetic parameters, predict the interspecies and fourth-tryptophan dependence of Cry4 magnetic sensitivity WITHOUT fitting to the experimental phenotype (Xu et al. 2021)?

## Independence protocol (anti-circularity, per ChatGPT critique)
- Parameter-generation stage: chain-edge distances from structures (ESM/AFDB, as in attempt-1 pipeline); lifetime tau from the M1.2 back-ET model locked 2026-09-25 (tetrad 1.1us / triad 0.14us anchors from literature kinetics, NOT from Xu's species comparison; beta=1.4 A^-1).
- Phenotype-comparison stage: simulate, then compare to Xu WITHOUT any refitting. Any post-hoc parameter change invalidates the attempt.

## Locked metric
Sensitivity curve A(B) = anisotropy amplitude (Eq. 4 of paper) vs field strength B. Lifetime-weighted sensitivity integral: S = integral of A(B) over B in [25, 65] uT (geomagnetic natural-variation window; justification: Earth's field varies 25-65 uT globally; locked before computation). Units: fractional yield anisotropy x uT.

## Locked predictions and gates (ALL must hold for PASS)
- P1 (WT ordering): S(erCry4) > S(ClCry4) AND S(erCry4) > S(GgCry4) — Xu's measured ordering.
- P2 (WDF mutant direction): virtual erCry4 W369F: larger short-time-window response (S_short over identical integration with triad lifetime) but the LONG-lived component vanishes: specifically, with triad-anchor tau, peak A increases relative to WT-with-triad-tau while lifetime-weighted S drops toward triad values. Direction must match Xu's "larger but transient" measured result.
- P3 (ablation, anti-circularity) - AMENDED 00:13 before any computation: the original (ii) was ill-defined in M1.2 (tau is derived from distances; no geometry-independent species tau exists). Replacement: (i) equalized-d_term control (all species assigned the grand-mean d_term): species S-differences must shrink by >=80% vs P1 — proving differences arise from structure-derived kinetics, not numerical artifacts; (ii) anchor-robustness control: P1 ordering must be invariant across tetrad-anchor tau0 in the literature range [1.0, 1.3] us (7 values) — proving the result is not an artifact of one anchor choice.
- P4 (robustness): P1 ordering survives in >=95% of 100 perturbation draws (chain distances +/-10% uniform, tau +/-20% uniform; sensitivity via locked tau-grid interpolation, no new solves beyond the grid).
- P5 (field window): primary curve B in 0-200 uT (justification: contains the full geomagnetic regime and the Hiscock response region); secondary wide curve 0-2 mT computed and reported.

## Compute plan (locked)
B grid: 12 points 0-200 uT + 8 points 0.2-2 mT; 65-orientation Lebedev grid as validated. Configs: 3 species WT + erCry4-WDF + 2 ablations = 6 primary + tau-grid (15 taus, 3 species) for P4. Est. ~40 min CPU serial. Results to results/h1_attempt4_curves.json. If ANY gate fails: recorded honestly; paper reports outcome; further redirection only via new user-directed ChatGPT consult.


## COMPUTE AMENDMENT (00:59, before any gate evaluation; logged)
Stage-1 at the locked resolution (65 orientations x 20 B x 25 taus) measured ~15 min/tau (~6h total) - intractable before the deadline. Numerical resolution amended: orientation grid 65 -> 25 (every 3rd Lebedev point, same sphere coverage), primary B grid 12 -> 9 points (same 0-200 uT range), parallelized across taus. GATES, METRIC, WINDOW, TAU GRID UNCHANGED. The 2 partially computed tau points (0.1, 0.2) are discarded; the curve is recomputed uniformly at the amended resolution. This changes compute granularity only, not any decision rule.


## COMPUTE AMENDMENT 2 (03:00, before gate evaluation; logged)
The sandbox has recycled processes repeatedly tonight (PIDs reused, daemons killed mid-run), giving ~2 tau points/hour effective throughput. Tau grid amended 25 -> 16 points (log-dense 0.1-10 us, same span; all gate taus - species WT/triad - fall inside). Engine made incremental (per-tau JSONL resume) + supervisor-restarted. Gates, metric, window, anchors UNCHANGED. Previously computed partial curve points are kept only if computed at the current amended resolution; otherwise recomputed.
