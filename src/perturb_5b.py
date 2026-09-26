#!/usr/bin/env python3
"""W02 EXPLORATORY 5b robustness (judge round-6 weakness 2): is the ~44x
persistence extension stable under stochastic lifetime perturbation?
Locked definition (PREREG_ATTEMPT5B): I = A_peak x tau x (1-exp(-T/tau)), T=10us.
A_peak from locked extended-grid values (tetrad 0.192125, triad 0.192123 -
near-identical, knife-edge disclosed). tau_tetrad=0.1495us, tau_triad=0.0034us.
Post-hoc exploratory; locked 5b gates unchanged."""
import json, numpy as np
T = 10.0
A_TET, A_TRI = 0.19212547395153934, 0.19212289889122622
TAU_TET, TAU_TRI = 0.1495, 0.0034
def I(A, tau): return A * tau * (1 - np.exp(-T / tau))
locked_ratio = I(A_TET, TAU_TET) / I(A_TRI, TAU_TRI)
rng = np.random.default_rng(20260926)
res = {'label': 'EXPLORATORY post-hoc (judge round-6); locked 5b gates unchanged',
       'locked_ratio': locked_ratio, 'T_us': T}
for sigma in (0.05, 0.10, 0.20, 0.30):
    tt = TAU_TET * np.exp(rng.normal(0, sigma, 20000))
    tr = TAU_TRI * np.exp(rng.normal(0, sigma, 20000))
    r = I(A_TET, tt) / I(A_TRI, tr)
    res[f'lognormal_sigma_{sigma}'] = {
        'ratio_median': round(float(np.median(r)), 1),
        'ratio_p5': round(float(np.percentile(r, 5)), 1),
        'ratio_p95': round(float(np.percentile(r, 95)), 1),
        'frac_ratio_above_1': round(float((r > 1).mean()), 4),
        'frac_ratio_above_locked_gate_1.25': round(float((r > 1.25).mean()), 4)}
# worst-case direction: tetrad tau low, triad tau high, at 30% each
worst = I(A_TET, TAU_TET*0.7) / I(A_TRI, TAU_TRI*1.3)
res['worst_case_30pct_opposed'] = round(float(worst), 1)
res['conclusion'] = ('persistence extension is structurally robust: ratio stays far above the locked 1.25 gate '
 'at all perturbation levels (extension is ~tau-ratio driven; selectivity preservation is the knife-edge, already disclosed)')
json.dump(res, open('results/h1_attempt5b_perturbation.json','w'), indent=1)
print(json.dumps(res, indent=1)[:900])
