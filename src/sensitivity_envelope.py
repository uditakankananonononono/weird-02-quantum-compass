#!/usr/bin/env python3
"""W02 post-hoc sensitivity envelope (judge round-5 weakness 2).
Labeled SENSITIVITY ANALYSIS - results for the frozen model were already seen;
this sweeps the lifetime-model parameter class over the LOCKED S(tau) curve.
Question: what fraction of the plausible geometry->lifetime parameter space
recovers the Xu ordering (robin erCry4 ABOVE both ClCry4 and GgCry4)?
S(tau) itself spans all lifetimes, so kS/kT variation is subsumed by the curve.
"""
import json, numpy as np

d = json.load(open('results/h1_attempt4_stage1.json'))
sp, curve = d['species'], d['tau_grid_S']
taus = sorted(float(t) for t in curve)
S = np.array([curve[f'{t:g}']['S'] if f'{t:g}' in curve else curve[str(t)]['S'] for t in taus])
def S_at(tau):
    return float(np.interp(tau, taus, S))

dterm = {s: sp[s]['d_term_WT'] for s in ('erCry4', 'ClCry4', 'GgCry4')}

def tau_of(d, beta, tau0, D0=13.0):
    # tau0 = lifetime (us) at reference distance D0 for the anchor
    K0 = 1e6 / tau0
    k_back = K0 * np.exp(-beta * max(d - D0, 0.0))
    return 1.0 / max(k_back * 1e-6, 1e-12)

rng = np.random.default_rng(20260926)
betas = np.linspace(0.3, 2.5, 45)          # /Angstrom, generous BET range
tau0s = np.linspace(0.5, 3.0, 51)          # us anchor, wider than locked [1.0,1.3]
jitters = [0.0, 0.25, 0.5, 1.0, 1.5]       # Angstrom structure uncertainty (fold xval 1.2-1.8A)
n_struct = 40                              # structure samples per jitter level

total, xu_ok = 0, 0
max_adv = 0.0         # max (S_er - max(S_cl,S_gg)) anywhere
worst_def = 0.0       # min (S_er - max(...)) i.e. strongest anti-prediction
per_jitter = {}
for j in jitters:
    t_j = ok_j = 0
    for _ in range(n_struct if j > 0 else 1):
        ds = {s: dterm[s] + (rng.normal(0, j) if j > 0 else 0.0) for s in dterm}
        for b in betas:
            for t0 in tau0s:
                Se = S_at(tau_of(ds['erCry4'], b, t0))
                Sc = S_at(tau_of(ds['ClCry4'], b, t0))
                Sg = S_at(tau_of(ds['GgCry4'], b, t0))
                adv = Se - max(Sc, Sg)
                max_adv = max(max_adv, adv); worst_def = min(worst_def, adv)
                t_j += 1; total += 1
                if Se > Sc and Se > Sg:
                    ok_j += 1; xu_ok += 1
    per_jitter[f'{j}A'] = {'n': t_j, 'xu_ordering_frac': round(ok_j / t_j, 6)}

res = {
 'label': 'SENSITIVITY ANALYSIS (post-hoc, frozen-model results already seen; judge round-5 weakness 2)',
 'grid': {'beta_per_A': [0.3, 2.5, 45], 'tau0_us': [0.5, 3.0, 51],
          'structure_jitter_A': jitters, 'structure_samples': n_struct,
          'D0_A': 13.0, 'S_tau_curve': 'locked attempt-4 16-point curve (kS/kT variation subsumed: curve spans all lifetimes)'},
 'total_parameter_points': total,
 'xu_ordering_fraction': round(xu_ok / total, 6),
 'max_robin_advantage': round(max_adv, 6),
 'max_robin_deficit': round(worst_def, 6),
 'per_jitter': per_jitter,
 'conclusion': ('geometry-lifetime model class constrained' if xu_ok == 0
                else f'Xu ordering recovered at {xu_ok/total:.4%} of parameter space'),
}
json.dump(res, open('results/h1_sensitivity_envelope.json', 'w'), indent=1)
print(json.dumps({k: res[k] for k in ('total_parameter_points','xu_ordering_fraction','max_robin_advantage','max_robin_deficit','conclusion')}, indent=1))
print(json.dumps(per_jitter, indent=1))
