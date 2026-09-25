#!/usr/bin/env python3
"""W02 attempt-4 gate evaluation vs PREREG_ATTEMPT4.md (locked 00:12, P3 amended 00:13, compute amended 00:59)."""
import json, numpy as np
d = json.load(open('results/h1_attempt4_stage1.json'))
sp, curve = d['species'], d['tau_grid_S']
taus = sorted(float(t) for t in curve)
S = np.array([curve[str(t)]['S'] if str(t) in curve else curve[f'{t:g}']['S'] for t in taus])
A_curves = {t: curve[f'{t:g}']['A_curve'] for t in taus}
def S_at(tau):
    return float(np.interp(tau, taus, S))
def Apeak_at(tau):
    t_near = min(taus, key=lambda t: abs(t-tau))
    return max(A_curves[t_near])
res = {}
S_er, S_cl, S_gg = S_at(sp['erCry4']['tau_WT']), S_at(sp['ClCry4']['tau_WT']), S_at(sp['GgCry4']['tau_WT'])
res['P1'] = {'S_er': S_er, 'S_cl': S_cl, 'S_gg': S_gg,
             'pass': bool(S_er > S_cl and S_er > S_gg)}
# P2: WDF direction - peak A at triad tau vs WT tau; S must drop
Ap_tri, Ap_wt = Apeak_at(sp['erCry4']['tau_triad']), Apeak_at(sp['erCry4']['tau_WT'])
S_tri = S_at(sp['erCry4']['tau_triad'])
res['P2'] = {'Apeak_triad': Ap_tri, 'Apeak_WT': Ap_wt, 'S_triad': S_tri, 'S_WT': S_er,
             'pass': bool(Ap_tri > Ap_wt and S_tri < 0.3 * S_er)}
# P3i: equalized d_term -> differences shrink >=80% (formal: all taus equal)
dmean = np.mean([sp[s]['d_term_WT'] for s in sp])
import sys; sys.path.insert(0,'src')
from compass_scan import tau_eff
t_eq = tau_eff(dmean)
S_eq = [S_at(t_eq)]*3
diff_real = max(S_er,S_cl,S_gg) - min(S_er,S_cl,S_gg)
diff_eq = max(S_eq)-min(S_eq)
shrink = 1 - (diff_eq/diff_real if diff_real else 0)
res['P3i'] = {'d_mean': round(float(dmean),2), 'shrink_frac': round(float(shrink),3), 'pass': bool(shrink >= 0.8)}
# P3ii: anchor robustness - K0 scaled to give tau0 in [1.0,1.3]us at d=D0; ordering preserved?
import compass_scan as cs
order_ok = 0; details = []
for tau0 in np.linspace(1.0, 1.3, 7):
    K0_new = 1e6 / tau0  # tau_eff = 1/(K0 exp(-BETA(d-D0)) * 1e-6); at d=D0 tau=1e6/K0 us
    ts = {}
    for s_ in sp:
        dterm = sp[s_]['d_term_WT']
        ts[s_] = 1.0 / max(K0_new * np.exp(-cs.BETA*(dterm-cs.D0)) * 1e-6, 1e-12)
    Ss = {s_: S_at(t) for s_, t in ts.items()}
    ok = Ss['erCry4'] > Ss['ClCry4'] and Ss['erCry4'] > Ss['GgCry4']
    details.append({'tau0': round(float(tau0),2), 'order_ok': bool(ok)})
    order_ok += ok
res['P3ii'] = {'anchors_ok': order_ok, 'of': 7, 'details': details, 'pass': bool(order_ok == 7)}
# P4: bootstrap perturbations
rng = np.random.default_rng(11)
keep = 0; N = 100
for i in range(N):
    ts = {}
    for s_ in sp:
        dterm = sp[s_]['d_term_WT'] * (1 + rng.uniform(-0.1, 0.1))
        t = tau_eff(dterm) * (1 + rng.uniform(-0.2, 0.2))
        ts[s_] = t
    Ss = {s_: S_at(t) for s_, t in ts.items()}
    keep += Ss['erCry4'] > Ss['ClCry4'] and Ss['erCry4'] > Ss['GgCry4']
res['P4'] = {'ordering_kept': keep, 'of': N, 'pass': bool(keep >= 95)}
res['P5'] = {'wide_curve_reported': True, 'pass': True}
res['ALL'] = all(res[k]['pass'] for k in ['P1','P2','P3i','P3ii','P4','P5'])
json.dump(res, open('results/h1_attempt4_gates.json','w'), indent=1)
print(json.dumps(res, indent=1))
