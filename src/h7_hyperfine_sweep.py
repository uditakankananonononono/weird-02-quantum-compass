#!/usr/bin/env python3
"""W02 queue #7 hyperfine-tensor sweep (amendment locked 16:22 IST, commit feb831f).
23 tensor sets at the three species' LOCKED WT taus (results/h1_attempt4_stage1.json).
B grid + window + 9-theta family identical to attempt-4 amended 08:02 spec. Report-only."""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from bench_common import completion_yield, N5, N10, Y45
from radical_pair import RadicalPair

B_PRIMARY = np.linspace(1e-6, 200e-6, 9)
WIN = (25e-6, 65e-6)
THETAS9 = np.linspace(0, np.pi, 37)[::4][:9]
TAUS = {'erCry4': 0.1495, 'ClCry4': 4.454, 'GgCry4': 4.2711}

def iso(T): return np.eye(3) * np.trace(T) / 3.0
def aniso(T): return T - iso(T)
def build(s_iso, s_an, per=None):
    ts = {'N5': N5.copy(), 'N10': N10.copy(), 'Y45': Y45.copy()}
    out = {}
    for k, T in ts.items():
        scale = per.get(k, 1.0) if per else 1.0
        out[k] = (iso(T) * s_iso + aniso(T) * s_an) * scale
    return out

SETS = []
for si in (0.70, 0.85, 1.00, 1.15, 1.30):
    for sa in (0.70, 1.00, 1.30):
        SETS.append((f'global_iso{si}_an{sa}', build(si, sa)))
for name in ('N5', 'N10', 'Y45'):
    for sc in (0.70, 1.30):
        SETS.append((f'{name}_x{sc}', build(1.0, 1.0, {name: sc})))
SETS.append(('control_isotropic', build(1.0, 0.0)))
SETS.append(('control_iso0', build(0.0, 1.0)))
assert len(SETS) == 23

def S_at(tau_us, tens):
    k = 1.0 / tau_us
    rp = RadicalPair(hyperfine_A=[tens['N5'], tens['N10']], hyperfine_B=[tens['Y45']], kS=k, kT=k)
    A = []
    for B in B_PRIMARY:
        phis = np.array([completion_yield(rp, float(B), float(t)) for t in THETAS9])
        A.append((phis.max() - phis.min()) / phis.mean())
    A = np.array(A)
    m = (B_PRIMARY >= WIN[0]) & (B_PRIMARY <= WIN[1])
    f = np.trapezoid if hasattr(np, 'trapezoid') else np.trapz
    return float(f(A[m], B_PRIMARY[m])), float(A.max() - A.min())

rows = []
for label, tens in SETS:
    r = {'set': label, 'S': {}, 'aniso_span': {}}
    for sp, tau in TAUS.items():
        s, span = S_at(tau, tens)
        r['S'][sp] = round(s, 6); r['aniso_span'][sp] = round(span, 6)
    adv = r['S']['erCry4'] - max(r['S']['ClCry4'], r['S']['GgCry4'])
    r['xu_ordering'] = bool(adv > 0); r['advantage'] = round(adv, 6)
    rows.append(r)
    print(label, 'xu' if r['xu_ordering'] else 'no', round(adv, 5), flush=True)

frac = sum(r['xu_ordering'] for r in rows) / len(rows)
out = {'amendment_commit': 'feb831f', 'n_sets': len(SETS), 'xu_ordering_fraction': round(frac, 4),
       'max_advantage': max(r['advantage'] for r in rows), 'max_deficit': min(r['advantage'] for r in rows),
       'taus_locked': TAUS, 'rows': rows}
json.dump(out, open('results/h7_hyperfine_sweep.json', 'w'), indent=1)
print('XU fraction', round(frac, 4), '->', 'results/h7_hyperfine_sweep.json')
