#!/usr/bin/env python3
"""W02 queue #15 MDE (amendment 2026-09-27 11:17 IST, locked pre-computation).
Replicates the locked 7f-b / 7g-b permutation streams (same seed, same statistic code via
exec-prefix / module import) to record null SDs; AUC MDE via Hanley-McNeil bisection."""
import json, random, sys
import numpy as np
sys.path.insert(0, 'src')
from scipy.stats import norm

SEED, NPERM = 260927, 10000
Z1 = norm.ppf(0.95) + norm.ppf(0.8)    # one-sided
Z2 = norm.ppf(0.975) + norm.ppf(0.8)   # two-sided
out = {'amendment': '2026-09-27 11:17 IST #15 MDE', 'power': 0.8, 'alpha_note': 'one-sided 0.05 (7f-b, AUC); two-sided 0.05 (7g-b)'}
sanity = {}

# ---- 7f-b helpers without running the module main body ----
src = open('src/attempt7f_b_gxp.py').read()
prefix = src.split('\nlab = labels(None)')[0]
ns = {}
exec(prefix, ns)
load_msa, labels, cry1_qcols, cry2_qcols, W_group = ns['load_msa'], ns['labels'], ns['cry1_qcols'], ns['cry2_qcols'], ns['W_group']
lab = labels(None)
msa1 = load_msa('data/cry1_panel_msa.fasta')
msa1 = {k: v for k, v in msa1.items() if not k.startswith('REF|')}
qc1 = cry1_qcols()
msa2, qc2 = cry2_qcols()
committed = json.load(open('results/attempt7f_b_gxp.json'))

def arm(msa, qc, lab_use, tag, committed_block):
    rows = [r for r in msa if (r.split('|')[1] if '|' in r else r) in lab_use]
    mig = [r for r in rows if lab_use[r.split('|')[1] if '|' in r else r] == 'migratory']
    sed = [r for r in rows if lab_use[r.split('|')[1] if '|' in r else r] == 'sedentary']
    Wm, Ws = W_group(msa, mig, qc), W_group(msa, sed, qc)
    dW = Wm - Ws
    if committed_block is not None:
        assert abs(Wm - committed_block['W_mig']) < 5e-6 and abs(Ws - committed_block['W_sed']) < 5e-6, (tag, Wm, Ws)
        sanity[tag] = 'observed W_mig/W_sed reproduce committed JSON'
    rng = random.Random(SEED)
    allrows = mig + sed; nm = len(mig); nulls = []
    for _ in range(NPERM):
        samp = set(rng.sample(allrows, nm))
        a = [r for r in allrows if r in samp]; b = [r for r in allrows if r not in samp]
        nulls.append(W_group(msa, a, qc) - W_group(msa, b, qc))
    sd = float(np.std(nulls))
    return {'n_mig': nm, 'n_sed': len(sed), 'DeltaW_obs': round(dW, 5), 'SD_null': round(sd, 6),
            'MDE_DeltaW': round(Z1 * sd, 6)}

sp1 = {k.split('|')[1] for k in msa1}
sp2 = set(msa2.keys())
inter = sp1 & sp2
lab_i = {sp: lab[sp] for sp in inter if sp in lab}
msa1i = {k: v for k, v in msa1.items() if k.split('|')[1] in lab_i}
msa2i = {k: v for k, v in msa2.items() if k in lab_i}
gxp = {}
gxp['CRY1'] = arm(msa1, qc1, lab, 'CRY1', committed['CRY1'])
print('CRY1 done', gxp['CRY1'], flush=True)
gxp['CRY2'] = arm(msa2, qc2, lab, 'CRY2', committed['CRY2'])
print('CRY2 done', gxp['CRY2'], flush=True)
gxp['CRY1_intersection'] = arm(msa1i, qc1, lab_i, 'CRY1_intersection', committed['intersection']['CRY1'])
print('CRY1_int done', gxp['CRY1_intersection'], flush=True)
gxp['CRY2_intersection'] = arm(msa2i, qc2, lab_i, 'CRY2_intersection', committed['intersection']['CRY2'])
print('CRY2_int done', gxp['CRY2_intersection'], flush=True)
out['gxp_7fb'] = gxp

# ---- 7g-b tail charge (module import re-runs its own locked tests, ~2s) ----
import attempt7g_ab as g
gb = json.load(open('results/attempt7g_ab.json'))['g7b']['results']
def tail_null_sd(feats, fi, tag):
    rows = [(sp, v[fi], v[2]) for sp, v in feats.items()]
    mig = [v for _, v, l in rows if l == 'migratory']; sed = [v for _, v, l in rows if l == 'sedentary']
    obs = (sum(mig)/len(mig)) - (sum(sed)/len(sed))
    assert abs(obs - gb[tag]['obs_delta']) < 1e-4, (tag, obs)
    sanity['7gb_' + tag] = 'observed delta reproduces committed JSON'
    allv = mig + sed; nm = len(mig); rng = random.Random(SEED); nulls = []
    for _ in range(NPERM):
        samp = set(rng.sample(range(len(allv)), nm))
        a = [allv[i] for i in range(len(allv)) if i in samp]; b = [allv[i] for i in range(len(allv)) if i not in samp]
        nulls.append((sum(a)/len(a)) - (sum(b)/len(b)))
    sd = float(np.std(nulls))
    return {'n_mig': nm, 'n_sed': len(sed), 'obs_delta': gb[tag]['obs_delta'], 'SD_null': round(sd, 6),
            'MDE_delta': round(Z2 * sd, 6)}
tail = {}
for fi, fname in [(0, 'net_charge'), (1, 'charged_density')]:
    tail['cry1_' + fname] = tail_null_sd(g.f1, fi, 'cry1_' + fname)
    tail['cry2_' + fname] = tail_null_sd(g.f2, fi, 'cry2_' + fname)
out['tail_7gb'] = tail
print('7g-b done', tail, flush=True)

# ---- classifier AUC (attempt 3a, 36/80, one-sided vs 0.5) ----
n1, n0 = 36, 80
def se_hm(A):
    Q1 = A / (2 - A); Q2 = 2 * A * A / (1 + A)
    return math.sqrt((A * (1 - A) + (n1 - 1) * (Q1 - A * A) + (n0 - 1) * (Q2 - A * A)) / (n1 * n0))
import math
lo, hi = 0.5, 0.999999
for _ in range(100):
    mid = (lo + hi) / 2
    if (mid - 0.5) / se_hm(mid) < Z1: lo = mid
    else: hi = mid
out['classifier_auc'] = {'test': 'attempt 3a LOO AUC vs 0.5, one-sided', 'n_mig': n1, 'n_sed': n0,
                         'observed_best_auc': 0.402, 'MDE_AUC': round((lo + hi) / 2, 4)}
out['sanity'] = sanity
out['framing'] = 'At the locked panel sizes the nulls exclude migration-specific effects >= MDE in each native unit; smaller effects cannot be excluded. Retrospective MDE, not observed power.'
json.dump(out, open('results/attempt8c_power.json', 'w'), indent=1)
print('WRITTEN results/attempt8c_power.json', flush=True)
