#!/usr/bin/env python3
"""W02 attempt 8b corrected (amendments 11:13 + 11:14 IST, locked pre-computation)."""
import json, random, sys
sys.path.insert(0, 'src')
import attempt8_redirect3 as r3

SEED, NPERM = 260927, 10000
table = r3.table
def pcons(t):
    s1, s2 = t['s1'], t['s2']
    if s1['dom'] == s2['dom'] and s1['dom_frac'] >= 0.95 and s2['dom_frac'] >= 0.95:
        return s1['dom']
    return None
def cry4state(t):
    s4 = t['s4']
    return (s4['dom_frac'] >= 10/11, s4['dom'])
keys = list(table)
def candidates(tab, states):
    return [r for r in keys if pcons(tab[r]) is not None and states[r][0] and states[r][1] != pcons(tab[r])]
states = {r: cry4state(table[r]) for r in keys}
cand = candidates(table, states)
rng = random.Random(SEED); cnt = 0
vals = list(states.values())
for _ in range(NPERM):
    perm = vals[:]; rng.shuffle(perm)
    st = dict(zip(keys, perm))
    if len(candidates(table, st)) >= len(cand): cnt += 1
p = (1 + cnt) / (NPERM + 1)
import numpy as np
def loc(r, t):
    if t['near']: return 'near-site'
    x = r3.g.ca6[r][1]
    return 'surface' if np.linalg.norm(x - r3.g.fad6.mean(axis=0)) > 12.0 else 'core'
detail = {str(r): {'paralog_consensus': pcons(table[r]), 'cry4': states[r][1],
                   'cry4_dom_frac': round(table[r]['s4']['dom_frac'], 3), 'loc': loc(r, table[r])} for r in cand}
out = {'amendments': ['2026-09-27 11:13 IST attempt8b', '2026-09-27 11:14 IST clarification'],
       'seed': SEED, 'nperm': NPERM, 'null': 'joint shuffle of (cry4 tight, symbol) over mapped residues',
       'n_mapped': len(keys), 'n_candidates': len(cand), 'enrichment_p_one_sided': round(p, 5),
       'candidates': sorted(cand), 'detail': detail,
       'framing': 'Cry4-distinguishing residues - family distinction, NOT magnetoreception causes'}
json.dump(out, open('results/attempt8b_specialization.json', 'w'), indent=1)
print('mapped:', len(keys), 'candidates:', len(cand), 'p:', round(p, 5))
for r in sorted(cand):
    d = detail[str(r)]
    print(r, d['paralog_consensus'], '->', d['cry4'], d['cry4_dom_frac'], d['loc'])
