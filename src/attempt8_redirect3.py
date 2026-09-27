#!/usr/bin/env python3
"""W02 attempt 8 / redirect 3 (amendment 2026-09-27 11:12 IST, locked pre-computation).
Cry4-specific constraint residues: cry4-conserved, paralog-variable (Cry1 AND Cry2), near-site (8.0A)."""
import json, random, sys
import numpy as np
sys.path.insert(0, 'src')
import attempt7g_ab as g  # locked geometry + panel machinery (deterministic re-run)

SEED, NPERM = 260927, 10000
AA = 'ACDEFGHIKLMNPQRSTVWY'

# cry4 alignment
seqs = {}
for line in open('data/mafft_ebi_crosscheck.fasta'):
    line = line.strip()
    if line.startswith('>'): k = line[1:]; seqs[k] = ''
    else: seqs[k] += line
ref4 = seqs['KX168611']

# pairwise cry4-ref -> 6PTZ sequence
from Bio.Align import PairwiseAligner
resnums6 = sorted(g.ca6); seq6 = ''.join(g.AA3.get(g.ca6[r][0], 'X') for r in resnums6)
al = PairwiseAligner().align(ref4.replace('-', ''), seq6)[0]
m4 = {}  # 6PTZ residue -> cry4 column
pos4 = 0
ungapped4 = ref4.replace('-', '')
# build cry4 column index per ungapped position
col_of_pos = []
for i, a in enumerate(ref4):
    if a != '-': col_of_pos.append(i)
for (s4, e4), (s6, e6) in zip(*al.aligned):
    for i in range(e4 - s4):
        m4[resnums6[s6 + i]] = col_of_pos[s4 + i]

# cry1 column per 6PTZ residue
res2col1 = {}
for c, r in g.col2res1.items(): res2col1.setdefault(r, c)
# cry2 column per 6PTZ residue (via 7e m6to2)
aa2 = {}
for line in open('data/folds/CRY2__REF__Columba_livia__XP_064920712.1.pdb'):
    if line.startswith('ATOM') and line[12:16].strip() == 'CA':
        aa2[int(line[22:26])] = line[17:20].strip()
m6to2 = {}
al2 = PairwiseAligner().align(seq6, ''.join(g.AA3.get(aa2[r], 'X') for r in sorted(aa2)))[0]
resnums2 = sorted(aa2)
for (s6, e6), (s2, e2) in zip(*al2.aligned):
    for i in range(e6 - s6): m6to2[resnums6[s6 + i]] = resnums2[s2 + i]
res2col2 = {}
for c, r in g.col2res2.items(): res2col2.setdefault(r, c)

rows1 = [v for k, v in g.msa1.items()]
rows2 = [v for k, v in g.msa2.items()]
rows4 = list(seqs.values())

def stats(rows, col):
    occ = [s[col] for s in rows if s[col] in AA]
    if len(occ) < 0.9 * len(rows) or not occ: return None
    from collections import Counter
    cnt = Counter(occ); dom, dn = cnt.most_common(1)[0]
    variant = len([1 for v in cnt.values() if v >= 0.05 * len(rows)]) >= 2
    return {'occ': len(occ) / len(rows), 'dom': dom, 'dom_frac': dn / len(occ), 'variant': variant}

fad_centroid = g.fad6.mean(axis=0)
chain = json.load(open('data/trp_chain_definition.json'))['chain_6ptz_numbering']
def near_site(r):
    x = g.ca6[r][1]
    if np.linalg.norm(x - fad_centroid) <= 8.0: return True
    return any(np.linalg.norm(x - g.ca6[t][1]) <= 8.0 for t in chain if t in g.ca6)

# per-residue table over mappable residues
table = {}
for r in resnums6:
    c4 = m4.get(r); c1 = res2col1.get(r); r2 = m6to2.get(r); c2 = res2col2.get(r2) if r2 else None
    if c4 is None or c1 is None or c2 is None: continue
    s4, s1, s2 = stats(rows4, c4), stats(rows1, c1), stats(rows2, c2)
    if not s4 or not s1 or not s2: continue
    table[r] = {'c4': c4, 'c1': c1, 'c2': c2, 's4': s4, 's1': s1, 's2': s2, 'near': near_site(r)}

def conserved4(s): return s['dom_frac'] >= 10 / 11
cand_both = [r for r, t in table.items() if conserved4(t['s4']) and t['s1']['variant'] and t['s2']['variant'] and t['near']]
cand_either = [r for r, t in table.items() if conserved4(t['s4']) and (t['s1']['variant'] or t['s2']['variant']) and t['near']]

# permutation null: preserve count of cry4-conserved residues among mapped table
keys = list(table)
n_cons = sum(1 for r in keys if conserved4(table[r]['s4']))
obs = len(cand_both)
rng = random.Random(SEED); cnt = 0
for _ in range(NPERM):
    draw = rng.sample(keys, n_cons)
    h = sum(1 for r in draw if table[r]['s1']['variant'] and table[r]['s2']['variant'] and table[r]['near'])
    if h >= obs: cnt += 1
p = (1 + cnt) / (NPERM + 1)

diag = {
 'near_site': sum(1 for t in table.values() if t['near']),
 'var_cry1': sum(1 for t in table.values() if t['s1']['variant']),
 'var_cry2': sum(1 for t in table.values() if t['s2']['variant']),
 'var_both': sum(1 for t in table.values() if t['s1']['variant'] and t['s2']['variant']),
 'near_and_var1': sum(1 for t in table.values() if t['near'] and t['s1']['variant']),
 'near_and_var2': sum(1 for t in table.values() if t['near'] and t['s2']['variant']),
 'near_and_var_either': sum(1 for t in table.values() if t['near'] and (t['s1']['variant'] or t['s2']['variant'])),
 'cry4_nonconserved': sum(1 for t in table.values() if not conserved4(t['s4'])),
}
print('DIAG:', diag)
out = {'amendment': '2026-09-27 11:12 IST attempt8 redirect3', 'seed': SEED, 'nperm': NPERM,
 'diagnostics': diag,
 'n_mapped_residues': len(table), 'n_cry4_conserved': n_cons,
 'candidates_both': cand_both, 'candidates_either_sensitivity': cand_either,
 'enrichment': {'obs': obs, 'p_one_sided': round(p, 5)},
 'detail': {str(r): {'cry4_cons': table[r]['s4']['dom'], 'cry1_cons': table[r]['s1']['dom'], 'cry2_cons': table[r]['s2']['dom'],
                     'cry1_var': table[r]['s1']['variant'], 'cry2_var': table[r]['s2']['variant']} for r in cand_both},
 'framing': 'magnetoreception-COMPATIBLE constraint residues - no causality claim (verdict anti-goal)'}
json.dump(out, open('results/attempt8_redirect3.json', 'w'), indent=1)
print('mapped residues:', len(table), 'cry4-conserved:', n_cons)
print('candidates (both):', cand_both)
print('candidates (either, sensitivity):', cand_either)
print('enrichment p:', p)
