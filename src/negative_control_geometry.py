#!/usr/bin/env python3
"""W02 EXPLORATORY negative control (judge round-6 weakness 1): does ANY
alternative geometric coordinate recover the Xu ordering under the same locked
lifetime model? Labeled exploratory (post-hoc); locked 5a/5b results unchanged.
Coordinates per species (same four chain Trps, same locked tau_eff + S(tau)):
  d_term       = locked terminal-Trp-ring->FAD min distance (control sanity: reproduces locked values)
  d_term_cent  = terminal-Trp ring-centroid -> FAD isoalloxazine-centroid distance (6PTZ frame)
  d_triad      = 3rd-position Trp ring -> FAD min distance (= locked triad coordinate, sanity)
  hop34        = terminal hop length (Trp3-centroid -> Trp4-centroid, species frame)
NOTE: erCry4 PDB (MN709784 model) contains only chain positions 3-4 (residues 347/398;
positions 1-2 absent), so full-tetrad coordinates (total path length, mean pairwise,
first-Trp->FAD) are NOT computable for erCry4 - disclosed data limitation.
Xu ordering test: S(tau(coord)) for erCry4 > both ClCry4 and GgCry4.
"""
import sys, json, numpy as np
sys.path.insert(0, 'src')
import compass_scan as cs
from compass_scan import tau_eff
from Bio.Align import PairwiseAligner

d = json.load(open('results/h1_attempt4_stage1.json'))
curve = d['tau_grid_S']
taus = sorted(float(t) for t in curve)
S = np.array([curve[f'{t:g}']['S'] if f'{t:g}' in curve else curve[str(t)]['S'] for t in taus])
def S_at(tau): return float(np.interp(tau, taus, S))

SPECIES = {
    'erCry4': ('CRY4__REFSEQ__Erithacus_rubecula__MN709784', 'MN709784'),
    'ClCry4': ('CRY4__REFSEQ__Columba_livia__KX168611', 'KX168611'),
    'GgCry4': ('CRY4__REFSEQ__Gallus_gallus__NM_001039596', 'NM_001039596'),
}
TRP_RING = cs.TRP_RING

def load_pdb(pdb):
    cas, rings = {}, {}
    for line in open(pdb):
        if line.startswith('ATOM'):
            r = int(line[22:26]); atom = line[12:16].strip()
            xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            if atom == 'CA': cas[r] = xyz
            if line[17:20].strip() == 'TRP' and atom in TRP_RING: rings.setdefault(r, []).append(xyz)
    return cas, {k: np.array(v) for k, v in rings.items()}

def seq_of(acc):
    for rec in open('data/cry4_subpanel.fasta'):
        pass
    from Bio import SeqIO
    for r in SeqIO.parse('data/cry4_subpanel.fasta', 'fasta'):
        if acc in r.id: return str(r.seq)
    raise KeyError(acc)

res = {}
_, msa_all = cs.msa_row('REF|tr|A0A386QUR4|A0A386QUR4_COLLI')
clseq = msa_all['REF|tr|A0A386QUR4|A0A386QUR4_COLLI'].replace('-', '')
defn = cs.load_def()
chain_cl = defn['ref_chain_residx']   # [395,372,318,369] FAD->surface

for sp, (fold, acc) in SPECIES.items():
    pdb = f'data/folds/{fold}.pdb'
    casS, ringS = load_pdb(pdb)
    aln = PairwiseAligner(); aln.mode='global'; aln.match_score=2; aln.mismatch_score=-1
    aln.open_gap_score=-5; aln.extend_gap_score=-0.5
    a = aln.align(seq_of(acc), clseq)[0]
    pairs = []; sp2cl = {}
    for (s0,s1),(c0,c1) in zip(*a.aligned):
        for i in range(s1-s0):
            pairs.append((s0+i+1, c0+i+1)); sp2cl[s0+i+1] = c0+i+1
    inv = {v:k for k,v in sp2cl.items()}
    chain_sp = [inv.get(c) for c in chain_cl]
    cents = {r: ringS[r].mean(0) for r in chain_sp if r is not None and r in ringS}
    hop34 = float(np.linalg.norm(cents[chain_sp[2]] - cents[chain_sp[3]]))
    d_term = cs.d_term_edge(pdb, chain_sp[3], pairs)       # locked coordinate (sanity)
    d_triad = cs.d_term_edge(pdb, chain_sp[2], pairs)      # locked triad coordinate (sanity)
    d_term_cent = cs.d_term_centroid(pdb, chain_sp[3], pairs)
    res[sp] = {'chain_sp': [int(x) for x in chain_sp], 'hop34': round(hop34,2),
               'd_term': round(d_term,2), 'd_triad': round(d_triad,2),
               'd_term_cent': round(d_term_cent,2),
               'd_term_locked': d['species'][sp]['d_term_WT'],
               'd_triad_locked': d['species'][sp]['d_term_triad']}

out = {'label': 'EXPLORATORY negative control (judge round-6); locked results unchanged',
       'species': res, 'ordering_tests': {}}
for coord in ('d_term', 'd_term_cent', 'd_triad', 'hop34'):
    row = {}
    for sp in SPECIES:
        tau = tau_eff(res[sp][coord])
        row[sp] = {'coord': res[sp][coord], 'tau': round(tau,4), 'S': S_at(tau)}
    ok = row['erCry4']['S'] > row['ClCry4']['S'] and row['erCry4']['S'] > row['GgCry4']['S']
    out['ordering_tests'][coord] = {'per_species': row, 'xu_ordering_recovered': bool(ok)}
json.dump(out, open('results/h1_negative_control_geometry.json','w'), indent=1)
for coord, t in out['ordering_tests'].items():
    print(coord, '-> Xu ordering recovered:', t['xu_ordering_recovered'],
          '| S:', {sp: f"{t['per_species'][sp]['S']:.2e}" for sp in SPECIES})
print('sanity d_term vs locked:', {sp: (res[sp]['d_term'], res[sp]['d_term_locked']) for sp in SPECIES})
