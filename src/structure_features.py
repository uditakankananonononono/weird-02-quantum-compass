#!/usr/bin/env python3
"""Structure features from ESM folds: Trp-chain geometry grounded in the erCry4 reference fold.
Chain definition is derived from the REFERENCE structure itself (no hand-typed numbering):
the erCry4 fold's Trp residues are clustered; the chain = the maximal spatially-sequential
Trp path through the domain core; its residue indices map to MSA columns for cross-species reads.
"""
import json, os, math, glob
import numpy as np

D='data/'
def parse_pdb(path):
    res = []
    seen = set()
    for line in open(path):
        if line.startswith('ATOM') and line[12:16].strip() == 'CA':
            idx = int(line[22:26])
            aa = line[17:20].strip()
            xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            b = float(line[60:66])
            res.append((idx, aa, xyz, b))
    return res

def dist(a,b): return float(np.linalg.norm(a-b))

# 1) reference erCry4 fold: find Trp chain
ref = parse_pdb(glob.glob(D+'folds/REF__*A0A2I4SZI9*.pdb')[0])
trps = [(i,aa,x,b) for i,aa,x,b in ref if aa=='TRP']
# chain = set of Trps where consecutive C-alpha distances are 8-30 A (ET-relevant hops)
best = []
for t in trps:
    chain = [t]
    cur = t
    pool = [u for u in trps if u[0] != t[0]]
    while pool:
        nxt = min(pool, key=lambda u: dist(cur[2], u[2]))
        d = dist(cur[2], nxt[2])
        if d <= 30.0:
            chain.append(nxt); pool.remove(nxt); cur = nxt
        else: break
    if len(chain) > len(best): best = chain
best.sort(key=lambda t: t[0])
chain_idx = [t[0] for t in best]
chain_gaps = [round(dist(a[2],b[2]),2) for a,b in zip(best, best[1:])]
print('erCry4 reference Trp chain (fold-derived):', chain_idx, 'gaps A:', chain_gaps)

# 2) map reference residue indices to MSA columns
msa = {}
name = None
for line in open(D+'cry1_panel_msa.fasta'):
    if line.startswith('>'): name = line[1:].strip(); msa[name]=''
    else: msa[name]+=line.strip()
refkey = [k for k in msa if 'A0A2I4SZI9' in k][0]
refaln = msa[refkey]
# residue index i (1-based) -> column
res2col = {}
ri = 0
for c, ch in enumerate(refaln):
    if ch != '-':
        ri += 1; res2col[ri] = c
cols = [res2col[i] for i in chain_idx if i in res2col]

AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
       'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
       'THR':'T','TRP':'W','TYR':'Y','VAL':'V'}

# 3) per-species structure features
out = {}
for pdb in sorted(glob.glob(D+'folds/*.pdb')):
    nm = os.path.basename(pdb)[:-4].replace('__','|')
    res = parse_pdb(pdb)
    plddt = np.mean([b for _,_,_,b in res])
    # Trp presence at chain columns: map species sequence to columns via its own alignment row
    key = nm if nm in msa else None
    if key is None:
        # REF rows and CRY4 rows named differently; try matching by accession substring
        acc = nm.split('|')[-1]
        key = next((k for k in msa if acc.split('.')[0] in k), None)
    at_chain = {}
    if key is not None:
        al = msa[key]
        for j, col in enumerate(cols):
            at_chain[f'chain{j+1}_aa'] = al[col] if col < len(al) else '-'
    # spatial: species' own Trp residues near reference chain centroid
    ref_centroid = np.mean([t[2] for t in best], axis=0)
    own_trps = [(i,x) for i,aa,x,b in res if aa=='TRP']
    near = [(i, round(dist(x, ref_centroid),2)) for i,x in own_trps if dist(x, ref_centroid) < 25]
    # consecutive Trp-Trp distances among domain Trps sorted by proximity graph
    pair_d = sorted(round(dist(a[1],b[1]),2) for i,a in enumerate(own_trps) for b in own_trps[i+1:] if dist(a[1],b[1])<30)
    out[nm] = {'plddt_mean': round(float(plddt),2),
               'n_trp_domain': len(own_trps),
               'chain_col_aa': at_chain,
               'chain_trp_present': sum(1 for v in at_chain.values() if v=='W'),
               'n_trp_near_refcentroid_25A': len(near),
               'trp_pair_dists_lt30A': pair_d[:12]}
json.dump(out, open(D+'cry1_structure_features_v1.json','w'), indent=1)
print('structure features:', len(out))
json.dump({'ref_chain_residx': chain_idx, 'ref_chain_gaps_A': chain_gaps, 'msa_cols': cols},
          open(D+'trp_chain_definition.json','w'), indent=1)
