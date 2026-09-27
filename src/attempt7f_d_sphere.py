#!/usr/bin/env python3
"""W02 7f-d (amendment 2026-09-27 09:03 IST, locked pre-evaluation): radical-pair-sphere localized canalization.
S_RP = locked Trp chain [395,372,318,369] (6PTZ numbering) UNION all residues with any heavy atom
within 5.0A of any FAD heavy atom or any chain-Trp heavy atom. Cry1: 6PTZ directly.
Cry2: AF fold superposed on 6PTZ CA (7e machinery), sphere computed in 6PTZ frame against 6PTZ FAD
+ mapped chain-Trp atoms. W_RP = mean max-freq over qualified S_RP columns (>=90% occupancy),
column bootstrap CI95 (10k, seed 260927, 7d arm-A scheme).
WIN: W_RP,Cry1 point > W_RP,Cry2 CI95 upper. Falsification: claim dies."""
import json, random, glob
import numpy as np
from collections import Counter
from Bio.Align import PairwiseAligner
from Bio.SVDSuperimposer import SVDSuperimposer
from pyfamsa import Aligner, Sequence

SEED, NBOOT, RCUT = 260927, 10000, 5.0
AA3={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','TRP':'W','TYR':'Y','VAL':'V','PRO':'P','SER':'S','THR':'T'}
CHAIN6 = json.load(open('data/trp_chain_definition.json'))['chain_6ptz_numbering']

def parse_pdb_all(path):
    res = {}   # resnum -> (aa3, array Nx3 heavy atoms)
    ca = {}
    fad = []
    for line in open(path):
        rec = line[:6].strip()
        if rec == 'ATOM':
            an = line[12:16].strip()
            if an == 'H' or an.startswith('H'): continue
            el = line[76:78].strip()
            if el == 'H': continue
            rn = int(line[22:26]); aa = line[17:20].strip()
            xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            res.setdefault(rn, (aa, []))[1].append(xyz)
            if an == 'CA': ca[rn] = xyz
        elif rec == 'HETATM' and line[17:20].strip() == 'FAD':
            fad.append(np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]))
    return res, ca, np.array(fad)

res6, ca6, fad6 = parse_pdb_all('/tmp/6ptz.pdb')
res2, ca2, _ = parse_pdb_all('data/folds/CRY2__REF__Columba_livia__XP_064920712.1.pdb')

def sphere(res, ref_atoms):
    out = set()
    for rn, (aa, atoms) in res.items():
        A = np.array(atoms)
        d = np.sqrt(((A[:, None, :] - ref_atoms[None, :, :])**2).sum(-1)).min()
        if d <= RCUT: out.add(rn)
    return out

chain_atoms6 = np.array([x for rn in CHAIN6 for x in res6[rn][1]])
ref6 = np.vstack([fad6, chain_atoms6])
SRP1 = set(CHAIN6) | sphere(res6, ref6)
print('CRY1 S_RP residues: %d' % len(SRP1), flush=True)

# superpose cry2 on 6PTZ via CA alignment (7e machinery)
resnums6 = sorted(ca6); seq6 = ''.join(AA3.get(res6[r][0],'X') for r in resnums6)
resnums2 = sorted(ca2); seq2 = ''.join(AA3.get(res2[r][0],'X') for r in resnums2)
al = PairwiseAligner().align(seq6, seq2)[0]
m6to2 = {}
for (s6,e6),(s2,e2) in zip(*al.aligned):
    for i in range(e6-s6): m6to2[resnums6[s6+i]] = resnums2[s2+i]
pairs = [(a,b) for a,b in m6to2.items()]
fixed = np.array([ca6[a] for a,_ in pairs]); moving = np.array([ca2[b] for _,b in pairs])
sup = SVDSuperimposer(); sup.set(fixed, moving); sup.run(); rot, tran = sup.get_rotran()
def X(x): return np.dot(x, rot) + tran
res2X = {rn: (aa, [X(x) for x in atoms]) for rn,(aa,atoms) in res2.items()}
chain2 = [m6to2[r6] for r6 in CHAIN6 if r6 in m6to2]
chain_atoms2X = np.array([x for rn in chain2 for x in res2X[rn][1]])
ref2 = np.vstack([fad6, chain_atoms2X])
SRP2 = set(chain2) | sphere(res2X, ref2)
print('CRY2 S_RP residues: %d (chain mapped: %s)' % (len(SRP2), chain2), flush=True)

# CRY1 restricted columns via geometry map
geo = json.load(open('data/attempt6_column_geometry.json'))
rp_cols1 = sorted(int(c) for c,g in geo.items() if g['resnum_6ptz'] in SRP1)
msa1 = {}
for line in open('data/cry1_panel_msa.fasta'):
    if line.startswith('>'): k=line[1:].strip(); msa1[k]=''
    else: msa1[k]+=line.strip()
rows1 = [k for k in msa1 if not k.startswith('REF|')]

# CRY2 MSA + restricted columns via REF row
files = sorted(glob.glob('data/cry2_panel/*.fa'))
refseq = ''.join(open('data/cry2_panel/Columba_livia.fa').read().split('\n')[1:]).strip()
seqs = [Sequence(fn.split('/')[-1][:-3].encode(), ''.join(open(fn).read().split('\n')[1:]).strip().encode()) for fn in files]
seqs.append(Sequence(b'REF', refseq.encode()))
msa2l = [(s.id.decode(), s.sequence.decode()) for s in Aligner().align(seqs)]
ref_row = [s for s in msa2l if s[0]=='REF'][0][1]
pos2col, pos = {}, 0
for i, a in enumerate(ref_row):
    if a != '-': pos += 1; pos2col[pos] = i
rp_cols2 = sorted(pos2col[r] for r in SRP2 if r in pos2col)
msa2 = {k: v for k, v in msa2l if k != 'REF'}
rows2 = list(msa2)

def restricted_W(msa, rows, cols):
    vals = []
    for c in cols:
        res = [msa[r][c] for r in rows if msa[r][c] in 'ACDEFGHIKLMNPQRSTVWY']
        if len(res) >= 0.9*len(rows):
            vals.append(Counter(res).most_common(1)[0][1]/len(res))
    return vals

v1 = restricted_W(msa1, rows1, rp_cols1)
v2 = restricted_W(msa2, rows2, rp_cols2)
W1 = sum(v1)/len(v1); W2 = sum(v2)/len(v2)
print('CRY1 W_RP=%.5f (n=%d cols)' % (W1, len(v1)), flush=True)
print('CRY2 W_RP=%.5f (n=%d cols)' % (W2, len(v2)), flush=True)

rng = random.Random(SEED)
def ci95(vals):
    n = len(vals); boots = sorted(sum(rng.choice(vals) for _ in range(n))/n for _ in range(NBOOT))
    return [boots[int(0.025*NBOOT)], boots[int(0.975*NBOOT)]]
ci1, ci2 = ci95(v1), ci95(v2)
print('CRY1 CI95 [%.5f, %.5f]' % tuple(ci1), flush=True)
print('CRY2 CI95 [%.5f, %.5f]' % tuple(ci2), flush=True)
win = W1 > ci2[1]
out = {'amendment':'2026-09-27 09:03 IST 7f-d','seed':SEED,'nboot':NBOOT,'rcut_A':RCUT,
 'S_RP':{'CRY1_n':len(SRP1),'CRY2_n':len(SRP2),'CRY1_residues_6ptz':sorted(SRP1),'CRY2_residues':sorted(SRP2),'chain2_mapped':chain2},
 'CRY1':{'W_RP':round(W1,5),'n_cols':len(v1),'CI95':[round(x,5) for x in ci1]},
 'CRY2':{'W_RP':round(W2,5),'n_cols':len(v2),'CI95':[round(x,5) for x in ci2]},
 'win':bool(win),'falsification':'W_RP,Cry1 <= W_RP,Cry2 CI95 upper -> localized-canalization claim dies'}
json.dump(out, open('results/attempt7f_d_sphere.json','w'), indent=1)
print('WIN:', win)
