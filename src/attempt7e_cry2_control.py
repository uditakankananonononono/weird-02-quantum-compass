#!/usr/bin/env python3
"""W02 attempt 7e (amendment 2026-09-27 03:08 IST, locked pre-outcome): Cry2 discrimination control.
Mirrors 7c/7d machinery exactly. Cry2 core = mapped Trp-triad + FAD-contact <=4.5A + Trp-contact
<=6.0A in the Cry2 ESMFold reference frame (FAD transferred from 6PTZ by CA superposition).
W over qualified non-core columns (>=90% occupancy), column bootstrap 10k seed 260927.
WIN: CI95_Cry1.lower > W_Cry2.point AND W_Cry1.point > CI95_Cry2.upper. Else falsification branch."""
import json, glob
import numpy as np
from collections import Counter
from pyfamsa import Aligner, Sequence
from Bio.Align import PairwiseAligner
from Bio.SVDSuperimposer import SVDSuperimposer
import random

SEED, NBOOT = 260927, 10000
AA3={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','TRP':'W','TYR':'Y','VAL':'V','PRO':'P','SER':'S','THR':'T'}

# 1. 6PTZ frame: CA coords + FAD atoms
ca6, fad = {}, []
for line in open('/tmp/6ptz.pdb'):
    if line.startswith('ATOM') and line[12:16].strip()=='CA':
        ca6[int(line[22:26])] = (line[17:20].strip(), np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])]))
    elif line.startswith('HETATM') and line[17:20].strip()=='FAD':
        fad.append(np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])]))
resnums6 = sorted(ca6)
seq6 = ''.join(AA3.get(ca6[r][0],'X') for r in resnums6)

# 2. Cry2 fold CA coords (residues 1..400 of Columba livia Cry2)
ca2 = {}
for line in open('data/folds/CRY2__REF__Columba_livia__XP_064920712.1.pdb'):
    if line.startswith('ATOM') and line[12:16].strip()=='CA':
        ca2[int(line[22:26])] = (line[17:20].strip(), np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])]))
resnums2 = sorted(ca2)
seq2 = ''.join(AA3.get(ca2[r][0],'X') for r in resnums2)

# 3. Align 6PTZ seq <-> Cry2 seq (same aligner as attempt6_geometry.py)
al = PairwiseAligner().align(seq6, seq2)[0]
p6, p2, pairs = [], [], []
for (s6,e6),(s2,e2) in zip(*al.aligned):
    for i in range(e6-s6):
        p6.append(resnums6[s6+i]); p2.append(resnums2[s2+i])
        pairs.append((resnums6[s6+i], resnums2[s2+i]))
m6to2 = dict(pairs)

# 4. Superpose Cry2 onto 6PTZ over aligned CA pairs; transform Cry2 CAs into 6PTZ frame
fixed = np.array([ca6[a][1] for a,_ in pairs])
moving = np.array([ca2[b][1] for _,b in pairs])
sup = SVDSuperimposer(); sup.set(fixed, moving); sup.run()
rot, tran = sup.get_rotran()
xyz2 = {r: np.dot(x, rot) + tran for r,(aa,x) in ca2.items()}  # Cry2 CAs in 6PTZ frame
print('aligned CA pairs:', len(pairs), 'RMSD: %.2f A' % sup.get_rms())

# 5. Triad mapping (locked chain def, 6PTZ numbering) -> Cry2 resnums; sanity: must be Trp
chain6 = json.load(open('data/trp_chain_definition.json'))['chain_6ptz_numbering']
triad2 = {r6: m6to2.get(r6) for r6 in chain6}
print('triad map 6PTZ->Cry2:', triad2, 'Cry2 AAs:', {r6: ca2[r2][0] for r6,r2 in triad2.items() if r2})
triad_xyz = [xyz2[r2] for r2 in triad2.values() if r2]

# 6. Core residues in Cry2 numbering (same thresholds as 7c)
core2 = set(r2 for r2 in triad2.values() if r2)
for r,(aa,x) in ca2.items():
    xt = xyz2[r]
    fd = min(np.linalg.norm(xt-f) for f in fad)
    td = min(np.linalg.norm(xt-t) for t in triad_xyz)
    if fd <= 4.5 or td <= 6.0: core2.add(r)
print('core residues (Cry2 numbering):', len(core2))

# 7. pyfamsa align Cry2 panel + REF (full-length Columba livia Cry2)
files = sorted(glob.glob('data/cry2_panel/*.fa'))
refseq = ''.join(open('data/cry2_panel/Columba_livia.fa').read().split('\n')[1:]).strip()
seqs = [Sequence(fn.split('/')[-1][:-3].encode(), ''.join(open(fn).read().split('\n')[1:]).strip().encode()) for fn in files]
seqs.append(Sequence(b'REF', refseq.encode()))
msa = [(s.id.decode(), s.sequence.decode()) for s in Aligner().align(seqs)]
ref_row = [s for s in msa if s[0]=='REF'][0][1]
rows = [k for k,_ in msa if k != 'REF']
print('MSA rows:', len(rows), 'length:', len(ref_row))

# 8. Map core residues -> MSA columns; geometry only exists for residues <= 400
pos2col, pos = {}, 0
for i, a in enumerate(ref_row):
    if a != '-':
        pos += 1
        pos2col[pos] = i
geocols = set(pos2col[r] for r in pos2col if r <= 400)
core_cols = set(pos2col[r] for r in core2 if r in pos2col)
qcols = [c for c in sorted(geocols) if c not in core_cols]
print('geometry cols:', len(geocols), 'core cols:', len(core_cols), 'qualified non-core:', len(qcols))

# 9. W statistic (identical rule to 7c)
rowseq = dict(msa)
vals = []
for c in qcols:
    res = [rowseq[r][c] for r in rows if rowseq[r][c] in 'ACDEFGHIKLMNPQRSTVWY']
    if len(res) >= 0.9*len(rows):
        vals.append(Counter(res).most_common(1)[0][1]/len(res))
W2 = sum(vals)/len(vals)
rng = random.Random(SEED)
boots = sorted(sum(rng.choice(vals) for _ in range(len(vals)))/len(vals) for _ in range(NBOOT))
ci2 = [boots[int(0.025*NBOOT)], boots[int(0.975*NBOOT)]]
print('CRY2 W=%.5f n=%d CI95 [%.5f, %.5f]' % (W2, len(vals), ci2[0], ci2[1]))

# 10. Locked win/falsification evaluation vs Cry1 (7d committed results)
r7d = json.load(open('results/attempt7d_robustness.json'))
ci1 = r7d['W_bootstrap_CI95']['CRY1']
W1 = r7d['W']['CRY1']
cond_a = ci1[0] > W2
cond_b = W1 > ci2[1]
win = cond_a and cond_b
print('Cry1 W=%.5f CI95 %s | cond_a (CI1.lower > W2): %s | cond_b (W1 > CI2.upper): %s | WIN: %s' % (W1, ci1, cond_a, cond_b, win))

out = {'amendment':'2026-09-27 03:08 IST 7e','seed':SEED,'nboot':NBOOT,
 'panel_species_ok':len(rows),'cry2_ref':'XP_064920712.1 Columba_livia (ESMFold, first 400aa, same crop convention)',
 'aligned_pairs':len(pairs),'superposition_rmsd':round(float(sup.get_rms()),3),
 'triad_map_6ptz_to_cry2':{str(k):v for k,v in triad2.items()},
 'n_core_residues':len(core2),'n_geometry_cols':len(geocols),'n_qualified_cols':len(vals),
 'W_cry2':round(W2,5),'ci95_cry2':[round(x,5) for x in ci2],
 'W_cry1':W1,'ci95_cry1':ci1,
 'cond_a_CI1lower_gt_W2':bool(cond_a),'cond_b_W1_gt_CI2upper':bool(cond_b),'win':bool(win),
 'falsification_branch':'if not win: canalization = generic photoreceptor purifying-selection pattern; downgrade 7c novelty claim; judge redirect for next discriminator'}
json.dump(out, open('results/attempt7e_cry2_control.json','w'), indent=1)
print('SAVED results/attempt7e_cry2_control.json')
