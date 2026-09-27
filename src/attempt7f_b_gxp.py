#!/usr/bin/env python3
"""W02 7f-b (amendment 2026-09-27 09:03 IST, locked pre-evaluation): gene x phenotype interaction.
DeltaW = W_migratory - W_sedentary per gene on its OWN panel MSA + qualified columns (7c rule,
>=90% occupancy within each subgroup). 10k label permutations per gene, seed 260927.
WIN: DeltaW_Cry1 > 0 at p<0.05 AND DeltaW_Cry2 n.s. Sensitivity: intersection species set."""
import json, random
from collections import Counter

SEED, NPERM = 260927, 10000

def load_msa(path):
    msa = {}
    for line in open(path):
        if line.startswith('>'): k=line[1:].strip(); msa[k]=''
        else: msa[k]+=line.strip()
    return msa

def labels(species):
    # labels from the Cry1 panel headers (locked panel list)
    lab = {}
    for k in load_msa('data/cry1_panel_msa.fasta'):
        parts = k.split('|')
        if len(parts) >= 3 and parts[0] in ('migratory','sedentary'):
            lab[parts[1]] = parts[0]
    return lab

def cry1_qcols():
    geo = json.load(open('data/attempt6_column_geometry.json'))
    chain6 = set(json.load(open('data/trp_chain_definition.json'))['chain_6ptz_numbering'])
    core = {int(c) for c,g in geo.items() if g['resnum_6ptz'] in chain6 or g['fad_d']<=4.5 or g['trp_d']<=6.0}
    return sorted(set(int(c) for c in geo) - core)

def cry2_qcols():
    # rebuilt from the 7e core residue set via the REF row mapping (identical to 7e)
    import glob
    import numpy as np
    from Bio.Align import PairwiseAligner
    from Bio.SVDSuperimposer import SVDSuperimposer
    AA3={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','TRP':'W','TYR':'Y','VAL':'V','PRO':'P','SER':'S','THR':'T'}
    ca6, fad = {}, []
    for line in open('/tmp/6ptz.pdb'):
        if line.startswith('ATOM') and line[12:16].strip()=='CA':
            ca6[int(line[22:26])] = (line[17:20].strip(), np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])]))
        elif line.startswith('HETATM') and line[17:20].strip()=='FAD':
            fad.append(np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])]))
    resnums6 = sorted(ca6); seq6 = ''.join(AA3.get(ca6[r][0],'X') for r in resnums6)
    ca2 = {}
    for line in open('data/folds/CRY2__REF__Columba_livia__XP_064920712.1.pdb'):
        if line.startswith('ATOM') and line[12:16].strip()=='CA':
            ca2[int(line[22:26])] = (line[17:20].strip(), np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])]))
    resnums2 = sorted(ca2); seq2 = ''.join(AA3.get(ca2[r][0],'X') for r in resnums2)
    al = PairwiseAligner().align(seq6, seq2)[0]
    m6to2 = {}
    for (s6,e6),(s2,e2) in zip(*al.aligned):
        for i in range(e6-s6): m6to2[resnums6[s6+i]] = resnums2[s2+i]
    pairs = [(a,b) for a,b in m6to2.items()]
    fixed = np.array([ca6[a][1] for a,_ in pairs]); moving = np.array([ca2[b][1] for _,b in pairs])
    sup = SVDSuperimposer(); sup.set(fixed, moving); sup.run(); rot, tran = sup.get_rotran()
    xyz2 = {r: np.dot(x, rot)+tran for r,(aa,x) in ca2.items()}
    chain6 = json.load(open('data/trp_chain_definition.json'))['chain_6ptz_numbering']
    triad2 = [m6to2[r6] for r6 in chain6 if r6 in m6to2]
    triad_xyz = [xyz2[r2] for r2 in triad2]
    core2 = set(triad2)
    for r in ca2:
        xt = xyz2[r]
        if min(np.linalg.norm(xt-f) for f in fad) <= 4.5 or min(np.linalg.norm(xt-t) for t in triad_xyz) <= 6.0:
            core2.add(r)
    # map to cry2 MSA columns via the REF row
    from pyfamsa import Aligner, Sequence
    import glob
    files = sorted(glob.glob('data/cry2_panel/*.fa'))
    refseq = ''.join(open('data/cry2_panel/Columba_livia.fa').read().split('\n')[1:]).strip()
    seqs = [Sequence(fn.split('/')[-1][:-3].encode(), ''.join(open(fn).read().split('\n')[1:]).strip().encode()) for fn in files]
    seqs.append(Sequence(b'REF', refseq.encode()))
    msa_list = [(s.id.decode(), s.sequence.decode()) for s in Aligner().align(seqs)]
    ref_row = [s for s in msa_list if s[0]=='REF'][0][1]
    pos2col, pos = {}, 0
    for i, a in enumerate(ref_row):
        if a != '-':
            pos += 1; pos2col[pos] = i
    geocols = set(pos2col[r] for r in pos2col if r <= 400)
    core_cols = set(pos2col[r] for r in core2 if r in pos2col)
    msa = {k: v for k, v in msa_list if k != 'REF'}
    return msa, sorted(geocols - core_cols)

def W_group(msa, rows, qcols):
    vals = []
    for c in qcols:
        res = [msa[r][c] for r in rows if msa[r][c] in 'ACDEFGHIKLMNPQRSTVWY']
        if len(res) >= 0.9*len(rows):
            vals.append(Counter(res).most_common(1)[0][1]/len(res))
    return sum(vals)/len(vals)

def run_gene(msa, qcols, lab, tag):
    # msa keys: cry1 = 'migratory|Species|ACC'; cry2 = species basename
    rows = list(msa.keys())
    def species_of(r): return r.split('|')[1] if '|' in r else r
    rows = [r for r in rows if species_of(r) in lab]
    mig = [r for r in rows if lab[species_of(r)]=='migratory']
    sed = [r for r in rows if lab[species_of(r)]=='sedentary']
    Wm, Ws = W_group(msa, mig, qcols), W_group(msa, sed, qcols)
    dW = Wm - Ws
    rng = random.Random(SEED)
    cnt = 0
    allrows = mig + sed; nm = len(mig)
    for _ in range(NPERM):
        samp = set(rng.sample(allrows, nm))
        a = [r for r in allrows if r in samp]; b = [r for r in allrows if r not in samp]
        if W_group(msa, a, qcols) - W_group(msa, b, qcols) >= dW: cnt += 1
    p = (1+cnt)/(NPERM+1)
    print('%s: n_mig=%d n_sed=%d W_mig=%.5f W_sed=%.5f DeltaW=%.5f p=%.5f' % (tag, len(mig), len(sed), Wm, Ws, dW, p), flush=True)
    return {'n_mig':len(mig),'n_sed':len(sed),'W_mig':round(Wm,5),'W_sed':round(Ws,5),'DeltaW':round(dW,5),'p_one_sided':round(p,5)}

lab = labels(None)
msa1 = load_msa('data/cry1_panel_msa.fasta')
msa1 = {k:v for k,v in msa1.items() if not k.startswith('REF|')}
q1 = cry1_qcols()
r1 = run_gene(msa1, q1, lab, 'CRY1')
msa2, q2 = cry2_qcols()
r2 = run_gene(msa2, q2, lab, 'CRY2')

# sensitivity: intersection species
sp1 = {k.split('|')[1] for k in msa1}
sp2 = set(msa2.keys())
inter = sp1 & sp2
lab_i = {s: lab[s] for s in inter if s in lab}
r1i = run_gene({k:v for k,v in msa1.items() if k.split('|')[1] in lab_i}, q1, lab_i, 'CRY1-intersect')
r2i = run_gene({k:v for k,v in msa2.items() if k in lab_i}, q2, lab_i, 'CRY2-intersect')

win = r1['p_one_sided'] < 0.05 and r1['DeltaW'] > 0 and r2['p_one_sided'] >= 0.05
out = {'amendment':'2026-09-27 09:03 IST 7f-b','seed':SEED,'nperm':NPERM,
 'CRY1':r1,'CRY2':r2,'intersection':{'CRY1':r1i,'CRY2':r2i,'n_species':len(lab_i)},
 'win':bool(win),
 'falsification':'DeltaW_Cry1 n.s. OR DeltaW_Cry2 also significant same direction'}
json.dump(out, open('results/attempt7f_b_gxp.json','w'), indent=1)
print('WIN:', win)
