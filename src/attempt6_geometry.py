#!/usr/bin/env python3
"""Attempt-6 input precompute: per-MSA-column distances in the 6PTZ reference frame.
Mechanical transform of already-grounded inputs (6PTZ + frozen MSA + locked chain def);
same alignment method as src/chain_define.py. No scoring, no outcomes."""
import json
import numpy as np
from Bio import SeqIO
from Bio.Align import PairwiseAligner

AA3={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'F','PHE':'F','TRP':'W','TYR':'Y','VAL':'V','PRO':'P','SER':'S','THR':'T','MET':'M'}
# (duplicate keys harmless)
ca = {}; fad = []
for line in open('/tmp/6ptz.pdb'):
    if line.startswith('ATOM') and line[12:16].strip()=='CA':
        r=int(line[22:26])
        ca[r]=(line[17:20].strip(), np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])]))
    elif line.startswith('HETATM') and line[17:20].strip()=='FAD':
        fad.append(np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])]))
resnums = sorted(ca)
seq6 = ''.join(AA3.get(ca[r][0],'X') for r in resnums)

refs={r.id.split('|')[1]: str(r.seq) for r in SeqIO.parse('data/ref_sequences.fasta','fasta')}
clcry4 = refs['A0A386QUR4']
al = PairwiseAligner().align(seq6, clcry4)[0]
p2u = {}
for (ps,pe),(us,ue) in zip(*al.aligned):
    for i in range(pe-ps): p2u[resnums[ps+i]] = us+i+1

msa={}
for line in open('data/cry1_panel_msa.fasta'):
    if line.startswith('>'): k=line[1:].strip(); msa[k]=''
    else: msa[k]+=line.strip()
refrow = next(k for k in msa if 'A0A386QUR4' in k)
u2col = {}
col = -1
for i,ch in enumerate(msa[refrow]):
    col += 1
    if ch != '-':
        pass
# build uniprot-index -> column by counting non-gap chars
u2col = {}
u = 0
for ci,ch in enumerate(msa[refrow]):
    if ch != '-':
        u += 1
        u2col[u] = ci

chain6 = json.load(open('data/trp_chain_definition.json'))['chain_6ptz_numbering']
chain_xyz = [ca[r][1] for r in chain6]
out = {}
for r,(aa,x) in ca.items():
    if r not in p2u or p2u[r] not in u2col: continue
    cidx = u2col[p2u[r]]
    out[str(cidx)] = {
        'aa6ptz': AA3.get(aa,'X'),
        'fad_d': round(float(min(np.linalg.norm(x-f) for f in fad)),3),
        'trp_d': round(float(min(np.linalg.norm(x-t) for t in chain_xyz)),3),
        'resnum_6ptz': r}
json.dump(out, open('data/attempt6_column_geometry.json','w'), indent=1)
print(json.dumps({'n_cols_mapped': len(out),
                  'chain_cols_present': [str(u2col.get(p2u.get(r,-1),-1)) in out for r in chain6]}))
