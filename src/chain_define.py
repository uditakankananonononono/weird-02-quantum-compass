#!/usr/bin/env python3
"""Build trp_chain_definition.json from PDB 6PTZ (grounded 20:35).
Chain (6PTZ numbering, FAD->surface): W395 -> W372 -> W318 -> W369.
Edge hops (A): FAD-395 3.63, 395-372 5.05, 372-318 3.82, 318-369 3.73.
Maps to MSA columns via the REF ClCry4 row and to erCry4 residue indices via its row.
"""
import json
from Bio import SeqIO
from Bio.Align import PairwiseAligner

CHAIN_6PTZ = [395, 372, 318, 369]
HOPS = {'FAD_to_W395': 3.63, 'W395_W372': 5.05, 'W372_W318': 3.82, 'W318_W369': 3.73}

AA3={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
seq6=''; seen=set()
for line in open('/tmp/6ptz.pdb'):
    if line.startswith('ATOM') and line[12:16].strip()=='CA':
        r=int(line[22:26])
        if r not in seen: seen.add(r); seq6+=AA3.get(line[17:20].strip(),'X')
resnums6 = sorted(seen)

refs={r.id.split('|')[1]: str(r.seq) for r in SeqIO.parse('data/ref_sequences.fasta','fasta')}
clcry4 = refs['A0A386QUR4']
al = PairwiseAligner().align(seq6, clcry4)[0]
p2u = {}
for (ps,pe),(us,ue) in zip(*al.aligned):
    for i in range(pe-ps): p2u[resnums6[ps+i]] = us+i+1
chain_uniprot = [int(p2u[r]) for r in CHAIN_6PTZ]
print('6PTZ chain -> ClCry4(A0A386QUR4) residue numbers:', chain_uniprot)
assert all(clcry4[u-1]=='W' for u in chain_uniprot), 'not all W in ClCry4!'

msa={}
for line in open('data/cry1_panel_msa.fasta'):
    if line.startswith('>'): k=line[1:].strip(); msa[k]=''
    else: msa[k]+=line.strip()
clkey=[k for k in msa if 'A0A386QUR4' in k][0]
erkey=[k for k in msa if 'A0A2I4SZI9' in k][0]
def res2col(row):
    m={}; ri=0
    for c,ch in enumerate(row):
        if ch!='-': ri+=1; m[ri]=c
    return m
cl_r2c=res2col(msa[clkey]); er_r2c=res2col(msa[erkey])
cols=[int(cl_r2c[u]) for u in chain_uniprot]
# erCry4 residue indices at those columns
col2er={c:r for r,c in er_r2c.items()}
er_idx=[int(col2er[c]) for c in cols]
er_aas=[msa[erkey][c] for c in cols]
cl_aas=[msa[clkey][c] for c in cols]
print('MSA cols:', cols, 'ClCry4 aas:', cl_aas, 'erCry4 residx:', er_idx, 'erCry4 aas:', er_aas)
# triad-vs-tetrad check: erCry1/ClCry1 at those columns
for rk in ['Q5IZC5','X5D0Q5']:
    k=[kk for kk in msa if rk in kk][0]
    print(rk, 'at chain cols:', [msa[k][c] for c in cols])
out={'pdb_source':'6PTZ (Xu et al. 2021, ClCry4 Y319D + FAD; COMPND curation error verified 99.6% ClCry4 identity)',
     'chain_6ptz_numbering': CHAIN_6PTZ, 'chain_order':'FAD->surface',
     'edge_hops_A_6PTZ': HOPS, 'chain_clcry4_uniprot': chain_uniprot,
     'msa_cols': cols, 'ref_chain_residx': er_idx, 'ref_note':'ref_chain_residx in erCry4 A0A2I4SZI9 numbering',
     'clcry4_aas': cl_aas, 'ercry4_aas': er_aas}
json.dump(out, open('data/trp_chain_definition.json','w'), indent=1)
print('definition written')
