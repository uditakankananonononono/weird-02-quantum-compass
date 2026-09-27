#!/usr/bin/env python3
"""W02 queue #6 (12:31 IST amendment): agreement-is-not-function control. Report-only.
Per-residue pLDDT (11 committed Cry4 folds) vs 6PTZ solved-coverage vs functional-region
membership for the 37 attempt8b Cry4-distinguishing residues. Point-biserial r + counts."""
import json, glob, numpy as np
from Bio.Align import PairwiseAligner
from scipy.stats import pointbiserialr

CORE37 = json.load(open('results/attempt8b_specialization.json'))['candidates']
AA3={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','TRP':'W','TYR':'Y','VAL':'V','PRO':'P','SER':'S','THR':'T'}

def load_msa(p):
    m={}
    for l in open(p):
        if l.startswith('>'): k=l[1:].strip(); m[k]=''
        else: m[k]+=l.strip()
    return m

# 6PTZ: CA map, FAD, SS records
ca6, fad, ss6 = {}, [], {}
cur=None
for l in open('/tmp/6ptz.pdb'):
    if l.startswith('HELIX'):
        for r in range(int(l[21:25]), int(l[33:37])+1): ss6[r]='H'
    elif l.startswith('SHEET'):
        for r in range(int(l[22:26]), int(l[33:37])+1): ss6[r]='E'
    elif l.startswith('ATOM') and l[12:16].strip()=='CA':
        ca6[int(l[22:26])]=(l[17:20].strip(), np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])]))
    elif l.startswith('HETATM') and l[17:20].strip()=='FAD':
        fad.append(np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])]))
resnums6=sorted(ca6); seq6=''.join(AA3.get(ca6[r][0],'X') for r in resnums6)
fad_cent=np.mean(fad,axis=0)
chain6=json.load(open('data/trp_chain_definition.json'))['chain_6ptz_numbering']
trp_ca=[ca6[r][1] for r in chain6 if r in ca6]

msa=load_msa('data/mafft_ebi_crosscheck.fasta')
kx=msa['KX168611']; kx_seq=kx.replace('-','')
al=PairwiseAligner().align(kx_seq, seq6)[0]
# kx position (1-based) -> 6PTZ resnum, and kx position -> MSA column
kxpos2col={}; pos=0
for i,a in enumerate(kx):
    if a!='-': pos+=1; kxpos2col[pos]=i
kxpos2r6={}
for (sk,ek),(s6,e6) in zip(*al.aligned):
    for i in range(ek-sk): kxpos2r6[sk+1+i]=resnums6[s6+i]
r6tokxpos={v:k for k,v in kxpos2r6.items()}

# folds
folds={}
for fp in sorted(glob.glob('data/folds/CRY4__*.pdb')):
    name=fp.split('/')[-1][:-4]; acc=name.split('__')[3]
    cas={}
    for l in open(fp):
        if l.startswith('ATOM') and l[12:16].strip()=='CA':
            cas[int(l[22:26])]=(l[17:20].strip(), np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])]), float(l[60:66]))
    bs=np.array([v[2] for k,v in sorted(cas.items())])
    scale01 = np.median(bs) <= 1.5
    folds[acc]={'name':name,'cas':cas,'scale':'0-1' if scale01 else '0-100'}

def fold_plddt(acc, col):
    row=msa.get(acc)
    if row is None: return None
    pos=0; idx=None
    for i,a in enumerate(row):
        if a!='-': pos+=1
        if i==col:
            idx=pos if a!='-' else None; break
    if idx is None: return None
    c=folds[acc]['cas'].get(idx)
    if c is None: return None
    v=c[2]*100.0 if folds[acc]['scale']=='0-1' else c[2]
    return v

rows=[]
for r6 in CORE37:
    kxpos=r6tokxpos.get(r6); col=kxpos2col.get(kxpos) if kxpos else None
    solved = r6 in ca6
    in_shell=False
    if solved:
        x=ca6[r6][1]
        in_shell = (np.linalg.norm(x-fad_cent)<=8.0) or any(np.linalg.norm(x-t)<=8.0 for t in trp_ca)
    is319 = (r6==319)
    vals={acc: fold_plddt(acc,col) for acc in folds} if col is not None else {}
    vals={k:v for k,v in vals.items() if v is not None}
    med=float(np.median(list(vals.values()))) if vals else None
    rows.append({'residue_6ptz':r6,'msa_col':col,'solved_6ptz':solved,
      'crystal_ss':ss6.get(r6,'C') if solved else None,
      'in_fad_trp_8A_shell':in_shell,'is_position_319':is319,
      'functional_member':in_shell or is319,
      'plddt_per_fold':{k:round(v,2) for k,v in vals.items()},
      'plddt_median':round(med,2) if med is not None else None,
      'n_folds_mapped':len(vals)})

# SS agreement: pigeon fold KX168611 CA i,i+3 vs crystal records, on solved core residues
pig=folds['KX168611']['cas']
agree=disagree=0; ss_detail={}
for r in rows:
    if not r['solved_6ptz'] or r['crystal_ss'] in (None,'C'): continue
    kxpos=r6tokxpos.get(r['residue_6ptz'])
    if kxpos is None or kxpos not in pig or kxpos+3 not in pig: continue
    d=np.linalg.norm(pig[kxpos][1]-pig[kxpos+3][1])
    call='H' if d<=7.0 else ('E' if d>=10.0 else 'C')
    ok = call==r['crystal_ss']
    agree+=ok; disagree+=(not ok)
    ss_detail[r['residue_6ptz']]={'crystal':r['crystal_ss'],'fold_call':call,'d_i_i3':round(float(d),2),'agree':ok}

y=[1 if r['functional_member'] else 0 for r in rows]
x=[r['plddt_median'] if r['plddt_median'] is not None else np.nan for r in rows]
mask=[i for i in range(len(rows)) if not np.isnan(x[i])]
r_pb,p_pb=pointbiserialr([y[i] for i in mask],[x[i] for i in mask])
out={
 'amendment':'2026-09-27 12:31 IST queue #6 (delegated executor, report-only)',
 'executor_note':'computed by delegated executor at ref 4845c5e; uncommitted, returned to lane-29 for review',
 'plddt_source':'11 committed Cry4 folds in data/folds (scale detection per file: median B<=1.5 -> 0-1 x100, else 0-100; 11:28 clarification)',
 'mapping':'KX168611 pairwise-aligned to 6PTZ CA sequence (PairwiseAligner, same as 8/8b); per-fold residue = non-gap walk of its MSA row to the mapped column',
 'n_residues':len(rows),
 'coverage_counts':{
   'solved_6ptz':sum(r['solved_6ptz'] for r in rows),
   'not_solved_6ptz':sum(not r['solved_6ptz'] for r in rows),
   'in_fad_trp_8A_shell':sum(r['in_fad_trp_8A_shell'] for r in rows),
   'position_319':1,
   'functional_members':sum(r['functional_member'] for r in rows),
   'with_plddt':sum(r['plddt_median'] is not None for r in rows),
   'census_entries_context':len(json.load(open('data/rcsb_cry_census_details.json')))},
 'census_limitation':'data/rcsb_cry_census_details.json carries method/resolution/title only, no per-structure residue ranges; per-residue solved coverage is therefore computed against 6PTZ (the numbering anchor, Y319D) only. Census is context (34 solved CRY structures), not per-residue coverage.',
 'ss_method':'6PTZ HELIX/SHEET records vs pigeon fold (KX168611) CA d(i,i+3) heuristic (<=7.0A helix-call, >=10.0A sheet-call); DSSP unavailable in executor env (pydssp pulls multi-GB torch, barred by protocol) - coarse agreement only',
 'ss_agreement':{'n_compared':agree+disagree,'agree':agree,'disagree':disagree,'detail':ss_detail},
 'point_biserial':{'r':round(float(r_pb),4),'p_two_sided':round(float(p_pb),6),'n':len(mask),
   'x':'plddt_median across folds (0-100)','y':'functional-region membership (8A shell or position 319)'},
 'per_residue':rows,
 'interpretation_lock':'REPORT-ONLY. High agreement/confidence does NOT imply function; this amendment can only LIMIT interpretation, never upgrade a residue to functional. No causal-residue claims (locked anti-goal).'}
json.dump(out, open('results/attempt8f_structure_agreement.json','w'), indent=1)
print('r=',out['point_biserial']['r'],'p=',out['point_biserial']['p_two_sided'],'n=',len(mask))
print('coverage:',out['coverage_counts'])
print('ss:',agree,'/',agree+disagree)
meds=[r['plddt_median'] for r in rows if r['plddt_median'] is not None]
print('plddt median range:',min(meds),max(meds))
