#!/usr/bin/env python3
"""W02 feature extractor v1: frozen Cry1 panel -> sequence features.
Panel rebuilt deterministically under FROZEN RULES R1/R2/R5 (panel_frozen.md).
"""
import json, re
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from pyfamsa import Aligner, Sequence

D = 'data/'
pl = json.load(open(D+'paralog_labels.json'))
spmap = json.load(open(D+'cry_aves_species.json'))
mig = json.load(open(D+'migration_labels_raw.json'))

def lab(v): return v if isinstance(v, str) else v.get('label')

# species -> list of (accession, length) for Cry1
seqs = {r.id: str(r.seq) for r in SeqIO.parse(D+'cry_aves_proteins.fasta','fasta')}
by_sp = defaultdict(list)
for acc, l in pl.items():
    if lab(l) != 'Cry1': continue
    key = acc.split('.')[0]
    info = spmap.get(key) or spmap.get(acc)
    if not info: continue
    sp = info['organism']
    # find matching seq record (headers may carry version)
    s = seqs.get(acc) or seqs.get(key) or next((v for k,v in seqs.items() if k.startswith(key)), None)
    if s: by_sp[sp].append((acc, len(s), s))

def rep(cands):  # R5: longest, tie-break lowest accession
    return sorted(cands, key=lambda t: (-t[1], t[0]))[0]

panel = {'migratory': {}, 'sedentary': {}, 'partial': {}}
for sp, cands in by_sp.items():
    m = mig.get(sp)
    if m == 3: panel['migratory'][sp] = rep(cands)
    elif m == 1: panel['sedentary'][sp] = rep(cands)
    elif m == 2: panel['partial'][sp] = rep(cands)

n = {k: len(v) for k,v in panel.items()}
print('panel counts (expect 36 mig / 80 sed per freeze):', n)
assert n['migratory'] == 36 and n['sedentary'] == 80, 'COUNT MISMATCH vs freeze'
json.dump({k: {sp: [a, L] for sp,(a,L,s) in v.items()} for k,v in panel.items()},
          open(D+'panel_frozen_members.json','w'), indent=1)

# --- MSA (FAMSA) over frozen primary contrast + references ---
refs = list(SeqIO.parse(D+'ref_sequences.fasta','fasta'))
items = []
for grp in ['migratory','sedentary']:
    for sp,(a,L,s) in panel[grp].items():
        items.append((f'{grp}|{sp.replace(" ","_")}|{a}', s))
for r in refs:
    items.append((f'REF|{r.id}', str(r.seq)))
al = Aligner()
seq_objs = [Sequence(bytes(i,'utf8'), bytes(s,'utf8')) for i,s in items]
msa = al.align(seq_objs)
with open(D+'cry1_panel_msa.fasta','w') as fh:
    for m in msa: fh.write(f'>{m.id.decode()}\n{m.sequence.decode()}\n')
print('MSA written:', len(items), 'seqs')

# --- features v1 (sequence-level; structure features join later) ---
import numpy as np
aln = {m.id.decode(): m.sequence.decode() for m in msa}
ref1 = aln['REF|sp|Q5IZC5|CRY1_ERIRU'] if 'REF|sp|Q5IZC5|CRY1_ERIRU' in aln else aln[[k for k in aln if 'Q5IZC5' in k][0]]
ref4 = aln[[k for k in aln if 'A0A2I4SZI9' in k][0]]

def tail_len(aligned, ref):
    # residues after the last aligned column where ref has a residue
    last = max(i for i,c in enumerate(ref) if c != '-')
    return sum(1 for c in aligned[last+1:] if c != '-')

def pairwise_ident(a, b):
    pairs = [(x,y) for x,y in zip(a,b) if x!='-' and y!='-']
    return sum(x==y for x,y in pairs)/len(pairs) if pairs else 0.0

def trp_signature(s):
    pos = [i for i,c in enumerate(s) if c=='W']
    gaps = [b-a for a,b in zip(pos,pos[1:])]
    return {'trp_total': len(pos), 'trp_cterm150': sum(1 for p in pos if p > len(s)-150),
            'trp_gaps': gaps}

def kd(s):
    scale = {'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,
             'H':-3.2,'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,
             'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}
    vals = [scale[c] for c in s if c in scale]
    return sum(vals)/len(vals) if vals else 0.0

feats = {}
for i,(name, raw) in enumerate(items):
    if name.startswith('REF|'): continue
    a = aln[name]
    s = raw
    pa = ProteinAnalysis(s)
    grp = name.split('|')[0]
    feats[name] = {'group': grp, 'length': len(s),
        'tail_after_ref': tail_len(a, ref1),
        'ident_erCry1': round(pairwise_ident(a, ref1),4),
        'ident_erCry4': round(pairwise_ident(a, ref4),4),
        'pI': round(pa.isoelectric_point(),3),
        'charge_pH7': round(pa.charge_at_pH(7.0),2),
        'aromatic_frac': round(sum(1 for c in s if c in 'FWY')/len(s),4),
        'kd_mean': round(kd(s),3),
        **trp_signature(s)}
json.dump(feats, open(D+'cry1_features_v1.json','w'), indent=1)
print('features:', len(feats))
# quick sanity: group means for headline features
for grp in ['migratory','sedentary']:
    sub = [f for f in feats.values() if f['group']==grp]
    print(grp, 'n=',len(sub),
          'mean_len', round(np.mean([f['length'] for f in sub]),1),
          'mean_tail', round(np.mean([f['tail_after_ref'] for f in sub]),1),
          'mean_trp_ct', round(np.mean([f['trp_cterm150'] for f in sub]),2),
          'mean_ident_erCry1', round(np.mean([f['ident_erCry1'] for f in sub]),3))
