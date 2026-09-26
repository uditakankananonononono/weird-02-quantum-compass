#!/usr/bin/env python3
"""W02 attempt-6 variant-effect panel scorer. LOCKED PRE-OUTCOME (2026-09-26 ~20:20 IST):
this file is committed BEFORE any evaluation run; the formula, weights, and null
definitions below are frozen by commit hash. No tuning after outcome inspection
(amendment item 4). Deliberately excludes conservation from the main score -
conservation-only ranking is one of the two locked null baselines.

Deliverable (amendment item 1): ranked table of natural panel variants with
predicted magnetosensitivity-relevant effect + structural rationale.
Benchmark (item 2): enrichment of locked truth-set positions in top-k=20 vs
(a) random ranking null, (b) conservation-only ranking null; 1000 perms, one-sided p<0.05.
"""
import json, math, glob
import numpy as np
from Bio import SeqIO

D = 'data/'
TOPK = 20            # locked
NPERM = 1000         # locked
GAPCAP = 0.90        # locked: column scored only if >=90% gap-free across non-REF rows
RNG_SEED = 260926    # locked

# ---- fixed physicochemical delta table (locked; charge/aromatic/volume) ----
CHARGE = {'D':-1,'E':-1,'K':1,'R':1,'H':0.5}
AROMATIC = {'F','W','Y','H'}
VOL = {'A':88.6,'R':173.4,'N':114.1,'D':111.1,'C':108.5,'Q':143.8,'E':138.4,'G':60.1,
       'H':153.2,'I':166.7,'L':166.7,'K':168.6,'M':162.9,'F':189.9,'P':112.7,'S':89.0,
       'T':116.1,'W':227.8,'Y':193.6,'V':140.0}
def phys_delta(a, b):
    dq = abs(CHARGE.get(a,0)-CHARGE.get(b,0))
    da = abs((a in AROMATIC)-(b in AROMATIC))
    dv = abs(VOL[a]-VOL[b])/170.0  # normalized to ~max
    return dq + 1.5*da + 0.5*dv    # weights locked

# ---- load frozen MSA ----
msa = {}
for line in open(D+'cry1_panel_msa.fasta'):
    if line.startswith('>'): k=line[1:].strip(); msa[k]=''
    else: msa[k]+=line.strip()
rows = [k for k in msa if not k.startswith('REF|')]
L = len(next(iter(msa.values())))

# ---- reference-fold geometry (precomputed by src/attempt6_geometry.py from 6PTZ; locked inputs) ----
geo = json.load(open(D+'attempt6_column_geometry.json'))  # {str(col): {fad_d, trp_d, resnum_6ptz, aa6ptz}}
chain6 = set(json.load(open(D+'trp_chain_definition.json'))['chain_6ptz_numbering'])
chain_cols = {int(c) for c,g in geo.items() if g['resnum_6ptz'] in chain6}
fad_d = {int(c): g['fad_d'] for c,g in geo.items()}
trp_d = {int(c): g['trp_d'] for c,g in geo.items()}

W_CHAIN_BREAK, W_CHAIN_GAIN = 3.0, 2.0   # locked
def prox_w(d): return math.exp(-max(d-4.5,0)/6.0)  # locked decay; FAD contact zone full weight

variants = []
for col in range(L):
    residues = [msa[r][col] for r in rows if msa[r][col] in VOL]  # canonical AA only
    if len(residues) < GAPCAP*len(rows): continue
    if col not in fad_d: continue  # outside folded domain map
    from collections import Counter
    cons, freq = Counter(residues).most_common(1)[0]
    for a in set(residues):
        if a == cons: continue
        carriers = [r for r in rows if msa[r][col]==a]
        pd = phys_delta(cons, a)
        s = prox_w(fad_d[col])*pd + math.exp(-trp_d[col]/6.0)*pd
        rationale = []
        if col in chain_cols:
            if cons=='W' and a!='W': s += W_CHAIN_BREAK; rationale.append('Trp-chain break')
            elif cons!='W' and a=='W': s += W_CHAIN_GAIN; rationale.append('novel chain-adjacent Trp')
        if fad_d[col] <= 4.5: rationale.append('FAD contact zone')
        if trp_d[col] <= 6.0: rationale.append('ET-pathway adjacent')
        variants.append({'col': col, 'from': cons, 'to': a, 'score': round(s,4),
                         'n_carriers': len(carriers), 'rationale': '; '.join(rationale) or 'peripheral'})

variants.sort(key=lambda v: -v['score'])
json.dump({'locked': {'topk': TOPK, 'nperm': NPERM, 'gapcap': GAPCAP, 'seed': RNG_SEED,
                      'weights': {'chain_break': W_CHAIN_BREAK, 'chain_gain': W_CHAIN_GAIN}},
           'n_variants': len(variants), 'ranking': variants},
          open('results/attempt6_ranked_panel.json','w'), indent=1)
print(json.dumps({'n_variants': len(variants), 'top3': [f"{v['from']}{v['col']}{v['to']}:{v['score']}" for v in variants[:3]]}))
