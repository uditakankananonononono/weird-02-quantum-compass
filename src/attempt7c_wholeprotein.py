#!/usr/bin/env python3
"""W02 7c PRIMARY TEST (amendment 21:25 IST): whole-protein canalization contrast.
W = mean conservation over qualified non-core columns. Win: Cry1 W ranks 1/6 AND
all 5 pairwise permutation tests (pooled background columns, label re-sampling,
10k perms, seed 260926) one-sided p<0.05."""
import json, random, glob
from collections import Counter
from pyfamsa import Aligner, Sequence

SEED, NPERM = 260926, 10000
refs = json.load(open('data/attempt7b_uniprot_refs.json'))['refs']

def bg_cons(gene):
    files = sorted(glob.glob(f'data/controls7b/*__{gene}.fa'))
    seqs = [Sequence(fn.split('/')[-1].encode(), ''.join(open(fn).read().split('\n')[1:]).encode()) for fn in files]
    r = refs[gene]
    seqs.append(Sequence(b'REF', r['seq'].encode()))
    msa = [(s.id.decode(), s.sequence.decode()) for s in Aligner().align(seqs)]
    ref_row = [s for s in msa if s[0]=='REF'][0][1]
    core_cols, pos = set(), 0
    for i, a in enumerate(ref_row):
        if a != '-':
            pos += 1
            for f in r['features']:
                if f['start'] <= pos <= f['end']: core_cols.add(i)
    seqrows = [t[1] for t in msa if t[0] != 'REF']
    vals = []
    for i in range(len(ref_row)):
        if i in core_cols: continue
        res = [s[i] for s in seqrows if s[i] != '-']
        if len(res) >= 0.9*len(seqrows):
            vals.append(Counter(res).most_common(1)[0][1]/len(res))
    return vals

def cry1_bg():
    msa = {}
    for line in open('data/cry1_panel_msa.fasta'):
        if line.startswith('>'): k=line[1:].strip(); msa[k]=''
        else: msa[k]+=line.strip()
    rows = [k for k in msa if not k.startswith('REF|')]
    geo = json.load(open('data/attempt6_column_geometry.json'))
    chain6 = set(json.load(open('data/trp_chain_definition.json'))['chain_6ptz_numbering'])
    core_c = {int(c) for c,g in geo.items() if g['resnum_6ptz'] in chain6 or g['fad_d']<=4.5 or g['trp_d']<=6.0}
    geocols = {int(c) for c in geo}
    vals = []
    L = len(next(iter(msa.values())))
    for i in range(L):
        if i in core_c or i not in geocols: continue
        res = [msa[r][i] for r in rows if msa[r][i] in 'ACDEFGHIKLMNPQRSTVWY']
        if len(res) >= 0.9*len(rows):
            vals.append(Counter(res).most_common(1)[0][1]/len(res))
    return vals

cry = cry1_bg()
W = {'CRY1': sum(cry)/len(cry)}
pools = {'CRY1': cry}
for g in ['PKM','LDHA','RHO','OPN4','GAPDH']:
    v = bg_cons(g); pools[g] = v; W[g] = sum(v)/len(v)
    print(g, 'W=%.4f n=%d' % (W[g], len(v)), flush=True)
print('CRY1 W=%.4f n=%d' % (W['CRY1'], len(cry)))

rng = random.Random(SEED)
pair_p = {}
for g in ['PKM','LDHA','RHO','OPN4','GAPDH']:
    obs = W['CRY1'] - W[g]
    a, b = pools['CRY1'], pools[g]
    pooled = a + b; na = len(a)
    cnt = 0
    for _ in range(NPERM):
        samp = rng.sample(pooled, na)
        wa = sum(samp)/na
        rest_sum = sum(pooled) - sum(samp)
        wb = rest_sum/(len(pooled)-na)
        if wa - wb >= obs: cnt += 1
    pair_p[g] = round((1+cnt)/(NPERM+1), 5)
    print('pair CRY1 vs', g, 'diff=%.4f p=%.5f' % (obs, pair_p[g]), flush=True)

rank = 1 + sum(1 for g in ['PKM','LDHA','RHO','OPN4','GAPDH'] if W[g] > W['CRY1'])
win = rank == 1 and all(p < 0.05 for p in pair_p.values())
out = {'amendment':'2026-09-26 21:25 IST 7c','nperm':NPERM,'seed':SEED,
 'W': {k: round(v,5) for k,v in W.items()}, 'pairwise_p': pair_p,
 'cry1_rank_of_6': rank, 'win': bool(win)}
json.dump(out, open('results/attempt7c_wholeprotein.json','w'), indent=1)
print('WIN:', win, 'rank:', rank)
