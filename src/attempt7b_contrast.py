#!/usr/bin/env python3
"""W02 7b PRIMARY TEST (amendment 2026-09-26 20:58 IST): cross-protein canalization contrast.
Per protein: align panel + UniProt reference (pyfamsa), core = columns hit by reference
FT functional residues; C = mean_cons(core) - mean_cons(rest); within-protein null:
10,000 random column sets of the same size, one-sided upper p. Cross-protein: Cry1 wins
if its C ranks 1/6 AND its within-protein p < 0.05. Seed 260926."""
import json, random, glob
from collections import Counter
from pyfamsa import Aligner, Sequence

SEED, NPERM = 260926, 10000
refs = json.load(open('data/attempt7b_uniprot_refs.json'))['refs']
ledger = json.load(open('data/attempt7b_controls_ledger.json'))

def cons_profile(msa_rows):
    seqs = [t[1] for t in msa_rows]
    L = len(seqs[0]); cons = []
    for i in range(L):
        res = [s[i] for s in seqs if len(s) > i and s[i] != '-']
        if not res: cons.append(None); continue
        cons.append(Counter(res).most_common(1)[0][1]/len(res))
    return cons

def eval_protein(gene):
    files = sorted(glob.glob(f'data/controls7b/*__{gene}.fa'))
    seqs = []
    for fn in files:
        lines = open(fn).read().split('\n')
        seqs.append(Sequence(fn.split('/')[-1].encode(), ''.join(lines[1:]).encode()))
    r = refs[gene]
    seqs.append(Sequence(b'REF', r['seq'].encode()))
    aln = Aligner().align(seqs)
    msa = [(s.id.decode(), s.sequence.decode()) for s in aln]
    ref_row = [s for s in msa if s[0]=='REF'][0][1]
    # reference FT residues -> MSA columns
    core_cols = set()
    pos = 0
    for i, a in enumerate(ref_row):
        if a != '-':
            pos += 1
            for f in r['features']:
                if f['start'] <= pos <= f['end']: core_cols.add(i)
    cons = cons_profile([s for s in msa if s[0] != 'REF'])
    valid = [i for i in range(len(cons)) if cons[i] is not None]
    core = [i for i in core_cols if i in set(valid)]
    rest = [i for i in valid if i not in core_cols]
    if len(core) < 2: return {'gene':gene,'error':'core<2'}
    C = sum(cons[i] for i in core)/len(core) - sum(cons[i] for i in rest)/len(rest)
    rng = random.Random(SEED)
    null = []
    for _ in range(NPERM):
        samp = rng.sample(valid, len(core))
        srest = [i for i in valid if i not in set(samp)]
        null.append(sum(cons[i] for i in samp)/len(samp) - sum(cons[i] for i in srest)/len(srest))
    null.sort()
    p = (1 + sum(x >= C for x in null))/(NPERM+1)
    return {'gene':gene,'n_seqs':len(files),'n_core_cols':len(core),
            'core_mean_cons':round(sum(cons[i] for i in core)/len(core),4),
            'rest_mean_cons':round(sum(cons[i] for i in rest)/len(rest),4),
            'C':round(C,4),'null_median':round(null[NPERM//2],4),'p_within':round(p,5)}

results = {}
for gene in ['PKM','LDHA','RHO','OPN4','GAPDH']:
    results[gene] = eval_protein(gene)
    print(gene, results[gene], flush=True)

# Cry1 contrast from the locked attempt-7 inputs (same MSA, same 36-col core)
import math
msa = {}
for line in open('data/cry1_panel_msa.fasta'):
    if line.startswith('>'): k=line[1:].strip(); msa[k]=''
    else: msa[k]+=line.strip()
rows = [k for k in msa if not k.startswith('REF|')]
geo = json.load(open('data/attempt6_column_geometry.json'))
chain6 = set(json.load(open('data/trp_chain_definition.json'))['chain_6ptz_numbering'])
core_c = {int(c) for c,g in geo.items() if g['resnum_6ptz'] in chain6 or g['fad_d']<=4.5 or g['trp_d']<=6.0}
cons = []
L = len(next(iter(msa.values())))
for i in range(L):
    res = [msa[r][i] for r in rows if msa[r][i] in 'ACDEFGHIKLMNPQRSTVWY']
    cons.append(Counter(res).most_common(1)[0][1]/len(res) if len(res)>=0.9*len(rows) else None)
valid = [i for i in range(L) if cons[i] is not None and i in {int(c) for c in geo}]
core = [i for i in valid if i in core_c]; rest = [i for i in valid if i not in core_c]
C = sum(cons[i] for i in core)/len(core) - sum(cons[i] for i in rest)/len(rest)
rng = random.Random(SEED)
null = []
for _ in range(NPERM):
    samp = set(rng.sample(valid, len(core)))
    srest = [i for i in valid if i not in samp]
    null.append(sum(cons[i] for i in samp)/len(samp) - sum(cons[i] for i in srest)/len(srest))
null.sort()
p = (1 + sum(x >= C for x in null))/(NPERM+1)
results['CRY1'] = {'gene':'CRY1','n_seqs':len(rows),'n_core_cols':len(core),
 'core_mean_cons':round(sum(cons[i] for i in core)/len(core),4),
 'rest_mean_cons':round(sum(cons[i] for i in rest)/len(rest),4),
 'C':round(C,4),'null_median':round(null[NPERM//2],4),'p_within':round(p,5)}
print('CRY1', results['CRY1'])

rank = 1 + sum(1 for g in ['PKM','LDHA','RHO','OPN4','GAPDH'] if results[g].get('C',-9) > results['CRY1']['C'])
out = {'amendment':'2026-09-26 20:58 IST 7b','nperm':NPERM,'seed':SEED,'results':results,
 'cry1_rank_of_6':rank, 'win': bool(rank==1 and results['CRY1']['p_within']<0.05)}
json.dump(out, open('results/attempt7b_canalization.json','w'), indent=1)
print('WIN:', out['win'], 'rank:', rank)
