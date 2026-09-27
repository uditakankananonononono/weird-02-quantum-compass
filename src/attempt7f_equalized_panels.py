#!/usr/bin/env python3
"""W02 queue #18 equalized panels (amendment 80f21d0, locked pre-computation).
Species intersection across CRY1 + 5 controls; identical 7b pipeline restricted.
Report-only sensitivity; original 7b stands. Seed 260926."""
import json, random, glob
from collections import Counter
from pyfamsa import Aligner, Sequence

SEED, NPERM = 260926, 10000
refs = json.load(open('data/attempt7b_uniprot_refs.json'))['refs']

# species sets
ctrl_species = {}
for gene in ['PKM','LDHA','RHO','OPN4','GAPDH']:
    ctrl_species[gene] = {fn.split('/')[-1].split('__')[0] for fn in glob.glob(f'data/controls7b/*__{gene}.fa')}
msa = {}
for line in open('data/cry1_panel_msa.fasta'):
    if line.startswith('>'): k=line[1:].strip(); msa[k]=''
    else: msa[k]+=line.strip()
cry1_species = {k.split('|')[1] for k in msa if not k.startswith('REF|')}
inter = set(cry1_species)
for g, s in ctrl_species.items(): inter &= s
inter = sorted(inter)

def cons_profile(rows):
    seqs = [t[1] for t in rows]
    L = len(seqs[0]); cons = []
    for i in range(L):
        res = [s[i] for s in seqs if len(s) > i and s[i] != '-']
        cons.append(None if not res else Counter(res).most_common(1)[0][1]/len(res))
    return cons

def eval_control(gene):
    files = sorted(fn for fn in glob.glob(f'data/controls7b/*__{gene}.fa')
                   if fn.split('/')[-1].split('__')[0] in set(inter))
    seqs = []
    for fn in files:
        lines = open(fn).read().split('\n')
        seqs.append(Sequence(fn.split('/')[-1].encode(), ''.join(lines[1:]).encode()))
    r = refs[gene]
    seqs.append(Sequence(b'REF', r['seq'].encode()))
    aln = Aligner().align(seqs)
    msa_rows = [(s.id.decode(), s.sequence.decode()) for s in aln]
    ref_row = [s for s in msa_rows if s[0]=='REF'][0][1]
    core_cols = set(); pos = 0
    for i, a in enumerate(ref_row):
        if a != '-':
            pos += 1
            for f in r['features']:
                if f['start'] <= pos <= f['end']: core_cols.add(i)
    cons = cons_profile([s for s in msa_rows if s[0] != 'REF'])
    valid = [i for i in range(len(cons)) if cons[i] is not None]
    core = [i for i in core_cols if i in set(valid)]
    rest = [i for i in valid if i not in core_cols]
    if len(core) < 2: return {'gene':gene,'error':'core<2'}
    C = sum(cons[i] for i in core)/len(core) - sum(cons[i] for i in rest)/len(rest)
    rng = random.Random(SEED)
    null = []
    for _ in range(NPERM):
        samp = set(rng.sample(valid, len(core)))
        srest = [i for i in valid if i not in samp]
        null.append(sum(cons[i] for i in samp)/len(samp) - sum(cons[i] for i in srest)/len(srest))
    null.sort()
    p = (1 + sum(x >= C for x in null))/(NPERM+1)
    return {'gene':gene,'n_seqs':len(files),'n_core_cols':len(core),
            'C':round(C,4),'null_median':round(null[NPERM//2],4),'p_within':round(p,5)}

results = {}
for gene in ['PKM','LDHA','RHO','OPN4','GAPDH']:
    results[gene] = eval_control(gene)
    print(gene, results[gene], flush=True)

# Cry1 restricted
rows = {k:v for k,v in msa.items() if not k.startswith('REF|') and k.split('|')[1] in set(inter)}
geo = json.load(open('data/attempt6_column_geometry.json'))
chain6 = set(json.load(open('data/trp_chain_definition.json'))['chain_6ptz_numbering'])
core_c = {int(c) for c,g in geo.items() if g['resnum_6ptz'] in chain6 or g['fad_d']<=4.5 or g['trp_d']<=6.0}
L = len(next(iter(rows.values())))
cons = []
for i in range(L):
    res = [v[i] for v in rows.values() if v[i] in 'ACDEFGHIKLMNPQRSTVWY']
    cons.append(None if not res else Counter(res).most_common(1)[0][1]/len(res))
valid = [i for i in range(L) if cons[i] is not None]
core = sorted(core_c & set(valid))
rest = [i for i in valid if i not in core_c]
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
                   'C':round(C,4),'null_median':round(null[NPERM//2],4),'p_within':round(p,5)}
print('CRY1', results['CRY1'], flush=True)

ranked = sorted(results.values(), key=lambda r: -r['C'])
for i, r in enumerate(ranked, 1): r['rank_equalized'] = i
out = {'amendment_commit':'80f21d0','rule':'report-only sensitivity; original 7b stands',
       'intersection_n_species': len(inter), 'intersection_species': inter,
       'results': results}
json.dump(out, open('results/attempt7f_equalized_panels.json','w'), indent=1)
print('intersection', len(inter), 'species; CRY1 rank', results['CRY1']['rank_equalized'], 'of 6')
