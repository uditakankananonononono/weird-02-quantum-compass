#!/usr/bin/env python3
"""Assemble Cry4 sub-panel: 4 RefSeq species (R5 rep) + 7 tree-verified mined orthologs."""
import json, gzip, subprocess, tempfile, os
from collections import defaultdict
from Bio import SeqIO
import pyhmmer

D='data/'
# 1) RefSeq Cry4 representatives (R5)
pl=json.load(open(D+'paralog_labels.json')); spmap=json.load(open(D+'cry_aves_species.json'))
seqs={r.id.split('.')[0]: (r.id, str(r.seq)) for r in SeqIO.parse(D+'cry_aves_proteins.fasta','fasta')}
by_sp=defaultdict(list)
for acc,l in pl.items():
    if (l if isinstance(l,str) else l.get('label'))!='Cry4': continue
    key=acc.split('.')[0]; info=spmap.get(key)
    if info and key in seqs: by_sp[info['organism']].append(seqs[key])
out={}
for sp,cands in by_sp.items():
    rid,s = sorted(cands, key=lambda t:(-len(t[1]), t[0]))[0]
    out[sp]=(rid,s)
print('refseq cry4:', {k:v[0] for k,v in out.items()})

# 2) mined orthologs: top hmmsearch hit per extra proteome (seed HMM = the one built from 12 classified Cry4)
with tempfile.NamedTemporaryFile('wb',suffix='.hmm',delete=False) as tf:
    tf.write(open(D+'cry4_seed.hmm','rb').read()); hmmpath=tf.name
alphabet=pyhmmer.easel.Alphabet.amino()
with pyhmmer.plan7.HMMFile(hmmpath) as hf:
    hmm=next(hf)
mined={}
for f in sorted(os.listdir(D+'proteomes_extra')):
    sp=f.replace('.pep.fa.gz','').replace('_',' ')
    path=os.path.join(D+'proteomes_extra',f)
    opener=gzip.open if f.endswith('.gz') else open
    with opener(path,'rt') as fh:
        recs=list(SeqIO.parse(fh,'fasta'))
    with tempfile.NamedTemporaryFile('w',suffix='.faa',delete=False) as tf:
        for r in recs: tf.write(f'>{r.id}\n{r.seq}\n')
        faapath=tf.name
    with pyhmmer.easel.SequenceFile(faapath, digital=True, alphabet=alphabet) as sf:
        ds=sf.read_block()
    res=list(pyhmmer.hmmsearch(hmm, ds, cpus=4))
    hits=list(res[0]) if res else []
    if hits:
        h=hits[0]
        sid=h.name if isinstance(h.name,str) else h.name.decode()
        seq=next(str(r.seq) for r in recs if r.id==sid)
        mined[sp]=(sid, seq, h.evalue)
    os.unlink(faapath)
os.unlink(hmmpath)
print('mined:', {k:(v[0],f'{v[2]:.1e}') for k,v in mined.items()})

with open(D+'cry4_subpanel.fasta','w') as fh:
    for sp,(rid,s) in sorted(out.items()): fh.write(f'>REFSEQ|{sp.replace(" ","_")}|{rid}\n{s}\n')
    for sp,(sid,s,ev) in sorted(mined.items()): fh.write(f'>MINED|{sp.replace(" ","_")}|{sid}\n{s}\n')
json.dump({'refseq':{k:v[0] for k,v in out.items()},'mined':{k:[v[0],v[2]] for k,v in mined.items()}},
          open(D+'cry4_subpanel_members.json','w'), indent=1)
print('total cry4 subpanel:', len(out)+len(mined))
