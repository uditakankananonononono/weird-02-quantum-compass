import json, subprocess
from collections import Counter, defaultdict
from Bio import Phylo
from pyfamsa import Aligner, Sequence
import os
os.chdir('/home/sandbox/weird10/projects/W02')

meta = json.load(open('data/cry_aves_species.json'))
ref = {}
for a, v in meta.items():
    t = v['title'].lower()
    if 'cryptochrome 4' in t: ref[a] = 'Cry4'
    elif 'cryptochrome 1' in t and 'like' not in t: ref[a] = 'Cry1'
    elif 'cryptochrome 2' in t and 'like' not in t: ref[a] = 'Cry2'
    elif 'dash' in t: ref[a] = 'DASH'
    elif 'photolyase' in t: ref[a] = 'photolyase'
print('reference seeds:', Counter(ref.values()), flush=True)

seqs = {}
name, buf = None, []
for line in open('data/cry_aves_proteins.fasta'):
    line = line.rstrip()
    if line.startswith('>'):
        if name: seqs[name] = ''.join(buf)
        name, buf = line[1:], []
    else: buf.append(line)
if name: seqs[name] = ''.join(buf)
keep = {a: s for a, s in seqs.items() if len(s) >= 300}
print('seqs >=300aa:', len(keep), flush=True)

aln = Aligner().align([Sequence(a.encode(), s.encode()) for a, s in keep.items()])
with open('/tmp/w02_cry.aln','w') as fh:
    for s in aln:
        nm = s.id.decode() if isinstance(s.id, bytes) else s.id
        sq = s.sequence.decode() if isinstance(s.sequence, bytes) else s.sequence
        fh.write(f'>{nm}\n{sq}\n')
print('aligned', flush=True)
subprocess.run(['/home/sandbox/bin/FastTree','-wag','/tmp/w02_cry.aln'], stdout=open('/tmp/w02_cry.nwk','w'), stderr=subprocess.DEVNULL)
print('tree done', flush=True)

tree = Phylo.read('/tmp/w02_cry.nwk','newick')
depth, paths = {}, {}
def walk(c, d, p):
    depth[c] = d; paths[c] = p
    for x in c.clades: walk(x, d + (x.branch_length or 0.0), p + [c])
walk(tree.root, 0.0, [])
tips = tree.get_terminals()
ref_tips = [t for t in tips if t.name.split('.')[0] in ref]
def lca(a,b):
    pa,pb = paths[a],paths[b]; i=0
    while i < min(len(pa),len(pb)) and pa[i] is pb[i]: i+=1
    return pa[i-1] if i else tree.root
def dist(a,b):
    c = lca(a,b); return depth[a]+depth[b]-2*depth[c]
pred = {}
for t in tips:
    acc = t.name.split('.')[0]
    ds = sorted(((dist(t,rt), rt) for rt in ref_tips), key=lambda x:x[0])[:3]
    votes = defaultdict(float)
    for d, rt in ds:
        votes[ref[rt.name.split('.')[0]]] += 1.0/(d+1e-6)
    pred[acc] = max(votes.items(), key=lambda x:x[1])[0]
print('classified paralogs:', Counter(pred.values()), flush=True)
json.dump(pred, open('data/paralog_labels.json','w'), indent=1)
sp4 = sorted({meta[a]['organism'] for a,v in pred.items() if v=='Cry4' and a in meta})
print('species with classified Cry4:', len(sp4), flush=True)
json.dump(sp4, open('data/species_with_cry4.json','w'), indent=1)
