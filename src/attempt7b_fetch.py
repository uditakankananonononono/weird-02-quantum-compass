#!/usr/bin/env python3
"""W02 7b: fetch control-protein orthologs (PKM, LDHA, RHO, OPN4, GAPDH) for the
avian panel via NCBI eutils. Resumable: skips species+gene already fetched.
Rule (locked 21:07 IST): esearch <species>[orgn] AND <gene>[gene] AND refseq[filter],
longest protein, one per species; misses logged, no substitution."""
import json, os, time, urllib.parse, urllib.request

BASE = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
OUT = 'data/controls7b'
os.makedirs(OUT, exist_ok=True)
spp = sorted({v['organism'] for v in json.load(open('data/cry_aves_species.json')).values()})
genes = ['PKM','LDHA','RHO','OPN4','GAPDH']
ledger_p = 'data/attempt7b_controls_ledger.json'
ledger = json.load(open(ledger_p)) if os.path.exists(ledger_p) else {}

def get(url):
    for t in range(3):
        try:
            with urllib.request.urlopen(url, timeout=30) as r: return r.read().decode()
        except Exception as e:
            time.sleep(2*t+1)
    return None

done = 0
for sp in spp:
    for g in genes:
        key = f"{sp}|{g}"
        if key in ledger and ledger[key].get('status') in ('ok','miss'): continue
        q = urllib.parse.quote(f'"{sp}"[orgn] AND {g}[gene] AND refseq[filter]')
        r = get(BASE+f'esearch.fcgi?db=protein&term={q}&retmax=5&retmode=json')
        ids = json.loads(r)['esearchresult'].get('idlist',[]) if r else []
        best = None
        if ids:
            f = get(BASE+f'efetch.fcgi?db=protein&id={",".join(ids)}&rettype=fasta&retmode=text')
            if f:
                recs = [x for x in f.split('>') if x.strip()]
                if recs:
                    best = max(recs, key=len)
        if best:
            hdr, seq = best.split('\n',1)
            acc = hdr.split()[0]
            fn = f"{OUT}/{sp.replace(' ','_')}__{g}.fa"
            open(fn,'w').write('>'+hdr.rstrip()+'\n'+''.join(seq.split())+'\n')
            ledger[key] = {'status':'ok','acc':acc,'len':len(''.join(seq.split())),'file':fn}
            done += 1
        else:
            ledger[key] = {'status':'miss'}
        time.sleep(0.37)
    json.dump(ledger, open(ledger_p,'w'), indent=1)
    print(sp, 'cumulative ok:', sum(1 for v in ledger.values() if v['status']=='ok'), flush=True)
print('DONE, fetched this run:', done)
