#!/usr/bin/env python3
"""W02 7e: fetch Cry2 protein orthologs for the 116-species avian panel (NCBI,
resumable, sha256 ledger; misses logged, no substitution). Amendment 03:08 IST."""
import json, subprocess, time, hashlib, os
os.makedirs('data/cry2_panel', exist_ok=True)
LEDGER = 'data/cry2_panel_ledger.json'
led = json.load(open(LEDGER)) if os.path.exists(LEDGER) else {}
species = [l.strip() for l in open('/tmp/w02_panel_species.txt') if l.strip()]
todo = [s for s in species if s not in led]
print(f'{len(todo)} to fetch, {len(led)} done', flush=True)
def run(cmd):
    return subprocess.run(cmd, capture_output=True, timeout=40)
n = 0
for sp in todo:
    q = sp.replace('_', '+')
    ok = False
    for attempt in range(3):
        try:
            r = run(['curl','-sf','--max-time','25',
                f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=protein&term=cry2%5Bgene%5D+AND+{q}%5Borganism%5D&retmax=3&retmode=json'])
            if r.returncode == 0:
                ids = json.loads(r.stdout)['esearchresult']['idlist']
                if ids:
                    f = run(['curl','-sf','--max-time','25',
                        f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id={ids[0]}&rettype=fasta&retmode=text'])
                    if f.returncode == 0 and f.stdout.startswith(b'>'):
                        seq = f.stdout
                        acc = seq.split(b'\n')[0].split(b' ')[0][1:].decode()
                        open(f'data/cry2_panel/{sp}.fa','wb').write(seq)
                        led[sp] = {'status':'ok','accession':acc,'sha256':hashlib.sha256(seq).hexdigest(),'bytes':len(seq)}
                        ok = True
                        break
        except Exception:
            pass
        time.sleep(2 + attempt*2)
    if not ok:
        led[sp] = {'status':'miss'}
    n += 1
    if n % 15 == 0:
        json.dump(led, open(LEDGER,'w'), indent=1)
        print(f'  {n} processed', flush=True)
    time.sleep(1.0)
json.dump(led, open(LEDGER,'w'), indent=1)
okn = sum(1 for v in led.values() if v['status']=='ok')
print(f'DONE: {okn} ok / {len(led)} species')
