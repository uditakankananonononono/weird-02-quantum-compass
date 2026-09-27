#!/usr/bin/env python3
"""W02 queue #11 (13:12 IST amendment + 13:15 IST addendum 37010e4): episodic-selection screen
on avian Cry4, 10-species 5v5 (Pavo cristatus excluded - no committed label exists anywhere).
ARM1 CDS ledger (sha256) -> ARM2 codon back-map onto committed Cry4 subpanel protein MSA ->
ARM3 pairwise Nei-Gojobori dN/dS, MWU cross-class vs within-class. Report-only screen."""
import json, hashlib, numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.codonalign import codonseq
from scipy.stats import mannwhitneyu

SPP=[("REFSEQ","Columba_livia","KX168611","sedentary"),("REFSEQ","Erithacus_rubecula","MN709784","migratory"),
("REFSEQ","Gallus_gallus","NM_001039596","sedentary"),("REFSEQ","Passer_domesticus","AY494987","sedentary"),
("MINED","Catharus_ustulatus","XP_032936297.1","migratory"),("MINED","Ficedula_albicollis","ENSFALP00000002532.1","migratory"),
("MINED","Hirundo_rustica","XP_039941661.1","migratory"),("MINED","Serinus_canaria","ENSSCAP00000009299.1","sedentary"),
("MINED","Struthio_camelus","XP_068774195.1","sedentary"),("MINED","Zonotrichia_albicollis","ENSZALP00000018676.1","migratory")]

def load_msa(p):
    m={}
    for l in open(p):
        if l.startswith('>'): k=l[1:].strip(); m[k]=''
        else: m[k]+=l.strip()
    return m
cry4_msa=load_msa('data/mafft_ebi_crosscheck.fasta')   # committed Cry4 subpanel MSA (accession-keyed)
cry1_msa=load_msa('data/cry1_panel_msa.fasta')          # for the identity cross-check only

# ARM 1: ledger + CDS load
ledger={}
cds={}
for src,sp,acc,lab in SPP:
    b=open(f'data/cds/{acc}.fna','rb').read()
    rec=SeqIO.read(f'data/cds/{acc}.fna','fasta')
    s=str(rec.seq).upper()
    ledger[acc]={"species":sp,"source_db":src,"sha256":hashlib.sha256(b).hexdigest(),"n_nt":len(s)}
    cds[sp]=(s,lab)
json.dump({"arm":"ARM1 CDS ledger","accessions":ledger,"failed":[],"n_ok":len(ledger)},
          open('results/attempt8h_cds_ledger.json','w'), indent=1)

def mask_ambiguous(cds_seq):
    """KX168611.fna carries 2 IUPAC ambiguity codes (R,Y). Mask non-ACGT bases to N;
    codons containing N become '---' after back-mapping so pairwise NG skips them.
    Disclosed: affects Columba_livia (KX168611) only, 2 codons."""
    return "".join(b if b in "ACGT" else "N" for b in cds_seq.upper())
for sp in cds: cds[sp]=(mask_ambiguous(cds[sp][0]), cds[sp][1])

# identity cross-check: CDS translation vs both MSAs (justifies back-map target)
ident={}
for src,sp,acc,lab in SPP:
    s,_=cds[sp]
    tr=str(Seq(s[:len(s)-len(s)%3]).translate(to_stop=False))
    row4=cry4_msa[acc].replace('-','')
    key1=[h for h in cry1_msa if sp in h][0]
    row1=cry1_msa[key1].replace('-','')
    def idfrac(a,b):
        n=min(len(a),len(b)); return sum(x==y for x,y in zip(a[:n],b[:n]))/max(len(a),len(b))
    ident[sp]={"vs_cry4_subpanel_msa":round(idfrac(tr,row4),4),"vs_cry1_panel_row":round(idfrac(tr,row1),4),
               "cds_aa":len(tr),"msa_row_aa":len(row4)}
    assert ident[sp]["vs_cry4_subpanel_msa"]>0.95, f"{sp} CDS does not match Cry4 subpanel MSA"

# ARM 2: back-map CDS onto the committed Cry4 subpanel protein MSA rows (10 kept species)
aln={}
for src,sp,acc,lab in SPP:
    row=cry4_msa[acc]; s,_=cds[sp]; out=[]; i=0
    for a in row:
        if a=='-': out.append('---')
        else:
            codon=s[i:i+3]; i+=3
            out.append(codon if 'N' not in codon else '---')
    aln[sp]=''.join(out)
L=len(aln[SPP[0][1]])
assert all(len(v)==L for v in aln.values())
labels={sp:lab for _,sp,_,lab in SPP}

# ARM 3: pairwise NG dN/dS over all 45 pairs (Biopython codonseq, method NG)
sp10=[sp for _,sp,_,_ in SPP]
pairs=[]
for i in range(len(sp10)):
    for j in range(i+1,len(sp10)):
        a,b=sp10[i],sp10[j]
        try:
            sa=codonseq.CodonSeq(aln[a]); sb=codonseq.CodonSeq(aln[b])
            dN,dS=codonseq.cal_dn_ds(sa, sb, method='NG86')
            dN=float(dN); dS=float(dS)
            if dN!=dN or dS!=dS or dS<=0: w=None
            else: w=dN/dS
        except Exception as e:
            dN=dS=w=None
        cls='cross' if labels[a]!=labels[b] else 'within'
        pairs.append({"a":a,"b":b,"class":cls,"la":labels[a],"lb":labels[b],
                      "dN":None if dN is None else round(dN,6),"dS":None if dS is None else round(dS,6),
                      "omega":None if w is None else round(w,6)})
cross=[p["omega"] for p in pairs if p["class"]=='cross' and p["omega"] is not None]
within=[p["omega"] for p in pairs if p["class"]=='within' and p["omega"] is not None]
U,p2=mannwhitneyu(cross,within,alternative='two-sided')
out={"amendment":"2026-09-27 13:12 IST queue #11 + 13:15 IST addendum 37010e4 (10-species 5v5; delegated executor, report-only)",
 "executor_note":"computed by delegated executor at ref 8c424cb; uncommitted, returned to lane-29 for review",
 "pavo_exclusion":"Pavo cristatus excluded: no committed label exists anywhere in the repo (h1_scoring.json, n_partial=0) and no protein-MSA row in cry1_panel_msa.fasta; exclusion disclosed verbatim per addendum",
 "tool_availability_verbatim":"which yn00 codeml paml hyphy HYPHYMP -> all empty; PAML/HyPhy NOT present. Locked fallback used: pairwise Nei-Gojobori (Biopython 1.88 Bio.codonalign codonseq.cal_dn_ds method='NG') + Mann-Whitney. Disclosed heuristic screen, not a branch-site test.",
 "arm1":{"accessions_ok":len(ledger),"failed":0,"ledger_file":"results/attempt8h_cds_ledger.json",
   "routes":"4 REFSEQ nuccore efetch fasta_cds_na (KX168611, MN709784, NM_001039596, AY494987); 3 XP protein efetch fasta_cds_na (XP_032936297.1, XP_039941661.1, XP_068774195.1); 3 ENS via Ensembl REST: lookup/id on Translation id to resolve Parent transcript, then sequence/id(transcript)?type=cds (direct sequence/id on the Translation id silently returns PROTEIN - alphabet-verified during retrieval)"},
 "arm2":{"backmap_target":"data/mafft_ebi_crosscheck.fasta (committed Cry4 subpanel protein MSA, accession-keyed)",
   "deviation_disclosure":"Amendment text named cry1_panel_msa.fasta; identity check (below) shows CDS translations match the Cry4 subpanel MSA (>0.95) and NOT the Cry1 panel rows (different accessions, Cry1-family proteins). The committed Cry4 subpanel MSA was used; executing against Cry1 rows would back-map Cry4 CDS onto Cry1 protein (invalid).",
   "identity_check":ident,"alignment_columns":L},
 "arm3":{"n_pairs":len(pairs),"n_cross":len(cross),"n_within":len(within),
   "omega_cross_median":round(float(np.median(cross)),6),"omega_within_median":round(float(np.median(within)),6),
   "omega_cross_mean":round(float(np.mean(cross)),6),"omega_within_mean":round(float(np.mean(within)),6),
   "mwu_U":float(U),"mwu_p_two_sided":float(p2),
   "mwu_power_note":"5v5 classes; minimum attainable two-sided MWU p at 25-vs-20 pairs is limited (lane-29: 0.00794 at label level) - weakly powered, reported as such",
   "pairs":pairs},
 "interpretation_lock":"REPORT-ONLY screen. A significant elevation is a SCREEN signal for discussion, never an adaptive-evolution claim; no causal residues, no migratory-prediction claim."}
json.dump(out, open('results/attempt8h_dnds.json','w'), indent=1)
print('identity ok; pairs:',len(pairs),'cross:',len(cross),'within:',len(within))
print('omega cross med',out['arm3']['omega_cross_median'],'within med',out['arm3']['omega_within_median'])
print('MWU U',U,'p',p2)
