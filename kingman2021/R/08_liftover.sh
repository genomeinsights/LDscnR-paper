#!/usr/bin/env bash
## Build the gasAcu1 <-> gasAcu1-4 liftOver chains and lift everything that has to
## cross assemblies.
##
## Why this exists: the 3sp dataset is in **gasAcu1** (Chr1..Chr21, arabic) while every
## Kingman product is in **gasAcu1-4** (chrI..chrXXI, roman). The published .chain files
## live on Dryad behind an anti-bot gate, so the chain is reconstructed here from the
## hub's bigChain tracks on FigShare, then reversed with UCSC chainSwap.
##
## Needs on PATH: bigBedToBed, liftOver, chainSwap  (http://hgdownload.soe.ucsc.edu/admin/exe/)
set -euo pipefail
OUT=$(cd "$(dirname "$0")/../data/liftover" && pwd)
TR=~/gitlab/LD-scaling-genome-scans/empirical_data/kingman2021/tracks
cd "$OUT"

## ---- 1. fetch the bigChain pair and rebuild a .chain --------------------------------
[[ -s g14ToG1.bigChain.bb      ]] || curl -sSL -o g14ToG1.bigChain.bb      https://ndownloader.figshare.com/files/39711067
[[ -s g14ToG1.bigChain.link.bb ]] || curl -sSL -o g14ToG1.bigChain.link.bb https://ndownloader.figshare.com/files/39711046
bigBedToBed g14ToG1.bigChain.bb      chain.bed
bigBedToBed g14ToG1.bigChain.link.bb link.bed

## bigChain: chrom start end id score strand tSize qName qSize qStart qEnd chainScore
## bigLink : chrom start end id qStart
python3 - <<'PY'
import collections
chains={}
for l in open("chain.bed"):
    f=l.rstrip("\n").split("\t")
    chains[f[3]]=dict(t=f[0],tStart=int(f[1]),tEnd=int(f[2]),strand=f[5],tSize=int(f[6]),
                      q=f[7],qSize=int(f[8]),qStart=int(f[9]),qEnd=int(f[10]),chainScore=f[11])
links=collections.defaultdict(list)
for l in open("link.bed"):
    f=l.rstrip("\n").split("\t"); links[f[3]].append((int(f[1]),int(f[2]),int(f[4])))
n=nb=0
with open("gasAcu1-4ToGasAcu1.chain","w") as out:
    for cid,c in chains.items():
        L=sorted(links[cid])
        if not L: continue
        out.write("chain %s %s %d + %d %d %s %d %s %d %d %s\n" % (
            c["chainScore"],c["t"],c["tSize"],c["tStart"],c["tEnd"],
            c["q"],c["qSize"],c["strand"],c["qStart"],c["qEnd"],cid))
        for i,(ts,te,qs) in enumerate(L):
            size=te-ts
            if i+1<len(L):
                out.write("%d\t%d\t%d\n"%(size,L[i+1][0]-te,L[i+1][2]-(qs+size)))
            else:
                out.write("%d\n"%size)
            nb+=1
        out.write("\n"); n+=1
print(f"reconstructed {n} chains / {nb} blocks")
PY
chainSwap gasAcu1-4ToGasAcu1.chain gasAcu1ToGasAcu1-4.chain

## ---- 2. lift the EcoPeak/TempoPeak sets into gasAcu1, keeping the published p-values --
for s in c155.specific c155.sensitive c150.specific c150.sensitive; do
  awk -F'\t' -v OFS='\t' '{n=split($4,a,"="); ps=(n>1?a[n]:"NA");
                           m=split($5,b,"="); pw=(m>1?b[m]:"NA"); print $1,$2,$3,ps"|"pw}' \
      "$TR/gasAcu1-4.$s.50kb.final.peaks.bed" > pin.bed
  liftOver -minMatch=0.5 pin.bed gasAcu1-4ToGasAcu1.chain "pv_$s.bed" "pv_un_$s.bed"
  echo "$s: $(wc -l < pin.bed) -> $(wc -l < pv_$s.bed)"
done
rm -f pin.bed

## ---- 3. lift the 3sp LFMM outlier regions the other way, into gasAcu1-4 ---------------
python3 - <<'PY'
import csv, os
ROMAN=["I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV",
       "XVI","XVII","XVIII","XIX","XX","XXI"]
src=os.path.expanduser("~/gitlab/LDscnR-paper/kingman2021/data/regions_tau0.05_lmin10_rho0.60.csv")
with open("lfmm_g1.bed","w") as f:
    for r in csv.DictReader(open(src)):
        if r["method"]!="LFMM": continue
        n=int(r["Chr"].replace("Chr",""))
        f.write(f"chr{ROMAN[n-1]}\t{int(r['start'])-1}\t{r['end']}\tLFMM_{r['region']}\n")
PY
liftOver -minMatch=0.5 lfmm_g1.bed gasAcu1ToGasAcu1-4.chain lfmm_g14.bed lfmm_un.bed
echo "3sp LFMM regions: $(wc -l < lfmm_g1.bed) -> $(wc -l < lfmm_g14.bed)"
