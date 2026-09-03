#!/usr/bin/env python3
"""Build a cross-based (pedigree) recombination map for threespine stickleback
from Glazer et al. 2015 G3 (doi:10.1534/g3.115.017905) supplementary File S1-S3.

Outputs, in BOTH gasAcu1 ("old") and Glazer revised gasAcu1-4 ("new") coordinates:
  glazer2015_markers.tsv    marker-level: cM in each of the two F2 crosses + bin coords
  glazer2015_rate_bins.tsv  per-bin local recombination rate (cM/Mb), per cross + mean
"""
import csv, hashlib, os, sys, openpyxl
from collections import defaultdict

# Preconditions on the published inputs. These are checked rather than documented
# because the natural way to re-fetch them FAILS SILENTLY: Wiley (and to a lesser
# extent the PMC file endpoint) sit behind Cloudflare, which answers a scripted
# curl/wget with HTTP 200 and a ~6 KB HTML interstitial. That lands on disk with a
# plausible name and parses as "something", so a rebuild can half-succeed and
# produce a map that looks fine. See README.md, "Re-fetching the source files".
SOURCES = {
    "source/supp_g3.115.017905_FileS1.xlsx":
        (1342078, "f8bd177572bf8c2b8388455357bc4f5faa05f35d1fe92b6180b6684a983888bb"),
    "source/supp_g3.115.017905_FileS2.xlsx":
        (1857861, "c695fca8dbe02dda1f430f0c5b390d59b18ae6c0d60562fe6eeb258a97aef07a"),
    "source/supp_g3.115.017905_FileS3.xlsx":
        (244412,  "d3cc2fb083fd1efcb068ad65394b4be57f4e2c6d627e9b03193cabbfde50e273"),
}


def check_sources():
    """Refuse to build on inputs that are absent, truncated, or not what they claim."""
    for path, (want_size, want_sha) in sorted(SOURCES.items()):
        if not os.path.exists(path):
            sys.exit(f"MISSING: {path}\n"
                     f"  source/ is gitignored, so a fresh clone will not have it.\n"
                     f"  README.md gives the DOI and says why a scripted fetch will not work.")
        size = os.path.getsize(path)
        with open(path, "rb") as fh:
            got = hashlib.sha256(fh.read()).hexdigest()
        if got != want_sha:
            hint = ("  That size is in the range of a Cloudflare interstitial -- you have "
                    "an HTML\n  block page, not the file. Fetch it through a browser."
                    if size < 50_000 else
                    "  Size matches but the hash does not: the file was revised upstream, or "
                    "edited\n  locally. Do not build on it silently.")
            sys.exit(f"WRONG CONTENT: {path}\n"
                     f"  expected {want_size:,} bytes  sha256 {want_sha[:16]}...\n"
                     f"  got      {size:,} bytes  sha256 {got[:16]}...\n{hint}")
    print(f"source files verified: {len(SOURCES)} of {len(SOURCES)}")


check_sources()

def genetic_map(path):
    wb = openpyxl.load_workbook(path, read_only=True)
    ws = wb["Genetic Map"]
    out = {}
    for i, row in enumerate(ws.iter_rows(values_only=True)):
        if i == 0 or row is None or row[0] is None or row[1] is None:
            continue
        out[str(row[1]).strip()] = (str(row[0]).strip(), float(row[2]))
    return out

ftc  = genetic_map("source/supp_g3.115.017905_FileS1.xlsx")   # Fishtrap Ck (FW) x Little Campbell R (marine)
bepa = genetic_map("source/supp_g3.115.017905_FileS2.xlsx")   # Bear Paw Lake (FW) x Little Campbell R (marine)

wb = openpyxl.load_workbook("source/supp_g3.115.017905_FileS3.xlsx", read_only=True)
ws = wb.worksheets[0]
bins = {}
for i, row in enumerate(ws.iter_rows(values_only=True)):
    if i == 0 or row is None or row[0] is None:
        continue
    m = str(row[0]).strip()
    bins[m] = dict(scaffold=row[1],
                   oldChr=str(row[2]).strip(), oldStart=int(row[3]), oldEnd=int(row[4]),
                   binLength=int(row[5]),
                   newChr=str(row[6]).strip(), newStart=int(row[7]), newEnd=int(row[8]),
                   newOrientation=row[9])

print(f"FTC map markers  : {len(ftc)}")
print(f"BEPA map markers : {len(bepa)}")
print(f"File S3 bins     : {len(bins)}")

# ---- marker table -----------------------------------------------------------
rows = []
for m, b in bins.items():
    f = ftc.get(m); p = bepa.get(m)
    if f is None and p is None:
        continue
    chrom = (f or p)[0]
    rows.append(dict(marker=m, chr=chrom, scaffold=b["scaffold"],
                     cM_FTC=f[1] if f else "", cM_BEPA=p[1] if p else "",
                     old_chr=b["oldChr"], old_start=b["oldStart"], old_end=b["oldEnd"],
                     old_mid=(b["oldStart"] + b["oldEnd"]) // 2,
                     new_chr=b["newChr"], new_start=min(b["newStart"], b["newEnd"]),
                     new_end=max(b["newStart"], b["newEnd"]),
                     new_mid=(b["newStart"] + b["newEnd"]) // 2,
                     bin_length=b["binLength"]))
print(f"markers with map position and coordinates: {len(rows)}")

def chrkey(c):
    try: return (0, int(c))
    except ValueError: return (1, c)

rows.sort(key=lambda r: (chrkey(r["chr"]), r["new_mid"]))
with open("glazer2015_markers.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(rows)

# ---- per-bin local rate -----------------------------------------------------
# For each cross: order that cross's markers along the REVISED assembly (the one whose
# scaffold order the map itself validated), then take the local slope of cM vs Mb over
# the two flanking markers. Rate is reported against BOTH coordinate systems, so a SNP
# in gasAcu1 coordinates can be looked up directly via old_chr/old_start/old_end.
def local_rates(mapdict, key="new_mid"):
    by_chr = defaultdict(list)
    for r in rows:
        if r["marker"] in mapdict:
            by_chr[r["chr"]].append(r)
    rate = {}
    tot_cM = tot_bp = 0.0
    for c, rs in by_chr.items():
        rs = sorted(rs, key=lambda r: r[key])
        cm = [mapdict[r["marker"]][1] for r in rs]
        bp = [r[key] for r in rs]
        for i, r in enumerate(rs):
            lo, hi = max(0, i - 1), min(len(rs) - 1, i + 1)
            dcm, dbp = abs(cm[hi] - cm[lo]), abs(bp[hi] - bp[lo])
            rate[r["marker"]] = round(dcm / (dbp / 1e6), 4) if dbp > 0 else ""
        tot_cM += max(cm) - min(cm)
        tot_bp += max(bp) - min(bp)
    return rate, tot_cM, tot_bp

r_ftc,  L_ftc,  bp_ftc  = local_rates(ftc)
r_bepa, L_bepa, bp_bepa = local_rates(bepa)
print(f"FTC  total map length {L_ftc:8.1f} cM over {bp_ftc/1e6:6.1f} Mb -> {L_ftc/(bp_ftc/1e6):.2f} cM/Mb")
print(f"BEPA total map length {L_bepa:8.1f} cM over {bp_bepa/1e6:6.1f} Mb -> {L_bepa/(bp_bepa/1e6):.2f} cM/Mb")

out = []
for r in rows:
    a, b = r_ftc.get(r["marker"], ""), r_bepa.get(r["marker"], "")
    vals = [v for v in (a, b) if v != ""]
    out.append({**{k: r[k] for k in ("marker", "chr", "old_chr", "old_start", "old_end",
                                     "new_chr", "new_start", "new_end", "bin_length")},
                "cM_FTC": r["cM_FTC"], "cM_BEPA": r["cM_BEPA"],
                "rate_FTC_cMperMb": a, "rate_BEPA_cMperMb": b,
                "rate_mean_cMperMb": round(sum(vals) / len(vals), 4) if vals else ""})
with open("glazer2015_rate_bins.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(out)

both = [o for o in out if o["rate_FTC_cMperMb"] != "" and o["rate_BEPA_cMperMb"] != ""]
print(f"bins with a rate in both crosses: {len(both)}")
import statistics
xs = [o["rate_FTC_cMperMb"] for o in both]; ys = [o["rate_BEPA_cMperMb"] for o in both]
mx, my = statistics.mean(xs), statistics.mean(ys)
num = sum((x-mx)*(y-my) for x, y in zip(xs, ys))
den = (sum((x-mx)**2 for x in xs) * sum((y-my)**2 for y in ys)) ** 0.5
print(f"cross-cross correlation of local rate: r = {num/den:.3f}")
print(f"bins anchored to a gasAcu1 chromosome (old_chr != Un): "
      f"{sum(1 for o in out if o['old_chr'] not in ('Un','None'))}")
