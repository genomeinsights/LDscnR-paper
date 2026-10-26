# `module_9sp/data/` — raw inputs, copied in

Two files, verified by SHA-256 against `PROVENANCE.csv` by `R/01_inputs.R`. Neither is tracked:
`9sp_data.rds` (290 MB) is over GitHub's 100 MB limit and `lfmm_F.rds` (13 MB) is kept beside it.

| file | contents |
|---|---|
| `9sp_data.rds` | list: `GT` (149 × 1,355,325 dosages), `map`, `GRM`, `pheno` (population, lineage, locality), `ecotype_bin` |
| `lfmm_F.rds` | LFMM F-values over the full pre-MAF map; producing code preserved in `R/legacy_source.R` |

## Getting the files

Both files are archived, byte for byte, in the LD-scaling Zenodo deposit
(doi:[10.5281/zenodo.18466057](https://doi.org/10.5281/zenodo.18466057), `9sp.zip`, CC-BY-4.0).
They derive from the nine-spined data of Fang et al. (2021, *Mol Biol Evol* 38:4205,
doi:10.1093/molbev/msab144; raw data doi:10.5281/zenodo.4722879). From the repository root:

```bash
Rscript R/00_fetch_inputs.R --install
```

downloads the archives, checks them against Zenodo's MD5, extracts the files to `rebuild_check/`,
checks each against the SHA-256 in `PROVENANCE.csv`, and copies verified files here (an existing
file is never overwritten). Use `--zip-dir=PATH` for local copies of `3sp.zip`/`9sp.zip`.
