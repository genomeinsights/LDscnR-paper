# `module_3sp/data/` — raw inputs, copied in

Two files, copied here so the pipeline is self-contained on this filesystem. Both were
verified byte-identical to their originals by SHA-256 at copy time.

**The originals are recorded in `PROVENANCE.csv` as history, not as a dependency.**
Nothing in the pipeline resolves a path outside this repository to find or check them.
That is deliberate: an earlier draft of `00_config.R` kept `origin_*` paths pointing at
`LD-scaling-genome-scans` so `01_inputs.R` could re-verify against them, which quietly
reintroduced the external dependency the copy existed to remove — the module would have
warned, or failed, whenever that repo was absent or moved.

`01_inputs.R` therefore checks each file against the SHA-256 written here. That catches
the thing worth catching — a local copy that has been corrupted or replaced — without
needing any other repository to exist.

| file | size | why it is not in git |
|---|---|---|
| `3sp_data.RData` | 159 MB | over GitHub's 100 MB hard limit; bzip2 only reaches 103 MB |
| `lfmm_F.rds` | 6.4 MB | trackable, but kept beside the file it belongs with |

`PROVENANCE.csv` and this README **are** tracked, so the repository carries proof of
which bytes were used without carrying the bytes.

`lfmm_F.rds` deserves one more line: it holds LFMM F-values over the full pre-MAF map,
dated September 2025, and **no script in any repository produces it**. Copying it here
makes it present; it does not make it reproducible. See `LFMM_SOURCE` in `00_config.R`.
