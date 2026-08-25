# Local patch to Nemo 2.4.2

`~/Nemo/nemo-release` carries one uncommitted change, needed for the two-phase
(burn-in → `source_pop`) design on macOS. Rebuild with:

```bash
cd ~/Nemo/nemo-release && make MAC_ARM=1 GSL_PATH=/opt/homebrew/opt/gsl/ -j8 && cp bin/nemo2.4.2-macARM ~/Nemo/
```

## `src/binarydataloader.cc` — chunk the population read

`BinaryDataLoader::extractPop` passed the full stored size to a single `read()`.
A stored metapopulation here is 18,432 individuals × 1e6 bitstring loci ≈ **4.6 GB**,
over the `INT_MAX` limit on a single `read()`.

* **Linux** silently caps the request at `0x7ffff000` and returns a short read, which
  the surrounding loop then continues — so this never appeared on the cluster.
* **macOS** fails the call outright with `EINVAL`, and Nemo reports it as
  `Binary file appears corrupted: ... read data Invalid argument (reading in
  4617591472 bytes, read -1 so far)` — misleading, since the file is fine.

The fix asks for at most one chunk per call and lets the existing loop do the rest.
The loop already handled partial reads, so nothing else changed.

Worth reporting upstream: any macOS user with a stored population over 2 GB hits
this, and the error message points at the wrong thing.
