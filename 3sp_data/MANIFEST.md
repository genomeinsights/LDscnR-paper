# `3sp_data/` — external inputs, not in git

These files are the panel's analysis inputs. They are **not tracked** (4.5 GB;
`3sp_LDscnR_data.rds` alone is 846 MB). This manifest *is* tracked, so a fresh
clone can tell which files it is missing and whether a fetched copy is the right
one — the gap that made these inputs invisible to anyone who did not already
have them on disk.

## Verifying

```bash
./verify_manifest.sh            # check every file against the table
./verify_manifest.sh --write    # regenerate the table after a deliberate change
```

It follows symlinks, so it checks the bytes a script would actually read.

## Where the bundle lives

`3sp_LDscnR_data.rds` was moved out of the repository on 2026-09-03 and a symlink
left in its place, so the ~114 scripts across four repositories that hardcode
`~/gitlab/LDscnR-paper/3sp_data/3sp_LDscnR_data.rds` continue to resolve. The
bytes now live at

    ~/gitlab/LD-scaling-genome-scans/empirical_data/3sp/3sp_LDscnR_data.rds

alongside the rest of the bulk 3sp data, in a directory that repository ignores
(`.gitignore:37`). Verified byte-identical across the move.

On a machine without that layout, put the file anywhere and symlink it to the
path above, or edit the readers. `verify_manifest.sh` will confirm you have the
right file either way.

## Provenance

`3sp_LDscnR_data.rds` is produced by
`module_sticklebacks_LDscnR/regen_3sp_data.R`. The `3sp_latent_basis*.rds` come
from `superseded/build_latent_basis_3sp.R`, which writes
`3sp_latent_basis%s.rds` via `sprintf()` — the basis tag is a command-line
argument, so the four variants differ by an invocation nothing else records.
Only `3sp_latent_basis.rds` is read by anything today.

`rec_maps/` is excluded here: it is tracked in git and carries its own README and
build scripts.

## Files

wrote table for 18 files
| file | bytes | sha256 | stored at |
|---|---|---|---|
| `3sp.gds` | 27513174 | `cd7f4030bf909e3cf7b799f0db2c13ef39962928b1320dc117092694eddf9250` | here |
| `3sp_LDscnR_data.rds` | 887106787 | `5b6f984ac558e08d05eef3a866e616139792f5ab47d47097f8a78631856a0867` | ~/gitlab/LD-scaling-genome-scans/empirical_data/3sp/3sp_LDscnR_data.rds |
| `3sp_latent_basis.rds` | 4754 | `59bf1fd13fc2bf735e0da378d206c8d10b5dee01066285521bc30bfaf268162e` | here |
| `3sp_latent_basis_pruned.rds` | 4755 | `2c8bc01e542dec43274a3b3c1524379339f5d7494d7b74090c0fc424732b83a2` | here |
| `3sp_latent_basis_thin250.rds` | 4748 | `94dfaf45d06c9362fed8aadcf3ea2a92f52d83da3f7509b2b36bb716cb1abe34` | here |
| `3sp_latent_basis_thin50.rds` | 4749 | `46b278c77e0be65a7f4fbabf81f490a5bac2060e590b2e42950d92d9a8ccd2cf` | here |
| `data_for_emmax_perm.RData` | 148057953 | `2f781b8d6fa1e969136f4009af6ee6d2633a9616244997dcadaf2ffa97bf49c2` | here |
| `emmax_perm_reginal.rds` | 596138178 | `5d12bee227486975fd831c240c4bdf477ff3fddf19b25faa9b4f9d4c69e81b39` | here |
| `grm_null.rds` | 183815 | `8964f675e72a9ec43d3df377d5cedf38e72a8fa715bbffe0fffc47bf0d21d36b` | here |
| `grm_phenotypes.rds` | 1287845 | `89d4ce4a2c81bd7b114feacedafdc7af41b40965c1084128f77f9fc423a0d3d3` | here |
| `ld_decay_3sp.rds` | 333548968 | `dd7ef7695237dbef3ed8503dd8b302e453b694b8186f14a2f7ddd992dc3dc7e1` | here |
| `ld_decay_3sp_MAF01.rds` | 318274812 | `2430fcfa82e02b682c754f3bc5be8b95e45b8718b7aec482756eb6d9290cd777` | here |
| `ld_ws_3sp_MAF01.rds` | 87167166 | `5cc27888979cd4aec4619fd02aa9d2f74a6a21fe64995be4172e77030676c728` | here |
| `out_emx_global.rds` | 596127758 | `4911a99cafb448980f83e26acfa0d59f91e0c8511fdea18948e15ebb04d073a6` | here |
| `out_emx_grm.rds` | 596121337 | `fca79f6dbd50114dbfbf1ed7c785e6c5d879e244694c422a93496b301d3a48c7` | here |
| `out_emx_grm_0.rds` | 596002904 | `50389c734fa1926f930d7d09e5a61da335b6ff9de06aef9348e57856c3ece49f` | here |
| `out_emx_grm_05.rds` | 596111965 | `bf2d43d9e790170801cc93f51387b42e21810536c6c90a24ea0b50319d18e438` | here |
| `run_3sp_perm_global.R` | 1366 | `97396bfa4a40adab895e754bdec03b91714c3342d41ae2e22b5d6f91002339bf` | here |
