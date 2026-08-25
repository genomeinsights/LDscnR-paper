## =====================================================================
## module_sticklebacks_LDscnR / build_latent_basis_3sp.R
##
## The SHARED latent structure basis for the symmetric 2x2 null design
## ({EMMAX, LFMM} x {K-MVN, latent}). LFMM corrects structure with a hard top-K
## latent-factor truncation, so its "home-field" specificity null draws phenotypes
## from that same subspace -- the exact analogue of the K-MVN null that is EMMAX's
## home field. Both engines must see the IDENTICAL basis for the cross to be fair,
## so we compute it ONCE here and both the EMMAX builder (emmax_latent_null_3sp.R)
## and the LFMM cluster runner read this file.
##
## Basis = top-K principal components of the centered/scaled genotype matrix
## (an unsupervised, phenotype-INDEPENDENT estimate of LFMM's latent space -- using
## the fitted, env-conditional U would be circular). Computed from the n x n Gram
## matrix so no full SVD of the 117 x 790k matrix is needed. K matches LFMM's K=5.
##
## The latent null phenotype is then  y = U %*% (d * rnorm(K)),  i.e. a draw from
## MVN(0, U diag(d^2) U^T) -- the top-K projection of the genotype covariance --
## residualised on the observed ecotype, exactly as structured_null() does for K-MVN.
##
## Run from the LDscnR-paper root:
##   Rscript module_sticklebacks_LDscnR/build_latent_basis_3sp.R
## Writes module_sticklebacks_LDscnR/data/3sp_latent_basis.rds  (U [n x K], d [K], K)
## =====================================================================
suppressMessages({ library(data.table) })

## Optional arg "pruned": build the basis from the LD-independent GRM markers only,
## so the top-K PCs reflect NEUTRAL structure rather than the adaptive LD (large
## inversions / Eda) that dominates PCs computed from all SNPs. This makes the latent
## null the clean PC analogue of the LD-pruned GRM (EMMAX's kinship).
## Modes (arg 1):
##   (none)            -- all SNPs (faithful to what LFMM's latent factors actually see)
##   pruned            -- the LD-independent GRM markers (w_j threshold); ~94% here, ~= raw
##   thin=<kb>         -- keep 1 SNP per <kb> physical window: collapses LD blocks
##                        (inversions/Eda) so the top-K PCs reflect NEUTRAL structure,
##                        not the adaptive LD that otherwise dominates them
##   noregions         -- drop the discovered candidate regions + Chr1 inversion, then PCA
a      <- commandArgs(trailingOnly = TRUE)
mode   <- if (length(a) >= 1) a[1] else "all"
suppressMessages(library(data.table))
BND <- "module_sticklebacks_LDscnR/data/3sp_LDscnR_data.rds"
tag <- if (mode == "all") "" else paste0("_", gsub("[^A-Za-z0-9]", "", mode))
OUT <- sprintf("module_sticklebacks_LDscnR/data/3sp_latent_basis%s.rds", tag)
K   <- 5L                                        # == K_LFMM

d <- readRDS(BND); G <- d$GTs; n <- nrow(G); map <- as.data.table(d$map)
if (mode == "pruned") {
  keep <- colnames(G) %in% d$grm_markers
} else if (grepl("^thin=", mode)) {
  win <- as.numeric(sub("thin=", "", mode)) * 1000            # window in bp
  mk  <- map[match(colnames(G), marker)]
  mk[, w := floor(Pos / win)]
  sel <- mk[, .I[1], by = .(Chr, w)]$V1                        # first SNP per (chr,window)
  keep <- seq_len(ncol(G)) %in% sel
} else {                                                        # "all" (default)
  keep <- rep(TRUE, ncol(G))
}
cat(sprintf("mode=%s: using %d of %d SNPs (%.1f%%)\n", mode, sum(keep), ncol(G), 100*mean(keep)))
G <- G[, keep, drop = FALSE]
stopifnot(!anyNA(G))
Gs   <- scale(G)                                 # center + unit-variance per SNP (Patterson-style)
Gs[!is.finite(Gs)] <- 0                          # guard any constant column (none expected)
Gram <- tcrossprod(Gs) / (ncol(Gs) - 1L)         # n x n genotype covariance among individuals
eg   <- eigen(Gram, symmetric = TRUE)
U    <- eg$vectors[, seq_len(K), drop = FALSE]    # top-K PCs (n x K)
dval <- sqrt(pmax(eg$values[seq_len(K)], 0))      # per-PC sd (spectrum of the retained subspace)

basis <- list(U = U, d = dval, K = K,
              var_explained = eg$values[seq_len(K)] / sum(eg$values))
saveRDS(basis, OUT)
cat(sprintf("wrote %s\n  U %dx%d ; top-%d var explained = %s (cum %.1f%%)\n",
            OUT, nrow(U), ncol(U), K,
            paste(sprintf("%.1f%%", 100 * basis$var_explained), collapse=" "),
            100 * sum(basis$var_explained)))
