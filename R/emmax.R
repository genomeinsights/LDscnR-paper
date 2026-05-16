## ======================================================================
## EMMAX (permutation-free) + EMMA REML helpers
## Package style note: avoid library() inside functions; call pkg::fun()
## ======================================================================

library(data.table)

#' EMMAX single-locus association scan (permutation-free)
#'
#' Performs an EMMAX-style association test for each SNP (each column of \code{X})
#' against phenotype \code{Y}, accounting for relatedness via a kinship matrix
#' \code{K}. Variance components are estimated once under the null model
#' (intercept + optional fixed-effect covariates). Each SNP is then tested by
#' comparing residual sums of squares (RSS) under the null vs. SNP-augmented model
#' after whitening by the fitted covariance.
#'
#' @param Y Numeric vector of length \code{n}, phenotype.
#' @param X Numeric matrix (\code{n x m}), genotypes (individuals x SNPs).
#' @param K Numeric matrix (\code{n x n}), kinship/relatedness matrix.
#' @param Covar Optional fixed-effect covariates. Can be a vector (length \code{n}),
#'   a matrix (\code{n x p}) or a data.frame. Included alongside an intercept.
#' @param nbchunks Integer >= 1. Number of SNP chunks processed sequentially to
#'   reduce peak memory during whitening. Default \code{2}.
#'
#' @return A list with elements:
#' \itemize{
#'   \item \code{F}: numeric vector (length \code{m}), F statistics.
#'   \item \code{pval}: numeric vector (length \code{m}), asymptotic p-values.
#'   \item \code{Rsq}: numeric vector (length \code{m}), pseudo-R² from RSS ratio.
#' }
#'
#' @export
emmax <- function(Y, X, K, Covar = NULL, nbchunks = 2L) {
  
  ## -----------------------------
  ## Basic checks
  ## -----------------------------
  n <- length(Y)
  stopifnot(is.matrix(X), nrow(X) == n)
  stopifnot(is.matrix(K), nrow(K) == n, ncol(K) == n)
  nbchunks <- as.integer(nbchunks)
  stopifnot(nbchunks >= 1L)
  
  m <- ncol(X)
  
  ## -----------------------------
  ## Fixed-effect design: intercept + covariates
  ## -----------------------------
  X0 <- matrix(1, n, 1)  # intercept
  
  if (!is.null(Covar)) {
    Cmat <- if (is.data.frame(Covar)) as.matrix(Covar) else as.matrix(Covar)
    stopifnot(nrow(Cmat) == n)
    # Drop covariate columns with zero variance (prevents singular X'X)
    keep <- apply(Cmat, 2, function(z) is.finite(stats::var(z)) && stats::var(z) > 0)
    if (any(keep)) X0 <- cbind(X0, Cmat[, keep, drop = FALSE])
  }
  
  ## -----------------------------
  ## Kinship normalization
  ## -----------------------------
  K_norm <- (n - 1) / sum((diag(n) - matrix(1, n, n) / n) * K) * K
  
  ## -----------------------------
  ## Fit null model once (variance components via REML)
  ## -----------------------------
  null <- emma.REMLE(Y, X0, K_norm)
  
  ## Whitening transform M:  M' M = (vg*K + ve*I)^-1
  V <- null$vg * K_norm + null$ve * diag(n)
  M <- solve(chol(V))
  
  Y_t  <- crossprod(M, Y)
  X0_t <- crossprod(M, X0)
  
  ## Null RSS: intercept + covariates
  RSS_H0 <- sum(stats::lsfit(X0_t, Y_t, intercept = FALSE)$residuals^2)
  
  ## Degrees of freedom:
  ## Reduced model has q0 = ncol(X0) parameters, full adds 1 SNP parameter
  q0 <- ncol(X0)
  df1 <- 1L
  df2 <- n - q0 - df1
  if (df2 <= 0) stop("Not enough degrees of freedom: check n and number of covariates.")
  
  ## -----------------------------
  ## SNP scan in chunks (memory-friendly)
  ## -----------------------------
  nbchunks <- min(nbchunks, m)
  cuts <- unique(round(seq(0, m, length.out = nbchunks + 1L)))
  
  RSSf <- numeric(m)
  
  for (j in seq_len(length(cuts) - 1L)) {
    idx <- (cuts[j] + 1L):cuts[j + 1L]
    if (!length(idx)) next
    
    X_t <- crossprod(M, X[, idx, drop = FALSE])
    
    RSSf[idx] <- apply(X_t, 2, function(x) {
      sum(stats::lsfit(cbind(X0_t, x), Y_t, intercept = FALSE)$residuals^2)
    })
  }
  
  ## -----------------------------
  ## Test statistics (same structure as your original)
  ## -----------------------------
  Rsq <- 1 - 1 / (RSS_H0 / RSSf)                          # pseudo-R²
  F   <- (RSS_H0 / RSSf - 1) * (df2 / df1)                # F statistic
  pval <- stats::pf(F, df1, df2, lower.tail = FALSE)      # asymptotic p-values
  
  list(F = F, pval = pval, Rsq = Rsq)
}

## ======================================================================
## EMMA REML helpers
## ======================================================================

#' Restricted maximum likelihood for EMMA (variance components)
#'
#' Fits variance components under mixed model: y ~ X (fixed) + u (random),
#' with Var(u) proportional to K and residual variance ve*I.
#'
#' @keywords internal
#' Restricted maximum likelihood for EMMA (no Z support)
#'
#' Estimates variance components (vg, ve) under a mixed model with covariance
#' vg*K + ve*I. This version only supports the common use case Z = NULL.
#'
#' @keywords internal
emma.REMLE <- function(y, X, K, ngrids = 100, llim = -10, ulim = 10,
                       esp = 1e-10, eig.R = NULL) {
  
  n <- length(y)
  q <- ncol(X)
  stopifnot(nrow(X) == n)
  stopifnot(nrow(K) == n, ncol(K) == n)
  
  ## Guard against singular fixed-effect design
  if (det(crossprod(X, X)) == 0) {
    warning("X is singular")
    return(list(REML = 0, delta = 0, ve = 0, vg = 0))
  }
  
  ## Eigen decomposition for REML problem (Z = NULL)
  if (is.null(eig.R)) eig.R <- emma.eigen.R.wo.Z(K, X)
  
  etas <- crossprod(eig.R$vectors, y)
  
  logdelta <- (0:ngrids) / ngrids * (ulim - llim) + llim
  m <- length(logdelta)
  delta <- exp(logdelta)
  
  Lambdas <- matrix(eig.R$values, n - q, m) + matrix(delta, n - q, m, byrow = TRUE)
  Etasq   <- matrix(etas * etas, n - q, m)
  
  ## REML log-likelihood and derivative wrt delta (on log-scale grid)
  LL  <- 0.5 * (
    (n - q) * (log((n - q) / (2 * pi)) - 1 - log(colSums(Etasq / Lambdas))) -
      colSums(log(Lambdas))
  )
  
  dLL <- 0.5 * delta * (
    (n - q) * colSums(Etasq / (Lambdas * Lambdas)) / colSums(Etasq / Lambdas) -
      colSums(1 / Lambdas)
  )
  
  ## Find candidate maxima by detecting sign changes in dLL
  optlogdelta <- numeric(0)
  optLL <- numeric(0)
  
  if (dLL[1] < esp) {
    optlogdelta <- c(optlogdelta, llim)
    optLL <- c(optLL, emma.delta.REML.LL.wo.Z(llim, eig.R$values, etas))
  }
  if (dLL[m - 1] > -esp) {
    optlogdelta <- c(optlogdelta, ulim)
    optLL <- c(optLL, emma.delta.REML.LL.wo.Z(ulim, eig.R$values, etas))
  }
  
  for (i in 1:(m - 1)) {
    if ((dLL[i] * dLL[i + 1] < -esp * esp) && (dLL[i] > 0) && (dLL[i + 1] < 0)) {
      r <- uniroot(emma.delta.REML.dLL.wo.Z,
                   lower = logdelta[i], upper = logdelta[i + 1],
                   lambda = eig.R$values, etas = etas)
      optlogdelta <- c(optlogdelta, r$root)
      optLL <- c(optLL, emma.delta.REML.LL.wo.Z(r$root, eig.R$values, etas))
    }
  }
  
  ## Pick best candidate
  maxdelta <- exp(optlogdelta[which.max(optLL)])
  maxLL <- max(optLL)
  
  ## vg (called maxva in original code) and ve
  vg <- sum(etas * etas / (eig.R$values + maxdelta)) / (n - q)
  ve <- vg * maxdelta
  
  list(REML = maxLL, delta = maxdelta, ve = ve, vg = vg)
}

#' Eigen decomposition for EMMA without Z
#' @keywords internal
emma.eigen.R.wo.Z <- function(K, X) {
  n <- nrow(X)
  q <- ncol(X)
  S <- diag(n) - X %*% solve(crossprod(X, X)) %*% t(X)
  eig <- eigen(S %*% (K + diag(1, n)) %*% S, symmetric = TRUE)
  stopifnot(!is.complex(eig$values))
  list(
    values  = eig$values[1:(n - q)] - 1,
    vectors = eig$vectors[, 1:(n - q), drop = FALSE]
  )
}

#' dLL for REML root finding (no Z)
#' @keywords internal
emma.delta.REML.dLL.wo.Z <- function(logdelta, lambda, etas) {
  nq <- length(etas)
  delta <- exp(logdelta)
  etasq <- etas * etas
  ldelta <- lambda + delta
  0.5 * (nq * sum(etasq / (ldelta * ldelta)) / sum(etasq / ldelta) - sum(1 / ldelta))
}

#' LL for REML objective (no Z)
#' @keywords internal
emma.delta.REML.LL.wo.Z <- function(logdelta, lambda, etas) {
  nq <- length(etas)
  delta <- exp(logdelta)
  0.5 * (nq * (log(nq / (2 * pi)) - 1 - log(sum(etas * etas / (lambda + delta)))) -
           sum(log(lambda + delta)))
}

#' Eigen decomposition for EMMA without Z
#' @keywords internal
emma.eigen.R.wo.Z <- function(K, X) {
  n <- nrow(X)
  q <- ncol(X)
  S <- diag(n) - X %*% solve(crossprod(X, X)) %*% t(X)
  eig <- eigen(S %*% (K + diag(1, n)) %*% S, symmetric = TRUE)
  stopifnot(!is.complex(eig$values))
  list(values = eig$values[1:(n - q)] - 1,
       vectors = eig$vectors[, 1:(n - q), drop = FALSE])
}

#' Eigen decomposition for EMMA with Z
#' @keywords internal
emma.eigen.R.w.Z <- function(Z, K, X, complete = TRUE) {
  if (!complete) {
    vids <- colSums(Z) > 0
    Z <- Z[, vids, drop = FALSE]
    K <- K[vids, vids, drop = FALSE]
  }
  n <- nrow(Z)
  t <- ncol(Z)
  q <- ncol(X)
  
  SZ <- Z - X %*% solve(crossprod(X, X)) %*% crossprod(X, Z)
  
  eig <- eigen(K %*% crossprod(Z, SZ), symmetric = FALSE, EISPACK = TRUE)
  if (is.complex(eig$values)) {
    eig$values  <- Re(eig$values)
    eig$vectors <- Re(eig$vectors)
  }
  
  qr.X <- qr.Q(qr(X))
  list(
    values  = eig$values[1:(t - q)],
    vectors = qr.Q(qr(cbind(SZ %*% eig$vectors[, 1:(t - q), drop = FALSE], qr.X)),
                   complete = TRUE)[, c(1:(t - q), (t + 1):n), drop = FALSE]
  )
}

#' dLL for REML root finding (no Z)
#' @keywords internal
emma.delta.REML.dLL.wo.Z <- function(logdelta, lambda, etas) {
  nq <- length(etas)
  delta <- exp(logdelta)
  etasq <- etas * etas
  ldelta <- lambda + delta
  0.5 * (nq * sum(etasq / (ldelta * ldelta)) / sum(etasq / ldelta) - sum(1 / ldelta))
}

#' LL for REML objective (no Z)
#' @keywords internal
emma.delta.REML.LL.wo.Z <- function(logdelta, lambda, etas) {
  nq <- length(etas)
  delta <- exp(logdelta)
  0.5 * (nq * (log(nq / (2 * pi)) - 1 - log(sum(etas * etas / (lambda + delta)))) -
           sum(log(lambda + delta)))
}
