## ---------------------------------------------------------------------------
##  Consensus clustering in tangent PCA coordinates: bootstrap k-means,
##  co-association matrix, spectral partition, and a real representative
##  shape per cluster.
## ---------------------------------------------------------------------------

# Tangent coordinates (n x km) at the Procrustes mean: each configuration is
# centered, scaled to unit size, rotated onto the mean without reflection and
# projected onto the tangent space at the mean.
.tangent_coords <- function(X) {
  k <- dim(X)[1L]; m <- dim(X)[2L]; n <- dim(X)[3L]
  nrm <- function(A) { A <- sweep(A, 2L, colMeans(A)); A / sqrt(sum(A^2)) }
  mu <- nrm(shape_mean(X))
  V <- matrix(0, n, k * m)
  for (i in seq_len(n)) {
    A <- nrm(X[, , i])
    s <- svd(crossprod(A, mu)); U <- s$u
    if (det(U %*% t(s$v)) < 0) U[, m] <- -U[, m]
    Ai <- A %*% (U %*% t(s$v))
    V[i, ] <- as.vector(Ai - mu * sum(Ai * mu))
  }
  V
}

# Tangent coordinates at the Procrustes mean, reduced by PCA.
.tangent_pca <- function(X, var_explained, max_dim, min_dim = 2L) {
  n <- dim(X)[3L]
  V <- .tangent_coords(X)
  p <- stats::prcomp(V, center = TRUE, scale. = FALSE,
                     rank. = min(max_dim, n - 1L))
  cv <- cumsum(p$sdev^2) / sum(p$sdev^2)
  dv <- which(cv >= var_explained)[1L]
  dv <- if (is.na(dv)) max_dim else min(max_dim, dv)
  dv <- min(max(min_dim, dv), ncol(p$x))
  list(scores = p$x, dim = dv, cumvar = cv)
}

# Normalized spectral partition of a co-association matrix.
.spectral_partition <- function(C, K) {
  dg <- rowSums(C)
  M <- C / sqrt(outer(dg, dg))
  U <- eigen(M, symmetric = TRUE)$vectors[, seq_len(K), drop = FALSE]
  U <- U / pmax(sqrt(rowSums(U^2)), 1e-12)
  stats::kmeans(U, K, nstart = 20L)$cluster
}

#' Consensus clustering of shapes in tangent PCA coordinates
#'
#' Clusters \eqn{k \times m \times n} landmark configurations with a
#' bootstrap ensemble of k-means fits in the leading principal components of
#' the tangent coordinates, a co-association matrix, and a spectral partition
#' of that matrix. Unlike \code{\link{bagclust1}}, no label alignment is
#' needed, and the result carries a real observed shape as representative of
#' each cluster.
#'
#' @param X Numeric array of dimension \eqn{k \times m \times n}, with
#'   \eqn{m} equal to 2 or 3.
#' @param K Number of clusters (>= 2).
#' @param B Number of bootstrap replicates (default 30).
#' @param var_explained Proportion of tangent variance that fixes the number
#'   of principal components (default 0.90).
#' @param max_dim Largest number of components kept (default 10). The number
#'   used is at least \code{max(2, K - 1)} and at most \code{max_dim}.
#' @param seed Optional integer seed. The caller's random-number stream is
#'   restored on exit. With \code{NULL} the current stream is used.
#' @param base Base clusterer of each replicate: \code{"kmeans"} (default,
#'   the original behavior) or \code{"vii"}, the spherical-variance-per-cluster
#'   Gaussian mixture of \code{\link{shape_mixture}} (EM with 3 starts),
#'   which assigns all \eqn{n} observations by maximum posterior. Use
#'   \code{"vii"} when clusters have unequal spread.
#' @param ... With \code{base = "kmeans"}, further arguments for
#'   \code{\link[stats]{kmeans}} in each replicate (defaults
#'   \code{nstart = 1}, \code{iter.max = 30}). Ignored for \code{"vii"}.
#'
#' @details
#' The tangent coordinates are those of the Procrustes mean shape
#' (\code{\link{shape_mean}}): each configuration is centered, scaled to unit
#' size, rotated onto the mean without reflection, and projected onto the
#' tangent space at the mean. A principal component analysis gives the scores;
#' the smallest number of components with cumulative variance of at least
#' \code{var_explained} is used, between \code{max(2, K - 1)} and \code{max_dim}:
#' K groups need at least K - 1 dimensions to be separated by k-means. With
#' K = 3 this is the earlier floor of 2.
#'
#' Each replicate draws a bootstrap sample of the scores, drops duplicate
#' rows, runs k-means with a single start, and assigns all \eqn{n}
#' observations to the nearest of the \eqn{K} centers. The co-association
#' matrix holds the share of replicates in which two observations get the
#' same label. The partition is the k-means (20 starts) of the row-normalized
#' first \eqn{K} eigenvectors of the symmetrically normalized matrix.
#'
#' The representative of a cluster is the member with the smallest Riemannian
#' distance \eqn{\rho} to the full Procrustes mean of that cluster.
#'
#' @return An object of class \code{"shape_ensemble"}, a list with
#'   \describe{
#'     \item{cluster}{Integer vector of length n, values in 1:K.}
#'     \item{coassoc}{The n x n co-association matrix.}
#'     \item{dim}{Number of principal components used.}
#'     \item{representative}{Integer vector of length K: index of the real
#'       observation closest to the Procrustes mean of each cluster.}
#'     \item{centers}{Array \eqn{k \times m \times K} with the cluster means.}
#'     \item{B, K, var_explained, base}{The options used.}
#'     \item{n_fallback}{With \code{base = "vii"}, the number of replicates
#'       whose mixture fit fell back to k-means.}
#'   }
#'
#' @examples
#' sim <- sim_shapes3d(n = 40, sigma = 1, k = 8, seed = 1)
#' fit <- shape_ensemble(sim$X, K = 2, B = 10, seed = 1)
#' fit$dim
#' fit$representative
#' table(sim$group, fit$cluster)
#'
#' @seealso \code{\link{bagclust1}}, \code{\link{shape_kmeans}}
#' @export
shape_ensemble <- function(X, K, B = 30L, var_explained = 0.90,
                           max_dim = 10L, seed = NULL,
                           base = c("kmeans", "vii"), ...) {
  base <- match.arg(base)
  X <- .check_array(X, "X", min_n = 3L)
  n <- dim(X)[3L]
  stopifnot(is.numeric(K), length(K) == 1L, K >= 2, K <= n,
            is.numeric(B), length(B) == 1L, B >= 1,
            is.numeric(var_explained), length(var_explained) == 1L,
            var_explained > 0, var_explained <= 1,
            is.numeric(max_dim), length(max_dim) == 1L, max_dim >= 2)
  K <- as.integer(K); B <- as.integer(B); max_dim <- as.integer(max_dim)
  .local_seed(seed)
  km_args <- list(nstart = 1L, iter.max = 30L)
  dots <- list(...)
  km_args[names(dots)] <- dots

  tp <- .tangent_pca(X, var_explained, max(max_dim, K - 1L), max(2L, K - 1L))
  S <- tp$scores[, seq_len(tp$dim), drop = FALSE]
  C <- matrix(0, n, n); n_fb <- 0L
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    Sb <- unique(S[idx, , drop = FALSE])
    if (base == "vii") {
      f <- .vii_em(Sb, K, nstart = 3L)
      n_fb <- n_fb + f$fallback
      lab <- max.col(.vii_post(S, f$mu, f$s2, f$pi_k)$z, ties.method = "first")
      C <- C + outer(lab, lab, "==")
      next
    }
    cen <- suppressWarnings(do.call(stats::kmeans,
                                    c(list(Sb, K), km_args)))$centers
    lab <- max.col(-as.matrix(stats::dist(rbind(cen, S)))[-seq_len(K),
                                                       seq_len(K), drop = FALSE],
                   ties.method = "first")
    C <- C + outer(lab, lab, "==")
  }
  C <- C / B
  cl <- .spectral_partition(C, K)

  centers <- array(0, c(dim(X)[1:2], K))
  rep_i <- integer(K)
  for (h in seq_len(K)) {
    idx <- which(cl == h)
    if (length(idx) == 0L) { rep_i[h] <- NA_integer_; next }
    mu <- if (length(idx) > 1L) shape_mean(X[, , idx, drop = FALSE])
          else X[, , idx]
    centers[, , h] <- mu
    d <- shape_distance(X[, , idx, drop = FALSE], "rho",
                        Y = array(mu, c(dim(X)[1:2], 1L)))
    rep_i[h] <- idx[which.min(d)]
  }
  structure(list(cluster = cl, coassoc = C, dim = tp$dim,
                 representative = rep_i, centers = centers,
                 B = B, K = K, var_explained = var_explained,
                 base = base, n_fallback = if (base == "vii") n_fb else NULL),
            class = "shape_ensemble")
}

#' @export
print.shape_ensemble <- function(x, ...) {
  cat("Shape ensemble (tangent PCA, co-association, spectral)\n")
  cat("  n =", length(x$cluster), " K =", x$K, " B =", x$B,
      " PCA dimension =", x$dim, "\n")
  cat("  cluster sizes:", paste(tabulate(x$cluster, x$K), collapse = " "), "\n")
  cat("  representatives:", paste(x$representative, collapse = " "), "\n")
  invisible(x)
}
