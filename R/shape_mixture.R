## ---------------------------------------------------------------------------
##  Gaussian mixture with a spherical variance per cluster (mclust "VII") in
##  tangent PCA coordinates, fitted by EM written here (no new dependency).
##  Motivation: groups with unequal spread, where the Voronoi boundary of
##  k-means is not the Bayes boundary.
## ---------------------------------------------------------------------------

# k-means++ seeding on the rows of S (n x p); returns K row indices.
.kmpp_seeds <- function(S, K) {
  n <- nrow(S)
  idx <- integer(K); idx[1L] <- sample.int(n, 1L)
  d2 <- rowSums((S - rep(S[idx[1L], ], each = n))^2)
  for (j in seq_len(K)[-1L]) {
    if (sum(d2) <= 0) { rest <- setdiff(seq_len(n), idx[seq_len(j - 1L)])
      idx[j] <- rest[sample.int(length(rest), 1L)] }
    else idx[j] <- sample.int(n, 1L, prob = d2)
    d2 <- pmin(d2, rowSums((S - rep(S[idx[j], ], each = n))^2))
  }
  idx
}

# Posterior matrix (n x K) of a VII model, by log-sum-exp.
.vii_post <- function(S, mu, s2, pi_k) {
  p <- ncol(S); K <- nrow(mu)
  lg <- matrix(0, nrow(S), K)
  for (h in seq_len(K)) {
    d2 <- rowSums((S - rep(mu[h, ], each = nrow(S)))^2)
    lg[, h] <- log(pi_k[h]) - 0.5 * p * log(2 * pi * s2[h]) - 0.5 * d2 / s2[h]
  }
  m <- apply(lg, 1L, max)
  lse <- m + log(rowSums(exp(lg - m)))
  list(z = exp(lg - lse), loglik = sum(lse))
}

# One EM run from hard labels `lab`. Returns NULL on numerical failure.
.vii_em_run <- function(S, K, lab, floor_var, iter.max, tol) {
  n <- nrow(S); p <- ncol(S)
  z <- matrix(0, n, K); z[cbind(seq_len(n), lab)] <- 1
  ll_old <- -Inf; conv <- FALSE; n_reset <- 0L; it <- 0L
  repeat {
    it <- it + 1L
    nk <- colSums(z)
    # empty or degenerate component (less than one observation of weight):
    # restart it at the observation worst explained by the current model
    bad <- which(nk < 1)
    if (length(bad) > 0L) {
      n_reset <- n_reset + length(bad)
      if (n_reset > 5L * K) return(NULL)
      for (h in bad) {
        if (it > 1L) {
          sc <- rowSums(sapply(seq_len(K), function(g)
            exp(log(pi_k[g]) - 0.5 * p * log(2 * pi * s2[g]) -
                0.5 * rowSums((S - rep(mu[g, ], each = n))^2) / s2[g])))
          j <- which.min(sc)
        } else j <- sample.int(n, 1L)
        z[, h] <- 0; z[j, ] <- 0; z[j, h] <- 1
      }
      nk <- colSums(z)
    }
    pi_k <- nk / sum(nk)
    mu <- crossprod(z, S) / nk
    s2 <- numeric(K)
    for (h in seq_len(K))
      s2[h] <- sum(z[, h] * rowSums((S - rep(mu[h, ], each = n))^2)) / (p * nk[h])
    s2 <- pmax(s2, floor_var)
    po <- .vii_post(S, mu, s2, pi_k)
    if (!is.finite(po$loglik)) return(NULL)
    z <- po$z
    if (abs(po$loglik - ll_old) <= tol * (1 + abs(po$loglik))) { conv <- TRUE; break }
    ll_old <- po$loglik
    if (it >= iter.max) break
  }
  list(mu = mu, s2 = s2, pi_k = pi_k, z = z, loglik = po$loglik,
       converged = conv, iter = it, at_floor = any(s2 <= floor_var * (1 + 1e-8)))
}

# EM for a VII mixture with multiple starts and a k-means fallback.
# Returns list(mu, s2, pi_k, z, loglik, converged, fallback).
.vii_em <- function(S, K, nstart = 10L, init = "kmeans++", iter.max = 200L,
                    tol = 1e-8, var_floor = 1e-6) {
  n <- nrow(S); p <- ncol(S)
  tv <- sum(apply(S, 2L, stats::var)) / p
  floor_var <- var_floor * tv
  best <- NULL
  if (is.finite(tv) && tv > 0 && n >= K) {
    for (st in seq_len(nstart)) {
      lab <- tryCatch({
        if (init == "random") sample(rep_len(seq_len(K), n))
        else {
          ctr <- S[.kmpp_seeds(S, K), , drop = FALSE]
          lb <- max.col(-as.matrix(stats::dist(rbind(ctr, S)))[-seq_len(K), seq_len(K),
                                                            drop = FALSE], ties.method = "first")
          if (init == "kmeans")
            lb <- suppressWarnings(stats::kmeans(S, ctr, iter.max = 20L))$cluster
          lb
        }
      }, error = function(e) NULL)
      if (is.null(lab)) next
      r <- tryCatch(.vii_em_run(S, K, lab, floor_var, iter.max, tol),
                    error = function(e) NULL)
      if (is.null(r)) next
      # a start with a variance at the floor is a likelihood spike: only keep
      # it if no regular start exists
      if (is.null(best) || (best$at_floor && !r$at_floor) ||
          (best$at_floor == r$at_floor && r$loglik > best$loglik)) best <- r
    }
  }
  if (!is.null(best) && !best$at_floor) {
    best$fallback <- FALSE; return(best)
  }
  # documented fallback: k-means partition with per-group spherical variances
  km <- tryCatch(suppressWarnings(stats::kmeans(S, K, nstart = 5L, iter.max = 50L)),
                 error = function(e) NULL)
  if (is.null(km)) {
    if (!is.null(best)) { best$fallback <- TRUE; return(best) }
    stop("mixture fit and k-means fallback both failed")
  }
  lab <- km$cluster
  nk <- tabulate(lab, K)
  mu <- km$centers
  s2 <- vapply(seq_len(K), function(h) {
    if (nk[h] == 0L) return(max(tv, floor_var))
    max(sum((S[lab == h, , drop = FALSE] - rep(mu[h, ], each = nk[h]))^2) /
          (p * nk[h]), floor_var)
  }, 0)
  pi_k <- pmax(nk, 1e-8) / sum(pmax(nk, 1e-8))
  z <- matrix(0, n, K); z[cbind(seq_len(n), lab)] <- 1
  ll <- if (is.finite(tv) && tv > 0) .vii_post(S, mu, s2, pi_k)$loglik else NA_real_
  list(mu = mu, s2 = s2, pi_k = pi_k, z = z, loglik = ll,
       converged = FALSE, iter = 0L, at_floor = FALSE, fallback = TRUE)
}

#' Gaussian mixture clustering of shapes with a variance per cluster
#'
#' Clusters \eqn{k \times m \times n} landmark configurations with a Gaussian
#' mixture with spherical covariance \eqn{\sigma_h^2 I} for each cluster
#' (model \code{"VII"} of \pkg{mclust}) in the leading principal components of
#' the tangent coordinates. The EM algorithm is implemented here, with
#' safeguards for duplicated rows, tiny groups and degenerate components.
#'
#' @param X Numeric array of dimension \eqn{k \times m \times n}, with
#'   \eqn{m} equal to 2 or 3.
#' @param K Number of clusters (>= 2).
#' @param var_explained Proportion of tangent variance that fixes the number
#'   of principal components (default 0.90).
#' @param max_dim Largest number of components kept (default 10). The number
#'   used is at least \code{max(2, K - 1)}, as in \code{\link{shape_ensemble}}.
#' @param model Covariance model. Only \code{"VII"} (spherical, unequal
#'   volume) is implemented.
#' @param nstart Number of EM starts (default 10); the best log-likelihood
#'   among starts without a variance at the floor is kept.
#' @param init Initialization of each start: \code{"kmeans++"} (default,
#'   k-means++ seeding, nearest-seed labels), \code{"kmeans"} (seeding
#'   followed by Lloyd iterations) or \code{"random"} (random labels).
#' @param iter.max Maximum EM iterations per start (default 200).
#' @param tol Relative log-likelihood tolerance for convergence (default 1e-8).
#' @param var_floor Variance floor, as a fraction of the mean per-dimension
#'   variance of the scores (default 1e-6).
#' @param seed Optional integer seed. The caller's random-number stream is
#'   restored on exit. With \code{NULL} the current stream is used.
#' @param ... Not used; reserved.
#'
#' @details
#' The scores are those of \code{shape_ensemble}: tangent coordinates at the
#' Procrustes mean, then PCA, with the smallest number of components reaching
#' \code{var_explained}, between \code{max(2, K - 1)} and \code{max_dim}.
#'
#' Why a variance per cluster: k-means (and so a k-means ensemble) assigns
#' to the nearest center, which is the Bayes rule only when all clusters have
#' the same spread. With unequal spread the Bayes boundary is a sphere around
#' the compact cluster and the bisector misassigns the tail of the diffuse
#' one. This model estimates that spread.
#'
#' Safeguards. Variances are floored at \code{var_floor} times the mean
#' per-dimension variance of the scores, so duplicated rows (bootstrap
#' resamples) or a group of identical points never make the likelihood
#' singular. A component with less than one observation of total weight is
#' restarted at the observation worst explained by the current model (at most
#' \code{5 * K} restarts per start, then that start is abandoned). Starts that
#' end with a variance at the floor are treated as likelihood spikes and are
#' used only if no regular start exists. If no start succeeds, the partition
#' is the k-means (5 starts) partition with per-group spherical variances,
#' \code{fallback} is \code{TRUE} and \code{converged} is \code{FALSE}.
#'
#' @return An object of class \code{"shape_mixture"}, a list with
#'   \describe{
#'     \item{cluster}{Integer vector of length n, maximum posterior label.}
#'     \item{posterior}{n x K matrix of posterior probabilities.}
#'     \item{means}{K x dim matrix of component means in score coordinates.}
#'     \item{variances}{Length-K vector of spherical variances.}
#'     \item{proportions}{Length-K vector of mixing proportions.}
#'     \item{loglik}{Log-likelihood of the fit.}
#'     \item{dim}{Number of principal components used.}
#'     \item{representative}{Integer vector of length K: index of the real
#'       observation closest (Riemannian \eqn{\rho}) to the Procrustes mean
#'       of each cluster.}
#'     \item{centers}{Array \eqn{k \times m \times K} with the cluster means.}
#'     \item{converged}{Whether EM met the tolerance.}
#'     \item{fallback}{\code{TRUE} if the k-means fallback was used.}
#'     \item{K, model, var_explained}{The options used.}
#'   }
#'
#' @examples
#' sim <- sim_shapes3d(n = 40, sigma = 1, k = 8, seed = 1)
#' fit <- shape_mixture(sim$X, K = 2, seed = 1)
#' fit$variances
#' table(sim$group, fit$cluster)
#'
#' @seealso \code{\link{shape_ensemble}}, \code{\link{shape_kmeans}}
#' @export
shape_mixture <- function(X, K, var_explained = 0.90, max_dim = 10L,
                          model = "VII", nstart = 10L,
                          init = c("kmeans++", "kmeans", "random"),
                          iter.max = 200L, tol = 1e-8, var_floor = 1e-6,
                          seed = NULL, ...) {
  X <- .check_array(X, "X", min_n = 3L)
  n <- dim(X)[3L]
  model <- match.arg(model, "VII"); init <- match.arg(init)
  if (length(list(...)) > 0L) warning("extra arguments are ignored")
  stopifnot(is.numeric(K), length(K) == 1L, K >= 2, K <= n,
            is.numeric(var_explained), length(var_explained) == 1L,
            var_explained > 0, var_explained <= 1,
            is.numeric(max_dim), length(max_dim) == 1L, max_dim >= 2,
            is.numeric(nstart), nstart >= 1, is.numeric(iter.max), iter.max >= 1,
            is.numeric(var_floor), var_floor > 0)
  K <- as.integer(K); max_dim <- as.integer(max_dim)
  .local_seed(seed)
  tp <- .tangent_pca(X, var_explained, max(max_dim, K - 1L), max(2L, K - 1L))
  S <- tp$scores[, seq_len(tp$dim), drop = FALSE]
  fit <- .vii_em(S, K, as.integer(nstart), init, as.integer(iter.max), tol, var_floor)
  cl <- max.col(fit$z, ties.method = "first")
  rc <- .cluster_centers(X, cl, K)
  structure(list(cluster = cl, posterior = fit$z, means = fit$mu,
                 variances = fit$s2, proportions = fit$pi_k,
                 loglik = fit$loglik, dim = tp$dim,
                 representative = rc$representative, centers = rc$centers,
                 converged = fit$converged, fallback = fit$fallback,
                 K = K, model = model, var_explained = var_explained),
            class = "shape_mixture")
}

# Procrustes mean and real representative of each cluster.
.cluster_centers <- function(X, cl, K) {
  centers <- array(0, c(dim(X)[1:2], K))
  rep_i <- integer(K)
  for (h in seq_len(K)) {
    idx <- which(cl == h)
    if (length(idx) == 0L) { rep_i[h] <- NA_integer_; next }
    mu <- if (length(idx) > 1L) shape_mean(X[, , idx, drop = FALSE]) else X[, , idx]
    centers[, , h] <- mu
    d <- shape_distance(X[, , idx, drop = FALSE], "rho",
                        Y = array(mu, c(dim(X)[1:2], 1L)))
    rep_i[h] <- idx[which.min(d)]
  }
  list(centers = centers, representative = rep_i)
}

#' @export
print.shape_mixture <- function(x, ...) {
  cat("Shape mixture (tangent PCA, spherical Gaussian per cluster, VII)\n")
  cat("  n =", length(x$cluster), " K =", x$K, " PCA dimension =", x$dim, "\n")
  cat("  cluster sizes:", paste(tabulate(x$cluster, x$K), collapse = " "), "\n")
  cat("  variances:", paste(signif(x$variances, 3), collapse = " "), "\n")
  cat("  loglik =", format(x$loglik, digits = 6), " converged =", x$converged,
      " fallback =", x$fallback, "\n")
  invisible(x)
}
