## ---------------------------------------------------------------------------
##  BagClust1: bagged clustering ensemble (Dudoit & Fridlyand 2003) adapted
##  to shape data. Given a base clusterer, draw B bootstrap samples of the
##  n observations, cluster each, align labels to a reference clustering
##  by an optimal permutation, and take the majority vote.
##
##  Label alignment uses the Hungarian algorithm (O(K^3)) via
##  clue::solve_LSAP when clue is available, or an exhaustive search over
##  S_K when K is small (K <= 8).
## ---------------------------------------------------------------------------

#' BagClust1 for shape clustering
#'
#' @param Z Complex pre-shape matrix of size (k-1) x n.
#' @param K Number of clusters.
#' @param base Base clusterer. One of "pam" (default) or "rankkmed".
#' @param B Number of bootstrap replicates (default 100).
#' @param distance One of "rho" (default), "dF".
#' @param seed Optional integer seed for reproducibility.
#' @param parallel Logical. If TRUE, use \code{parallel::mclapply}
#'   (Unix/macOS) or \code{parallel::parLapply} (Windows).
#' @param n_cores Number of cores to use if \code{parallel = TRUE}.
#'   Default = \code{max(1, parallel::detectCores() - 1)}.
#' @param ... Additional arguments passed to the base clusterer
#'   (e.g., \code{m = 7} for rankkmed).
#'
#' @return A list with components:
#'   \describe{
#'     \item{cluster}{Integer vector of length n, the majority-vote
#'       consensus labels.}
#'     \item{votes}{n x K integer matrix of vote counts.}
#'     \item{reference}{Integer vector, the reference (non-bagged)
#'       clustering used for label alignment.}
#'     \item{B}{Number of bootstrap replicates.}
#'     \item{base}{Base clusterer name.}
#'   }
#'
#' @examples
#' set.seed(1)
#' k <- 4; n <- 60
#' # 3 clusters via mean shift
#' X <- array(0, dim = c(k, 2, n))
#' for (i in seq_len(n)) {
#'   g <- ((i - 1) %% 3) + 1
#'   mu <- matrix(c(0,0, 1,0, 0,1, 1,1), 4, 2, byrow = TRUE)
#'   X[,,i] <- mu + matrix(rnorm(k*2, sd = 0.1 * g), k, 2)
#' }
#' Z <- landmarks_to_preshape(X)
#' fit <- bagclust1(Z, K = 3, base = "pam", B = 30, seed = 1)
#' table(fit$cluster)
#'
#' @export
bagclust1 <- function(Z, K,
                      base = c("pam", "rankkmed"),
                      B = 100L,
                      distance = c("rho", "dF"),
                      seed = NULL,
                      parallel = FALSE,
                      n_cores = NULL,
                      ...) {
  base <- match.arg(base)
  distance <- match.arg(distance)

  if (!is.null(seed)) set.seed(seed)
  n <- ncol(Z)
  base_fun <- switch(base,
                     pam      = shape_pam,
                     rankkmed = shape_rankkmed)

  # Precompute the full-sample distance matrix so each bootstrap avoids
  # recomputing on-the-fly by sub-indexing.
  D_full <- switch(distance,
                   rho = rho_matrix(Z),
                   dF  = dF_matrix(Z))

  # Reference clustering on the full sample
  ref_fit <- base_fun(D_full, K = K, ...)
  ref_labels <- ref_fit$cluster

  # Bootstrap replicates
  one_replicate <- function(b) {
    idx <- sample.int(n, n, replace = TRUE)
    D_b <- D_full[idx, idx]
    fit <- base_fun(D_b, K = K, ...)
    # Align labels of the bootstrap fit to the reference on shared indices.
    aligned <- align_labels(fit$cluster, ref_labels[idx], K)
    # Push aligned labels back to the *original* index space, voting only
    # for observations that appeared in this bootstrap.
    out <- integer(n)
    for (i in seq_along(idx)) out[idx[i]] <- aligned[i]
    out
  }

  results <- if (parallel && requireNamespace("parallel", quietly = TRUE)) {
    if (is.null(n_cores)) {
      n_cores <- max(1L, parallel::detectCores(logical = TRUE) - 1L)
    }
    if (.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(n_cores)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      parallel::clusterExport(cl, c("D_full", "ref_labels", "base_fun",
                                    "align_labels", "K", "n"),
                              envir = environment())
      parallel::parLapply(cl, seq_len(B), one_replicate)
    } else {
      parallel::mclapply(seq_len(B), one_replicate,
                         mc.cores = n_cores)
    }
  } else {
    lapply(seq_len(B), one_replicate)
  }

  # Tally votes (zero votes for observations absent in a given bootstrap)
  votes <- matrix(0L, nrow = n, ncol = K)
  for (r in results) {
    for (i in seq_len(n)) {
      if (r[i] > 0L) votes[i, r[i]] <- votes[i, r[i]] + 1L
    }
  }
  cluster <- max.col(votes, ties.method = "first")

  structure(
    list(cluster   = as.integer(cluster),
         votes     = votes,
         reference = as.integer(ref_labels),
         B         = as.integer(B),
         base      = base),
    class = "bagclust1"
  )
}

#' @export
print.bagclust1 <- function(x, ...) {
  cat("BagClust1 fit\n")
  cat("  base clusterer:", x$base, "\n")
  cat("  B =", x$B, "bootstrap replicates\n")
  cat("  cluster sizes:", paste(table(x$cluster), collapse = " / "), "\n")
  invisible(x)
}

#' Align cluster labels via optimal permutation
#'
#' Given two label vectors on the same set of n objects, returns the
#' relabelling of \code{a} to \code{b} that maximises the number of
#' matches.
#'
#' Uses the Hungarian algorithm via \code{clue::solve_LSAP} if available;
#' otherwise falls back to exhaustive enumeration of \eqn{S_K} (feasible
#' for K <= 8).
#'
#' @param a Integer vector of labels in 1:K.
#' @param b Integer vector of labels in 1:K, same length as \code{a}.
#' @param K Number of clusters.
#'
#' @return Integer vector of the relabelled \code{a}.
#'
#' @export
align_labels <- function(a, b, K) {
  stopifnot(length(a) == length(b), K >= 1L)
  # Overlap matrix M[i, j] = |{a == i and b == j}|
  M <- matrix(0L, nrow = K, ncol = K)
  for (t in seq_along(a)) {
    if (a[t] >= 1L && a[t] <= K && b[t] >= 1L && b[t] <= K) {
      M[a[t], b[t]] <- M[a[t], b[t]] + 1L
    }
  }
  if (requireNamespace("clue", quietly = TRUE)) {
    # Convert to cost matrix (maximise -> minimise)
    C <- max(M) - M
    perm <- as.integer(clue::solve_LSAP(C))
  } else if (K <= 8L) {
    # Exhaustive search
    perms <- gtools_permutations(K)
    best <- 0L
    perm <- seq_len(K)
    for (p in seq_len(nrow(perms))) {
      pp <- perms[p, ]
      score <- sum(M[cbind(seq_len(K), pp)])
      if (score > best) { best <- score; perm <- pp }
    }
  } else {
    stop("K > 8 and 'clue' package not available. Install 'clue' for the Hungarian algorithm.")
  }
  # Relabel a
  perm[a]
}

# Minimal permutation generator (avoids gtools dependency)
gtools_permutations <- function(n) {
  if (n == 1L) return(matrix(1L, 1, 1))
  smaller <- gtools_permutations(n - 1L)
  out <- matrix(0L, nrow = factorial(n), ncol = n)
  r <- 1L
  for (i in seq_len(n)) {
    for (j in seq_len(nrow(smaller))) {
      # Insert i at each of n positions
      row <- integer(n)
      row[i] <- n
      others <- setdiff(seq_len(n), i)
      row[others] <- smaller[j, ]
      out[r, ] <- row
      r <- r + 1L
    }
  }
  out
}
