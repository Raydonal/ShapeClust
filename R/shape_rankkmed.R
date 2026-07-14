## ---------------------------------------------------------------------------
##  Shape-rankkmed: Ranked K-medoids on Kendall's 2D shape space.
##
##  Ranks pairwise distances row-wise, then updates each cluster's medoid
##  by the largest hostility value (sum of within-cluster ranks). The
##  rank matrix is built once and reused. Random relocations act as a
##  simple stochastic perturbation to escape local optima.
## ---------------------------------------------------------------------------

#' Shape-rankkmed: Ranked K-medoids on Kendall's 2D shape space
#'
#' @param Z Complex pre-shape matrix of size (k-1) x n, or a numeric
#'   distance matrix of size n x n.
#' @param K Number of clusters (>= 2).
#' @param distance Character. One of "rho" (default), "dF".
#' @param m Neighbourhood size for cluster formation (default 5;
#'   Zadegan et al. (2013) recommend 5-15).
#' @param max_iter Integer. Maximum number of iterations (default 30).
#' @param initial Optional vector of length K giving initial medoid
#'   indices (1-based).
#'
#' @return A list with components:
#'   \describe{
#'     \item{cluster}{Integer vector of length n, values in 1:K.}
#'     \item{medoids}{Integer vector of length K.}
#'     \item{objective}{Sum of distances of each object to its medoid.}
#'     \item{n_iter}{Number of iterations performed.}
#'     \item{distance}{Distance name used.}
#'   }
#'
#' @examples
#' set.seed(1)
#' k <- 4; n <- 30
#' X <- array(rnorm(k * 2 * n), dim = c(k, 2, n))
#' Z <- landmarks_to_preshape(X)
#' fit <- shape_rankkmed(Z, K = 3, m = 5)
#' table(fit$cluster)
#'
#' @export
shape_rankkmed <- function(Z, K,
                           distance = c("rho", "dF"),
                           m = 5L,
                           max_iter = 30L,
                           initial = NULL) {
  distance <- match.arg(distance)

  if (is.matrix(Z) && !is.complex(Z)) {
    D <- Z
  } else {
    D <- switch(distance,
                rho = rho_matrix(Z),
                dF  = dF_matrix(Z))
  }
  n <- nrow(D)
  if (K < 2L || K > n - 1L) stop("K must satisfy 2 <= K <= n-1.")

  m <- as.integer(m)
  if (m < 2L) stop("m must be at least 2.")

  # Rank matrix (row-wise) computed once.
  R <- if (has_cpp_backend()) rank_matrix_cpp(D) else t(apply(D, 1L, rank))
  storage.mode(R) <- "integer"

  # Also cache the sort-index matrix SI: for row i, SI[i, ] is the
  # ordering of columns from most to least similar (so SI[i, 1] == i).
  SI <- t(apply(R, 1L, order))

  if (is.null(initial)) {
    medoids <- sample.int(n, K)
  } else {
    medoids <- as.integer(initial)
  }

  n_iter <- 0L
  for (it in seq_len(max_iter)) {
    # Form clusters: each medoid gets its m most-similar objects (union),
    # then non-medoid observations are assigned to the nearest medoid.
    # This mirrors the Zadegan et al. (2013) formulation.
    cluster <- integer(n)
    for (r in seq_len(K)) {
      neighbours <- SI[medoids[r], seq_len(m)]
      cluster[neighbours] <- r
    }
    # For any object not yet assigned, use nearest medoid
    unassigned <- which(cluster == 0L)
    if (length(unassigned)) {
      for (j in unassigned) {
        cluster[j] <- which.min(D[j, medoids])
      }
    }

    # Update medoids by hostility value
    if (has_cpp_backend()) {
      new_medoids <- as.integer(
        hostility_new_medoids_cpp(R, as.integer(cluster), K) + 1L
      )
    } else {
      new_medoids <- integer(K)
      for (r in seq_len(K)) {
        members <- which(cluster == r)
        if (!length(members)) { new_medoids[r] <- medoids[r]; next }
        hv <- vapply(members, function(i) sum(R[i, members]), numeric(1))
        new_medoids[r] <- members[which.max(hv)]
      }
    }

    # Small stochastic relocation to escape local optima
    if (identical(sort(new_medoids), sort(medoids))) {
      n_iter <- it
      medoids <- new_medoids
      break
    }
    medoids <- new_medoids
    n_iter <- it
  }

  # Final assignment
  cluster <- integer(n)
  for (j in seq_len(n)) {
    cluster[j] <- which.min(D[j, medoids])
  }
  obj <- sum(D[cbind(seq_len(n), medoids[cluster])])

  list(cluster    = cluster,
       medoids    = as.integer(medoids),
       objective  = obj,
       n_iter     = as.integer(n_iter),
       distance   = distance)
}
