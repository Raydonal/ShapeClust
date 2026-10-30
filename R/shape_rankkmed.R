## ---------------------------------------------------------------------------
##  Shape-rankkmed: Ranked K-medoids on Kendall's 2D shape space.
##
##  Follows Zadegan, Mirzaie and Sadoughi (2013) as implemented in
##  kmed::rankkmed (Budiaji 2019). Distances are ranked row-wise once. At
##  each iteration every medoid gathers its m nearest objects; within that
##  group the object with the largest hostility value (row sum of the
##  group's ranks) becomes the new medoid. Medoids that collide are merged
##  and the missing ones are relocated at random to objects outside every
##  group, which is the stochastic step of the algorithm.
## ---------------------------------------------------------------------------

#' Shape-rankkmed: Ranked K-medoids on Kendall's 2D shape space
#'
#' @param Z Complex pre-shape matrix of size (k-1) x n, or a numeric
#'   distance matrix of size n x n.
#' @param K Number of clusters (>= 2).
#' @param distance Character. One of "rho" (default), "dF".
#' @param m Neighbourhood size for cluster formation (default 5;
#'   Zadegan et al. (2013) recommend 5-15).
#' @param max_iter Integer. Maximum number of iterations (default 30). The
#'   loop stops earlier at a fixed point (medoids unchanged and no
#'   relocation), from which further iterations cannot move.
#' @param initial Optional vector of length K giving initial medoid
#'   indices (1-based).
#'
#' @details
#' With the same distance matrix, initial medoids and random-number
#' stream, the partition equals that of \code{kmed::rankkmed} run with
#' \code{iterate = max_iter}.
#'
#' @references
#' Zadegan, S. M. R., Mirzaie, M. and Sadoughi, F. (2013). Ranked
#' k-medoids: a fast and accurate rank-based partitioning algorithm for
#' clustering large datasets. \emph{Knowledge-Based Systems}, 39, 133--143.
#' \doi{10.1016/j.knosys.2012.10.012}
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
  if (m < 2L || m > n) stop("m must satisfy 2 <= m <= n.")

  # Row-wise ranks (ties averaged, as in base::rank) and the matching
  # index matrix SI: SI[i, ] lists objects from nearest to farthest of i.
  R  <- rank_matrix_cpp(D)
  SI <- order_matrix_cpp(D)

  if (is.null(initial)) {
    medoids <- sort(sample.int(n, K))
  } else {
    stopifnot(length(unique(initial)) == K, all(initial %in% seq_len(n)))
    medoids <- as.integer(initial)
  }

  n_iter <- 0L
  for (it in seq_len(max_iter)) {
    n_iter <- it
    groups <- lapply(medoids, function(md) SI[md, seq_len(m)])
    # Hostility: row sums of the group's ranks; the largest wins.
    cand <- vapply(groups, function(g)
      g[which.max(rowSums(R[g, g, drop = FALSE]))], integer(1))
    in_groups <- unlist(groups)
    hits <- vapply(cand, function(x) sum(in_groups == x), integer(1))
    kept <- c(cand[hits == 1L], unique(cand[hits > 1L]))
    relocated <- length(kept) != K
    if (relocated) {
      pool <- setdiff(seq_len(n), in_groups)
      if (length(pool) < K - length(kept)) pool <- setdiff(seq_len(n), kept)
      kept <- c(kept, pool[sample.int(length(pool), K - length(kept))])
    }
    new_medoids <- sort(as.integer(kept))
    done <- !relocated && identical(new_medoids, as.integer(medoids))
    medoids <- new_medoids
    if (done) break
  }

  cluster <- max.col(-D[, medoids, drop = FALSE], ties.method = "first")
  obj <- sum(D[cbind(seq_len(n), medoids[cluster])])

  list(cluster    = as.integer(cluster),
       medoids    = as.integer(medoids),
       objective  = obj,
       n_iter     = as.integer(n_iter),
       distance   = distance)
}
