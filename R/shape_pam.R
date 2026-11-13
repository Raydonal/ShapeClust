## ---------------------------------------------------------------------------
##  Shape-PAM: Partitioning Around Medoids on Kendall's 2D shape space.
##
##  The algorithm is the standard PAM of Rousseeuw and Kaufman (1990) with
##  the four swap-cost cases of Ng and Han (2002), computed with the
##  Riemannian distance rho on the pre-shape sphere. The C++ backend
##  performs the swap-cost evaluation and the nearest-medoid assignment,
##  which are the two hot paths.
## ---------------------------------------------------------------------------

#' Shape-PAM: PAM on Kendall's 2D shape space
#'
#' @param Z Complex pre-shape matrix of size (k-1) x n produced by
#'   \code{\link{landmarks_to_preshape}}, OR a numeric distance matrix
#'   of size n x n, OR a real landmark array of dimension k x m x n with
#'   m = 2 or 3 (distances from \code{\link{shape_distance}}). If a
#'   distance matrix is supplied, it is used directly.
#' @param K Number of clusters (>= 2).
#' @param distance Character. One of "rho" (Riemannian; default),
#'   "dF" (full Procrustes). Ignored if \code{Z} is already a distance
#'   matrix.
#' @param initial Optional vector of length K giving initial medoid
#'   indices (1-based). If NULL, medoids are sampled uniformly at
#'   random from 1:n.
#' @param max_iter Integer. Maximum number of swap iterations
#'   (default 100).
#'
#' @return A list with components:
#'   \describe{
#'     \item{cluster}{Integer vector of length n, values in 1:K.}
#'     \item{medoids}{Integer vector of length K.}
#'     \item{objective}{Sum of distances of each object to its medoid.}
#'     \item{n_iter}{Number of accepted swaps performed.}
#'     \item{distance}{Distance name used.}
#'   }
#'
#' @examples
#' set.seed(1)
#' k <- 4; n <- 30
#' X <- array(rnorm(k * 2 * n), dim = c(k, 2, n))
#' Z <- landmarks_to_preshape(X)
#' fit <- shape_pam(Z, K = 3)
#' table(fit$cluster)
#'
#' @export
shape_pam <- function(Z, K,
                      distance = c("rho", "dF"),
                      initial = NULL,
                      max_iter = 100L) {
  distance <- match.arg(distance)

  # Build distance matrix if needed
  if (is.matrix(Z) && !is.complex(Z)) {
    D <- Z  # user passed a distance matrix directly
  } else if (is.array(Z) && length(dim(Z)) == 3L && !is.complex(Z)) {
    D <- shape_distance(Z, distance)  # real k x m x n landmark array (m = 2, 3)
  } else {
    D <- switch(distance,
                rho = rho_matrix(Z),
                dF  = dF_matrix(Z))
  }
  n <- nrow(D)
  if (K < 2L || K > n - 1L) stop("K must satisfy 2 <= K <= n-1.")

  # Initial medoids (1-based in R)
  if (is.null(initial)) {
    medoids <- sort(sample.int(n, K))
  } else {
    stopifnot(length(initial) == K, all(initial %in% 1:n))
    medoids <- as.integer(initial)
  }

  # Full sweep in C++ (assignment + best-swap search), falling back to the
  # pure-R implementation when the compiled backend is unavailable.
  if (has_cpp_backend()) {
    out <- shape_pam_cpp(D, as.integer(medoids - 1L), as.integer(max_iter))
    medoids <- as.integer(out$medoids) + 1L
    cluster <- as.integer(out$cluster) + 1L
    return(list(cluster   = cluster,
                medoids   = medoids,
                objective = out$objective,
                n_iter    = as.integer(out$n_iter),
                distance  = distance))
  }

  n_iter <- 0L
  repeat {
    NS <- .nearest_two(D, medoids)
    nearest <- NS[, 1L]; second <- NS[, 2L]
    best_TC <- 0; best_s <- 0L; best_h <- 0L
    for (s in seq_len(K)) {
      for (h in seq_len(n)) {
        if (h %in% medoids) next
        TCsh <- .swap_cost(D, medoids, nearest, second, s - 1L, h - 1L)
        if (TCsh < best_TC) { best_TC <- TCsh; best_s <- s; best_h <- h }
      }
    }
    if (best_TC < 0 && n_iter < max_iter) {
      medoids[best_s] <- best_h
      n_iter <- n_iter + 1L
    } else break
  }
  NS <- .nearest_two(D, medoids)
  cluster <- NS[, 1L] + 1L
  obj <- sum(D[cbind(seq_len(n), medoids[cluster])])

  list(cluster    = cluster,
       medoids    = as.integer(medoids),
       objective  = obj,
       n_iter     = n_iter,
       distance   = distance)
}

# Internal wrappers over the C++ backend, with pure-R fallback for CI or
# environments without a working compiler.

.nearest_two <- function(D, medoids) {
  if (has_cpp_backend()) {
    return(pam_nearest_two_cpp(D, as.integer(medoids - 1L)))
  }
  n <- nrow(D)
  NS <- matrix(0L, nrow = n, ncol = 2L)
  for (j in seq_len(n)) {
    d <- D[j, medoids]
    ord <- order(d)
    NS[j, 1L] <- ord[1L] - 1L
    NS[j, 2L] <- ord[2L] - 1L
  }
  NS
}

.swap_cost <- function(D, medoids, nearest, second, s, h) {
  if (has_cpp_backend()) {
    return(pam_swap_cost_cpp(D,
                             as.integer(medoids - 1L),
                             as.integer(nearest),
                             as.integer(second),
                             as.integer(s),
                             as.integer(h)))
  }
  # Pure-R fallback: mirror the C++ logic exactly
  n <- nrow(D)
  os_idx <- medoids[s + 1L]  # 1-based medoid
  TCsh <- 0
  for (j in seq_len(n)) {
    j0 <- j - 1L
    if (j0 == os_idx - 1L || j0 == h) next
    oi_pos <- nearest[j]
    d_j_h  <- D[j, h + 1L]
    d_j_os <- D[j, os_idx]
    if (oi_pos == s) {
      oi2 <- medoids[second[j] + 1L]
      d_j_oi <- D[j, oi2]
      TCsh <- TCsh + if (d_j_h >= d_j_oi) d_j_oi - d_j_os else d_j_h - d_j_os
    } else {
      oi_med <- medoids[oi_pos + 1L]
      d_j_oi <- D[j, oi_med]
      if (d_j_h < d_j_oi) TCsh <- TCsh + (d_j_h - d_j_oi)
    }
  }
  TCsh
}
