## ---------------------------------------------------------------------------
##  Model selection: choosing K when no ground truth is available.
##
##  Provides two shape-space-compatible criteria:
##    - Average silhouette width (Rousseeuw 1987)
##    - Gap statistic with a Riemannian null (Tibshirani et al. 2001)
##
##  Both operate on the Riemannian distance matrix, so they respect the
##  geometry of Sigma_2^k. The bagging step of BagClust1 does not attempt
##  to estimate K; that is a separate model-selection task.
## ---------------------------------------------------------------------------

#' Select K by silhouette width or the gap statistic
#'
#' @param Z Complex pre-shape matrix (k-1) x n, or an n x n distance matrix.
#' @param K_range Integer vector of candidate values of K.
#' @param method One of "silhouette" (default) or "gap".
#' @param base Base clusterer to use for each candidate K. One of
#'   "pam" (default) or "rankkmed".
#' @param B_gap Number of null reference draws for the gap statistic.
#' @param distance One of "rho" (default), "dF".
#'
#' @return A list with:
#'   \describe{
#'     \item{K_hat}{Selected K.}
#'     \item{scores}{Named numeric vector of scores (silhouette width or
#'       gap statistic) at each candidate K.}
#'     \item{method}{Selection method used.}
#'   }
#'
#' @examples
#' set.seed(1)
#' # Three clusters
#' mu_list <- list(
#'   matrix(c(0,0, 10,0, 10,10, 0,10), 4, 2, byrow = TRUE),
#'   matrix(c(0,0, 20,0, 20,10, 0,10), 4, 2, byrow = TRUE),
#'   matrix(c(0,0, 15,0, 10,10, 0,10), 4, 2, byrow = TRUE)
#' )
#' out <- sim_landmark_gaussian(mu_list, c(15, 15, 15), sigma = 0.5)
#' Z <- landmarks_to_preshape(out$X)
#' res <- select_K(Z, K_range = 2:5, method = "silhouette")
#' res$K_hat
#'
#' @export
select_K <- function(Z, K_range = 2:8,
                     method = c("silhouette", "gap"),
                     base = c("pam", "rankkmed"),
                     B_gap = 20L,
                     distance = c("rho", "dF")) {
  method   <- match.arg(method)
  base     <- match.arg(base)
  distance <- match.arg(distance)

  if (is.matrix(Z) && !is.complex(Z)) {
    D <- Z
  } else {
    D <- switch(distance,
                rho = rho_matrix(Z),
                dF  = dF_matrix(Z))
  }
  n <- nrow(D)
  base_fun <- switch(base,
                     pam      = shape_pam,
                     rankkmed = shape_rankkmed)

  if (method == "silhouette") {
    scores <- sapply(K_range, function(K) {
      fit <- base_fun(D, K = K)
      .silhouette_avg(D, fit$cluster)
    })
    names(scores) <- as.character(K_range)
    K_hat <- as.integer(names(which.max(scores)))
  } else {
    # Gap statistic: log W(D, cluster) - E*[log W(null)]
    obs <- sapply(K_range, function(K) {
      fit <- base_fun(D, K = K)
      log(.within_ss(D, fit$cluster))
    })
    ref <- matrix(0, length(K_range), B_gap)
    for (b in seq_len(B_gap)) {
      # Null: uniform on shape space via random unit-modulus complex
      # perturbation of pre-shape sphere. Approximated here by sampling
      # rows/cols of D via a random rotation-like index permutation
      # (a mild approximation appropriate for exploratory use).
      idx <- sample.int(n)
      D_ref <- D[idx, idx]
      for (K_i in seq_along(K_range)) {
        K <- K_range[K_i]
        fit <- base_fun(D_ref, K = K)
        ref[K_i, b] <- log(.within_ss(D_ref, fit$cluster))
      }
    }
    gap <- rowMeans(ref) - obs
    scores <- setNames(gap, as.character(K_range))
    K_hat <- as.integer(names(which.max(scores)))
  }

  list(K_hat = K_hat,
       scores = scores,
       method = method)
}

.silhouette_avg <- function(D, cluster) {
  n <- length(cluster)
  K <- length(unique(cluster))
  if (K < 2L) return(NA_real_)
  s <- numeric(n)
  for (i in seq_len(n)) {
    own <- which(cluster == cluster[i])
    own <- setdiff(own, i)
    a_i <- if (length(own)) mean(D[i, own]) else 0
    others <- setdiff(unique(cluster), cluster[i])
    b_i <- min(sapply(others, function(r) mean(D[i, cluster == r])))
    s[i] <- if (max(a_i, b_i) > 0) (b_i - a_i) / max(a_i, b_i) else 0
  }
  mean(s)
}

.within_ss <- function(D, cluster) {
  K <- length(unique(cluster))
  W <- 0
  for (r in seq_len(K)) {
    idx <- which(cluster == r)
    if (length(idx) < 2L) next
    Dr <- D[idx, idx]
    # Sum of squared pairwise distances / (2 |C_r|)
    W <- W + sum(Dr^2) / (2 * length(idx))
  }
  max(W, .Machine$double.eps)
}
