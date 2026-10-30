## ---------------------------------------------------------------------------
##  Model selection: choosing K when no ground truth is available.
##
##  Provides two shape-space-compatible criteria:
##    - Average silhouette width (Rousseeuw 1987)
##    - Gap statistic (Tibshirani et al. 2001) with the uniform distribution
##      on the shape space as the null reference
##
##  Both operate on the Riemannian distance matrix, so they respect the
##  geometry of Sigma_2^k. The bagging step of BagClust1 does not attempt
##  to estimate K; that is a separate model-selection task.
## ---------------------------------------------------------------------------

#' Select K by silhouette width or the gap statistic
#'
#' @param Z Complex pre-shape matrix (k-1) x n, or an n x n distance matrix
#'   (silhouette only: the gap statistic needs the pre-shapes to draw its
#'   null reference).
#' @param K_range Integer vector of candidate values of K.
#' @param method One of "silhouette" (default) or "gap".
#' @param base Base clusterer to use for each candidate K. One of
#'   "pam" (default) or "rankkmed".
#' @param B_gap Number of null reference samples for the gap statistic.
#' @param distance One of "rho" (default), "dF".
#'
#' @return A list with:
#'   \describe{
#'     \item{K_hat}{Selected K.}
#'     \item{scores}{Named numeric vector of scores (silhouette width or
#'       gap statistic) at each candidate K.}
#'     \item{se}{For \code{"gap"}, the simulation standard errors
#'       \eqn{s_K}; \code{NULL} otherwise.}
#'     \item{method}{Selection method used.}
#'   }
#'
#' @details
#' The gap statistic compares \eqn{\log W_K} with its mean under
#' \code{B_gap} samples of the same size drawn uniformly on the shape
#' space (uniform pre-shapes, i.e. a Complex Bingham with zero
#' concentration). The selected K is the smallest one with
#' \eqn{\mathrm{Gap}(K) \ge \mathrm{Gap}(K') - s_{K'}}, where \eqn{K'} is
#' the next candidate (Tibshirani, Walther and Hastie, 2001).
#'
#' @references
#' Rousseeuw, P. J. (1987). Silhouettes: a graphical aid to the
#' interpretation and validation of cluster analysis. \emph{Journal of
#' Computational and Applied Mathematics}, 20, 53--65.
#' \doi{10.1016/0377-0427(87)90125-7}
#'
#' Tibshirani, R., Walther, G. and Hastie, T. (2001). Estimating the number
#' of clusters in a data set via the gap statistic. \emph{Journal of the
#' Royal Statistical Society B}, 63(2), 411--423.
#' \doi{10.1111/1467-9868.00293}
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

  dist_fun <- switch(distance, rho = rho_matrix, dF = dF_matrix)
  if (is.matrix(Z) && !is.complex(Z)) {
    if (method == "gap")
      stop("The gap statistic needs the pre-shape matrix Z, not distances.")
    D <- Z
  } else {
    D <- dist_fun(Z)
  }
  n <- nrow(D)
  base_fun <- switch(base,
                     pam      = shape_pam,
                     rankkmed = shape_rankkmed)
  se <- NULL

  if (method == "silhouette") {
    scores <- vapply(K_range, function(K) {
      .silhouette_avg(D, base_fun(D, K = K)$cluster)
    }, numeric(1))
    names(scores) <- as.character(K_range)
    K_hat <- as.integer(names(which.max(scores)))
  } else {
    logW <- function(D) vapply(K_range, function(K) {
      log(.within_ss(D, base_fun(D, K = K)$cluster))
    }, numeric(1))
    obs <- logW(D)
    d <- nrow(Z)
    ref <- vapply(seq_len(B_gap), function(b) {
      Zr <- matrix(complex(real = stats::rnorm(d * n),
                           imaginary = stats::rnorm(d * n)), d, n)
      Zr <- sweep(Zr, 2L, sqrt(colSums(Mod(Zr)^2)), FUN = "/")
      logW(dist_fun(Zr))
    }, numeric(length(K_range)))
    ref <- matrix(ref, nrow = length(K_range))
    scores <- setNames(rowMeans(ref) - obs, as.character(K_range))
    se <- setNames(apply(ref, 1L, stats::sd) * sqrt(1 + 1 / B_gap),
                   as.character(K_range))
    nK <- length(K_range)
    ok <- if (nK > 1L) which(scores[-nK] >= scores[-1L] - se[-1L]) else 1L
    K_hat <- if (length(ok)) K_range[ok[1L]] else K_range[nK]
  }

  list(K_hat = as.integer(K_hat),
       scores = scores,
       se = se,
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
  W <- 0
  for (r in unique(cluster)) {
    idx <- which(cluster == r)
    if (length(idx) < 2L) next
    Dr <- D[idx, idx]
    # Sum of squared pairwise distances / (2 |C_r|)
    W <- W + sum(Dr^2) / (2 * length(idx))
  }
  max(W, .Machine$double.eps)
}
