## ---------------------------------------------------------------------------
##  Validation indices: Rand Index (RI), Fowlkes-Mallows Index (FMI),
##  Relative Gain, and shape-space between/within variances (VEG, VIG).
## ---------------------------------------------------------------------------

#' Rand Index
#'
#' Rand (1971).
#'
#' @param truth Integer vector of ground-truth labels.
#' @param pred  Integer vector of predicted labels, same length.
#'
#' @return A numeric in [0, 1].
#'
#' @examples
#' rand_index(c(1,1,2,2,3,3), c(1,1,2,2,3,3))
#' rand_index(c(1,1,2,2,3,3), c(1,2,2,1,3,3))
#'
#' @export
rand_index <- function(truth, pred) {
  stopifnot(length(truth) == length(pred))
  n <- length(truth)
  # Compute contingency table
  tab <- table(truth, pred)
  a <- sum(choose(tab, 2))
  b <- choose(sum(tab), 2) - sum(choose(rowSums(tab), 2)) -
       sum(choose(colSums(tab), 2)) + a
  (a + b) / choose(n, 2)
}

#' Fowlkes-Mallows Index
#'
#' Fowlkes and Mallows (1983).
#'
#' @inheritParams rand_index
#'
#' @return A numeric in [0, 1].
#'
#' @examples
#' fmi_index(c(1,1,2,2,3,3), c(1,1,2,2,3,3))
#' fmi_index(c(1,1,2,2,3,3), c(1,2,2,1,3,3))
#'
#' @export
fmi_index <- function(truth, pred) {
  stopifnot(length(truth) == length(pred))
  tab <- table(truth, pred)
  TP <- sum(choose(tab, 2))
  Prow <- sum(choose(rowSums(tab), 2))
  Pcol <- sum(choose(colSums(tab), 2))
  if (Prow == 0 || Pcol == 0) return(0)
  TP / sqrt(Prow * Pcol)
}

#' Relative Gain
#'
#' Percentage increase of a validation measure from a standard to a
#' bagged variant.
#'
#' @param V_std Numeric. Value of the validation index for the standard variant.
#' @param V_bag Numeric. Value for the bagged variant.
#'
#' @return Numeric percentage.
#'
#' @examples
#' relative_gain(0.80, 0.88)  # +10 percent
#'
#' @export
relative_gain <- function(V_std, V_bag) {
  stopifnot(is.numeric(V_std), is.numeric(V_bag))
  if (any(V_std == 0)) {
    return(ifelse(V_bag == 0, 0, Inf))
  }
  100 * (V_bag - V_std) / V_std
}

#' Riemannian between- and within-group variances (VEG, VIG)
#'
#' Estimates \code{VEG = sum_i rho(mu_i, mu)^2 / (K-1)} and
#' \code{VIG = sum_i sum_j rho(X_ij, mu_i)^2 / (n-K)} using the
#' Riemannian distance rho on Kendall's shape space.
#'
#' The group and overall means are full Procrustes means: the dominant
#' eigenvector of the complex sum of squares \eqn{\sum_j z_j z_j^*} of the
#' pre-shapes (Kent, 1994). This is the Frechet mean under the full
#' Procrustes distance.
#'
#' @param X Landmark array (k x 2 x n).
#' @param cluster Vector of group labels, length n.
#'
#' @return A list with VEG, VIG, and the ratio VEG/VIG.
#'
#' @references
#' Kent, J. T. (1994). The complex Bingham distribution and shape
#' analysis. \emph{Journal of the Royal Statistical Society B}, 56(2),
#' 285--299. \doi{10.1111/j.2517-6161.1994.tb01978.x}
#'
#' @examples
#' mu_list <- list(
#'   matrix(c(0,0, 10,0, 10,10, 0,10), 4, 2, byrow = TRUE),
#'   matrix(c(0,0, 20,0, 20,10, 0,10), 4, 2, byrow = TRUE)
#' )
#' out <- sim_landmark_gaussian(mu_list, c(15, 15), sigma = 1, seed = 1)
#' riemannian_variances(out$X, out$group)
#'
#' @export
riemannian_variances <- function(X, cluster) {
  stopifnot(length(dim(X)) == 3L, dim(X)[2L] == 2L)
  n <- dim(X)[3L]
  stopifnot(length(cluster) == n)
  groups <- unique(cluster)
  K <- length(groups)

  Z <- landmarks_to_preshape(X)
  mu <- .procrustes_mean(Z)
  mu_g <- lapply(groups, function(r) .procrustes_mean(Z[, cluster == r, drop = FALSE]))

  VEG <- sum(vapply(mu_g, function(m) rho_pair(m, mu)^2, numeric(1))) /
    max(1L, K - 1L)
  VIG <- 0
  for (r in seq_len(K)) {
    for (j in which(cluster == groups[r])) {
      VIG <- VIG + rho_pair(Z[, j], mu_g[[r]])^2
    }
  }
  VIG <- VIG / max(1L, n - K)

  list(VEG = VEG, VIG = VIG, ratio = VEG / max(VIG, .Machine$double.eps))
}

# Full Procrustes mean of the pre-shape columns of Z (Kent 1994)
.procrustes_mean <- function(Z) {
  S <- Z %*% Conj(t(Z))
  eigen(S, symmetric = TRUE)$vectors[, 1L]
}

# Internal utility: rho between two single pre-shape columns
rho_pair <- function(z1, z2) {
  cval <- min(1, Mod(as.complex(sum(Conj(z1) * z2))))
  acos(cval)
}

#' Paired one-sided Wilcoxon signed-rank test
#'
#' Convenience wrapper for the paired one-sided test of
#' \code{median(V_bag - V_std) > 0}.
#'
#' @param V_std Numeric vector of standard-variant values.
#' @param V_bag Numeric vector of bagged-variant values, same length.
#'
#' @return A numeric p-value.
#'
#' @examples
#' set.seed(1)
#' v_std <- runif(30, 0.5, 0.8)
#' v_bag <- v_std + rnorm(30, 0.05, 0.02)
#' wilcoxon_paired_gain(v_std, v_bag)
#'
#' @export
wilcoxon_paired_gain <- function(V_std, V_bag) {
  stopifnot(length(V_std) == length(V_bag))
  res <- suppressWarnings(
    stats::wilcox.test(V_bag, V_std,
                       alternative = "greater",
                       paired = TRUE, exact = FALSE)
  )
  as.numeric(res$p.value)
}
