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
#' Riemannian distance rho on Sigma_2^k. Group means are Frechet means
#' computed via the shapes package if available, otherwise approximated
#' by the observation with smallest sum of squared rho-distances to the
#' others within its group.
#'
#' @param X Landmark array (k x 2 x n).
#' @param cluster Integer vector of labels in 1..K, length n.
#'
#' @return A list with VEG, VIG, and the ratio VEG/VIG.
#'
#' @export
riemannian_variances <- function(X, cluster) {
  stopifnot(length(dim(X)) == 3L, dim(X)[2L] == 2L)
  n <- dim(X)[3L]
  stopifnot(length(cluster) == n)
  K <- length(unique(cluster))

  Z <- landmarks_to_preshape(X)
  # Overall Frechet mean approximation: minimiser over the sample
  D_all <- rho_matrix(Z)
  overall_idx <- which.min(rowSums(D_all^2))
  z_overall <- Z[, overall_idx, drop = FALSE]

  # Per-group Frechet mean approximations
  z_groups <- vector("list", K)
  for (r in seq_len(K)) {
    idx <- which(cluster == r)
    Zg <- Z[, idx, drop = FALSE]
    Dg <- rho_matrix(Zg)
    z_groups[[r]] <- Zg[, which.min(rowSums(Dg^2)), drop = FALSE]
  }

  # VEG: sum_i rho(mu_i, mu)^2 / (K - 1)
  VEG <- 0
  for (r in seq_len(K)) {
    VEG <- VEG + rho_pair(z_groups[[r]], z_overall)^2
  }
  VEG <- VEG / max(1L, K - 1L)

  # VIG: sum_i sum_j rho(X_ij, mu_i)^2 / (n - K)
  VIG <- 0
  for (r in seq_len(K)) {
    idx <- which(cluster == r)
    for (j in idx) {
      VIG <- VIG + rho_pair(Z[, j, drop = FALSE], z_groups[[r]])^2
    }
  }
  VIG <- VIG / max(1L, n - K)

  list(VEG = VEG, VIG = VIG, ratio = VEG / max(VIG, .Machine$double.eps))
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
