## ---------------------------------------------------------------------------
##  Simulators for the Monte Carlo study
##  - Complex Bingham (Kent 1994; sampler of Kent, Constable, Er 2004)
##  - Isotropic and anisotropic Gaussian perturbations of a template
##    landmark configuration (Preston & Wood 2010; Amaral et al. 2007)
## ---------------------------------------------------------------------------

#' Simulate a single Complex Bingham sample
#'
#' Uses the Simplex-truncation method of Kent, Constable and Er (2004).
#'
#' @param n Sample size.
#' @param lambda Numeric vector of concentration parameters
#'   \eqn{(\lambda_1, \dots, \lambda_{k-1})}; note the last eigenvalue
#'   is fixed at 0 by convention.
#' @param seed Optional integer seed; the caller's random-number stream
#'   is left as it was.
#'
#' @return A complex matrix of size (k-1) x n whose columns are
#'   unit-norm Complex Bingham draws.
#'
#' @examples
#' Z <- sim_complex_bingham(n = 5, lambda = rep(9, 3))
#' round(colSums(Mod(Z)^2), 12)  # unit norm
#'
#' @export
sim_complex_bingham <- function(n, lambda, seed = NULL) {
  .local_seed(seed)
  ## A pre-shape lives on the unit sphere of C^d with d = k - 1. The Complex
  ## Bingham has d eigenvalues, the largest of which is fixed at 0 by
  ## convention and marks the *pole* of the distribution. The caller supplies
  ## the d - 1 free concentrations in `lambda`, so the sampler returns
  ## d = length(lambda) + 1 complex coordinates, with the pole in the last
  ## position. The pole coordinate carries the concentration and must NOT be
  ## dropped.
  d <- length(lambda) + 1L
  Z <- matrix(complex(0), nrow = d, ncol = n)
  for (i in seq_len(n)) Z[, i] <- .cbingham_one(lambda)
  norms <- sqrt(colSums(Mod(Z)^2))
  norms[norms == 0] <- 1
  sweep(Z, 2L, norms, FUN = "/")
}

.cbingham_one <- function(lambda) {
  ## The squared moduli S = (|z_1|^2, ..., |z_k|^2) lie on the unit simplex
  ## with density proportional to exp(-sum_{j<k} lambda_j S_j); the kth
  ## eigenvalue is 0 by convention, so the kth coordinate is the pole and is
  ## obtained as the remainder S_k = 1 - sum_{j<k} S_j.
  ##
  ## Two exact samplers are combined, since neither is efficient over the
  ## whole parameter range:
  ##
  ##  (a) Simplex truncation (Kent, Constable and Er, 2004): draw
  ##      S_j ~ Texp(lambda_j) on [0, 1], j < k, and accept if their sum is
  ##      below 1. Efficient when the lambda_j are large (mass near the
  ##      pole), i.e. when E[sum_j S_j] is well below 1.
  ##
  ##  (b) Uniform-simplex rejection: draw S uniformly on the simplex and
  ##      accept with probability exp(-sum_j lambda_j S_j), which is bounded
  ##      by 1 and hence valid. Efficient when the lambda_j are small, which
  ##      is exactly where (a) rejects almost surely.
  ##
  ## Both target the same density, so the choice affects only cost.
  d <- length(lambda)          # number of non-pole coordinates, = k - 1
  k <- d + 1L

  texp_mean <- function(l) {   # mean of Texp(l) truncated to [0, 1]
    if (abs(l) < 1e-12) return(0.5)
    1 / l - exp(-l) / (1 - exp(-l))
  }
  use_kent <- sum(vapply(lambda, texp_mean, numeric(1))) < 0.8

  S <- NULL
  if (use_kent) {
    for (try in seq_len(20000L)) {
      U <- runif(d)
      Sp <- vapply(seq_len(d), function(j) {
        lj <- lambda[j]
        if (abs(lj) < 1e-12) U[j]
        else -(1 / lj) * log(1 - U[j] * (1 - exp(-lj)))
      }, numeric(1))
      if (sum(Sp) < 1) { S <- c(Sp, 1 - sum(Sp)); break }
    }
  }
  if (is.null(S)) {
    lam_full <- c(lambda, 0)
    repeat {                                  # batched for efficiency
      m <- 2000L
      E <- matrix(rexp(m * k), nrow = m)
      Sp <- E / rowSums(E)                    # uniform on the simplex
      acc <- exp(-as.vector(Sp %*% lam_full))
      hit <- which(runif(m) < acc)
      if (length(hit)) { S <- Sp[hit[1L], ]; break }
    }
  }
  S <- pmax(S, 0); S <- S / sum(S)
  theta <- runif(k, 0, 2 * pi)
  sqrt(S) * complex(real = cos(theta), imaginary = sin(theta))
}

#' Simulate landmark configurations from a Gaussian perturbation
#'
#' Perturbs a template configuration by a multivariate Gaussian in
#' landmark space, following Amaral et al. (2007) and Preston & Wood
#' (2010).
#'
#' @param mu_list List of template configurations, each a k x 2 matrix.
#'   Length = number of groups G.
#' @param n_per_group Integer vector of length G, sample sizes per group.
#' @param sigma Standard deviation of the perturbation.
#' @param regime One of "isotropic" (default) or "anisotropic".
#' @param gamma Anisotropy parameter (default 4).
#' @param seed Optional integer seed; the caller's random-number stream
#'   is left as it was.
#'
#' @return A list with:
#'   \describe{
#'     \item{X}{Numeric array of shape (k, 2, n_total).}
#'     \item{group}{Integer vector of group labels in 1..G.}
#'   }
#'
#' @examples
#' mu_list <- list(
#'   matrix(c(0,0, 10,0, 10,10, 0,10), 4, 2, byrow = TRUE),
#'   matrix(c(0,0, 20,0, 20,10, 0,10), 4, 2, byrow = TRUE),
#'   matrix(c(0,0, 15,0, 10,10, 0,10), 4, 2, byrow = TRUE)
#' )
#' out <- sim_landmark_gaussian(mu_list, c(10, 10, 10), sigma = 1)
#' dim(out$X); table(out$group)
#'
#' @export
sim_landmark_gaussian <- function(mu_list, n_per_group,
                                  sigma = 1,
                                  regime = c("isotropic", "anisotropic"),
                                  gamma = 4,
                                  seed = NULL) {
  regime <- match.arg(regime)
  .local_seed(seed)
  G <- length(mu_list)
  stopifnot(length(n_per_group) == G, all(sapply(mu_list, ncol) == 2L))
  k <- nrow(mu_list[[1L]])
  m <- 2L  # 2D
  n_total <- sum(n_per_group)
  X <- array(0, dim = c(k, 2, n_total))
  group <- integer(n_total)

  if (regime == "isotropic") {
    Sigma <- diag(sigma^2, k * m)
  } else {
    # Kronecker-product covariance (Preston & Wood 2010)
    A_k <- matrix(1, k, k) + (gamma - 1) * diag(k)
    A_m <- matrix(1, m, m) + (gamma - 1) * diag(m)
    Sigma <- (sigma^2 / gamma^2) * kronecker(A_k, A_m)
  }

  # Cholesky once
  L <- chol(Sigma)
  offset <- 0L
  for (g in seq_len(G)) {
    ng <- n_per_group[g]
    mu_vec <- as.vector(t(mu_list[[g]]))  # row-wise = (x1,y1,x2,y2,...)
    for (i in seq_len(ng)) {
      eps <- t(L) %*% rnorm(k * m)
      vec <- mu_vec + as.numeric(eps)
      # Back to k x 2 matrix
      Xg <- matrix(vec, ncol = 2, byrow = TRUE)
      X[, , offset + i] <- Xg
      group[offset + i] <- g
    }
    offset <- offset + ng
  }
  list(X = X, group = group)
}
