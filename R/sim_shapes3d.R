## ---------------------------------------------------------------------------
##  Simulator for the 3D study: noisy cubes and parallelepipeds. Port of
##  cubos.generator.iso/aniso() and parallep.generator.iso/aniso() from the
##  3D bagging scripts of Nascimento et al.
## ---------------------------------------------------------------------------

# Template landmark configurations. The 8-landmark cube and parallelepiped
# (box with one edge twice as long) are the vertices; the 34-landmark ones
# add edge and face points. Coordinates are those of the data sets
# cube8landm, parallelep8landm, cube34landm and parallelep34landm of the
# 'Anthropometry' package, which the scripts use.
.shape_template <- function(name, k) {
  v <- switch(paste(name, k),
    "cube 8" = c(0, 0, 0, 0, 10, 10, 10, 10,
                 0, 10, 0, 10, 10, 10, 0, 0,
                 10, 10, 0, 0, 10, 0, 10, 0),
    "parallelepiped 8" = c(0, 10, 0, 10, 10, 10, 0, 0,
                           0, 0, 0, 0, 20, 20, 20, 20,
                           10, 10, 0, 0, 10, 0, 10, 0),
    "cube 34" = c(10, 5, 10, 10, 10, 10, 10, 10, 10, 10, 10, 5, 0, 5, 0, 0, 0,
                  0, 0, 10, 5, 5, 5, 5, 0, 0, 0, 0, 0, 0, 5, 5, 5, 10,
                  0, 0, 0, 6.7, 10, 0, 3.3, 6.7, 3.3, 6.7, 10, 10, 10, 10, 10,
                  0, 6.7, 3.3, 0, 3.3, 10, 6.7, 3.3, 0, 6.7, 3.3, 10, 6.7, 3.3,
                  0, 6.7, 3.3, 0, 10,
                  10, 5, 0, 10, 10, 5, 5, 5, 0, 0, 5, 5, 5, 0, 10, 10, 5, 5, 5,
                  10, 10, 10, 10, 10, 10, 10, 0, 0, 0, 0, 0, 0, 0, 0),
    "parallelepiped 34" = c(10, 5, 10, 10, 10, 10, 10, 10, 10, 10, 10, 5, 0, 5,
                            0, 0, 0, 0, 0, 10, 5, 5, 5, 5, 0, 0, 0, 0, 0, 0, 5,
                            5, 5, 10,
                            0, 0, 0, 13, 20, 0, 6.7, 13, 6.7, 13, 20, 20, 20,
                            20, 20, 0, 13, 6.7, 0, 6.7, 20, 13, 6.7, 0, 13, 6.7,
                            20, 13, 6.7, 0, 13, 6.7, 0, 20,
                            10, 5, 0, 10, 10, 5, 5, 5, 0, 0, 5, 5, 5, 0, 10, 10,
                            5, 5, 5, 10, 10, 10, 10, 10, 10, 10, 0, 0, 0, 0, 0,
                            0, 0, 0),
    stop("no built-in template for '", name, "' with k = ", k, ".", call. = FALSE))
  matrix(v, nrow = k, ncol = 3L)
}

# Multivariate t draws, scale matrix `sigma`, centered at 0: the same
# computation, in the same order of random draws, as
# mvtnorm::rmvt(n, sigma, df) of mvtnorm >= 0.9-9994 (default method
# "eigen", type "shifted").
.rmvt_shifted <- function(n, sigma, df) {
  ev <- eigen(sigma, symmetric = TRUE)
  R <- t(ev$vectors %*% (t(ev$vectors) * sqrt(pmax(ev$values, 0))))
  z <- matrix(rnorm(n * ncol(sigma)), nrow = n, byrow = TRUE) %*% R
  z / sqrt(stats::rchisq(n, df) / df)
}

# von Mises draws by the algorithm of Best and Fisher (1979), as CircStats::rvm.
.rvm <- function(n, mean, k) {
  vm <- numeric(n)
  a <- 1 + (1 + 4 * (k^2))^0.5
  b <- (a - (2 * a)^0.5) / (2 * k)
  r <- (1 + b^2) / (2 * b)
  obs <- 1L
  while (obs <= n) {
    U1 <- runif(1, 0, 1)
    z <- cos(pi * U1)
    f <- (1 + r * z) / (r + z)
    c <- k * (r - f)
    U2 <- runif(1, 0, 1)
    if (c * (2 - c) - U2 > 0) {
      U3 <- runif(1, 0, 1)
      vm[obs] <- (sign(U3 - 0.5) * acos(f) + mean) %% (2 * pi)
      obs <- obs + 1L
    } else if (log(c / U2) + 1 - c >= 0) {
      U3 <- runif(1, 0, 1)
      vm[obs] <- (sign(U3 - 0.5) * acos(f) + mean) %% (2 * pi)
      obs <- obs + 1L
    }
  }
  vm
}

# One group, as one call of cubos.generator.iso/aniso() in the scripts.
.sim_group3d <- function(n, M, sigma, type, gamma, df, rotate) {
  k <- nrow(M)
  p <- 3L * k
  param1 <- seq(-pi, pi)
  param2 <- param1[sample.int(length(param1), 1L)]
  param3 <- sample.int(2049^3, 1L)
  Sigma <- if (type == "iso") {
    diag(sigma^2, p)
  } else {
    kronecker(sigma^2 * (matrix(1, k, k) + (gamma - 1) * diag(k)),
              (matrix(1, 3L, 3L) + (gamma - 1) * diag(3L)) / gamma^2)
  }
  sims <- .rmvt_shifted(n, Sigma, df)
  sims <- sims[, c(seq_len(k) * 3L - 2L, seq_len(k) * 3L - 1L, seq_len(k) * 3L),
               drop = FALSE]
  arr <- array(as.vector(M) + t(sims), dim = c(k, 3L, n))
  if (!rotate) return(list(X = arr, X_norot = arr))
  angle <- .rvm(1L, param3, param2)
  rot <- arr
  rot[, 1L, ] <- arr[, 1L, ] * cos(angle) - arr[, 2L, ] * sin(angle)
  rot[, 2L, ] <- arr[, 1L, ] * sin(angle) + arr[, 2L, ] * cos(angle)
  list(X = rot, X_norot = arr)
}

#' Simulate noisy cubes and parallelepipeds in 3D
#'
#' Simulator of the 3D clustering study: each group is a sample of noisy
#' copies of a template configuration (a cube, or a box with one edge twice
#' as long) with 8 or 34 landmarks, with isotropic or anisotropic
#' perturbation, rotated by a common random angle about the z axis.
#'
#' @param n Number of configurations per group: a single number, or a vector
#'   of length \code{K}.
#' @param sigma Dispersion (default 3). The covariance scale is
#'   \code{sigma^2}.
#' @param type \code{"iso"} (independent Gaussian-like noise on every
#'   coordinate, default) or \code{"aniso"} (landmark and coordinate
#'   covariance of Kronecker form, see Details).
#' @param k Number of landmarks of the built-in templates, 8 (default) or
#'   34. Ignored if \code{templates} is given.
#' @param groups Character vector naming the template of each group, from
#'   \code{"cube"} and \code{"parallelepiped"} (default: both, so K = 2).
#'   Names may be repeated or reordered.
#' @param templates Optional list of \eqn{k \times 3} matrices used instead
#'   of the built-in ones, one per group.
#' @param gamma Anisotropy parameter of \code{type = "aniso"} (default 4, as
#'   in the scripts).
#' @param df Degrees of freedom of the multivariate t noise (default 99, as
#'   in the scripts).
#' @param rotate Logical. If \code{TRUE} (default) each group is rotated by
#'   one random angle about the z axis, as in the scripts. Shape distances
#'   do not depend on it.
#' @param independent Logical. The scripts reseed the generator with the
#'   same seed before each group, so that all groups share the same noise
#'   realization and the same rotation angle. This is reproduced when
#'   \code{independent = FALSE} (default). With \code{TRUE} the seed is set
#'   once and the groups get independent noise.
#' @param seed Optional integer seed; the caller's random-number stream is
#'   left as it was. In the scripts, replicate number \code{c} used
#'   \code{semente = c}. If \code{NULL} and \code{independent = FALSE} one
#'   seed is drawn from the current stream.
#'
#' @return A list with components:
#'   \describe{
#'     \item{X}{Numeric array \eqn{k \times 3 \times N} of the rotated
#'       configurations, groups stacked in order.}
#'     \item{group}{Integer vector of length \eqn{N}: the true group
#'       (1 to K).}
#'     \item{X_norot}{The configurations before rotation.}
#'     \item{K}{Number of groups.}
#'     \item{type, k, sigma, groups}{The options used.}
#'   }
#'
#' @details
#' For group \eqn{g} with template \eqn{M_g} (a \eqn{k \times 3} matrix),
#' configurations are \eqn{M_g + E}, where \eqn{\mathrm{vec}(E)} is drawn
#' from a multivariate t with \code{df} degrees of freedom and scale matrix
#' \eqn{\Sigma}: \eqn{\Sigma = \sigma^2 I_{3k}} for \code{"iso"}, and
#' \eqn{\Sigma = \sigma^2 (J_k + (\gamma - 1) I_k) \otimes
#' (J_3 + (\gamma - 1) I_3)/\gamma^2} for \code{"aniso"}, with \eqn{J} a
#' matrix of ones, the landmark-major ordering \eqn{(x_1, y_1, z_1, x_2,
#' \dots)}. The random draws follow the scripts: first a concentration
#' (one of \code{seq(-pi, pi)}) and a mean angle for the von Mises draw,
#' then the noise, then the angle. The noise reproduces
#' \code{mvtnorm::rmvt()} (version 0.9-9994 or later, whose random numbers
#' fill the matrix by rows) and the angle reproduces
#' \code{CircStats::rvm()}, without depending on those packages. Note that
#' some concentrations in \code{seq(-pi, pi)} are negative and the mean
#' angle is a large integer; this is kept as it is in the scripts.
#'
#' Built-in 8-landmark templates are the vertices of a cube
#' \eqn{[0, 10]^3} and of the box \eqn{[0, 10] \times [0, 20] \times
#' [0, 10]} with a vertex labelling that is the mirror image of the cube's
#' (reflections matter, see \code{\link{shape_distance}}); the 34-landmark
#' templates add edge and face points. They are the coordinates of
#' \code{cube8landm}, \code{parallelep8landm}, \code{cube34landm} and
#' \code{parallelep34landm} in the 'Anthropometry' package.
#'
#' @seealso \code{\link{shape_kmeans}}, \code{\link{shape_distance}},
#'   \code{\link{sim_landmark_gaussian}} for 2D data.
#'
#' @examples
#' sim <- sim_shapes3d(n = 5, sigma = 1, type = "aniso", k = 8, seed = 7)
#' dim(sim$X)
#' sim$group
#' # Reproducible, and the caller's random-number stream is not touched
#' identical(sim$X, sim_shapes3d(n = 5, sigma = 1, type = "aniso", seed = 7)$X)
#'
#' @export
sim_shapes3d <- function(n, sigma = 3,
                         type = c("iso", "aniso"),
                         k = 8L,
                         groups = c("cube", "parallelepiped"),
                         templates = NULL,
                         gamma = 4, df = 99,
                         rotate = TRUE,
                         independent = FALSE,
                         seed = NULL) {
  type <- match.arg(type)
  if (!is.numeric(sigma) || length(sigma) != 1L || sigma < 0) {
    stop("sigma must be a non-negative number.", call. = FALSE)
  }
  if (is.null(templates)) {
    groups <- match.arg(groups, c("cube", "parallelepiped"), several.ok = TRUE)
    if (!k %in% c(8, 34)) stop("k must be 8 or 34.", call. = FALSE)
    templates <- lapply(groups, .shape_template, k = as.integer(k))
  } else {
    if (!is.list(templates) || !all(vapply(templates, function(M)
      is.matrix(M) && is.numeric(M) && ncol(M) == 3L, logical(1)))) {
      stop("'templates' must be a list of k x 3 numeric matrices.", call. = FALSE)
    }
    if (length(unique(vapply(templates, nrow, integer(1)))) != 1L) {
      stop("all templates need the same number of landmarks.", call. = FALSE)
    }
    groups <- if (is.null(names(templates))) paste0("template", seq_along(templates)) else names(templates)
  }
  K <- length(templates)
  n <- rep_len(as.integer(n), K)
  if (anyNA(n) || any(n < 1L)) stop("n must be positive.", call. = FALSE)

  .local_seed(seed)
  if (!independent && is.null(seed)) {
    seed <- sample.int(.Machine$integer.max, 1L)
    # The loop below reseeds the stream; hand the caller back the state it
    # has after drawing the seed.
    after <- .rng_get()
    on.exit(.rng_set(after), add = TRUE)
  }

  Xs <- vector("list", K)
  X0 <- vector("list", K)
  for (g in seq_len(K)) {
    if (!independent) set.seed(seed)
    grp <- .sim_group3d(n[g], templates[[g]], sigma, type, gamma, df, rotate)
    Xs[[g]] <- grp$X
    X0[[g]] <- grp$X_norot
  }
  k_used <- nrow(templates[[1L]])
  bind <- function(L) array(unlist(L), dim = c(k_used, 3L, sum(n)))
  list(X = bind(Xs),
       group = rep(seq_len(K), times = n),
       X_norot = bind(X0),
       K = K, type = type, k = k_used, sigma = sigma, groups = groups)
}
