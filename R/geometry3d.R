## ---------------------------------------------------------------------------
##  Geometry for real landmark arrays (k x m x n, m = 2 or 3).
##
##  The complex representation in geometry.R is specific to m = 2. Here a
##  configuration is a real k x m matrix; pre-shapes are formed with the
##  Helmert sub-matrix and unit norm, and distances come from the singular
##  values of the cross-product of two pre-shapes. The 2D complex API is not
##  touched: for m = 2 the real routines return the same distances as
##  rho_matrix() and dF_matrix() on landmarks_to_preshape().
## ---------------------------------------------------------------------------

# Validate a landmark array and return it as a double array.
.check_array <- function(X, what = "X", min_n = 1L) {
  if (!is.numeric(X) || length(dim(X)) != 3L) {
    stop("'", what, "' must be a numeric array of dimension k x m x n.",
         call. = FALSE)
  }
  d <- dim(X)
  if (!d[2L] %in% c(2L, 3L)) {
    stop("'", what, "' must have m = 2 or m = 3 coordinates per landmark.",
         call. = FALSE)
  }
  if (d[1L] < d[2L] + 1L) {
    stop("'", what, "' needs at least k = m + 1 landmarks.", call. = FALSE)
  }
  if (d[3L] < min_n) {
    stop("'", what, "' needs at least ", min_n, " configuration(s).", call. = FALSE)
  }
  if (anyNA(X) || any(!is.finite(X))) {
    stop("'", what, "' contains missing or non-finite values.", call. = FALSE)
  }
  storage.mode(X) <- "double"
  dimnames(X) <- NULL
  X
}

.distance_code <- function(distance) {
  match(distance, c("rho", "dF", "dP")) - 1L
}

#' Shape distances for real landmark arrays
#'
#' Pairwise distances in Kendall's shape space for configurations of
#' \eqn{k} landmarks in \eqn{m = 2} or \eqn{3} dimensions, computed in C++.
#'
#' @param X Numeric array of dimension \eqn{k \times m \times n}, with
#'   \eqn{m} equal to 2 or 3 and \eqn{k \ge m + 1}.
#' @param distance One of \code{"rho"} (Riemannian, default), \code{"dF"}
#'   (full Procrustes) or \code{"dP"} (partial Procrustes).
#' @param reflection Logical. If \code{FALSE} (default) reflections are not
#'   removed, so a configuration and its mirror image are different shapes.
#'   If \code{TRUE}, shapes are compared up to reflection too.
#' @param Y Optional second array of dimension \eqn{k \times m \times n_y}.
#'   If given, the cross-distances between \code{X} and \code{Y} are
#'   returned instead.
#'
#' @return A numeric matrix: the symmetric \eqn{n \times n} matrix of
#'   distances between the configurations of \code{X}, or the
#'   \eqn{n \times n_y} matrix between those of \code{X} and \code{Y}.
#'
#' @details
#' Each configuration is centered with the Helmert sub-matrix and scaled to
#' unit norm, giving a \eqn{(k-1) \times m} pre-shape \eqn{z}. If
#' \eqn{s_1 \ge \dots \ge s_m} are the singular values of \eqn{z^\top w},
#' the cosine of the Riemannian distance is
#' \eqn{c = s_1 + \dots + s_{m-1} \pm s_m}, with the sign of \eqn{s_m} equal to
#' the sign of \eqn{\det(z^\top w)} unless \code{reflection = TRUE}, in
#' which case all signs are positive. Then
#' \deqn{\rho = \arccos c, \quad d_F = \sqrt{1 - c^2}, \quad
#'   d_P = \sqrt{2(1 - c)},}
#' with \eqn{c} clamped to \eqn{[0, 1]}. For \eqn{m = 2} without
#' reflection \eqn{c} equals \eqn{|z^* w|}, so the result agrees with
#' \code{\link{rho_matrix}} and \code{\link{dF_matrix}}. Two configurations
#' that are exactly equal get distance 0.
#'
#' @references
#' Dryden, I. L. and Mardia, K. V. (2016). \emph{Statistical Shape
#' Analysis, with Applications in R}, 2nd ed. Wiley, sections 4.1 and 4.2.
#'
#' @seealso \code{\link{shape_mean}}, \code{\link{shape_kmeans}}
#'
#' @examples
#' set.seed(1)
#' X <- array(rnorm(8 * 3 * 6), dim = c(8, 3, 6))
#' D <- shape_distance(X)
#' round(D[1:3, 1:3], 3)
#' all.equal(shape_distance(X, "dF"), sin(D))
#'
#' # A mirror image is a different shape unless reflections are allowed
#' Xm <- X; Xm[, 1, ] <- -Xm[, 1, ]
#' shape_distance(X[, , 1, drop = FALSE], Y = Xm[, , 1, drop = FALSE])
#' shape_distance(X[, , 1, drop = FALSE], Y = Xm[, , 1, drop = FALSE],
#'                reflection = TRUE)
#'
#' @export
shape_distance <- function(X, distance = c("rho", "dF", "dP"),
                           reflection = FALSE, Y = NULL) {
  distance <- match.arg(distance)
  stopifnot(is.logical(reflection), length(reflection) == 1L, !is.na(reflection))
  X <- .check_array(X, "X")
  if (is.null(Y)) {
    return(shape_dist_cpp(X, X, .distance_code(distance), reflection, TRUE))
  }
  Y <- .check_array(Y, "Y")
  if (!identical(dim(X)[1:2], dim(Y)[1:2])) {
    stop("'X' and 'Y' must have the same k and m.", call. = FALSE)
  }
  shape_dist_cpp(X, Y, .distance_code(distance), reflection, FALSE)
}

# Centre and scale each configuration (slice) to unit Frobenius norm.
.normalize_configs <- function(X) {
  for (i in seq_len(dim(X)[3L])) {
    M <- X[, , i]
    M <- sweep(M, 2L, colMeans(M))
    nM <- sqrt(sum(M^2))
    X[, , i] <- if (nM > 0) M / nM else M
  }
  X
}

#' Mean shape of a sample of 3D or 2D landmark configurations
#'
#' Two notions of mean shape for real landmark arrays: the full Procrustes
#' mean, and the MDS (or \eqn{\Phi}) mean, an extrinsic mean obtained from
#' the Veronese-Whitney embedding.
#'
#' @param X Numeric array of dimension \eqn{k \times m \times n}, with
#'   \eqn{m} equal to 2 or 3.
#' @param type \code{"procrustes"} (default) or \code{"phi"}.
#' @param tol,max_iter Convergence tolerance (relative change of the
#'   objective) and iteration cap of the Procrustes iteration. Used for
#'   \code{type = "procrustes"} only: the \code{"phi"} mean has an explicit
#'   formula and needs the Procrustes mean only to choose between a shape
#'   and its mirror image, for which a looser fixed tolerance is used.
#'
#' @return A \eqn{k \times m} matrix: a centered configuration of unit size
#'   (Frobenius norm 1) representing the mean shape. A warning is given if
#'   the Procrustes iteration does not reach \code{tol} within
#'   \code{max_iter} steps.
#'
#' @details
#' \strong{Procrustes.} The full Procrustes mean minimizes
#' \eqn{\sum_i d_F^2(x_i, \mu)}. It is computed by the fixed-point iteration
#' \eqn{\mu \leftarrow \mathrm{normalize}(\sum_i \cos\rho_i \, z_i R_i)},
#' where \eqn{R_i \in SO(m)} aligns the pre-shape \eqn{z_i} to \eqn{\mu},
#' started at the first configuration. Each step does not decrease
#' \eqn{\sum_i \cos^2\rho_i}. For \eqn{m = 2} it agrees with the dominant
#' eigenvector of \eqn{\sum_i z_i z_i^*}. It is the mean that
#' \code{shapes::procGPA()} returns in \code{$mshape}, up to the looser
#' convergence tolerance of that function.
#'
#' \strong{phi.} The MDS mean takes the \eqn{m} leading eigenvectors
#' \eqn{v_j} (eigenvalues \eqn{\lambda_j}) of
#' \eqn{n^{-1}\sum_i z_i z_i^\top}, forms the pre-shape
#' \eqn{[\sqrt{\lambda_1} v_1, \dots, \sqrt{\lambda_m} v_m] /
#' \sqrt{\lambda_1 + \dots + \lambda_m}} and lifts it with the transposed
#' Helmert sub-matrix. This is \code{shapes::MDSshape(x)} with its defaults
#' (\code{alpha = 1}, \code{projalpha = 1/2}). Because
#' \eqn{z z^\top} does not change under a reflection of \eqn{z}, the
#' last axis is flipped if that brings the mean closer to the Procrustes
#' mean, exactly as in \code{shapes::MDSshape}.
#'
#' @references
#' Kent, J. T. (1994). The complex Bingham distribution and shape analysis.
#' \emph{Journal of the Royal Statistical Society B}, 56, 285--299.
#' \doi{10.1111/j.2517-6161.1994.tb01978.x}
#'
#' Dryden, I. L. and Mardia, K. V. (2016). \emph{Statistical Shape
#' Analysis, with Applications in R}, 2nd ed. Wiley, chapters 6 and 7.
#'
#' @seealso \code{\link{shape_kmeans}}, which updates its centers with
#'   these means.
#'
#' @examples
#' set.seed(2)
#' base <- matrix(rnorm(8 * 3), 8, 3)
#' X <- array(rep(base, 20), dim = c(8, 3, 20)) +
#'   array(rnorm(8 * 3 * 20, sd = 0.05), dim = c(8, 3, 20))
#' mu <- shape_mean(X)
#' nu <- shape_mean(X, type = "phi")
#' # Both means are close to the template
#' shape_distance(array(c(mu, nu, base), dim = c(8, 3, 3)))
#'
#' @export
shape_mean <- function(X, type = c("procrustes", "phi"),
                       tol = 1e-12, max_iter = 1000L) {
  type <- match.arg(type)
  X <- .check_array(X, "X", min_n = 1L)
  if (type == "phi") {
    # The Procrustes mean only fixes the orientation of the phi mean, so a
    # loose tolerance is enough and non-convergence is not an error.
    mu <- procrustes_mean_cpp(X, 1e-8, 200L)$mean
    return(phi_mean_cpp(X, mu))
  }
  res <- procrustes_mean_cpp(X, tol, as.integer(max_iter))
  if (!res$converged) {
    warning("the Procrustes iteration did not converge in ", max_iter,
            " steps.", call. = FALSE)
  }
  res$mean
}
