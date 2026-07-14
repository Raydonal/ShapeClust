## ---------------------------------------------------------------------------
##  Geometry: from raw landmark configurations to complex pre-shape vectors,
##  and pairwise shape-space distances.
##
##  These functions provide the geometric foundation of the package. They
##  operate on 2D landmark configurations (arrays k x 2 x n) and produce
##  complex pre-shape matrices Z of size (k-1) x n whose columns lie on
##  the unit complex sphere S^{2k-3}_C.
##
##  All heavy computation is delegated to the C++ backend in src/. The
##  R functions here are thin, well-documented wrappers.
## ---------------------------------------------------------------------------

#' Helmert sub-matrix
#'
#' Returns the (k-1) x k Helmert sub-matrix used to remove translation
#' from a landmark configuration.
#'
#' @param k Integer. Number of landmarks (>= 2).
#'
#' @return A real matrix of size (k-1) x k.
#'
#' @details
#' The full Helmert matrix \eqn{H^F} is the k x k orthogonal matrix with
#' first row \eqn{(1/\sqrt{k}, \dots, 1/\sqrt{k})} and, for
#' \eqn{j \ge 2}, jth row \eqn{(-1, \dots, -1, j-1, 0, \dots, 0)/\sqrt{j(j-1)}}
#' where the entry \eqn{j-1} occurs at column \eqn{j}. Dropping the first
#' row yields the Helmert sub-matrix \eqn{H}.
#'
#' Applying \eqn{H} to a configuration removes translation, and
#' \eqn{H^\top H = I_k - k^{-1} \mathbf{1}_k \mathbf{1}_k^\top} is the
#' centring matrix -- so \eqn{H} is the orthonormal representation of
#' centring.
#'
#' @examples
#' H <- helmert_submatrix(4)
#' # Rows are orthonormal:
#' round(H %*% t(H), 12)
#'
#' @export
helmert_submatrix <- function(k) {
  k <- as.integer(k)
  if (k < 2L) stop("k must be at least 2.")
  H <- matrix(0.0, nrow = k - 1L, ncol = k)
  for (j in 2:k) {
    denom <- sqrt(j * (j - 1))
    H[j - 1L, 1:(j - 1L)] <- -1 / denom
    H[j - 1L, j]          <- (j - 1) / denom
  }
  H
}

#' Convert 2D landmark configurations to complex pre-shape vectors
#'
#' Given a 3D array of shape (k, 2, n) of Cartesian coordinates, returns
#' the (k-1) x n complex matrix Z whose columns are the corresponding
#' unit-norm complex pre-shapes.
#'
#' @param X Numeric array of shape (k, 2, n).
#' @param center If TRUE (default), remove translation via the Helmert
#'   sub-matrix. If FALSE, assume the configurations are already centred.
#'
#' @return A complex matrix of size (k-1) x n.
#'
#' @examples
#' set.seed(1)
#' k <- 4; n <- 5
#' X <- array(rnorm(k * 2 * n), dim = c(k, 2, n))
#' Z <- landmarks_to_preshape(X)
#' # Columns have unit modulus:
#' round(colSums(Mod(Z)^2), 12)
#'
#' @export
landmarks_to_preshape <- function(X, center = TRUE) {
  stopifnot(length(dim(X)) == 3L, dim(X)[2L] == 2L)
  k <- dim(X)[1L]
  n <- dim(X)[3L]

  # Vectorise to complex: each landmark becomes x + i*y
  Z0 <- matrix(complex(real      = as.numeric(X[, 1L, , drop = TRUE]),
                       imaginary = as.numeric(X[, 2L, , drop = TRUE])),
               nrow = k, ncol = n)

  if (center) {
    H <- helmert_submatrix(k)
    Z_H <- H %*% Z0
  } else {
    # If already centred, drop mean row equivalent: use first (k-1)
    # coordinates after Helmert projection anyway for consistency.
    H <- helmert_submatrix(k)
    Z_H <- H %*% Z0
  }
  # Normalise each column to unit complex norm
  norms <- sqrt(colSums(Mod(Z_H)^2))
  # Avoid division by zero
  norms[norms == 0] <- 1
  sweep(Z_H, 2L, norms, FUN = "/")
}


## ---------------------------------------------------------------------------
##  Distance matrix wrappers
##
##  Cost: O(n^2 (k-1)) via BLAS-backed Hermitian Gram matrix.
##  If the C++ backend is unavailable (e.g. package not compiled), we fall
##  back to a pure-R implementation that is slower but numerically
##  identical.
## ---------------------------------------------------------------------------

#' Riemannian geodesic distance matrix on Kendall's 2D shape space
#'
#' Returns the n x n symmetric matrix of intrinsic (Riemannian) distances
#' between the columns of a complex pre-shape matrix Z, in radians in
#' \eqn{[0, \pi/2]}.
#'
#' @param Z Complex matrix of size (k-1) x n, unit-norm columns.
#'
#' @return A symmetric numeric matrix.
#'
#' @details
#' \deqn{\rho(z_1, z_2) = \arccos |z_1^* z_2|}
#'
#' Cost: \eqn{\mathcal{O}(n^2 (k-1))} via a single BLAS-backed Hermitian
#' Gram matrix.
#'
#' @examples
#' k <- 4; n <- 5
#' set.seed(1)
#' X <- array(rnorm(k * 2 * n), dim = c(k, 2, n))
#' Z <- landmarks_to_preshape(X)
#' D <- rho_matrix(Z)
#' isSymmetric(D)
#'
#' @export
rho_matrix <- function(Z) {
  stopifnot(is.complex(Z), is.matrix(Z))
  if (has_cpp_backend()) {
    return(rho_matrix_cpp(Z))
  }
  # Pure R fallback
  G <- Conj(t(Z)) %*% Z
  M <- pmin(Mod(G), 1)
  D <- acos(M)
  diag(D) <- 0
  D
}

#' Full Procrustes distance matrix
#'
#' \deqn{d_F(z_1, z_2) = \sqrt{1 - |z_1^* z_2|^2}}
#'
#' @param Z Complex pre-shape matrix of size (k-1) x n.
#'
#' @return A symmetric numeric matrix.
#'
#' @export
dF_matrix <- function(Z) {
  stopifnot(is.complex(Z), is.matrix(Z))
  if (has_cpp_backend()) {
    return(dF_matrix_cpp(Z))
  }
  G <- Conj(t(Z)) %*% Z
  c2 <- pmin(Mod(G)^2, 1)
  D <- sqrt(1 - c2)
  diag(D) <- 0
  D
}

# Internal utility: is the compiled backend available?
has_cpp_backend <- function() {
  exists("rho_matrix_cpp", mode = "function",
         envir = asNamespace("ShapeClust"), inherits = FALSE)
}
