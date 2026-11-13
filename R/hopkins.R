## Hopkins statistic of clustering tendency for landmark shapes.

# Procrustes-align an array to a mean shape; rows are configurations.
.align_coords <- function(X) {
  k <- dim(X)[1L]; d <- dim(X)[2L]; n <- dim(X)[3L]
  mu <- shape_mean(X)
  out <- matrix(0, n, k * d)
  for (i in seq_len(n)) {
    M <- X[, , i]
    M <- sweep(M, 2L, colMeans(M))
    M <- M / sqrt(sum(M^2))
    s <- svd(crossprod(mu, M))
    D <- diag(c(rep(1, d - 1L), sign(det(s$u %*% t(s$v)))), d)
    out[i, ] <- as.vector(M %*% (s$v %*% D %*% t(s$u)))
  }
  out
}

.hopkins_ratio <- function(u, w, p) {
  su <- sum(u^p); su / (su + sum(w^p))
}

# Draw m points from the Gaussian with the mean and covariance of A.
.tangent_draw <- function(A, m) {
  n <- nrow(A)
  mu <- colMeans(A)
  s <- svd(sweep(A, 2L, mu), nu = 0L)
  r <- sum(s$d > 1e-8 * s$d[1L])
  Z <- matrix(rnorm(m * r), m, r) %*% (s$d[seq_len(r)] / sqrt(n - 1)
       * t(s$v[, seq_len(r), drop = FALSE]))
  sweep(Z, 2L, mu, "+")
}

# H on tangent coordinates; A is n x q.
.hop_tangent <- function(A, m, p) {
  n <- nrow(A)
  idx <- sample.int(n, m)
  G <- .tangent_draw(A, m)
  all_d <- as.matrix(stats::dist(rbind(G, A)))
  Du <- all_d[seq_len(m), m + seq_len(n), drop = FALSE]
  D <- all_d[m + seq_len(n), m + seq_len(n)]
  diag(D) <- Inf
  .hopkins_ratio(apply(Du, 1L, min), apply(D[idx, , drop = FALSE], 1L, min), p)
}

.hop_uniform <- function(X, m, p, distance) {
  n <- dim(X)[3L]
  idx <- sample.int(n, m)
  U <- array(rnorm(prod(dim(X)[1:2]) * m), c(dim(X)[1:2], m))
  Du <- shape_distance(U, Y = X, distance = distance)
  D <- shape_distance(X, distance = distance)
  diag(D) <- Inf
  .hopkins_ratio(apply(Du, 1L, min), apply(D[idx, , drop = FALSE], 1L, min), p)
}

#' Hopkins test of clustering tendency for landmark shapes
#'
#' Hopkins statistic with a Monte Carlo reference distribution, for a sample
#' of shapes given as a real \eqn{k \times d \times n} landmark array
#' (\eqn{d = 2} or \eqn{3}).
#'
#' @param X Numeric array \eqn{k \times d \times n}, \eqn{d \in \{2, 3\}}.
#' @param null \code{"tangent"} (default) or \code{"uniform"}; see Details.
#' @param m Number of artificial and of sampled points. Default
#'   \code{round(0.1 * n)}, at least 1.
#' @param power Exponent \eqn{p}; default 1.
#' @param B Number of Monte Carlo reference draws.
#' @param distance Shape distance (\code{"rho"}, \code{"dF"}, \code{"dP"}) for
#'   \code{null = "uniform"}. The tangent null uses Euclidean distance in the
#'   aligned coordinates.
#' @param seed Optional integer seed; the caller's random-number stream is
#'   restored on exit.
#'
#' @details
#' The statistic is \eqn{H = \sum u^p / (\sum u^p + \sum w^p)}, where
#' \eqn{u} are the nearest-neighbor distances from \eqn{m} artificial points
#' to the data and \eqn{w} those from \eqn{m} sampled data points to the other
#' data. Large values indicate clustering.
#'
#' With \code{null = "tangent"}, the artificial points come from a unimodal
#' Gaussian with the mean and covariance of the Procrustes-aligned data
#' coordinates (centered, unit-size configurations rotated onto the Procrustes
#' mean), and distances are Euclidean in those coordinates. This tests
#' unimodality against clustering. With \code{null = "uniform"}, the
#' artificial points are uniform on shape space (Gaussian landmarks) and
#' distances are \code{distance}. \strong{The uniform null rejects for any
#' concentrated sample, clustered or not, so it says nothing about clusters.}
#'
#' The p-value does not assume \eqn{H = 0.5}. For each of \code{B} draws a
#' full data set of \eqn{n} objects is simulated from the null (for the
#' tangent null, the Gaussian fitted to the data; for the uniform null,
#' \eqn{n} uniform shapes), \eqn{H} is computed on it by the same procedure,
#' and \eqn{p = (1 + \#\{H_b \ge H\}) / (B + 1)}.
#'
#' The exponent usually recommended is the dimension of the space. With a
#' Monte Carlo reference the choice affects power, not validity.
#'
#' @return An object of class \code{"hopkins_shape"}: a list with
#'   \code{statistic}, \code{p.value}, \code{null}, \code{power}, \code{m},
#'   \code{B}, and \code{reference} (the \code{B} reference values).
#'
#' @seealso \code{\link{shape_distance}}, \code{\link{shape_mean}}
#' @examples
#' sim <- sim_shapes3d(n = 30, sigma = 1, k = 8, seed = 1)
#' hopkins_shape(sim$X, B = 19, seed = 2)
#' @export
hopkins_shape <- function(X, null = c("tangent", "uniform"), m = NULL,
                          power = NULL, B = 199, distance = "rho",
                          seed = NULL) {
  null <- match.arg(null)
  distance <- match.arg(distance, c("rho", "dF", "dP"))
  X <- .check_array(X, "X", min_n = 3L)
  n <- dim(X)[3L]
  if (is.null(m)) m <- max(1L, round(0.1 * n))
  if (is.null(power)) power <- 1
  if (!is.numeric(m) || length(m) != 1L || m < 1 || m >= n)
    stop("'m' must be an integer in [1, n - 1].", call. = FALSE)
  if (!is.numeric(power) || length(power) != 1L || power <= 0)
    stop("'power' must be a positive number.", call. = FALSE)
  if (!is.numeric(B) || length(B) != 1L || B < 1)
    stop("'B' must be a positive integer.", call. = FALSE)
  m <- as.integer(m); B <- as.integer(B)
  .local_seed(seed)

  if (null == "tangent") {
    A <- .align_coords(X)
    H <- .hop_tangent(A, m, power)
    ref <- vapply(seq_len(B), function(b)
      .hop_tangent(.tangent_draw(A, n), m, power), numeric(1))
  } else {
    H <- .hop_uniform(X, m, power, distance)
    ref <- vapply(seq_len(B), function(b) {
      Xb <- array(rnorm(prod(dim(X))), dim(X))
      .hop_uniform(Xb, m, power, distance)
    }, numeric(1))
  }
  structure(list(statistic = H, p.value = (1 + sum(ref >= H)) / (B + 1),
                 null = null, power = power, m = m, B = B, reference = ref),
            class = "hopkins_shape")
}

#' @export
print.hopkins_shape <- function(x, ...) {
  cat("Hopkins test of clustering tendency (", x$null, " null)\n", sep = "")
  cat("H = ", format(x$statistic, digits = 4), ", p = ",
      format(x$p.value, digits = 3), "  (m = ", x$m, ", power = ", x$power,
      ", B = ", x$B, ")\n", sep = "")
  invisible(x)
}
