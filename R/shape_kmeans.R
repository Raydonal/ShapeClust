## ---------------------------------------------------------------------------
##  Shape k-means: Lloyd's algorithm on the shape space with a Procrustes
##  (or MDS/phi) mean as the cluster center. Port of LloydShapes() and
##  LloydShapes.phi() from the 3D bagging scripts (Nascimento et al.),
##  with the distance and mean computations in C++.
## ---------------------------------------------------------------------------

#' Shape k-means for real landmark arrays
#'
#' Lloyd's algorithm for clustering configurations of \eqn{k} landmarks in
#' two or three dimensions. Observations are assigned to the nearest center
#' in shape space and each center is then replaced by the mean shape of its
#' cluster, with several random restarts.
#'
#' @param X Numeric array of dimension \eqn{k \times m \times n}, with
#'   \eqn{m} equal to 2 or 3.
#' @param K Number of clusters (>= 2).
#' @param mean Mean shape used to update the centers: \code{"procrustes"}
#'   (full Procrustes mean, as \code{shapes::procGPA()}; default) or
#'   \code{"phi"} (MDS mean, as \code{shapes::MDSshape()}). See
#'   \code{\link{shape_mean}}.
#' @param alg_steps Maximum number of Lloyd steps per restart (default 10).
#' @param niter Number of random restarts (default 10).
#' @param tol Stopping rule: a restart stops when the relative decrease of
#'   the objective between two steps is below \code{tol} (default 1e-4).
#' @param distance Distance in which the objective is measured, one of
#'   \code{"rho"} (default), \code{"dF"} or \code{"dP"}. All three are
#'   increasing functions of the same cosine, so assignments do not depend
#'   on this choice; the objective value and the stopping rule do.
#'
#' @return A list with components:
#'   \describe{
#'     \item{cluster}{Integer vector of length n, values in 1:K.}
#'     \item{centers}{Numeric array \eqn{k \times m \times K} with the
#'       cluster means, each centered with unit norm.}
#'     \item{objective}{Mean squared distance between each observation and
#'       the center it was assigned to, at the step that gave the returned
#'       clustering (the centers in force at the start of that step).}
#'     \item{n_iter}{Index of the Lloyd step, within the restart that gave
#'       the returned clustering, at which the objective was lowest.}
#'     \item{initials}{List with the K observations drawn as initial centers
#'       in each restart.}
#'     \item{mean, distance}{The options used.}
#'   }
#'
#' @details
#' The algorithm follows \code{LloydShapes()} of the scripts used for the
#' 3D bagging study. Initial centers are \code{K} observations drawn with
#' \code{sample.int()} at the start of each restart. In each step every
#' observation goes to its nearest center (\code{max.col()} with random tie
#' breaking, as in the scripts), the centers are updated, and the objective
#' is evaluated with the distances of the assignment step, i.e. before the
#' update. The clustering and centers kept are those of the step with the
#' lowest objective over all steps of all restarts. A cluster with a single
#' member has that member as its center. If a cluster becomes empty its
#' center is left unchanged (the script stops with an error in that case).
#' Reflections are not removed (see \code{\link{shape_distance}}).
#'
#' The caller's random-number stream is advanced; use
#' \code{set.seed()} for reproducibility.
#'
#' @references
#' Dryden, I. L. and Mardia, K. V. (2016). \emph{Statistical Shape
#' Analysis, with Applications in R}, 2nd ed. Wiley.
#'
#' @seealso \code{\link{shape_mean}}, \code{\link{bagclust1}} with
#'   \code{base = "kmeans"}.
#'
#' @examples
#' sim <- sim_shapes3d(n = 10, sigma = 0.5, seed = 1)
#' set.seed(1)
#' fit <- shape_kmeans(sim$X, K = 2)
#' table(fit$cluster, sim$group)
#' fit$objective
#'
#' @export
shape_kmeans <- function(X, K,
                         mean = c("procrustes", "phi"),
                         alg_steps = 10L,
                         niter = 10L,
                         tol = 1e-4,
                         distance = c("rho", "dF", "dP")) {
  mean <- match.arg(mean)
  distance <- match.arg(distance)
  X <- .check_array(X, "X", min_n = 3L)
  n <- dim(X)[3L]
  K <- as.integer(K)
  if (K < 2L || K > n - 1L) stop("K must satisfy 2 <= K <= n-1.", call. = FALSE)
  alg_steps <- as.integer(alg_steps)
  niter <- as.integer(niter)
  if (alg_steps < 1L || niter < 1L) stop("alg_steps and niter must be >= 1.")
  type <- .distance_code(distance)

  update_center <- if (mean == "procrustes") {
    function(Xh) procrustes_mean_cpp(Xh, 1e-10, 1000L)$mean
  } else {
    function(Xh) shape_mean(Xh, type = "phi")
  }

  vopt <- 1e8
  copt <- NULL
  asig_opt <- NULL
  iter_opt <- NA_integer_
  initials <- vector("list", niter)

  for (iter in seq_len(niter)) {
    initials[[iter]] <- sample.int(n, K)
    centers <- X[, , initials[[iter]], drop = FALSE]
    obj <- numeric(0)
    asig_step <- list()
    cent_step <- list()
    for (step in seq_len(alg_steps)) {
      dist <- shape_dist_cpp(X, centers, type, FALSE, FALSE)
      asig <- max.col(-dist)
      for (h in seq_len(K)) {
        idx <- which(asig == h)
        if (length(idx) == 1L) {
          centers[, , h] <- X[, , idx]
        } else if (length(idx) > 1L) {
          centers[, , h] <- update_center(X[, , idx, drop = FALSE])
        }
      }
      obj[step] <- sum(dist[cbind(seq_len(n), asig)]^2) / n
      asig_step[[step]] <- asig
      cent_step[[step]] <- centers
      if (step > 1L) {
        aux <- obj[step]
        aux1 <- obj[step - 1L]
        if (aux1 <= 0 || (aux1 - aux) / aux1 < tol) break
      }
    }
    best <- which.min(obj)
    if (obj[best] < vopt) {
      vopt <- obj[best]
      copt <- cent_step[[best]]
      asig_opt <- asig_step[[best]]
      iter_opt <- best
    }
  }

  list(cluster   = as.integer(asig_opt),
       centers   = .normalize_configs(copt),
       objective = vopt,
       n_iter    = as.integer(iter_opt),
       initials  = initials,
       mean      = mean,
       distance  = distance)
}
