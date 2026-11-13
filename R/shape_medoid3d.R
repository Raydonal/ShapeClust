## ---------------------------------------------------------------------------
##  Randomized medoid searches on a shape distance matrix: CLARANS and the
##  hill-climbing search of the 3D bagging scripts. Both accept a distance
##  matrix, as shape_pam() does, or the landmarks themselves.
## ---------------------------------------------------------------------------

# Turn the first argument of a clusterer into an n x n distance matrix:
# a real k x m x n array and a complex pre-shape matrix are converted with
# `distance`, and a real matrix is taken to be a distance matrix already.
.as_distance_matrix <- function(D, distance, reflection = FALSE) {
  if (is.array(D) && length(dim(D)) == 3L && !is.complex(D)) {
    return(shape_distance(D, distance, reflection))
  }
  if (is.complex(D)) {
    if (reflection) stop("'reflection' needs a real landmark array.", call. = FALSE)
    return(switch(distance,
                  rho = rho_matrix(D),
                  dF  = dF_matrix(D),
                  dP  = sqrt(2 * (1 - cos(rho_matrix(D))))))
  }
  if (!is.matrix(D) || !is.numeric(D) || nrow(D) != ncol(D)) {
    stop("'D' must be a square distance matrix, a k x m x n landmark array ",
         "or a complex pre-shape matrix.", call. = FALSE)
  }
  if (anyNA(D)) stop("'D' contains missing values.", call. = FALSE)
  storage.mode(D) <- "double"
  D
}

#' CLARANS on a shape distance matrix
#'
#' Clustering Large Applications based on RANdomized Search (Ng and Han,
#' 2002), in the faster variant of Schubert and Rousseeuw (2021): random
#' restarts of a randomized medoid-swap search. Own implementation in C++
#' that uses R's random-number stream.
#'
#' @param D Numeric \eqn{n \times n} distance matrix, or a real landmark
#'   array \eqn{k \times m \times n}, or a complex pre-shape matrix, from
#'   which distances are computed with \code{distance}.
#' @param K Number of clusters (>= 2).
#' @param numlocal Number of random restarts (default 2).
#' @param maxneighbor Search effort of each restart. A value below 1 is a
#'   rate and means \code{r = ceiling(maxneighbor * (n - K))}; a value of 1
#'   or more is taken as \code{r} itself. A restart ends after \code{r - 1}
#'   consecutive unsuccessful swap trials, so \code{r = 1} tries no swap at
#'   all. Default 0.0125.
#' @param distance Distance used if \code{D} is not a distance matrix: one
#'   of \code{"rho"} (default), \code{"dF"}, \code{"dP"}.
#' @param reflection Passed to \code{\link{shape_distance}} when \code{D}
#'   is a landmark array.
#'
#' @return A list with components:
#'   \describe{
#'     \item{cluster}{Integer vector of length n, values in 1:K: the index
#'       of the nearest medoid.}
#'     \item{medoids}{Integer vector of length K, 1-based object indices.}
#'     \item{objective}{Sum of the distances of the objects to their
#'       medoid.}
#'     \item{n_iter}{Number of accepted swaps in the restart that gave the
#'       returned medoids.}
#'     \item{distance}{Distance name, or \code{NA} if \code{D} was a
#'       distance matrix.}
#'   }
#'
#' @details
#' A restart draws K distinct random medoids. It then repeatedly draws a
#' random non-medoid \eqn{h}, evaluates the change of the total cost for
#' swapping \eqn{h} with each of the K medoids at once, and performs the
#' best swap if it lowers the cost; the restart ends after \eqn{r - 1}
#' consecutive trials without a swap, with \eqn{r} defined at
#' \code{maxneighbor}. The best restart is returned. The defaults are those
#' of the 3D bagging scripts, which called \code{fastkmedoids::fastclarans()}
#' with \code{numlocal = 2} and \code{maxneighbor = 0.0125}. This function
#' follows the same algorithm but draws from R's generator, so for the same
#' \code{set.seed()} it does not return the partition of that package, and
#' nothing beyond base R is needed.
#'
#' With the default \code{maxneighbor}, \eqn{r = 1} whenever
#' \eqn{n - K \le 80}, so no swap is ever tried and the result is the better
#' of \code{numlocal} random medoid sets; for \eqn{n = 100} and \eqn{K = 2}
#' (\eqn{r = 2}) one failed trial ends a restart. Raise \code{maxneighbor}
#' (for instance to 0.3, or to a count such as 100) for a real search.
#'
#' @references
#' Ng, R. T. and Han, J. (2002). CLARANS: a method for clustering objects
#' for spatial data mining. \emph{IEEE Transactions on Knowledge and Data
#' Engineering}, 14(5), 1003--1016. \doi{10.1109/TKDE.2002.1033770}
#'
#' Schubert, E. and Rousseeuw, P. J. (2021). Fast and eager k-medoids
#' clustering: O(k) runtime improvement of the PAM, CLARA, and CLARANS
#' algorithms. \emph{Information Systems}, 101, 101804.
#' \doi{10.1016/j.is.2021.101804}
#'
#' @seealso \code{\link{shape_pam}}, \code{\link{shape_hill}},
#'   \code{\link{bagclust1}}
#'
#' @examples
#' sim <- sim_shapes3d(n = 10, sigma = 0.5, seed = 1)
#' D <- shape_distance(sim$X)
#' set.seed(1)
#' fit <- shape_clarans(D, K = 2)
#' table(fit$cluster, sim$group)
#'
#' @export
shape_clarans <- function(D, K, numlocal = 2L, maxneighbor = 0.0125,
                          distance = c("rho", "dF", "dP"),
                          reflection = FALSE) {
  distance <- match.arg(distance)
  given_D <- is.matrix(D) && !is.complex(D)
  D <- .as_distance_matrix(D, distance, reflection)
  n <- nrow(D)
  K <- as.integer(K)
  if (K < 2L || K > n - 1L) stop("K must satisfy 2 <= K <= n-1.", call. = FALSE)
  numlocal <- as.integer(numlocal)
  if (numlocal < 1L) stop("numlocal must be >= 1.", call. = FALSE)
  if (!is.numeric(maxneighbor) || length(maxneighbor) != 1L || maxneighbor <= 0) {
    stop("maxneighbor must be a positive number.", call. = FALSE)
  }
  out <- shape_clarans_cpp(D, K, numlocal, as.numeric(maxneighbor))
  list(cluster   = as.integer(out$cluster),
       medoids   = as.integer(out$medoids),
       objective = out$objective,
       n_iter    = as.integer(out$n_iter),
       distance  = if (given_D) NA_character_ else distance)
}

#' Hill-climbing medoid search on a shape distance matrix
#'
#' Random-neighbor hill climbing for K medoids, ported from
#' \code{hill_climbing_shapes()} of the 3D bagging scripts.
#'
#' @param D Numeric \eqn{n \times n} distance matrix, or a real landmark
#'   array \eqn{k \times m \times n}, or a complex pre-shape matrix, from
#'   which distances are computed with \code{distance}.
#' @param K Number of clusters (>= 2).
#' @param max_iter Number of neighbors tried (default 300, as in the
#'   scripts).
#' @param from_best Logical. If \code{FALSE} (default), faithful to the
#'   script, each neighbor is a modification of the \emph{initial} medoid
#'   set. If \code{TRUE}, each neighbor modifies the best set found so far,
#'   which is hill climbing in the usual sense.
#' @param initial Optional vector of length K with initial medoid indices
#'   (1-based). If \code{NULL} they are drawn with \code{sample.int()}.
#' @param distance Distance used if \code{D} is not a distance matrix: one
#'   of \code{"rho"} (default), \code{"dF"}, \code{"dP"}.
#' @param reflection Passed to \code{\link{shape_distance}} when \code{D}
#'   is a landmark array.
#'
#' @return A list with components:
#'   \describe{
#'     \item{cluster}{Integer vector of length n, values in 1:K: the index
#'       of the nearest medoid.}
#'     \item{medoids}{Integer vector of length K, 1-based object indices.}
#'     \item{objective}{Sum of the distances of the objects to their
#'       nearest medoid.}
#'     \item{n_iter}{Number of neighbors that improved the objective.}
#'     \item{distance}{Distance name, or \code{NA} if \code{D} was a
#'       distance matrix.}
#'   }
#'
#' @details
#' Starting from \code{initial}, each of the \code{max_iter} iterations
#' replaces one randomly chosen medoid by a randomly chosen object,
#' redrawing (cumulatively, on the same candidate set) while the set
#' contains a repeated object, and keeps the candidate set if its objective
#' is strictly lower. Objects are assigned with \code{max.col()}-style
#' random tie breaking, as in the script, using R's generator, so for the
#' same seed and distances the result equals that of the script.
#'
#' In the script every neighbor is built from the initial set, never from
#' the best set found, so it is a random search around the starting point
#' and not a hill climb in the strict sense. This is the default so that
#' results can be compared with the published scripts;
#' \code{from_best = TRUE} gives the usual hill climbing.
#'
#' @seealso \code{\link{shape_pam}}, \code{\link{shape_clarans}},
#'   \code{\link{bagclust1}}
#'
#' @examples
#' sim <- sim_shapes3d(n = 10, sigma = 0.5, seed = 1)
#' D <- shape_distance(sim$X)
#' set.seed(1)
#' fit <- shape_hill(D, K = 2, max_iter = 100)
#' table(fit$cluster, sim$group)
#'
#' @export
shape_hill <- function(D, K, max_iter = 300L, from_best = FALSE,
                       initial = NULL,
                       distance = c("rho", "dF", "dP"),
                       reflection = FALSE) {
  distance <- match.arg(distance)
  given_D <- is.matrix(D) && !is.complex(D)
  D <- .as_distance_matrix(D, distance, reflection)
  n <- nrow(D)
  K <- as.integer(K)
  if (K < 2L || K > n - 1L) stop("K must satisfy 2 <= K <= n-1.", call. = FALSE)
  max_iter <- as.integer(max_iter)
  if (max_iter < 0L) stop("max_iter must be >= 0.", call. = FALSE)
  if (is.null(initial)) {
    initial <- sample.int(n, K)
  } else {
    if (length(initial) != K || anyDuplicated(initial) || !all(initial %in% seq_len(n))) {
      stop("'initial' must be K distinct indices in 1:n.", call. = FALSE)
    }
    initial <- as.integer(initial)
  }
  out <- shape_hill_cpp(D, as.integer(initial - 1L), max_iter, isTRUE(from_best))
  list(cluster   = as.integer(out$cluster),
       medoids   = as.integer(out$medoids),
       objective = out$objective,
       n_iter    = as.integer(out$n_iter),
       distance  = if (given_D) NA_character_ else distance)
}
