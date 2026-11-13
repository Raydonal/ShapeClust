## ---------------------------------------------------------------------------
##  Self-attention smoothing of shapes before clustering.
##
##  Each observation attends to the others with softmax weights of its
##  squared tangent distances, as in a transformer layer with identity
##  query, key and value maps (nothing is trained). The per-row temperature
##  is set so that the attention row has a fixed perplexity, as in t-SNE.
##  The layer replaces each point by (1 - alpha) x_i + alpha sum_j a_ij x_j;
##  stacking layers gives the self-attention particle dynamics, which is
##  blurring mean shift. The smoothed points are then clustered.
## ---------------------------------------------------------------------------

# Row-stochastic attention matrix with zero diagonal and fixed perplexity.
# D2: n x n squared distances. Returns the matrix and the inverse
# temperatures. Bisection on log(beta) for all rows at once.
.attention_matrix <- function(D2, perplexity, tol = 1e-5, max_iter = 100L) {
  n <- nrow(D2)
  diag(D2) <- Inf
  D0 <- D2 - apply(D2, 1L, min)            # row minimum is 0, diagonal Inf
  scale <- apply(D0, 1L, function(r) mean(r[is.finite(r)]))
  scale[!(scale > 0)] <- 1
  target <- log(perplexity)
  lo <- rep(log(1e-6), n) - log(scale)
  hi <- rep(log(1e6), n) - log(scale)
  rowent <- function(lb) {
    W <- exp(-exp(lb) * D0)
    W <- W / rowSums(W)
    L <- W * log(W); L[W == 0] <- 0
    list(W = W, H = -rowSums(L))
  }
  for (it in seq_len(max_iter)) {
    mid <- (lo + hi) / 2
    H <- rowent(mid)$H
    big <- H > target                      # too flat: raise beta
    lo[big] <- mid[big]; hi[!big] <- mid[!big]
    if (max(hi - lo) < tol) break
  }
  lb <- (lo + hi) / 2
  list(A = rowent(lb)$W, beta = exp(lb))
}

# One or more attention layers on the rows of H.
.attention_smooth <- function(H, perplexity, layers = 1L, alpha = 1) {
  A1 <- NULL
  for (l in seq_len(layers)) {
    G <- tcrossprod(H)
    sq <- diag(G)
    D2 <- pmax(outer(sq, sq, "+") - 2 * G, 0)
    am <- .attention_matrix(D2, perplexity)
    if (l == 1L) A1 <- am
    H <- (1 - alpha) * H + alpha * (am$A %*% H)
  }
  list(H = H, A = A1$A, beta = A1$beta)
}

# Default perplexity: half the expected cluster size, between 5 and 30.
# With alpha = 0.5 and one layer this was the least harmful setting in the
# Euclidean pilot, where no setting
# beat k-means with 10 starts.
.default_perplexity <- function(n, K) {
  as.numeric(min(30, max(5, floor(n / (2 * K)))))
}

#' Self-attention smoothing and clustering of shapes
#'
#' Replaces each configuration by an attention-weighted average of the
#' others in tangent coordinates, then clusters the smoothed points. The
#' attention is the scaled dot-product attention of transformers with no
#' trained projections: queries, keys and values are the tangent coordinates
#' themselves, and the weights are a softmax of minus the squared distances.
#'
#' @param X Numeric array of dimension \eqn{k \times m \times n}, with
#'   \eqn{m} equal to 2 or 3.
#' @param K Number of clusters (>= 2).
#' @param perplexity Effective number of neighbors each observation attends
#'   to. Each row's temperature is set by bisection so that
#'   \eqn{\exp(-\sum_j a_{ij}\log a_{ij})} equals this value. Default
#'   \code{min(30, max(5, floor(n / (2 K))))}.
#' @param layers Number of attention layers (default 1). Weights are
#'   recomputed from the smoothed points at each layer.
#' @param alpha Residual weight in \eqn{(0, 1]}: each layer returns
#'   \eqn{(1-\alpha) x_i + \alpha \sum_j a_{ij} x_j}. Default 0.5.
#' @param reduce \code{"pca"} (default) works in the leading principal
#'   components of the tangent coordinates, chosen as in
#'   \code{\link{shape_ensemble}}; \code{"none"} uses all tangent
#'   coordinates.
#' @param var_explained,max_dim PCA rule, as in \code{\link{shape_ensemble}}.
#' @param base Clusterer applied to the smoothed points: \code{"kmeans"}
#'   (default, \code{\link[stats]{kmeans}} with \code{nstart} starts) or
#'   \code{"pam"} (\code{\link{shape_pam}} on their Euclidean distances).
#' @param nstart Number of k-means starts (default 10).
#' @param seed Optional integer seed. The caller's random-number stream is
#'   restored on exit.
#'
#' @details
#' The attention matrix has zero diagonal, so each point is replaced by an
#' average of other points; with \code{alpha < 1} part of the point is kept,
#' as in the residual connection of a transformer. Stacking layers with
#' recomputed weights is the blurring mean shift of Cheng (1995), and it is
#' also the self-attention dynamics studied by Geshkovski et al. (2023),
#' under which points merge into clusters.
#'
#' The smoothing helps when the attention mass of a point falls mostly on
#' its own group. The returned \code{self_mass} is the average share of
#' attention that falls inside the point's final cluster, a label-free proxy
#' for that condition. Values near \eqn{1/K} mean the attention did not
#' separate the groups.
#'
#' Under a mixture of independent observations no function of the whole
#' sample classifies a point better than the Bayes rule applied to that
#' point alone, so the smoothing cannot remove overlap between groups; it can
#' only reduce the loss due to estimating centers and searching partitions.
#' In Gaussian simulations it did not beat k-means with 10 starts in any
#' setting tried. The function is experimental.
#'
#' @return An object of class \code{"shape_attention"}, a list with
#'   \describe{
#'     \item{cluster}{Integer vector of length n, values in 1:K.}
#'     \item{smoothed}{n x d matrix of smoothed coordinates.}
#'     \item{attention}{n x n attention matrix of the first layer.}
#'     \item{beta}{Inverse temperature of each row in the first layer.}
#'     \item{self_mass}{Mean first-layer attention inside the own cluster.}
#'     \item{dim}{Number of coordinates used.}
#'     \item{perplexity, layers, alpha, base, K}{The options used.}
#'   }
#'
#' @references
#' Vaswani, A. et al. (2017). Attention is all you need. NeurIPS 30.
#'
#' Cheng, Y. (1995). Mean shift, mode seeking, and clustering. IEEE TPAMI
#' 17(8), 790--799.
#'
#' Geshkovski, B., Letrouit, C., Polyanskiy, Y. and Rigollet, P. (2023).
#' The emergence of clusters in self-attention dynamics. NeurIPS 36.
#'
#' @examples
#' sim <- sim_shapes3d(n = 60, sigma = 1, k = 8, seed = 1)
#' fit <- shape_attention(sim$X, K = 2, seed = 1)
#' fit
#' table(sim$group, fit$cluster)
#'
#' @seealso \code{\link{shape_ensemble}}, \code{\link{shape_kmeans}}
#' @export
shape_attention <- function(X, K, perplexity = NULL, layers = 1L, alpha = 0.5,
                            reduce = c("pca", "none"), var_explained = 0.90,
                            max_dim = 10L, base = c("kmeans", "pam"),
                            nstart = 10L, seed = NULL) {
  X <- .check_array(X, "X", min_n = 4L)
  reduce <- match.arg(reduce); base <- match.arg(base)
  n <- dim(X)[3L]
  stopifnot(is.numeric(K), length(K) == 1L, K >= 2, K < n,
            is.numeric(layers), length(layers) == 1L, layers >= 1,
            is.numeric(alpha), length(alpha) == 1L, alpha > 0, alpha <= 1,
            is.numeric(nstart), length(nstart) == 1L, nstart >= 1)
  K <- as.integer(K); layers <- as.integer(layers)
  if (is.null(perplexity)) perplexity <- .default_perplexity(n, K)
  stopifnot(is.numeric(perplexity), length(perplexity) == 1L,
            perplexity >= 1, perplexity < n - 1)
  .local_seed(seed)

  if (reduce == "pca") {
    tp <- .tangent_pca(X, var_explained, max(max_dim, K - 1L), max(2L, K - 1L))
    H0 <- tp$scores[, seq_len(tp$dim), drop = FALSE]
  } else {
    H0 <- scale(.tangent_coords(X), center = TRUE, scale = FALSE)
    attr(H0, "scaled:center") <- NULL
  }
  sm <- .attention_smooth(H0, perplexity, layers, alpha)

  cl <- if (base == "kmeans") {
    stats::kmeans(sm$H, K, nstart = as.integer(nstart), iter.max = 30L)$cluster
  } else {
    shape_pam(as.matrix(stats::dist(sm$H)), K)$cluster
  }
  self_mass <- mean(rowSums(sm$A * outer(cl, cl, "==")))

  structure(list(cluster = as.integer(cl), smoothed = sm$H,
                 attention = sm$A, beta = sm$beta, self_mass = self_mass,
                 dim = ncol(H0), perplexity = perplexity, layers = layers,
                 alpha = alpha, base = base, K = K),
            class = "shape_attention")
}

#' @export
print.shape_attention <- function(x, ...) {
  cat("Shape attention clustering\n")
  cat("  n =", length(x$cluster), " K =", x$K, " dim =", x$dim,
      " perplexity =", x$perplexity, " layers =", x$layers,
      " alpha =", x$alpha, " base =", x$base, "\n")
  cat("  cluster sizes:", paste(tabulate(x$cluster, x$K), collapse = " "), "\n")
  cat(sprintf("  attention mass inside own cluster: %.3f (1/K = %.3f)\n",
              x$self_mass, 1 / x$K))
  invisible(x)
}
