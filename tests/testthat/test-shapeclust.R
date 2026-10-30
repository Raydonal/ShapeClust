library(testthat)
library(ShapeClust)

test_that("Helmert sub-matrix has orthonormal rows and centring identity", {
  for (k in c(3L, 4L, 8L, 24L)) {
    H <- helmert_submatrix(k)
    expect_equal(nrow(H), k - 1L)
    expect_equal(ncol(H), k)
    # Rows orthonormal
    HH <- H %*% t(H)
    expect_equal(HH, diag(k - 1L), tolerance = 1e-10)
    # Centring identity: H^T H = I_k - 1/k * 1 1^T
    C <- diag(k) - matrix(1 / k, k, k)
    expect_equal(t(H) %*% H, C, tolerance = 1e-10)
  }
})

test_that("landmarks_to_preshape produces unit-norm complex columns", {
  set.seed(1)
  k <- 4; n <- 20
  X <- array(rnorm(k * 2 * n), dim = c(k, 2, n))
  Z <- landmarks_to_preshape(X)
  expect_equal(nrow(Z), k - 1L)
  expect_equal(ncol(Z), n)
  expect_true(is.complex(Z))
  norms <- colSums(Mod(Z)^2)
  expect_equal(norms, rep(1, n), tolerance = 1e-10)
})

test_that("landmarks_to_preshape is translation-invariant", {
  set.seed(1)
  k <- 4; n <- 5
  X <- array(rnorm(k * 2 * n), dim = c(k, 2, n))
  # Translate every landmark by the same vector per observation
  shift <- matrix(c(3, -2), nrow = k, ncol = 2, byrow = TRUE)
  X_shift <- X
  for (i in seq_len(n)) X_shift[,,i] <- X[,,i] + shift
  Z1 <- landmarks_to_preshape(X)
  Z2 <- landmarks_to_preshape(X_shift)
  expect_equal(Mod(Z1 - Z2), matrix(0, k - 1L, n), tolerance = 1e-10)
})

test_that("rho_matrix returns a symmetric matrix in [0, pi/2]", {
  set.seed(2)
  k <- 4; n <- 10
  X <- array(rnorm(k * 2 * n), dim = c(k, 2, n))
  Z <- landmarks_to_preshape(X)
  D <- rho_matrix(Z)
  expect_true(isSymmetric(D))
  expect_true(all(diag(D) == 0))
  expect_true(all(D >= 0))
  expect_true(all(D <= pi / 2 + 1e-10))
})

test_that("dF equals sin(rho) numerically", {
  set.seed(3)
  k <- 5; n <- 8
  X <- array(rnorm(k * 2 * n), dim = c(k, 2, n))
  Z <- landmarks_to_preshape(X)
  D <- rho_matrix(Z)
  D_F <- dF_matrix(Z)
  expect_equal(D_F, sin(D), tolerance = 1e-8)
})

test_that("shape_pam returns a valid partition on well-separated data", {
  set.seed(4)
  # Three well-separated clusters
  mu_list <- list(
    matrix(c(0,0, 10,0, 10,10, 0,10), 4, 2, byrow = TRUE),
    matrix(c(0,0, 20,0, 20,10, 0,10), 4, 2, byrow = TRUE),
    matrix(c(0,0, 15,0, 10,10, 0,10), 4, 2, byrow = TRUE)
  )
  sim <- sim_landmark_gaussian(mu_list, c(20, 20, 20),
                               sigma = 0.5, regime = "isotropic", seed = 4)
  Z <- landmarks_to_preshape(sim$X)
  fit <- shape_pam(Z, K = 3)
  expect_length(fit$cluster, 60L)
  expect_true(all(fit$cluster %in% 1:3))
  # Should have decent recovery
  expect_gt(rand_index(sim$group, fit$cluster), 0.8)
})

test_that("shape_rankkmed returns a valid partition on well-separated data", {
  set.seed(5)
  mu_list <- list(
    matrix(c(0,0, 10,0, 10,10, 0,10), 4, 2, byrow = TRUE),
    matrix(c(0,0, 20,0, 20,10, 0,10), 4, 2, byrow = TRUE),
    matrix(c(0,0, 15,0, 10,10, 0,10), 4, 2, byrow = TRUE)
  )
  sim <- sim_landmark_gaussian(mu_list, c(20, 20, 20),
                               sigma = 0.5, regime = "isotropic", seed = 5)
  Z <- landmarks_to_preshape(sim$X)
  fit <- shape_rankkmed(Z, K = 3, m = 5)
  expect_length(fit$cluster, 60L)
  expect_true(all(fit$cluster %in% 1:3))
})

test_that("bagclust1 produces a valid partition and stabilises", {
  set.seed(6)
  mu_list <- list(
    matrix(c(0,0, 10,0, 10,10, 0,10), 4, 2, byrow = TRUE),
    matrix(c(0,0, 20,0, 20,10, 0,10), 4, 2, byrow = TRUE),
    matrix(c(0,0, 15,0, 10,10, 0,10), 4, 2, byrow = TRUE)
  )
  sim <- sim_landmark_gaussian(mu_list, c(20, 20, 20),
                               sigma = 1.5, regime = "isotropic", seed = 6)
  Z <- landmarks_to_preshape(sim$X)
  fit_bag <- bagclust1(Z, K = 3, base = "pam", B = 20, seed = 6)
  expect_s3_class(fit_bag, "bagclust1")
  expect_length(fit_bag$cluster, 60L)
  expect_true(all(fit_bag$cluster %in% 1:3))
  expect_true(all(rowSums(fit_bag$votes) > 0L))
})

test_that("align_labels finds the optimal permutation", {
  set.seed(7)
  truth <- rep(1:3, each = 10)
  # Permuted version
  perm_map <- c(3L, 1L, 2L)
  observed <- perm_map[truth]
  aligned <- align_labels(observed, truth, 3L)
  expect_equal(aligned, truth)
})

test_that("rand_index and fmi_index return 1 for perfect agreement", {
  x <- rep(1:3, each = 10)
  expect_equal(rand_index(x, x), 1)
  expect_equal(fmi_index(x, x), 1)
})

test_that("Complex Bingham draws lie on the pre-shape sphere", {
  Z <- sim_complex_bingham(n = 20, lambda = c(9, 9, 9), seed = 8)
  norms <- colSums(Mod(Z)^2)
  expect_equal(norms, rep(1, 20), tolerance = 1e-10)
})

test_that("select_K recovers K on strongly separated data", {
  set.seed(9)
  mu_list <- list(
    matrix(c(0,0, 10,0, 10,10, 0,10), 4, 2, byrow = TRUE),
    matrix(c(0,0, 20,0, 20,10, 0,10), 4, 2, byrow = TRUE),
    matrix(c(0,0, 15,0, 10,10, 0,10), 4, 2, byrow = TRUE)
  )
  sim <- sim_landmark_gaussian(mu_list, c(30, 30, 30),
                               sigma = 0.3, regime = "isotropic", seed = 9)
  Z <- landmarks_to_preshape(sim$X)
  res <- select_K(Z, K_range = 2:5, method = "silhouette")
  expect_equal(res$K_hat, 3L)
})

test_that("Complex Bingham keeps its concentration pole (groups are separable)", {
  # With lambda = 200 the three rotated Bingham groups must be recoverable.
  # This guards against dropping the pole coordinate, which would turn the
  # sampler into near-uniform noise on the pre-shape sphere.
  set.seed(11)
  d <- 3L                                  # k = 4 landmarks
  ed <- c(0, 0, 1)
  Zs <- lapply(1:3, function(g) {
    Z0 <- sim_complex_bingham(40, lambda = rep(200, d - 1L))
    v <- rep(0, d); v[g] <- 1
    u <- v - ed
    U <- if (sqrt(sum(u^2)) < 1e-12) diag(d) else {
      u <- u / sqrt(sum(u^2)); diag(d) - 2 * outer(u, u)
    }
    U %*% Z0
  })
  Z <- do.call(cbind, Zs)
  Z <- sweep(Z, 2L, sqrt(colSums(Mod(Z)^2)), "/")
  truth <- rep(1:3, each = 40)

  # Sampler must return d = length(lambda) + 1 coordinates
  expect_equal(nrow(Z), d)

  fit <- shape_pam(Z, K = 3)
  expect_gt(rand_index(truth, fit$cluster), 0.95)
})

## ---- Regression tests for the 0.3.0 fixes --------------------------------

three_groups <- function(n_per = 20, sigma = 1, seed = 1) {
  mu_list <- list(
    matrix(c(0,0, 10,0, 10,10, 0,10), 4, 2, byrow = TRUE),
    matrix(c(0,0, 20,0, 20,10, 0,10), 4, 2, byrow = TRUE),
    matrix(c(0,0, 15,0, 10,10, 0,10), 4, 2, byrow = TRUE)
  )
  sim_landmark_gaussian(mu_list, rep(n_per, 3), sigma = sigma, seed = seed)
}

test_that("shape_rankkmed matches kmed::rankkmed (Zadegan et al. 2013)", {
  skip_if_not_installed("kmed")
  sim <- three_groups(sigma = 2, seed = 21)
  D <- rho_matrix(landmarks_to_preshape(sim$X))
  for (s in 1:10) {
    set.seed(s); ours <- shape_rankkmed(D, K = 3, m = 5, max_iter = 10)
    set.seed(s); ref  <- kmed::rankkmed(D, ncluster = 3, m = 5, iterate = 10)
    expect_equal(ours$medoids, as.integer(ref$medoid))
    expect_equal(ours$cluster, as.integer(ref$cluster))
  }
})

test_that("rank_matrix_cpp averages ties like base::rank", {
  D <- as.matrix(dist(c(0, 0, 1, 3, 3, 3)))
  expect_equal(ShapeClust:::rank_matrix_cpp(D), t(apply(D, 1L, rank)),
               ignore_attr = TRUE)
  expect_equal(ShapeClust:::order_matrix_cpp(D), t(apply(D, 1L, order)),
               ignore_attr = TRUE)
})

test_that("bagclust1 votes for every observation in every replicate", {
  Z <- landmarks_to_preshape(three_groups(sigma = 1.5, seed = 22)$X)
  fit <- bagclust1(Z, K = 3, B = 15, seed = 3)
  expect_true(all(rowSums(fit$votes) == 15L))
})

test_that("bagclust1 gives the same result serially and in parallel", {
  Z <- landmarks_to_preshape(three_groups(sigma = 2, seed = 23)$X)
  a <- bagclust1(Z, K = 3, base = "rankkmed", B = 10, seed = 4)
  b <- bagclust1(Z, K = 3, base = "rankkmed", B = 10, seed = 4,
                 parallel = TRUE, n_cores = 2L)
  expect_identical(a$votes, b$votes)
})

test_that("functions with a seed leave the caller's stream untouched", {
  Z <- landmarks_to_preshape(three_groups(seed = 24)$X)
  set.seed(99); x0 <- runif(1)
  set.seed(99); invisible(bagclust1(Z, K = 3, B = 3, seed = 1)); x1 <- runif(1)
  set.seed(99); invisible(sim_complex_bingham(3, rep(9, 3), seed = 1)); x2 <- runif(1)
  expect_identical(c(x1, x2), c(x0, x0))
})

test_that("gap statistic separates structured from uniform data", {
  Z <- landmarks_to_preshape(three_groups(n_per = 25, sigma = 0.3, seed = 25)$X)
  set.seed(10)
  res <- select_K(Z, K_range = 2:5, method = "gap", B_gap = 10)
  expect_gt(res$scores[["3"]], res$scores[["2"]])
  expect_true(all(res$scores > 0))
  expect_error(select_K(rho_matrix(Z), method = "gap"), "pre-shape")
})

test_that("riemannian_variances uses Procrustes means and ranks separation", {
  near <- three_groups(sigma = 3, seed = 26)
  far  <- three_groups(sigma = 0.3, seed = 26)
  v_near <- riemannian_variances(near$X, near$group)
  v_far  <- riemannian_variances(far$X, far$group)
  expect_gt(v_far$ratio, v_near$ratio)
  # The mean of identical pre-shapes is that pre-shape (up to rotation).
  Z1 <- landmarks_to_preshape(far$X[, , c(1, 1, 1)])
  mu <- ShapeClust:::.procrustes_mean(Z1)
  expect_equal(ShapeClust:::rho_pair(mu, Z1[, 1]), 0, tolerance = 1e-6)
})
