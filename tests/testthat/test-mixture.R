adj_rand_for_test <- function(a, b) {
  t <- table(a, b); c2 <- function(x) x * (x - 1) / 2
  s <- sum(c2(t)); ra <- sum(c2(rowSums(t))); cb <- sum(c2(colSums(t)))
  e <- ra * cb / c2(length(a)); (s - e) / ((ra + cb) / 2 - e)
}

# Unequal-spread groups: a compact cluster and a diffuse one close to it, in
# a few dimensions, embedded in a landmark array (tangent PCA recovers them).
.unequal_array <- function(seed, n1 = 60, n2 = 90, k = 8, s1 = 0.02, s2 = 0.12) {
  set.seed(seed)
  base <- matrix(rnorm(k * 3), k, 3)
  dir <- matrix(rnorm(k * 3), k, 3); dir <- dir / sqrt(sum(dir^2))
  mk <- function(n, s, shift) {
    out <- array(0, c(k, 3, n))
    for (i in seq_len(n)) out[, , i] <- base + shift * dir + s * matrix(rnorm(k * 3), k, 3)
    out
  }
  X <- array(c(mk(n1, s1, 0), mk(n2, s2, 0.12)), c(k, 3, n1 + n2))
  list(X = X, group = rep(1:2, c(n1, n2)))
}

test_that("shape_mixture recovers unequal-spread groups where k-means fails", {
  ari_mix <- ari_km <- numeric(0)
  for (s in 1:5) {
    d <- .unequal_array(s)
    m <- shape_mixture(d$X, 2, seed = 1)
    km <- shape_kmeans(d$X, 2, niter = 10)
    ari_mix[s] <- adj_rand_for_test(d$group, m$cluster)
    ari_km[s] <- adj_rand_for_test(d$group, km$cluster)
    expect_false(m$fallback)
  }
  expect_gt(mean(ari_mix), mean(ari_km) + 0.1)
  expect_gt(mean(ari_mix), 0.8)
})

test_that("shape_mixture returns valid structure and is reproducible", {
  d <- .unequal_array(3)
  set.seed(9); s0 <- .Random.seed
  a <- shape_mixture(d$X, 2, seed = 5)
  expect_identical(.Random.seed, s0)
  b <- shape_mixture(d$X, 2, seed = 5)
  expect_identical(a$cluster, b$cluster)
  expect_equal(a$loglik, b$loglik)
  expect_s3_class(a, "shape_mixture")
  n <- length(d$group)
  expect_equal(dim(a$posterior), c(n, 2)); expect_equal(rowSums(a$posterior), rep(1, n))
  expect_equal(sum(a$proportions), 1)
  expect_true(all(a$variances > 0))
  expect_equal(dim(a$means), c(2, a$dim))
  expect_gte(a$dim, 2); expect_lte(a$dim, 10)
  for (h in 1:2) expect_equal(a$cluster[a$representative[h]], h)
  expect_true(a$converged)
  expect_output(print(a), "Shape mixture")
  expect_error(shape_mixture(d$X, 1)); expect_error(shape_mixture(d$X, 1000))
  expect_error(shape_mixture(d$X, 2, model = "VVV"))
})

test_that("shape_mixture handles duplicated rows and small n", {
  d <- .unequal_array(2, n1 = 10, n2 = 12)
  idx <- sample.int(22, 22, replace = TRUE); set.seed(1)
  Xd <- d$X[, , c(idx, idx[1:5], rep(1L, 4))]
  expect_error(fit <- shape_mixture(Xd, 2, seed = 1), NA)
  expect_true(all(is.finite(fit$variances)) && all(fit$variances > 0))
  # n = 6 with K = 3 and a fully duplicated group
  X6 <- d$X[, , c(1, 1, 1, 11, 12, 13)]
  expect_error(f6 <- shape_mixture(X6, 3, seed = 2), NA)
  expect_length(f6$cluster, 6)
  X4 <- d$X[, , c(1, 2, 11, 12)]
  expect_error(shape_mixture(X4, 2, seed = 3), NA)
  # all identical rows: fallback is flagged, no error
  Xi <- d$X[, , c(1, 1, 1, 1, 2, 2)]
  expect_error(fi <- shape_mixture(Xi, 2, seed = 1), NA)
  expect_type(fi$fallback, "logical")
})

test_that("internal EM agrees closely with mclust VII on a clean example", {
  skip_if_not_installed("mclust")
  suppressMessages(library(mclust))  # Mclust() needs the package attached
  set.seed(11)
  S <- rbind(matrix(rnorm(120, 0, 0.3), 60, 2),
             matrix(rnorm(160, 0, 1.2), 80, 2) + rep(c(3, 3), each = 80))
  e <- ShapeClust:::.vii_em(S, 2, nstart = 10, tol = 1e-10)
  m <- suppressWarnings(mclust::Mclust(S, G = 2, modelNames = "VII", verbose = FALSE))
  expect_gt(adj_rand_for_test(m$classification, max.col(e$z)), 0.95)
  expect_equal(e$loglik, m$loglik, tolerance = 1e-3)
  expect_equal(sort(e$s2), sort(unname(m$parameters$variance$sigmasq)), tolerance = 1e-2)
})

test_that("shape_ensemble base = 'vii' works, and 'kmeans' stays the default", {
  d <- .unequal_array(4)
  a <- shape_ensemble(d$X, 2, B = 8, seed = 3)
  b <- shape_ensemble(d$X, 2, B = 8, seed = 3, base = "kmeans")
  expect_identical(a$cluster, b$cluster); expect_equal(a$coassoc, b$coassoc)
  v <- shape_ensemble(d$X, 2, B = 8, seed = 3, base = "vii")
  v2 <- shape_ensemble(d$X, 2, B = 8, seed = 3, base = "vii")
  expect_identical(v$cluster, v2$cluster)
  expect_equal(v$base, "vii"); expect_gte(v$n_fallback, 0)
  expect_true(isSymmetric(v$coassoc)); expect_true(all(diag(v$coassoc) == 1))
  expect_gt(adj_rand_for_test(d$group, v$cluster), 0.8)
  expect_error(shape_ensemble(d$X, 2, base = "pam"))
})
