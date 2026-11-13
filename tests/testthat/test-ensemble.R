# Reference implementation from the authors' 3D simulation study,
# written as in metodos/confirm/round2.R, hybrid/arms.R and cov/f2_arms.R.
.e1_ref <- function(X, K, B = 30L) {
  k <- dim(X)[1]; N <- dim(X)[3]
  nrm <- function(A) { A <- sweep(A, 2, colMeans(A)); A / sqrt(sum(A^2)) }
  mu <- nrm(shape_mean(X)); V <- matrix(0, N, 3 * k)
  for (i in seq_len(N)) {
    A <- nrm(X[, , i]); s <- svd(crossprod(A, mu)); U <- s$u
    if (det(U %*% t(s$v)) < 0) U[, 3] <- -U[, 3]
    Ai <- A %*% (U %*% t(s$v)); V[i, ] <- as.vector(Ai - mu * sum(Ai * mu))
  }
  p <- prcomp(V, center = TRUE, scale. = FALSE, rank. = min(10L, N - 1L))
  cv <- cumsum(p$sdev^2) / sum(p$sdev^2)
  dv <- max(2L, min(10L, which(cv >= 0.90)[1]))
  S <- p$x[, seq_len(dv), drop = FALSE]
  raw <- matrix(0L, N, B)
  for (b in seq_len(B)) {
    idx <- sample.int(N, N, replace = TRUE)
    Sb <- unique(S[idx, , drop = FALSE])
    cen <- suppressWarnings(kmeans(Sb, K, nstart = 1L, iter.max = 30L))$centers
    raw[, b] <- max.col(-as.matrix(dist(rbind(cen, S)))[-(1:K), 1:K, drop = FALSE],
                        ties.method = "first")
  }
  C <- matrix(0, N, N)
  for (m in seq_len(B)) { l <- raw[, m]; C <- C + outer(l, l, "==") }
  C <- C / B
  dg <- rowSums(C); M <- C / sqrt(outer(dg, dg))
  U <- eigen(M, symmetric = TRUE)$vectors[, 1:K, drop = FALSE]
  U <- U / pmax(sqrt(rowSums(U^2)), 1e-12)
  list(labels = kmeans(U, K, nstart = 20)$cluster, dim = dv, C = C)
}

test_that("shape_ensemble reproduces the E1 arm for the same RNG state", {
  for (s in 1:3) {
    sim <- sim_shapes3d(n = 40, sigma = 1 + s / 2, k = 8, seed = s)
    set.seed(100 + s); ref <- .e1_ref(sim$X, 2L)
    set.seed(100 + s); fit <- shape_ensemble(sim$X, 2L)
    expect_identical(as.integer(fit$cluster), as.integer(ref$labels))
    expect_identical(fit$dim, as.integer(ref$dim))
    expect_equal(fit$coassoc, ref$C)
  }
})

test_that("shape_ensemble returns valid structure, representatives and seed handling", {
  sim <- sim_shapes3d(n = 20, sigma = 1, k = 8, seed = 4)
  set.seed(9); s0 <- .Random.seed
  a <- shape_ensemble(sim$X, 2, B = 10, seed = 2)
  expect_identical(.Random.seed, s0)
  b <- shape_ensemble(sim$X, 2, B = 10, seed = 2)
  expect_identical(a$cluster, b$cluster)
  expect_s3_class(a, "shape_ensemble")
  expect_equal(dim(a$coassoc), rep(length(sim$group), 2))
  expect_true(isSymmetric(a$coassoc))
  expect_true(all(diag(a$coassoc) == 1))
  expect_length(a$representative, 2)
  # each representative belongs to its own cluster and is the nearest member
  for (h in 1:2) {
    expect_equal(a$cluster[a$representative[h]], h)
    idx <- which(a$cluster == h)
    d <- shape_distance(sim$X[, , idx], "rho",
                        Y = array(a$centers[, , h], c(8, 3, 1)))
    expect_equal(a$representative[h], idx[which.min(d)])
  }
  expect_output(print(a), "Shape ensemble")
  expect_gte(a$dim, 2); expect_lte(a$dim, 10)
})

test_that("shape_ensemble recovers separated groups and checks arguments", {
  sim <- sim_shapes3d(n = 40, sigma = 0.5, k = 8, seed = 1)
  fit <- shape_ensemble(sim$X, 2, seed = 1)
  expect_gt(rand_index(sim$group, fit$cluster), 0.9)
  expect_error(shape_ensemble(sim$X, 1))
  expect_error(shape_ensemble(sim$X, 100))
  expect_error(shape_ensemble(sim$X[, 1, , drop = FALSE], 2))
  expect_equal(shape_ensemble(sim$X, 2, max_dim = 3, seed = 1)$dim <= 3, TRUE)
})

test_that("shape_ensemble accepts 2D arrays", {
  set.seed(8)
  X <- array(rep(matrix(rnorm(14), 7, 2), 30) + rnorm(7 * 2 * 30, sd = 0.2), c(7, 2, 30))
  fit <- shape_ensemble(X, 2, B = 8, seed = 1)
  expect_length(fit$cluster, 30)
})
