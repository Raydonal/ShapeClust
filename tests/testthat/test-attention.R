test_that("attention rows are stochastic, zero-diagonal and hit the perplexity", {
  set.seed(1)
  H <- matrix(rnorm(60 * 4), 60, 4)
  G <- tcrossprod(H); sq <- diag(G)
  D2 <- pmax(outer(sq, sq, "+") - 2 * G, 0)
  for (P in c(3, 10, 25)) {
    am <- ShapeClust:::.attention_matrix(D2, P)
    expect_equal(rowSums(am$A), rep(1, 60), tolerance = 1e-12)
    expect_true(all(diag(am$A) == 0))
    L <- am$A * log(am$A); L[am$A == 0] <- 0
    expect_equal(exp(-rowSums(L)), rep(P, 60), tolerance = 1e-3)
  }
})

test_that("one layer with alpha = 1 is the attention average; alpha mixes", {
  set.seed(2)
  H <- matrix(rnorm(30 * 3), 30, 3)
  s1 <- ShapeClust:::.attention_smooth(H, 5, 1L, 1)
  expect_equal(s1$H, s1$A %*% H)
  s2 <- ShapeClust:::.attention_smooth(H, 5, 1L, 0.3)
  expect_equal(s2$H, 0.7 * H + 0.3 * (s2$A %*% H))
  # the attention average stays in the convex hull: shrinks the spread
  expect_lt(sum(scale(s1$H, scale = FALSE)^2), sum(scale(H, scale = FALSE)^2))
})

test_that("shape_attention returns valid structure and handles the seed", {
  sim <- sim_shapes3d(n = 20, sigma = 1, k = 8, seed = 3)   # 2 x 20 = 40 shapes
  set.seed(9); s0 <- .Random.seed
  a <- shape_attention(sim$X, 2, seed = 5)
  expect_identical(.Random.seed, s0)
  b <- shape_attention(sim$X, 2, seed = 5)
  expect_identical(a$cluster, b$cluster)
  expect_s3_class(a, "shape_attention")
  expect_length(a$cluster, 40)
  expect_true(all(a$cluster %in% 1:2))
  expect_equal(dim(a$attention), c(40, 40))
  expect_gte(a$self_mass, 0); expect_lte(a$self_mass, 1)
  expect_equal(a$perplexity, ShapeClust:::.default_perplexity(40, 2))
  expect_output(print(a), "Shape attention")
  expect_equal(ncol(shape_attention(sim$X, 2, reduce = "none", seed = 1)$smoothed), 24)
  expect_length(shape_attention(sim$X, 2, base = "pam", layers = 2,
                                alpha = 0.5, seed = 1)$cluster, 40)
})

test_that("shape_attention recovers separated groups and checks arguments", {
  sim <- sim_shapes3d(n = 30, sigma = 0.5, k = 8, seed = 1)  # 60 shapes
  fit <- shape_attention(sim$X, 2, seed = 1)
  expect_gt(rand_index(sim$group, fit$cluster), 0.9)
  expect_gt(fit$self_mass, 0.9)
  expect_error(shape_attention(sim$X, 1))
  expect_error(shape_attention(sim$X, 2, alpha = 0))
  expect_error(shape_attention(sim$X, 2, perplexity = 59))
})

test_that("the tangent-coordinate refactor leaves shape_ensemble unchanged", {
  sim <- sim_shapes3d(n = 15, sigma = 1, k = 8, seed = 2)
  tp <- ShapeClust:::.tangent_pca(sim$X, 0.9, 10L)
  V <- ShapeClust:::.tangent_coords(sim$X)
  p <- prcomp(V, center = TRUE, scale. = FALSE, rank. = 10L)
  expect_equal(tp$scores, p$x)
})
