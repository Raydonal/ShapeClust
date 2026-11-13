test_that("hopkins_shape rejects for separated groups", {
  sim <- sim_shapes3d(n = 30, sigma = 1, k = 8, seed = 1)
  h <- hopkins_shape(sim$X, B = 99, seed = 2)
  expect_s3_class(h, "hopkins_shape")
  expect_length(h$reference, 99)
  expect_lt(h$p.value, 0.1)
  expect_output(print(h), "Hopkins")
})

test_that("hopkins_shape type I error is not grossly inflated", {
  skip_on_cran()
  set.seed(11)
  mu <- matrix(rnorm(16), 8, 2 + 1)[, 1:3]
  p <- vapply(1:40, function(i) {
    X <- array(rep(mu, 60) + rnorm(8 * 3 * 60, sd = 0.3), c(8, 3, 60))
    hopkins_shape(X, B = 49)$p.value
  }, numeric(1))
  expect_lt(mean(p <= 0.1), 0.25)
})

test_that("hopkins_shape is reproducible and keeps the RNG stream", {
  X <- sim_shapes3d(n = 20, sigma = 2, k = 8, seed = 3)$X
  set.seed(5); s0 <- .Random.seed
  a <- hopkins_shape(X, B = 9, seed = 4)
  expect_identical(.Random.seed, s0)
  b <- hopkins_shape(X, B = 9, seed = 4)
  expect_identical(a$reference, b$reference)
  expect_identical(a$statistic, b$statistic)
  u <- hopkins_shape(X, null = "uniform", B = 9, seed = 4)
  expect_true(u$p.value > 0 && u$p.value <= 1)
})

test_that("hopkins_shape works for 2D input and checks arguments", {
  set.seed(8)
  X <- array(rep(matrix(rnorm(14), 7, 2), 25) + rnorm(7 * 2 * 25, sd = 0.2),
             c(7, 2, 25))
  h <- hopkins_shape(X, B = 19, seed = 1)
  expect_true(h$statistic > 0 && h$statistic < 1)
  expect_error(hopkins_shape(X, m = 30))
})
