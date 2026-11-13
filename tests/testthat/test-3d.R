## Tests for the real-array (m = 2 or 3) geometry, clusterers, simulator and
## the generalized bagclust1().

# ---- helpers ---------------------------------------------------------------

rand_rotation <- function(m = 3L) {
  Q <- qr.Q(qr(matrix(rnorm(m * m), m, m)))
  if (det(Q) < 0) Q[, 1L] <- -Q[, 1L]
  Q
}

# Apply an independent random similarity (rotation, scale, translation) to
# every configuration.
random_similarity <- function(X) {
  m <- dim(X)[2L]
  for (i in seq_len(dim(X)[3L])) {
    M <- X[, , i] %*% rand_rotation(m) * runif(1, 0.3, 4)
    X[, , i] <- sweep(M, 2L, rnorm(m, sd = 5), "+")
  }
  X
}

mirror <- function(X) { X[, 1L, ] <- -X[, 1L, ]; X }

# Distance matrices of shapes, entry by entry
shapes_matrix <- function(X, f) {
  n <- dim(X)[3L]
  D <- matrix(0, n, n)
  for (i in seq_len(n)) for (j in seq_len(n)) D[i, j] <- f(X[, , i], X[, , j])
  D
}

# A transcription, on a distance matrix, of hill_climbing_shapes() of the 3D
# scripts: all random draws go through sample() and max.col().
ref_hill_D <- function(D, K, max_iterations, from_best = FALSE) {
  n <- nrow(D)
  crit <- function(cand) {
    md <- D[, cand, drop = FALSE]
    list(belong = max.col(-md), est = sum(apply(md, 1L, min)), candidates = cand)
  }
  initial <- sample(1:n, K)
  best <- crit(initial)
  for (iteration in 1:max_iterations) {
    new_solution <- if (from_best) best$candidates else initial
    random_index <- sample(1:K, 1)
    new_neighbor <- sample(1:n, 1)
    new_solution[random_index] <- new_neighbor
    while (any(duplicated(new_solution))) {
      random_index <- sample(1:K, 1)
      new_neighbor <- sample(1:n, 1)
      new_solution[random_index] <- new_neighbor
    }
    cand <- crit(new_solution)
    if (cand$est < best$est) best <- cand
  }
  list(medoids = best$candidates, cluster = best$belong, objective = best$est)
}

two_groups <- function(n = 10, sigma = 0.5, seed = 1, ...) {
  sim_shapes3d(n = n, sigma = sigma, seed = seed, ...)
}

# ---- distances -------------------------------------------------------------

test_that("shape_distance matches shapes::riemdist and shapes::procdist in 3D", {
  skip_if_not_installed("shapes")
  set.seed(101)
  for (k in c(4L, 8L)) {
    X <- array(rnorm(k * 3 * 6), dim = c(k, 3, 6))
    rho <- shapes_matrix(X, function(a, b) shapes::riemdist(a, b))
    expect_equal(shape_distance(X, "rho"), rho, tolerance = 1e-7)
    expect_equal(shape_distance(X, "dF"),
                 shapes_matrix(X, function(a, b) shapes::procdist(a, b)),
                 tolerance = 1e-7)
    expect_equal(shape_distance(X, "dP"),
                 shapes_matrix(X, function(a, b)
                   shapes::procdist(a, b, type = "partial")),
                 tolerance = 1e-7)
    expect_equal(shape_distance(X, "rho", reflection = TRUE),
                 shapes_matrix(X, function(a, b)
                   shapes::riemdist(a, b, reflect = TRUE)),
                 tolerance = 1e-7)
  }
})

test_that("shape_distance returns a valid distance matrix", {
  set.seed(102)
  X <- array(rnorm(8 * 3 * 7), dim = c(8, 3, 7))
  for (d in c("rho", "dF", "dP")) {
    D <- shape_distance(X, d)
    expect_true(isSymmetric(D))
    expect_true(all(diag(D) == 0))
    expect_true(all(D >= 0))
  }
  D <- shape_distance(X)
  expect_true(all(D <= pi / 2 + 1e-12))
  expect_equal(shape_distance(X, "dF"), sin(D), tolerance = 1e-12)
  # d_P = sqrt(2 (1 - cos rho))
  expect_equal(shape_distance(X, "dP"), sqrt(2 * (1 - cos(D))), tolerance = 1e-12)
  # cross-distances are blocks of the full matrix
  expect_equal(shape_distance(X[, , 1:3], Y = X[, , 4:7]), D[1:3, 4:7],
               tolerance = 1e-12)
  expect_equal(dim(shape_distance(X[, , 2, drop = FALSE], Y = X)), c(1L, 7L))
})

test_that("shape_distance in 2D equals the complex rho and dF matrices", {
  set.seed(103)
  X <- array(rnorm(6 * 2 * 9), dim = c(6, 2, 9))
  Z <- landmarks_to_preshape(X)
  expect_equal(shape_distance(X, "rho"), rho_matrix(Z), tolerance = 1e-10)
  expect_equal(shape_distance(X, "dF"), dF_matrix(Z), tolerance = 1e-10)
})

test_that("shape_distance is invariant to rotation, translation and scale", {
  set.seed(104)
  for (m in 2:3) {
    X <- array(rnorm(7 * m * 6), dim = c(7, m, 6))
    Xs <- random_similarity(X)
    for (d in c("rho", "dF", "dP")) {
      expect_equal(shape_distance(X, d), shape_distance(Xs, d), tolerance = 1e-7)
    }
    # each configuration is at distance ~0 from its own transformed copy
    own <- diag(shape_distance(X, Y = Xs))
    expect_true(all(own < 1e-6))
  }
})

test_that("reflections are kept unless reflection = TRUE", {
  set.seed(105)
  for (m in 2:3) {
    X <- array(rnorm(8 * m * 5), dim = c(8, m, 5))
    Xm <- mirror(X)
    own <- diag(shape_distance(X, Y = Xm))
    expect_true(all(own > 1e-3))
    own_r <- diag(shape_distance(X, Y = Xm, reflection = TRUE))
    expect_true(all(own_r < 1e-6))
    # allowing reflections can only reduce distances
    expect_true(all(shape_distance(X, reflection = TRUE) <=
                      shape_distance(X) + 1e-12))
    # a planar configuration is its own mirror image up to rotation
  }
  P <- array(rnorm(6 * 3 * 3), dim = c(6, 3, 3)); P[, 3L, ] <- 0
  expect_true(all(diag(shape_distance(P, Y = mirror(P))) < 1e-6))
})

test_that("shape_distance validates its input", {
  expect_error(shape_distance(matrix(1, 4, 3)), "dimension")
  expect_error(shape_distance(array(rnorm(40), c(5, 4, 2))), "m = 2 or m = 3")
  X <- array(rnorm(5 * 3 * 2), c(5, 3, 2)); X[1] <- NA
  expect_error(shape_distance(X), "missing")
  expect_error(shape_distance(array(rnorm(30), c(5, 3, 2)),
                              Y = array(rnorm(24), c(4, 3, 2))), "same k")
})

# ---- means -----------------------------------------------------------------

noisy_sample <- function(m = 3L, k = 8L, n = 25L, sd = 0.1, seed = 1L) {
  set.seed(seed)
  base <- matrix(rnorm(k * m), k, m)
  array(rep(base, n), dim = c(k, m, n)) + array(rnorm(k * m * n, sd = sd), c(k, m, n))
}

test_that("shape_mean('procrustes') agrees with shapes::procGPA", {
  skip_if_not_installed("shapes")
  for (m in 2:3) {
    X <- noisy_sample(m = m, sd = 0.25, seed = 7L + m)
    mu <- shape_mean(X)
    expect_equal(dim(mu), c(8L, m))
    expect_equal(sqrt(sum(mu^2)), 1, tolerance = 1e-12)
    expect_equal(colMeans(mu), rep(0, m), tolerance = 1e-12, ignore_attr = TRUE)
    g <- shapes::procGPA(X)$mshape
    expect_lt(shapes::riemdist(mu, g), 1e-3)
    # same minimizer of the sum of squared full Procrustes distances
    f <- function(M) sum(apply(X, 3L, function(x) sin(shapes::riemdist(x, M))^2))
    expect_lte(f(mu), f(g) + 1e-9)
  }
})

test_that("shape_mean('procrustes') in 2D is the dominant eigenvector", {
  X <- noisy_sample(m = 2L, k = 6L, sd = 0.3, seed = 3L)
  Z <- landmarks_to_preshape(X)
  ev <- eigen(Z %*% Conj(t(Z)))$vectors[, 1L]
  mu <- shape_mean(X)
  zmu <- landmarks_to_preshape(array(mu, dim = c(6, 2, 1)))
  expect_equal(Mod(sum(Conj(zmu[, 1L]) * ev)), 1, tolerance = 1e-8)
})

test_that("shape_mean('phi') agrees with shapes::MDSshape", {
  skip_if_not_installed("shapes")
  for (m in 2:3) {
    X <- noisy_sample(m = m, sd = 0.25, seed = 20L + m)
    nu <- shape_mean(X, type = "phi")
    md <- shapes::MDSshape(X)
    expect_equal(sqrt(sum(nu^2)), 1, tolerance = 1e-12)
    # same shape, same orientation class (shapes are compared without reflection)
    expect_lt(shapes::riemdist(nu, md), 1e-6)
    expect_lt(shapes::riemdist(nu, md, reflect = TRUE), 1e-6)
  }
})

test_that("shape_mean is a function of the shapes only", {
  X <- noisy_sample(sd = 0.2, seed = 5L)
  set.seed(5)
  Xs <- random_similarity(X)
  expect_lt(shape_distance(array(c(shape_mean(X), shape_mean(Xs)), c(8, 3, 2)))[1, 2], 1e-6)
  expect_lt(shape_distance(array(c(shape_mean(X, "phi"), shape_mean(Xs, "phi")),
                                 c(8, 3, 2)))[1, 2], 1e-6)
  # one configuration is its own mean
  one <- shape_mean(X[, , 1, drop = FALSE])
  expect_lt(shape_distance(array(c(one, X[, , 1]), c(8, 3, 2)))[1, 2], 1e-6)
})

# ---- simulator -------------------------------------------------------------

test_that("sim_shapes3d returns the documented structure", {
  for (k in c(8L, 34L)) for (type in c("iso", "aniso")) {
    s <- sim_shapes3d(n = 4, sigma = 1, type = type, k = k, seed = 1)
    expect_equal(dim(s$X), c(k, 3L, 8L))
    expect_equal(dim(s$X_norot), c(k, 3L, 8L))
    expect_equal(s$group, rep(1:2, each = 4))
    expect_equal(s$K, 2L)
  }
  s <- sim_shapes3d(n = c(3, 5), seed = 1)
  expect_equal(dim(s$X)[3L], 8L)
  expect_equal(s$group, c(rep(1L, 3), rep(2L, 5)))
  s1 <- sim_shapes3d(n = 4, groups = "cube", seed = 1)
  expect_equal(s1$K, 1L)
  expect_error(sim_shapes3d(n = 4, k = 10), "8 or 34")
})

test_that("sim_shapes3d is reproducible and leaves the caller's stream alone", {
  a <- sim_shapes3d(n = 5, sigma = 2, type = "aniso", seed = 12)
  b <- sim_shapes3d(n = 5, sigma = 2, type = "aniso", seed = 12)
  expect_identical(a$X, b$X)
  expect_false(identical(a$X, sim_shapes3d(n = 5, sigma = 2, type = "aniso", seed = 13)$X))
  set.seed(77); x0 <- runif(1)
  set.seed(77); invisible(sim_shapes3d(n = 3, seed = 1)); x1 <- runif(1)
  set.seed(77); invisible(sim_shapes3d(n = 3, seed = NULL)); x2 <- runif(1)
  set.seed(77); invisible(sim_shapes3d(n = 3, seed = 1, independent = TRUE)); x3 <- runif(1)
  expect_identical(c(x1, x3), c(x0, x0))
  # without a seed the stream advances by exactly the seed draw
  set.seed(77); invisible(sample.int(.Machine$integer.max, 1L)); x4 <- runif(1)
  expect_identical(x2, x4)
})

test_that("sim_shapes3d reproduces the scripts: shared noise across groups", {
  # The scripts reseed before each group, so both groups get the same noise.
  s <- sim_shapes3d(n = 6, sigma = 2, type = "iso", k = 8, seed = 5, rotate = FALSE)
  e1 <- s$X_norot[, , 1:6] - array(ShapeClust:::.shape_template("cube", 8), c(8, 3, 6))
  e2 <- s$X_norot[, , 7:12] - array(ShapeClust:::.shape_template("parallelepiped", 8), c(8, 3, 6))
  expect_equal(e1, e2, tolerance = 1e-12)
  # independent = TRUE gives different noise
  s2 <- sim_shapes3d(n = 6, sigma = 2, type = "iso", k = 8, seed = 5,
                     rotate = FALSE, independent = TRUE)
  f1 <- s2$X_norot[, , 1:6] - array(ShapeClust:::.shape_template("cube", 8), c(8, 3, 6))
  f2 <- s2$X_norot[, , 7:12] - array(ShapeClust:::.shape_template("parallelepiped", 8), c(8, 3, 6))
  expect_gt(max(abs(f1 - f2)), 0.1)
})

test_that("the multivariate t and von Mises draws reproduce mvtnorm and CircStats", {
  skip_if_not_installed("mvtnorm")
  k <- 4L
  S <- kronecker(2 * (matrix(1, k, k) + 3 * diag(k)),
                 (matrix(1, 3, 3) + 3 * diag(3)) / 16)
  for (s in 1:3) {
    set.seed(s); a <- mvtnorm::rmvt(6, S, df = 99)
    set.seed(s); b <- ShapeClust:::.rmvt_shifted(6, S, 99)
    expect_equal(a, b, tolerance = 1e-10)
  }
  # reference values computed with CircStats::rvm 0.2-7
  set.seed(3); expect_equal(ShapeClust:::.rvm(1, 1234567, -0.141592653589793),
                            2.16586641299739, tolerance = 1e-8)
  set.seed(3); expect_equal(ShapeClust:::.rvm(1, 1234567, 2.14159265358979),
                            2.58560859349757, tolerance = 1e-8)
  set.seed(10); expect_equal(ShapeClust:::.rvm(1, 5, 2), 4.31029589589889,
                             tolerance = 1e-8)
  set.seed(11); expect_equal(ShapeClust:::.rvm(1, 987654321, -0.1415927),
                             2.57368090358796, tolerance = 1e-8)
})

test_that("built-in templates have the documented vertices", {
  cube <- ShapeClust:::.shape_template("cube", 8)
  expect_true(all(cube %in% c(0, 10)))
  expect_equal(nrow(unique(cube)), 8L)
  par <- ShapeClust:::.shape_template("parallelepiped", 8)
  expect_equal(sort(unique(as.vector(par))), c(0, 10, 20))
  expect_equal(dim(ShapeClust:::.shape_template("cube", 34)), c(34L, 3L))
  expect_equal(dim(ShapeClust:::.shape_template("parallelepiped", 34)), c(34L, 3L))
  # the cube and the box are different shapes
  expect_gt(shape_distance(array(c(cube, par), c(8, 3, 2)))[1, 2], 0.1)
})

# ---- clusterers ------------------------------------------------------------

test_that("shape_kmeans recovers well-separated groups (both means)", {
  sim <- two_groups(n = 10, sigma = 0.5, seed = 1)
  for (mean in c("procrustes", "phi")) {
    set.seed(3)
    fit <- shape_kmeans(sim$X, K = 2, mean = mean)
    expect_length(fit$cluster, 20L)
    expect_equal(dim(fit$centers), c(8L, 3L, 2L))
    expect_equal(rand_index(sim$group, fit$cluster), 1)
    expect_true(is.finite(fit$objective) && fit$objective >= 0)
    expect_equal(length(fit$initials), 10L)
  }
  # 34 landmarks, anisotropic noise, and a 2D example
  sim <- two_groups(n = 8, sigma = 1, seed = 2, type = "aniso", k = 34)
  set.seed(4)
  expect_equal(rand_index(sim$group, shape_kmeans(sim$X, 2)$cluster), 1)
  mu_list <- list(matrix(c(0,0, 10,0, 10,10, 0,10), 4, 2, byrow = TRUE),
                  matrix(c(0,0, 20,0, 15,10, 0,10), 4, 2, byrow = TRUE))
  g <- sim_landmark_gaussian(mu_list, c(15, 15), sigma = 0.3,
                             regime = "isotropic", seed = 1)
  set.seed(1)
  expect_gt(rand_index(g$group, shape_kmeans(g$X, 2)$cluster), 0.95)
})

test_that("shape_kmeans is deterministic given the seed and its objective is consistent", {
  sim <- two_groups(n = 8, sigma = 2, seed = 3)
  set.seed(9); a <- shape_kmeans(sim$X, 3)
  set.seed(9); b <- shape_kmeans(sim$X, 3)
  expect_identical(a, b)
  expect_true(all(a$cluster %in% 1:3))
  # centers are centered and of unit size
  for (h in 1:3) {
    expect_equal(sqrt(sum(a$centers[, , h]^2)), 1, tolerance = 1e-12)
    expect_equal(colMeans(a$centers[, , h]), rep(0, 3), tolerance = 1e-12)
  }
  expect_error(shape_kmeans(sim$X, 1), "K must")
  expect_error(shape_kmeans(sim$X, 16), "K must")
})

test_that("shape_kmeans follows the LloydShapes logic (R transcription, shapes)", {
  skip_if_not_installed("shapes")
  sim <- two_groups(n = 6, sigma = 3, seed = 4, type = "aniso")
  X <- sim$X
  n <- dim(X)[3L]
  ref_lloyd <- function(X, K, algSteps = 10, niter = 4, stopCr = 1e-4) {
    n <- dim(X)[3L]; vopt <- 1e8; dist <- matrix(0, n, K)
    for (iter in 1:niter) {
      obj <- list(); mean_sh <- list(); asig_step <- list()
      init <- sample(1:n, K, replace = FALSE)
      ms <- X[, , init]
      for (step in 1:algSteps) {
        for (h in 1:K) dist[, h] <- apply(X, 3, shapes::riemdist, y = ms[, , h])
        asig <- max.col(-dist)
        for (h in 1:K) {
          if (table(asig == h)[2] == 1) ms[, , h] <- X[, , asig == h]
          else ms[, , h] <- shapes::procGPA(X[, , asig == h], distances = TRUE,
                                            pcaoutput = TRUE)$mshape
          mean_sh[[step]] <- ms
        }
        obj[[step]] <- sum(dist[cbind(1:n, asig)]^2) / n
        asig_step[[step]] <- asig
        if (step > 1 && ((obj[[step - 1]] - obj[[step]]) / obj[[step - 1]]) < stopCr) break
      }
      if (min(unlist(obj)) < vopt) {
        vopt <- min(unlist(obj)); b <- which.min(unlist(obj))
        copt <- mean_sh[[b]]; asig_opt <- asig_step[[b]]
      }
    }
    list(asig = asig_opt, cases = copt, vopt = vopt)
  }
  set.seed(8); ref <- ref_lloyd(X, 2, niter = 4)
  set.seed(8); new <- shape_kmeans(X, 2, niter = 4)
  expect_identical(new$cluster, as.integer(ref$asig))
  expect_equal(new$objective, ref$vopt, tolerance = 1e-6)
  for (h in 1:2) expect_lt(shapes::riemdist(new$centers[, , h], ref$cases[, , h]), 1e-3)
})

test_that("shape_hill reproduces the transcription of the script, including the RNG stream", {
  set.seed(99)
  n <- 40
  D <- matrix(sample(1:3, n * n, TRUE), n)   # many ties on purpose
  D <- pmin(D, t(D)); diag(D) <- 0
  for (s in 1:5) {
    set.seed(s); r <- ref_hill_D(D, 3, 60); next_r <- runif(1)
    set.seed(s); m <- shape_hill(D, 3, 60);  next_m <- runif(1)
    expect_identical(m$medoids, as.integer(r$medoids))
    expect_identical(m$cluster, as.integer(r$cluster))
    expect_identical(m$objective, r$objective)
    expect_identical(next_m, next_r)
    set.seed(s); r <- ref_hill_D(D, 3, 60, from_best = TRUE)
    set.seed(s); m <- shape_hill(D, 3, 60, from_best = TRUE)
    expect_identical(m$medoids, as.integer(r$medoids))
    expect_identical(m$objective, r$objective)
  }
})

test_that("shape_hill and shape_clarans recover well-separated groups", {
  sim <- two_groups(n = 10, sigma = 0.5, seed = 1)
  D <- shape_distance(sim$X)
  set.seed(5); h <- shape_hill(D, 2)
  expect_equal(rand_index(sim$group, h$cluster), 1)
  set.seed(5); h2 <- shape_hill(sim$X, 2, max_iter = 50, from_best = TRUE)
  expect_equal(rand_index(sim$group, h2$cluster), 1)
  set.seed(5); cl <- shape_clarans(D, 2)
  expect_equal(rand_index(sim$group, cl$cluster), 1)
  # a landmark array gives the same as its distance matrix
  set.seed(6); a <- shape_clarans(sim$X, 2)
  set.seed(6); b <- shape_clarans(D, 2)
  expect_identical(a$medoids, b$medoids)
  expect_identical(a$distance, "rho")
  # consistency of the returned quantities
  for (fit in list(h, cl)) {
    expect_equal(length(unique(fit$medoids)), 2L)
    expect_equal(fit$objective, sum(D[cbind(1:20, fit$medoids[fit$cluster])]),
                 tolerance = 1e-12)
    expect_true(all(D[cbind(1:20, fit$medoids[fit$cluster])] ==
                      apply(D[, fit$medoids], 1, min)))
  }
  expect_error(shape_clarans(D, 20), "K must")
  expect_error(shape_hill(D, 2, initial = c(1, 1)), "distinct")
  expect_error(shape_clarans(D, 2, maxneighbor = -1), "positive")
})

test_that("shape_clarans draws its restarts from R's generator", {
  set.seed(31)
  D <- shape_distance(two_groups(n = 15, sigma = 4, seed = 6)$X)
  n <- nrow(D); K <- 3L
  # maxneighbor = 1 allows no swap trial: the result is the best of the
  # random initial medoid sets, which can be rebuilt with sample.int().
  draw_set <- function() {
    got <- integer()
    while (length(got) < K) { r <- sample.int(n, 1L); if (!r %in% got) got <- c(got, r) }
    got
  }
  for (s in 1:3) {
    set.seed(s); sets <- list(draw_set(), draw_set())
    cost <- vapply(sets, function(m) sum(apply(D[, m], 1, min)), 0)
    set.seed(s); fit <- shape_clarans(D, K, numlocal = 2, maxneighbor = 1)
    expect_identical(fit$medoids, as.integer(sets[[which.min(cost)]]))
    expect_equal(fit$objective, min(cost), tolerance = 1e-12)
    expect_identical(fit$n_iter, 0L)
  }
  # from the same initial medoids a search is never worse than no search
  for (s in 1:5) {
    set.seed(s); thorough <- shape_clarans(D, K, numlocal = 1, maxneighbor = 0.9)
    set.seed(s); none <- shape_clarans(D, K, numlocal = 1, maxneighbor = 1)
    expect_lte(thorough$objective, none$objective + 1e-12)
    expect_gte(thorough$n_iter, 0L)
  }
  # reproducible
  set.seed(2); a <- shape_clarans(D, K); set.seed(2); b <- shape_clarans(D, K)
  expect_identical(a, b)
})

# ---- bagclust1 -------------------------------------------------------------

test_that("bagclust1 on 2D complex input is unchanged from version 0.3.0", {
  fx <- readRDS(testthat::test_path("ref-bagclust1-2d.rds"))
  f1 <- bagclust1(fx$Z, K = 3, base = "pam", B = 20, seed = 11)
  expect_identical(f1$cluster, fx$pam$cluster)
  expect_identical(f1$votes, fx$pam$votes)
  expect_identical(f1$reference, fx$pam$reference)
  f2 <- bagclust1(fx$Z, K = 3, base = "rankkmed", B = 20, distance = "dF",
                  seed = 5, m = 7)
  expect_identical(f2$cluster, fx$rankkmed_dF$cluster)
  expect_identical(f2$votes, fx$rankkmed_dF$votes)
  expect_identical(f2$reference, fx$rankkmed_dF$reference)
})

test_that("bagclust1 accepts a real 2D array and agrees with the complex input", {
  set.seed(41)
  mu_list <- list(matrix(c(0,0, 10,0, 10,10, 0,10), 4, 2, byrow = TRUE),
                  matrix(c(0,0, 20,0, 15,10, 0,10), 4, 2, byrow = TRUE))
  g <- sim_landmark_gaussian(mu_list, c(15, 15), sigma = 0.6,
                             regime = "isotropic", seed = 2)
  a <- bagclust1(landmarks_to_preshape(g$X), K = 2, B = 12, seed = 1)
  b <- bagclust1(g$X, K = 2, B = 12, seed = 1)
  expect_identical(a$cluster, b$cluster)
  expect_equal(rand_index(g$group, b$cluster), 1)
})

test_that("bagclust1 works for every base clusterer on 3D arrays", {
  sim <- two_groups(n = 8, sigma = 0.5, seed = 1)
  for (base in c("pam", "rankkmed", "kmeans", "clarans", "hill")) {
    args <- switch(base, rankkmed = list(m = 4), hill = list(max_iter = 40), list())
    fit <- do.call(bagclust1, c(list(sim$X, K = 2, base = base, B = 8, seed = 1), args))
    expect_s3_class(fit, "bagclust1")
    expect_equal(fit$base, base)
    expect_true(all(rowSums(fit$votes) == 8L))
    expect_equal(rand_index(sim$group, fit$cluster), 1, info = base)
    expect_equal(dim(fit$votes), c(16L, 2L))
  }
  fit <- bagclust1(sim$X, K = 2, base = "kmeans", B = 6, seed = 2, mean = "phi")
  expect_equal(rand_index(sim$group, fit$cluster), 1)
  fit <- bagclust1(sim$X, K = 2, base = "pam", B = 6, seed = 2, distance = "dP",
                   reflection = TRUE)
  expect_equal(rand_index(sim$group, fit$cluster), 1)
  Z <- landmarks_to_preshape(array(rnorm(5 * 2 * 12), c(5, 2, 12)))
  expect_s3_class(bagclust1(Z, K = 2, B = 4, seed = 1, distance = "dP"), "bagclust1")
  expect_error(bagclust1(Z, K = 2, base = "kmeans", B = 2), "real landmark array")
  expect_error(bagclust1(sim$X, K = 2, base = "kmeans", B = 2, reflection = TRUE),
               "reflection")
  expect_error(bagclust1(Z, K = 2, B = 2, reflection = TRUE), "reflection")
})

test_that("bagclust1 on 3D arrays is reproducible with a seed and keeps the caller's stream", {
  sim <- two_groups(n = 8, sigma = 3, seed = 2)
  for (base in c("kmeans", "clarans", "hill")) {
    a <- bagclust1(sim$X, K = 2, base = base, B = 6, seed = 7)
    b <- bagclust1(sim$X, K = 2, base = base, B = 6, seed = 7)
    expect_identical(a$votes, b$votes)
    expect_identical(a$reference, b$reference)
    set.seed(5); x0 <- runif(1)
    set.seed(5); invisible(bagclust1(sim$X, K = 2, base = base, B = 3, seed = 1)); x1 <- runif(1)
    expect_identical(x0, x1)
  }
})

test_that("bagclust1 on 3D arrays gives the same result serially and in parallel", {
  sim <- two_groups(n = 8, sigma = 3, seed = 3)
  for (base in c("kmeans", "clarans", "hill")) {
    a <- bagclust1(sim$X, K = 2, base = base, B = 6, seed = 4)
    b <- bagclust1(sim$X, K = 2, base = base, B = 6, seed = 4,
                   parallel = TRUE, n_cores = 2L)
    expect_identical(a$votes, b$votes)
  }
})

test_that("shape_pam and shape_rankkmed accept a real landmark array", {
  sim <- two_groups(n = 10, sigma = 3, seed = 2)
  for (dist in c("rho", "dF")) {
    D <- shape_distance(sim$X, dist)
    set.seed(3); a <- shape_pam(sim$X, K = 2, distance = dist)
    set.seed(3); b <- shape_pam(D, K = 2, distance = dist)
    expect_identical(a, b)
    set.seed(3); a <- shape_rankkmed(sim$X, K = 2, distance = dist)
    set.seed(3); b <- shape_rankkmed(D, K = 2, distance = dist)
    expect_identical(a, b)
  }
})

test_that("an m = 2 real array gives the same partition as the complex pre-shapes", {
  set.seed(8)
  X <- array(rnorm(5 * 2 * 20), dim = c(5, 2, 20))
  Z <- landmarks_to_preshape(X)
  ini <- c(2L, 9L, 15L)
  expect_equal(shape_pam(X, K = 3, initial = ini)$cluster,
               shape_pam(Z, K = 3, initial = ini)$cluster)
  expect_equal(shape_rankkmed(X, K = 3, initial = ini)$cluster,
               shape_rankkmed(Z, K = 3, initial = ini)$cluster)
})

test_that("shape_pam runs on the 3D macaque skulls of 'shapes'", {
  skip_if_not_installed("shapes")
  data(macaques, package = "shapes", envir = environment())
  X <- macaques$x
  expect_identical(dim(X), c(7L, 3L, 18L))
  set.seed(1)
  fit <- shape_pam(X, K = 2)
  expect_length(fit$cluster, 18L)
  expect_setequal(fit$cluster, 1:2)
})
