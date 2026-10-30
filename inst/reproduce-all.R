## =============================================================================
##  reproduce-all.R
##  Runs every experiment in the paper and the Supplementary Material and
##  writes every figure to paper/figures/.
##
##  Design (as reported in the paper):
##    * R_MC   = 100 Monte Carlo replications per cell
##    * B      = 100 bootstrap replicates for BagClust1
##    * K      = 3 groups (simulations), K = 2 (real data)
##    * Regimes: Complex Bingham pre-shapes; isotropic landmarks;
##               anisotropic landmarks
##    * k in {4, 24} (Bingham) and {4, 34} (landmark regimes)
##    * n in {30, 150}
##    * Scenario 1 = well separated; Scenario 2 = overlapping
##
##  Outputs: 24 biplots + 48 violin plots + 2 real-data biplots.
## =============================================================================

suppressPackageStartupMessages({
  library(ShapeClust)
  library(ggplot2)
  library(shapes)
})

FIGDIR <- "paper/figures"
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

R_MC <- 100L
B    <- 100L
K    <- 3L

set.seed(20260713)

## ---------------------------------------------------------------- templates
## Mean configurations (square, rectangle, trapezoid); Table S3 of the
## Supplementary Material gives the k = 4 coordinates.
base_square    <- matrix(c(0,0,  0,10, 10,10, 10,0), 4, 2, byrow = TRUE)
base_rectangle <- matrix(c(0,0, 10,0,   0,20, 10,20), 4, 2, byrow = TRUE)
base_trapezoid <- matrix(c(0,0, 40,0,   0,30, 20,30), 4, 2, byrow = TRUE)

## Densify a polygon to k landmarks by interpolating along its perimeter.
densify <- function(P, k) {
  if (k == nrow(P)) return(P)
  cl <- rbind(P, P[1, ])
  seg <- sqrt(rowSums(diff(cl)^2))
  cum <- c(0, cumsum(seg))
  tot <- cum[length(cum)]
  s <- seq(0, tot, length.out = k + 1)[1:k]
  out <- t(vapply(s, function(si) {
    j <- max(which(cum <= si))
    j <- min(j, nrow(cl) - 1)
    w <- (si - cum[j]) / max(seg[j], 1e-12)
    cl[j, ] + w * (cl[j + 1, ] - cl[j, ])
  }, numeric(2)))
  out
}

templates <- function(k) {
  list(densify(base_square, k),
       densify(base_rectangle, k),
       densify(base_trapezoid, k))
}

## ------------------------------------------------- Complex Bingham generator
## Three groups obtained by rotating the simulated pre-shapes towards
## group-specific reference directions (Amaral et al. 2007).
## Concentration of the Complex Bingham by scenario.
## The concentration needed to hold the between-group overlap fixed grows with
## the dimension d = k - 1 of the pre-shape sphere: with lambda = 9 the mass is
## already almost uniform over 23 complex coordinates, so the three groups would
## be indistinguishable by construction. We therefore calibrate the overlapping
## scenario so that the two landmark counts pose a comparable difficulty
## (lambda = 9 for k = 4, lambda = 30 for k = 24, both giving a standard-PAM FMI
## of about 0.8-0.9); the well-separated scenario uses lambda = 200 in both.
bingham_lambda <- function(k, scenario) {
  if (scenario == 1L) return(200)
  if (k <= 4L) 9 else 30
}

sim_bingham_groups <- function(n_per, k, lambda) {
  d <- k - 1L                      # complex dimension of the pre-shape sphere
  ed <- rep(0, d); ed[d] <- 1      # the sampler puts the pole in coordinate d
  Zs <- lapply(1:3, function(g) {
    Z0 <- sim_complex_bingham(n_per, lambda = rep(lambda, d - 1L))
    ## Group-specific mean shape: rotate the pole e_d onto a direction v_g by a
    ## real Householder reflection (an isometry of the pre-shape sphere, so it
    ## maps the Complex Bingham onto a Complex Bingham with the same
    ## concentration but a different modal shape).
    v <- rep(0, d); v[((g - 1L) %% d) + 1L] <- 1
    u <- v - ed
    U <- if (sqrt(sum(u^2)) < 1e-12) diag(d) else {
      u <- u / sqrt(sum(u^2)); diag(d) - 2 * outer(u, u)
    }
    U %*% Z0
  })
  Z <- do.call(cbind, Zs)
  norms <- sqrt(colSums(Mod(Z)^2)); norms[norms == 0] <- 1
  list(Z = sweep(Z, 2L, norms, "/"), group = rep(1:3, each = n_per))
}

## ------------------------------------------------------------- one MC cell
## Returns a data frame with the four validation numbers per replication.
run_cell <- function(regime, k, n, scenario, distance = "rho") {
  n_per <- n %/% 3L
  reps <- vector("list", R_MC)

  for (r in seq_len(R_MC)) {
    if (regime == "bing") {
      lambda <- bingham_lambda(k, scenario)
      sim <- sim_bingham_groups(n_per, k, lambda)
      Z <- sim$Z; truth <- sim$group
    } else {
      sigma  <- if (scenario == 1L) 3 else 9
      rg     <- if (regime == "iso") "isotropic" else "anisotropic"
      sim <- sim_landmark_gaussian(templates(k), rep(n_per, 3),
                                   sigma = sigma, regime = rg)
      Z <- landmarks_to_preshape(sim$X); truth <- sim$group
    }

    D <- if (distance == "rho") rho_matrix(Z) else dF_matrix(Z)

    p_std <- shape_pam(D, K = K)$cluster
    p_bag <- bagclust1(Z, K = K, base = "pam", B = B, distance = distance)$cluster
    r_std <- shape_rankkmed(D, K = K, m = 5)$cluster
    r_bag <- bagclust1(Z, K = K, base = "rankkmed", B = B, m = 5,
                       distance = distance)$cluster

    reps[[r]] <- data.frame(
      rep = r,
      pam_IR_std   = rand_index(truth, p_std),
      pam_IR_bag   = rand_index(truth, p_bag),
      pam_IFM_std  = fmi_index(truth, p_std),
      pam_IFM_bag  = fmi_index(truth, p_bag),
      rkm_IR_std   = rand_index(truth, r_std),
      rkm_IR_bag   = rand_index(truth, r_bag),
      rkm_IFM_std  = fmi_index(truth, r_std),
      rkm_IFM_bag  = fmi_index(truth, r_bag)
    )
  }
  out <- do.call(rbind, reps)
  out$regime <- regime; out$k <- k; out$n <- n; out$scenario <- scenario
  out
}

## ------------------------------------------------------------------ biplots
## PCA of the pre-shapes (real and imaginary parts stacked), coloured by group.
biplot_cell <- function(regime, k, n, scenario, file, title) {
  n_per <- n %/% 3L
  if (regime == "bing") {
    lambda <- bingham_lambda(k, scenario)
    sim <- sim_bingham_groups(n_per, k, lambda)
    Z <- sim$Z; truth <- sim$group
    lab <- c("Group 1", "Group 2", "Group 3")
  } else {
    sigma <- if (scenario == 1L) 3 else 9
    rg    <- if (regime == "iso") "isotropic" else "anisotropic"
    sim <- sim_landmark_gaussian(templates(k), rep(n_per, 3),
                                 sigma = sigma, regime = rg)
    Z <- landmarks_to_preshape(sim$X); truth <- sim$group
    lab <- c("Squares", "Rectangles", "Trapezoids")
  }
  M <- rbind(Re(Z), Im(Z))              # real coordinates of the pre-shape
  pca <- prcomp(t(M), center = TRUE, scale. = FALSE)
  df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                   Group = factor(lab[truth], levels = lab))

  p <- ggplot(df, aes(PC1, PC2, colour = Group)) +
    geom_point(size = 1.9, alpha = 0.85) +
    scale_colour_manual(values = c("#D55E00", "#009E73", "#0072B2")) +
    labs(title = title, x = "PC1", y = "PC2") +
    theme_bw(base_size = 9) +
    theme(legend.position = "top",
          legend.title = element_text(size = 8),
          plot.title = element_text(size = 9, hjust = 0.5))
  ggsave(file.path(FIGDIR, file), p, width = 4.2, height = 3.4, device = "pdf")
}

## ------------------------------------------------------------- violin plots
## Paired violins: x = scenario, fill = approach (standard vs bagged).
violin_plot <- function(dat, algo, index, k, n, file) {
  std <- dat[[paste0(algo, "_", index, "_std")]]
  bag <- dat[[paste0(algo, "_", index, "_bag")]]
  df <- data.frame(
    value    = c(std, bag),
    approach = factor(rep(c("Standard", "B = 100"), each = nrow(dat)),
                      levels = c("B = 100", "Standard")),
    scenario = factor(rep(dat$scenario, 2),
                      levels = 1:2,
                      labels = c("Scenario 1", "Scenario 2"))
  )
  ylab <- if (index == "IR") "RI" else "FMI"
  ttl  <- sprintf("%s and n = %d", ylab, n)

  p <- ggplot(df, aes(scenario, value, fill = approach)) +
    geom_violin(position = position_dodge(width = 0.8),
                alpha = 0.45, colour = "grey30", linewidth = 0.3,
                trim = FALSE) +
    geom_boxplot(position = position_dodge(width = 0.8),
                 width = 0.18, outlier.size = 0.5, linewidth = 0.3,
                 alpha = 0.9) +
    scale_fill_manual(values = c("B = 100" = "#1B9E77", "Standard" = "#D95F02"),
                      name = "Approach") +
    labs(title = ttl, x = "Scenarios", y = ylab) +
    theme_bw(base_size = 9) +
    theme(legend.position = "top",
          plot.title = element_text(size = 9, hjust = 0.5))
  ggsave(file.path(FIGDIR, file), p, width = 4.2, height = 3.4, device = "pdf")
}

## =============================================================================
##  MAIN GRID
## =============================================================================
regimes <- list(
  list(name = "bing",  ks = c(4L, 24L), tag = "Concentration"),
  list(name = "iso",   ks = c(4L, 34L), tag = "Scenario"),
  list(name = "aniso", ks = c(4L, 34L), tag = "Scenario")
)

CACHE <- "paper/cache"
dir.create(CACHE, showWarnings = FALSE)

cell_cached <- function(regime, k, n, sc) {
  f <- file.path(CACHE, sprintf("%s_k%d_n%d_s%d.rds", regime, k, n, sc))
  if (file.exists(f)) return(readRDS(f))
  out <- run_cell(regime, k, n, sc)
  saveRDS(out, f)
  out
}

all_results <- list()
t0 <- Sys.time()

for (rg in regimes) {
  for (k in rg$ks) {
    for (n in c(30L, 150L)) {
      cells <- lapply(1:2, function(sc) cell_cached(rg$name, k, n, sc))
      dat <- do.call(rbind, cells)
      key <- sprintf("%s_k%d_n%d", rg$name, k, n)
      all_results[[key]] <- dat

      # two violin figures per (algorithm, index)
      for (algo in c("pam", "rkm")) {
        prefix <- if (algo == "pam") "pam" else "rankkmed"
        for (index in c("IR", "IFM")) {
          violin_plot(dat, algo, index, k, n,
                      sprintf("%s_violin_%s_%s_k%d_n%d.pdf",
                              prefix, index, rg$name, k, n))
        }
      }
      cat(sprintf("[%s] done  (%.1f min elapsed)\n", key,
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
      flush.console()
    }
  }
}

saveRDS(all_results, "paper/results_all_cells.rds")

## ------------------------------------------------------------------ biplots
## Naming follows the paper: Bingham uses suffixes A,B,E,F; the landmark
## regimes use A,B,C,D. Both mean (scenario 1 / scenario 2) x (n=30 / n=150).
bing_suffix <- c("A", "B", "E", "F")
land_suffix <- c("A", "B", "C", "D")

for (rg in regimes) {
  sfx <- if (rg$name == "bing") bing_suffix else land_suffix
  for (k in rg$ks) {
    combos <- list(c(30, 1), c(30, 2), c(150, 1), c(150, 2))
    for (i in seq_along(combos)) {
      n <- combos[[i]][1]; sc <- combos[[i]][2]
      ttl <- sprintf("%s %d, n = %d and k = %d", rg$tag, sc, n, k)
      biplot_cell(rg$name, k, n, sc,
                  sprintf("Fig.%s_k%d_%s.pdf", rg$name, k, sfx[i]), ttl)
    }
  }
}
cat("biplots done\n")

## =============================================================================
##  REAL DATA
## =============================================================================
real_biplot <- function(X, truth, labels, title, file) {
  Z <- landmarks_to_preshape(X)
  M <- rbind(Re(Z), Im(Z))
  pca <- prcomp(t(M), center = TRUE, scale. = FALSE)
  df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                   Group = factor(labels[truth], levels = labels))
  p <- ggplot(df, aes(PC1, PC2, colour = Group)) +
    geom_point(size = 2.2, alpha = 0.9) +
    scale_colour_manual(values = c("#0072B2", "#D55E00")) +
    labs(title = title, x = "PC1", y = "PC2") +
    theme_bw(base_size = 10) +
    theme(legend.position = "top",
          plot.title = element_text(size = 10, hjust = 0.5))
  ggsave(file.path(FIGDIR, file), p, width = 4.6, height = 3.8, device = "pdf")
}

real_table <- function(X, truth, name) {
  Z <- landmarks_to_preshape(X)
  D <- rho_matrix(Z)
  p_std <- shape_pam(D, K = 2)$cluster
  p_bag <- bagclust1(Z, K = 2, base = "pam", B = B, seed = 1)$cluster
  r_std <- shape_rankkmed(D, K = 2, m = 5)$cluster
  r_bag <- bagclust1(Z, K = 2, base = "rankkmed", B = B, seed = 1, m = 5)$cluster
  v <- riemannian_variances(X, truth)

  f <- function(cl) c(RI = rand_index(truth, cl), FMI = fmi_index(truth, cl))
  tab <- rbind(
    data.frame(algorithm = "Shape-PAM",      approach = "Standard", t(f(p_std))),
    data.frame(algorithm = "Shape-PAM",      approach = "Bagging",  t(f(p_bag))),
    data.frame(algorithm = "Shape-rankkmed", approach = "Standard", t(f(r_std))),
    data.frame(algorithm = "Shape-rankkmed", approach = "Bagging",  t(f(r_bag)))
  )
  cat(sprintf("\n=== %s  (VEG=%.3f VIG=%.3f ratio=%.2f) ===\n",
              name, v$VEG, v$VIG, v$ratio))
  print(tab, row.names = FALSE, digits = 3)
  cat(sprintf("  PAM       RG-RI = %+6.2f%%   RG-FMI = %+6.2f%%\n",
              relative_gain(tab$RI[1], tab$RI[2]),
              relative_gain(tab$FMI[1], tab$FMI[2])))
  cat(sprintf("  rankkmed  RG-RI = %+6.2f%%   RG-FMI = %+6.2f%%\n",
              relative_gain(tab$RI[3], tab$RI[4]),
              relative_gain(tab$FMI[3], tab$FMI[4])))
  attr(tab, "var") <- v
  tab
}

## Sand grains: 24 sea + 25 river, k = 50
data(sand)
X_sand <- sand$x
truth_sand <- as.integer(sand$group)          # 1 = river, 2 = sea
real_biplot(X_sand, truth_sand, c("River grains", "Sea grains"),
            "Sand grains", "sand_pca.pdf")
tab_sand <- real_table(X_sand, truth_sand, "SAND GRAINS (k=50, n=49, K=2)")

## Orang-utan crania: O'Higgins (1989), 24 females + 30 males, k = 8
data(pongof.dat); data(pongom.dat)
X_orang <- array(0, dim = c(8, 2, 54))
X_orang[, , 1:24]  <- pongof.dat
X_orang[, , 25:54] <- pongom.dat
truth_orang <- c(rep(1L, 24), rep(2L, 30))    # 1 = female, 2 = male
real_biplot(X_orang, truth_orang, c("Female", "Male"),
            "Orang-utan crania", "orang_pca.pdf")
tab_orang <- real_table(X_orang, truth_orang, "ORANG-UTAN CRANIA (k=8, n=54, K=2)")

saveRDS(list(sand = tab_sand, orang = tab_orang),
        "paper/results_real.rds")

cat(sprintf("\nTOTAL TIME: %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat("figures written:", length(list.files(FIGDIR, pattern = "\\.pdf$")), "\n")
