## =============================================================================
##  missing-figures.R — figures present in the original manuscript that were
##  not yet reproduced in v6:
##    (i)   landmark positions of a Complex Bingham pre-shape before and after
##          alignment to each of the three group reference directions (k = 4)
##    (ii)  the three predefined template objects (square, rectangle, trapezoid)
##    (iii) the sand-grain outlines, by provenance
##    (iv)  the eight cranial landmarks of the orang-utan
##  Plus the VEG/VIG variability table for the simulated pre-shape datasets.
## =============================================================================

suppressPackageStartupMessages({
  library(ShapeClust); library(ggplot2); library(shapes)
})
FIGDIR <- "paper/figures"
set.seed(20260714)

## helper: draw a k x 2 landmark configuration as a closed polygon
plot_config <- function(P, title, file, lab = TRUE, w = 3.2, h = 3.0) {
  df <- data.frame(x = P[, 1], y = P[, 2], id = seq_len(nrow(P)))
  cl <- rbind(df, df[1, ])
  p <- ggplot(df, aes(x, y)) +
    geom_path(data = cl, colour = "grey55", linewidth = 0.4) +
    geom_point(colour = "#D55E00", size = 2.6) +
    coord_equal() +
    labs(title = title, x = NULL, y = NULL) +
    theme_bw(base_size = 9) +
    theme(plot.title = element_text(size = 9, hjust = 0.5))
  if (lab) p <- p + geom_text(aes(label = id), vjust = -0.9, size = 2.6,
                              colour = "grey25")
  ggsave(file.path(FIGDIR, file), p, width = w, height = h, device = "pdf")
}

## -----------------------------------------------------------------------
## (i) Pre-shape alignment: original object and the three rotated groups
##     A complex pre-shape z in C^{k-1} is mapped back to landmark space by
##     inverting the Helmert projection, so the effect of the group rotation
##     can be seen directly on the landmarks.
## -----------------------------------------------------------------------
preshape_to_landmarks <- function(z, k) {
  H <- helmert_submatrix(k)
  w <- as.vector(t(H) %*% z)              # back to centred complex k-vector
  cbind(Re(w), Im(w))
}

d  <- 3L                                   # k = 4 landmarks -> d = k - 1
ed <- c(0, 0, 1)                           # the sampler's concentration pole
z0 <- sim_complex_bingham(1, lambda = rep(200, d - 1L))[, 1]

plot_config(preshape_to_landmarks(z0, 4L),
            expression(paste("Original object ", v[0])),
            "preshape_g0_k4_c1.pdf")

for (g in 1:3) {
  v <- rep(0, d); v[g] <- 1
  u <- v - ed
  U <- if (sqrt(sum(u^2)) < 1e-12) diag(d) else {
    u <- u / sqrt(sum(u^2)); diag(d) - 2 * outer(u, u)
  }
  zg <- as.vector(U %*% z0)
  plot_config(preshape_to_landmarks(zg, 4L),
              sprintf("Group %d reference direction", g),
              sprintf("preshape_g%d_k4_c1.pdf", g))
}

## -----------------------------------------------------------------------
## (ii) The three predefined mean configurations, k = 4
## -----------------------------------------------------------------------
sq <- matrix(c(0,0,  0,10, 10,10, 10,0), 4, 2, byrow = TRUE)
rc <- matrix(c(0,0, 10,0,   0,20, 10,20), 4, 2, byrow = TRUE)
tz <- matrix(c(0,0, 40,0,   0,30, 20,30), 4, 2, byrow = TRUE)
plot_config(sq, "Square",    "template_square.pdf")
plot_config(rc[c(1,2,4,3), ], "Rectangle", "template_rectangle.pdf")
plot_config(tz[c(1,2,4,3), ], "Trapezoid", "template_trapezoid.pdf")

## -----------------------------------------------------------------------
## (iii) Sand-grain outlines by provenance (50 landmarks each)
## -----------------------------------------------------------------------
data(sand)
X <- sand$x; grp <- as.integer(sand$group)      # 1 = river, 2 = sea
outline_panel <- function(idx, title, file) {
  dfs <- do.call(rbind, lapply(idx, function(i) {
    P <- X[, , i]
    P <- sweep(P, 2, colMeans(P))               # centre
    P <- P / sqrt(sum(P^2))                     # scale
    P <- rbind(P, P[1, ])
    data.frame(x = P[, 1], y = P[, 2], id = i)
  }))
  p <- ggplot(dfs, aes(x, y, group = id)) +
    geom_path(colour = "#0072B2", alpha = 0.45, linewidth = 0.3) +
    coord_equal() + labs(title = title, x = NULL, y = NULL) +
    theme_bw(base_size = 9) +
    theme(plot.title = element_text(size = 9, hjust = 0.5),
          axis.text = element_blank(), axis.ticks = element_blank())
  ggsave(file.path(FIGDIR, file), p, width = 3.4, height = 3.2, device = "pdf")
}
outline_panel(which(grp == 2), "Baltic Sea grains (n = 24)", "sand_sea.pdf")
outline_panel(which(grp == 1), "Selenchuk River grains (n = 25)", "sand_river.pdf")

## -----------------------------------------------------------------------
## (iv) Orang-utan cranial landmarks (mean configuration per sex)
## -----------------------------------------------------------------------
data(pongof.dat); data(pongom.dat)
lm_names <- c("nariale", "prosthion", "staphylion", "basion",
              "nasion", "bregma", "lambda", "opisthion")
mean_cfg <- function(A) {
  M <- apply(A, c(1, 2), mean)
  sweep(M, 2, colMeans(M))
}
Mf <- mean_cfg(pongof.dat); Mm <- mean_cfg(pongom.dat)
df <- rbind(
  data.frame(x = Mf[, 1], y = Mf[, 2], lab = lm_names, Sex = "Female"),
  data.frame(x = Mm[, 1], y = Mm[, 2], lab = lm_names, Sex = "Male"))
cl <- do.call(rbind, lapply(split(df, df$Sex), function(d) rbind(d, d[1, ])))
p <- ggplot(df, aes(x, y, colour = Sex)) +
  geom_path(data = cl, aes(group = Sex), linewidth = 0.35, alpha = 0.6) +
  geom_point(size = 2.4) +
  geom_text(data = subset(df, Sex == "Female"), aes(label = lab),
            vjust = -0.8, size = 2.3, colour = "grey25") +
  scale_colour_manual(values = c("#0072B2", "#D55E00")) +
  coord_equal() +
  labs(title = "Mean cranial configuration, 8 midline landmarks",
       x = NULL, y = NULL) +
  theme_bw(base_size = 9) +
  theme(legend.position = "top",
        plot.title = element_text(size = 9, hjust = 0.5))
ggsave(file.path(FIGDIR, "orang_landmarks.pdf"), p,
       width = 4.6, height = 3.6, device = "pdf")

## -----------------------------------------------------------------------
## (v) VEG / VIG variability table for the simulated pre-shape datasets
## -----------------------------------------------------------------------
src <- readLines(system.file("reproduce-all.R", package = "ShapeClust"))
cut <- grep("MAIN GRID", src)[1] - 2
eval(parse(text = paste(src[1:cut], collapse = "\n")))

rows <- list()
for (k in c(4L, 24L)) for (n in c(30L, 150L)) for (sc in 1:2) {
  lam <- bingham_lambda(k, sc)
  v <- replicate(20, {
    s <- sim_bingham_groups(n %/% 3L, k, lam)
    H <- helmert_submatrix(k)
    Xs <- array(0, dim = c(k, 2, ncol(s$Z)))
    for (i in seq_len(ncol(s$Z))) Xs[, , i] <- preshape_to_landmarks(s$Z[, i], k)
    rv <- riemannian_variances(Xs, s$group)
    c(rv$VEG, rv$VIG, rv$ratio)
  })
  rows[[length(rows) + 1]] <- data.frame(
    k = k, n = n, scenario = sc, lambda = lam,
    VEG = mean(v[1, ]), VIG = mean(v[2, ]), ratio = mean(v[3, ]))
}
tabv <- do.call(rbind, rows)
saveRDS(tabv, "paper/results_variability.rds")
print(tabv, row.names = FALSE, digits = 3)

cat("\nfigures now in figures/:", length(list.files(FIGDIR, "\\.pdf$")), "\n")
