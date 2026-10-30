## Wall-clock timings for Supplementary Sections S.1, S.6 and S.7.
## Run from the repository root on an otherwise idle machine:
##   Rscript inst/timing.R [paper_dir] [cores]
## Writes <paper_dir>/results_timing.rds and <paper_dir>/tab_scal.tex,
## tab_time.tex, tab_time_real.tex. Times are medians over REPS runs.
suppressPackageStartupMessages({library(ShapeClust); library(ggplot2); library(shapes)})
args <- commandArgs(trailingOnly = TRUE)
PAPER <- if (length(args) >= 1) args[1] else "paper"
CORES <- if (length(args) >= 2) as.integer(args[2]) else 3L
src <- readLines(system.file("reproduce-all.R", package = "ShapeClust"))
cut <- grep("MAIN GRID", src)[1] - 2
eval(parse(text = paste(src[1:cut], collapse = "\n")))
REPS <- 3L
tm <- function(expr) {
  f <- substitute(expr); e <- parent.frame()
  median(vapply(seq_len(REPS), function(r) system.time(eval(f, e))[["elapsed"]], 0))
}
fmt <- function(x, d = 3) sprintf(paste0("$%.", d, "f$"), x)
set.seed(20261024)

## ---- S.1: scalability on an (n, k) grid and in B ----------------------------
scal <- do.call(rbind, lapply(c(4L, 24L, 50L), function(k) do.call(rbind, lapply(
  c(50L, 100L, 200L, 500L, 1000L), function(n) {
    Z <- sim_complex_bingham(n, lambda = rep(9, k - 2L))
    D <- rho_matrix(Z)
    data.frame(n = n, k = k, dist = tm(rho_matrix(Z)),
               pam = tm(shape_pam(D, K = 3)), rkm = tm(shape_rankkmed(D, K = 3, m = 5)))
  }))))
Z200 <- sim_complex_bingham(200, lambda = rep(9, 22))
scalB <- do.call(rbind, lapply(c(10L, 25L, 50L, 100L), function(B)
  data.frame(B = B, time = tm(bagclust1(Z200, K = 3, base = "pam", B = B)))))
scalB$per <- scalB$time / scalB$B

## ---- S.6: BagClust1 (B = 100) serially and on CORES cores, every scenario ---
draw <- function(regime, k, n, sc) {
  if (regime == "bing") return(sim_bingham_groups(n %/% 3L, k, bingham_lambda(k, sc))$Z)
  rg <- if (regime == "iso") "isotropic" else "anisotropic"
  landmarks_to_preshape(sim_landmark_gaussian(templates(k), rep(n %/% 3L, 3),
                        sigma = if (sc == 1L) 3 else 9, regime = rg)$X)
}
combos <- list(bing = c(4L, 24L), iso = c(4L, 34L), aniso = c(4L, 34L))
grid <- do.call(rbind, lapply(names(combos), function(rg) expand.grid(
  regime = rg, k = combos[[rg]], n = c(30L, 150L), sc = 1:2,
  base = c("pam", "rankkmed"), stringsAsFactors = FALSE)))
grid$serial <- grid$par <- NA_real_
for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]; Z <- draw(g$regime, g$k, g$n, g$sc)
  grid$serial[i] <- tm(bagclust1(Z, K = 3, base = g$base, B = 100))
  grid$par[i]    <- tm(bagclust1(Z, K = 3, base = g$base, B = 100,
                                 parallel = TRUE, n_cores = CORES))
}

## ---- S.7: real data -----------------------------------------------------------
data(sand); data(pongof.dat); data(pongom.dat)
Xo <- array(0, dim = c(8, 2, 54)); Xo[, , 1:24] <- pongof.dat; Xo[, , 25:54] <- pongom.dat
real <- list(sand = landmarks_to_preshape(sand$x), orang = landmarks_to_preshape(Xo))
treal <- do.call(rbind, lapply(names(real), function(d) do.call(rbind, lapply(
  c("pam", "rankkmed"), function(b) data.frame(
    data = d, base = b,
    serial = tm(bagclust1(real[[d]], K = 2, base = b, B = 100)),
    par = tm(bagclust1(real[[d]], K = 2, base = b, B = 100, parallel = TRUE, n_cores = CORES)))))))

info <- list(R = R.version.string, cpu = tryCatch(sub(".*: ", "", grep("model name",
             readLines("/proc/cpuinfo"), value = TRUE)[1]), error = function(e) NA),
             cores = CORES, version = as.character(packageVersion("ShapeClust")), reps = REPS)
saveRDS(list(scal = scal, scalB = scalB, grid = grid, real = treal, info = info),
        file.path(PAPER, "results_timing.rds"))

## ---- LaTeX ---------------------------------------------------------------------
L <- c("\\begin{table}[h]", "\\centering",
  sprintf("\\caption{Median wall-clock time (seconds) over %d runs for the three main routines of \\texttt{ShapeClust} on a factorial grid of sample size $n$ and landmark count $k$ ($K=3$, $m=5$). Values below $1$ ms are reported as $0.000$. Referenced from Section~4 of the main paper.}", REPS),
  "\\label{tab:scal-nk}", "\\small", "\\begin{tabular}{cc rrr}", "\\toprule",
  "$n$ & $k$ & Distance matrix & Shape-PAM & Shape-rankkmed \\\\")
for (k in unique(scal$k)) {
  L <- c(L, "\\midrule")
  for (r in which(scal$k == k)) L <- c(L, sprintf("%4d & %2d & %s & %s & %s \\\\",
    scal$n[r], k, fmt(scal$dist[r]), fmt(scal$pam[r]), fmt(scal$rkm[r])))
}
L <- c(L, "\\bottomrule", "\\end{tabular}", "\\end{table}", "",
  "\\begin{table}[h]", "\\centering",
  "\\caption{Wall-clock time of BagClust1 (Shape-PAM base, $n=200$, $k=24$, $K=3$) as a function of the number of bootstrap replicates $B$. Referenced from Section~4 of the main paper.}",
  "\\label{tab:scal-B}", "\\small", "\\begin{tabular}{rrr}", "\\toprule",
  "$B$ & Time (s) & Time per replicate (s) \\\\", "\\midrule",
  sprintf("%3d & %s & %s \\\\", scalB$B, fmt(scalB$time, 2), fmt(scalB$per, 3)),
  "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(L, file.path(PAPER, "tab_scal.tex"))

regname <- c(bing = "Bingham pre-shape", iso = "isotropic landmark", aniso = "anisotropic landmark")
T <- c()
for (rg in names(combos)) for (b in c("pam", "rankkmed")) {
  ks <- combos[[rg]]; lab <- if (rg == "bing") "Conc.\\ " else "Scen.\\ "
  T <- c(T, "\\begin{table}[h]", "\\centering",
    sprintf("\\caption{Bagging-%s execution time (seconds, $B=100$), %s regime: serial and on %d cores.}",
            if (b == "pam") "PAM" else "rankkmed", regname[[rg]], CORES),
    sprintf("\\label{tab:time-%s-%s}", if (b == "pam") "pam" else "rkm", rg),
    "\\small", "\\begin{tabular}{cc cc cc}", "\\toprule",
    sprintf("& & \\multicolumn{2}{c}{$k=%d$} & \\multicolumn{2}{c}{$k=%d$} \\\\", ks[1], ks[2]),
    "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}",
    sprintf("$n$ & %s & Serial & Parallel & Serial & Parallel \\\\", lab), "\\midrule")
  for (n in c(30L, 150L)) for (sc in 1:2) {
    r1 <- grid[grid$regime == rg & grid$base == b & grid$n == n & grid$sc == sc & grid$k == ks[1], ]
    r2 <- grid[grid$regime == rg & grid$base == b & grid$n == n & grid$sc == sc & grid$k == ks[2], ]
    T <- c(T, sprintf("%s & %d & %s & %s & %s & %s \\\\", if (sc == 1L) sprintf("$%d$", n) else "",
                      sc, fmt(r1$serial, 2), fmt(r1$par, 2), fmt(r2$serial, 2), fmt(r2$par, 2)))
  }
  T <- c(T, "\\bottomrule", "\\end{tabular}", "\\end{table}", "")
}
writeLines(T, file.path(PAPER, "tab_time.tex"))

R <- c("\\begin{table}[h]", "\\centering",
  sprintf("\\caption{Execution time (seconds, $B=100$) of the bagging procedure on the two real datasets, serial and on %d cores.}", CORES),
  "\\label{tab:time-real}", "\\small", "\\begin{tabular}{lcc cc}", "\\toprule",
  "& \\multicolumn{2}{c}{Sand grains} & \\multicolumn{2}{c}{Orang-utan crania} \\\\",
  "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}",
  "Algorithm & Serial & Parallel & Serial & Parallel \\\\", "\\midrule")
for (b in c("pam", "rankkmed")) {
  s <- treal[treal$data == "sand" & treal$base == b, ]; o <- treal[treal$data == "orang" & treal$base == b, ]
  R <- c(R, sprintf("%s & %s & %s & %s & %s \\\\", if (b == "pam") "Shape-PAM" else "Shape-rankkmed",
                    fmt(s$serial, 2), fmt(s$par, 2), fmt(o$serial, 2), fmt(o$par, 2)))
}
R <- c(R, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(R, file.path(PAPER, "tab_time_real.tex"))
sp <- 100 * (1 - grid$par / grid$serial)
cat(sprintf("timing done. parallel time saving: %.0f%% to %.0f%% (median %.0f%%)\n",
            min(sp), max(sp), median(sp)))
