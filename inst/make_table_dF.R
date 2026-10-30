## Table S.8: conclusions about bagging under rho versus d_F.
## Needs <paper>/cache (inst/driver.R) and <paper>/cache_dF (inst/driver_dF.R).
## Usage: Rscript inst/make_table_dF.R [paper_dir]
suppressPackageStartupMessages(library(ShapeClust))
args <- commandArgs(trailingOnly = TRUE)
PAPER <- if (length(args) >= 1) args[1] else "paper"
load_cell <- function(dir, rg, k) do.call(rbind, lapply(
  file.path(PAPER, dir, sprintf("%s_k%d_n150_s2_c%d.rds", rg, k, 1:4)), readRDS))
cells <- list(list("bing", 24L, "Bingham, $k=24$"), list("iso", 34L, "Isotropic, $k=34$"),
              list("aniso", 34L, "Anisotropic, $k=34$"))
fp <- function(p) if (p < 0.001) "$<\\!0.001$" else sprintf("$%.3f$", p)
L <- c(); rows <- list()
for (c in cells) {
  r <- load_cell("cache", c[[1]], c[[2]]); f <- load_cell("cache_dF", c[[1]], c[[2]])
  stopifnot(identical(r$rep, f$rep))
  block <- c()
  for (a in c("pam", "rkm")) {
    s1 <- r[[paste0(a, "_IFM_std")]]; b1 <- r[[paste0(a, "_IFM_bag")]]
    s2 <- f[[paste0(a, "_IFM_std")]]; b2 <- f[[paste0(a, "_IFM_bag")]]
    p1 <- wilcoxon_paired_gain(s1, b1); p2 <- wilcoxon_paired_gain(s2, b2)
    same <- sign(median(b1) - median(s1)) == sign(median(b2) - median(s2)) &&
            (p1 < 0.05) == (p2 < 0.05)
    rows[[length(rows) + 1]] <- data.frame(cell = c[[1]], algo = a, rho_std = median(s1),
      rho_bag = median(b1), p_rho = p1, dF_std = median(s2), dF_bag = median(b2), p_dF = p2, same = same)
    block <- c(block, sprintf(" & %s & $%.3f\\to%.3f$ & %s & $%.3f\\to%.3f$ & %s & %s \\\\",
      if (a == "pam") "Shape-PAM     " else "Shape-rankkmed", median(s1), median(b1), fp(p1),
      median(s2), median(b2), fp(p2), if (same) "yes" else "\\textbf{no}"))
  }
  L <- c(L, sprintf("\\multirow{2}{*}{%s}", c[[3]]), block, "\\midrule")
}
L <- L[-length(L)]
writeLines(L, file.path(PAPER, "tab_dF_rows.tex"))
res <- do.call(rbind, rows); saveRDS(res, file.path(PAPER, "results_dF.rds"))
print(res, digits = 3, row.names = FALSE)
