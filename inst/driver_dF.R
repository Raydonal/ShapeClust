## Sensitivity to the shape-space distance (Supplementary Section S.8).
## Re-runs the three overlapping n = 150 scenarios of largest k with the full
## Procrustes distance d_F, using the same seeds as inst/driver.R, so each
## replication sees the same data under rho and under d_F.
## Run from the repository root: Rscript inst/driver_dF.R [paper_dir] [cores]
suppressPackageStartupMessages({library(ShapeClust); library(ggplot2); library(shapes); library(parallel)})
args <- commandArgs(trailingOnly = TRUE)
PAPER <- if (length(args) >= 1) args[1] else "paper"
CORES <- if (length(args) >= 2) as.integer(args[2]) else 1L
src <- readLines(system.file("reproduce-all.R", package = "ShapeClust"))
cut <- grep("MAIN GRID", src)[1] - 2
eval(parse(text = paste(src[1:cut], collapse = "\n")))
CACHE <- file.path(PAPER, "cache_dF"); dir.create(CACHE, showWarnings = FALSE)
cells <- list(c("bing", 24), c("iso", 34), c("aniso", 34))
jobs <- expand.grid(cell = seq_along(cells), ch = 1:4)
R_MC <- 25L; n <- 150L; sc <- 2L
invisible(mclapply(seq_len(nrow(jobs)), function(i) {
  cb <- cells[[jobs$cell[i]]]; rg <- cb[1]; k <- as.integer(cb[2]); ch <- jobs$ch[i]
  f <- file.path(CACHE, sprintf("%s_k%d_n%d_s%d_c%d.rds", rg, k, n, sc, ch))
  if (file.exists(f)) return(NULL)
  set.seed(20260713L + 100000L*ch + 1000L*sc + 10L*k + n)
  d <- run_cell(rg, k, n, sc, distance = "dF")
  d$rep <- d$rep + (ch - 1L)*25L
  saveRDS(d, f); cat(sprintf("+ dF %s k=%d c=%d\n", rg, k, ch)); NULL
}, mc.cores = CORES, mc.preschedule = FALSE))
cat("dF chunks:", length(list.files(CACHE)), "of 12\n")
