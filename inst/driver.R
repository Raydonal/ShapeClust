suppressPackageStartupMessages({library(ShapeClust); library(ggplot2); library(shapes)})
src <- readLines(system.file("reproduce-all.R", package = "ShapeClust"))
cut <- grep("MAIN GRID", src)[1] - 2
eval(parse(text = paste(src[1:cut], collapse = "\n")))
CACHE <- "paper/cache"; dir.create(CACHE, showWarnings = FALSE)
BUDGET <- 235                                  # seconds of work per invocation
grid <- expand.grid(ch = 1:4, sc = 1:2, n = c(30L,150L),
                    idx = 1:6, stringsAsFactors = FALSE)
combos <- list(c("bing",4), c("bing",24), c("iso",4), c("iso",34),
               c("aniso",4), c("aniso",34))
t0 <- Sys.time(); ndone <- 0L
for (i in seq_len(nrow(grid))) {
  g <- grid[i,]; cb <- combos[[g$idx]]; rg <- cb[1]; k <- as.integer(cb[2])
  f  <- file.path(CACHE, sprintf("%s_k%d_n%d_s%d_c%d.rds", rg, k, g$n, g$sc, g$ch))
  f0 <- file.path(CACHE, sprintf("%s_k%d_n%d_s%d_c0.rds", rg, k, g$n, g$sc))
  if (file.exists(f) || file.exists(f0)) next
  if (as.numeric(difftime(Sys.time(), t0, units="secs")) > BUDGET) break
  R_MC <<- 25L
  set.seed(20260713L + 100000L*g$ch + 1000L*g$sc + 10L*k + g$n)
  d <- run_cell(rg, k, as.integer(g$n), as.integer(g$sc))
  d$rep <- d$rep + (g$ch - 1L)*25L
  saveRDS(d, f); ndone <- ndone + 1L
  cat(sprintf("  + %s k=%d n=%d s=%d c=%d\n", rg, k, g$n, g$sc, g$ch))
}
tot <- length(list.files(CACHE, "\\.rds$"))
cat(sprintf("chunks done this call: %d | cache files: %d\n", ndone, tot))
