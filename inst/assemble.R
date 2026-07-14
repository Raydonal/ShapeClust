suppressPackageStartupMessages({library(ShapeClust); library(ggplot2); library(shapes)})
src <- readLines(system.file("reproduce-all.R", package = "ShapeClust"))
cut <- grep("MAIN GRID", src)[1] - 2
eval(parse(text = paste(src[1:cut], collapse="\n")))
CACHE <- "paper/cache"

load_cell <- function(rg, k, n, sc) {
  f0 <- file.path(CACHE, sprintf("%s_k%d_n%d_s%d_c0.rds", rg, k, n, sc))
  if (file.exists(f0)) return(readRDS(f0))
  fs <- file.path(CACHE, sprintf("%s_k%d_n%d_s%d_c%d.rds", rg, k, n, sc, 1:4))
  do.call(rbind, lapply(fs, readRDS))
}

regimes <- list(list(name="bing", ks=c(4L,24L)),
                list(name="iso",  ks=c(4L,34L)),
                list(name="aniso",ks=c(4L,34L)))
all_results <- list()
for (rg in regimes) for (k in rg$ks) for (n in c(30L,150L)) {
  dat <- rbind(load_cell(rg$name,k,n,1L), load_cell(rg$name,k,n,2L))
  stopifnot(nrow(dat) == 200)                 # 100 reps x 2 scenarios
  all_results[[sprintf("%s_k%d_n%d", rg$name,k,n)]] <- dat
  for (algo in c("pam","rkm")) {
    prefix <- if (algo=="pam") "pam" else "rankkmed"
    for (index in c("IR","IFM"))
      violin_plot(dat, algo, index, k, n,
                  sprintf("%s_violin_%s_%s_k%d_n%d.pdf", prefix, index, rg$name, k, n))
  }
}
saveRDS(all_results, "paper/results_all_cells.rds")
cat("violins:", length(list.files(FIGDIR, "violin.*pdf")), "\n")

# Summary table of medians + Wilcoxon p-values (all 24 cells)
rows <- list()
for (key in names(all_results)) {
  d <- all_results[[key]]
  for (sc in 1:2) {
    dd <- d[d$scenario==sc,]
    for (algo in c("pam","rkm")) for (idx in c("IR","IFM")) {
      s <- dd[[paste0(algo,"_",idx,"_std")]]; b <- dd[[paste0(algo,"_",idx,"_bag")]]
      rows[[length(rows)+1]] <- data.frame(cell=key, scenario=sc, algo=algo, index=idx,
        med_std=median(s), sd_std=sd(s), med_bag=median(b), sd_bag=sd(b),
        p=wilcoxon_paired_gain(s,b))
    }
  }
}
summ <- do.call(rbind, rows)
saveRDS(summ, "paper/results_summary.rds")
write.csv(summ, "paper/results_summary.csv", row.names=FALSE)
cat("summary rows:", nrow(summ), "\n")
print(head(summ[summ$index=="IFM" & summ$scenario==2, ], 12), row.names=FALSE, digits=3)
