suppressPackageStartupMessages({library(ShapeClust); library(ggplot2); library(shapes)})
src <- readLines(system.file("reproduce-all.R", package = "ShapeClust"))
cut <- grep("MAIN GRID", src)[1] - 2
eval(parse(text = paste(src[1:cut], collapse="\n")))
set.seed(99)
regimes <- list(list(name="bing", ks=c(4L,24L), tag="Concentration"),
                list(name="iso",  ks=c(4L,34L), tag="Scenario"),
                list(name="aniso",ks=c(4L,34L), tag="Scenario"))
bs <- c("A","B","E","F"); ls_ <- c("A","B","C","D")
for (rg in regimes) {
  sfx <- if (rg$name=="bing") bs else ls_
  for (k in rg$ks) {
    combos <- list(c(30,1), c(30,2), c(150,1), c(150,2))
    for (i in seq_along(combos)) {
      n <- combos[[i]][1]; sc <- combos[[i]][2]
      biplot_cell(rg$name, k, n, sc,
                  sprintf("Fig.%s_k%d_%s.pdf", rg$name, k, sfx[i]),
                  sprintf("%s %d, n = %d and k = %d", rg$tag, sc, n, k))
    }
  }
}
cat("biplots:", length(list.files(FIGDIR, "^Fig\\.")), "\n")
