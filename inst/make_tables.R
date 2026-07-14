s <- readRDS("paper/results_summary.rds")
s$cell <- as.character(s$cell)
parse_cell <- function(x) {
  p <- strsplit(x, "_")[[1]]
  list(regime=p[1], k=as.integer(sub("k","",p[2])), n=as.integer(sub("n","",p[3])))
}
fmtp <- function(p) if (p < 0.001) "$\\mathbf{<\\!0.001}$" else
  if (p < 0.05) sprintf("$\\mathbf{%.3f}$", p) else sprintf("$%.3f$", p)
regnames <- c(bing="Bingham", iso="Isotropic", aniso="Anisotropic")

## ---- Table 1 (main paper): FMI, scenario 2 (hard) + scenario 1, all cells
lines <- c()
for (rg in c("bing","iso","aniso")) {
  cells <- unique(s$cell[startsWith(s$cell, paste0(rg,"_"))])
  cells <- cells[order(sapply(cells, function(c) parse_cell(c)$k),
                       sapply(cells, function(c) parse_cell(c)$n))]
  first <- TRUE
  for (cl in cells) {
    pc <- parse_cell(cl)
    for (sc in 1:2) {
      pa <- s[s$cell==cl & s$scenario==sc & s$algo=="pam" & s$index=="IFM",]
      rk <- s[s$cell==cl & s$scenario==sc & s$algo=="rkm" & s$index=="IFM",]
      lab <- if (first) sprintf("\\multirow{8}{*}{%s}", regnames[[rg]]) else ""
      first <- FALSE
      lines <- c(lines, sprintf(
        " %s & $k\\!=\\!%d,\\,n\\!=\\!%d$ & %d & $%.3f\\to %.3f$ & %s & $%.3f\\to %.3f$ & %s \\\\",
        lab, pc$k, pc$n, sc, pa$med_std, pa$med_bag, fmtp(pa$p),
        rk$med_std, rk$med_bag, fmtp(rk$p)))
    }
  }
  lines <- c(lines, "\\midrule")
}
lines <- lines[-length(lines)]
writeLines(lines, "paper/tab_main.tex")
cat("main table rows:", length(lines), "\n")

## ---- Supplement: full per-cell tables (RI and FMI, medians + sd + p)
sup <- c()
for (rg in c("bing","iso","aniso")) {
  for (algo in c("pam","rkm")) {
    ks <- sort(unique(sapply(unique(s$cell[startsWith(s$cell, paste0(rg,"_"))]),
                             function(c) parse_cell(c)$k)))
    for (k in ks) {
      algname <- if (algo=="pam") "Shape-PAM" else "Shape-rankkmed"
      conc <- if (rg=="bing") "Conc." else "Scen."
      sup <- c(sup,
        "\\begin{table}[h]", "\\centering",
        sprintf("\\caption{%s versus its bagged variant, %s regime, $k=%d$. Medians ($\\tilde x$), standard deviations ($s$) and one-sided paired Wilcoxon $p$-values over $100$ Monte Carlo replications.}",
                algname, regnames[[rg]], k),
        sprintf("\\label{tab:%s-%s-k%d}", algo, rg, k),
        "\\footnotesize",
        "\\begin{tabular}{cc c cc cc}", "\\toprule",
        "& & & \\multicolumn{2}{c}{\\textbf{RI}} & \\multicolumn{2}{c}{\\textbf{FMI}} \\\\",
        "\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}",
        sprintf("$n$ & %s & Approach & $\\tilde x$ ($s$) & $p$-val.\\ & $\\tilde x$ ($s$) & $p$-val.\\ \\\\", conc),
        "\\midrule")
      for (n in c(30,150)) for (sc in 1:2) {
        cl <- sprintf("%s_k%d_n%d", rg, k, n)
        ri <- s[s$cell==cl & s$scenario==sc & s$algo==algo & s$index=="IR",]
        fm <- s[s$cell==cl & s$scenario==sc & s$algo==algo & s$index=="IFM",]
        sup <- c(sup,
          sprintf("$%d$ & %d & Std.\\ & $%.3f\\,(%.3f)$ & \\multirow{2}{*}{%s} & $%.3f\\,(%.3f)$ & \\multirow{2}{*}{%s} \\\\",
                  n, sc, ri$med_std, ri$sd_std, fmtp(ri$p), fm$med_std, fm$sd_std, fmtp(fm$p)),
          sprintf("     &   & $B{=}100$ & $%.3f\\,(%.3f)$ & & $%.3f\\,(%.3f)$ & \\\\",
                  ri$med_bag, ri$sd_bag, fm$med_bag, fm$sd_bag))
      }
      sup <- c(sup, "\\bottomrule", "\\end{tabular}", "\\end{table}", "")
    }
  }
}
writeLines(sup, "paper/tab_supp.tex")
cat("supplement tables:", sum(grepl("begin.table", sup)), "\n")

## ---- Real data tables
r <- readRDS("paper/results_real_mc.rds")
mk <- function(x, name, k, n) {
  M <- x$M
  g <- function(a,i,w) median(M[,paste0(a,"_",i,"_",w)])
  sdv <- function(a,i,w) sd(M[,paste0(a,"_",i,"_",w)])
  pv <- function(a,i) suppressWarnings(wilcox.test(M[,paste0(a,"_",i,"_bag")],
        M[,paste0(a,"_",i,"_std")], alternative="greater", paired=TRUE, exact=FALSE)$p.value)
  rg <- function(a,i) 100*(g(a,i,"bag")-g(a,i,"std"))/g(a,i,"std")
  c(sprintf("\\begin{table}[t]"), "\\centering",
    sprintf("\\caption{Real data: %s ($k=%d$, $n=%d$, $K=2$). Median (standard deviation) of each index over $100$ random initializations; RG is the relative gain of the median; $p$ from the one-sided paired Wilcoxon test. Riemannian diagnostics: $\\rmVEG=%.4f$, $\\rmVIG=%.4f$, ratio $=%.2f$.}",
            name, k, n, x$var$VEG, x$var$VIG, x$var$ratio),
    sprintf("\\label{tab:%s}", tolower(gsub(" .*","",name))),
    "\\small", "\\begin{tabular}{lc cc cc c}", "\\toprule",
    "Algorithm & Approach & RI & RG-RI (\\%) & FMI & RG-FMI (\\%) & $p$ (FMI) \\\\", "\\midrule",
    sprintf("\\multirow{2}{*}{Shape-PAM} & Standard & $%.3f\\,(%.3f)$ & \\multirow{2}{*}{$%+.2f$} & $%.3f\\,(%.3f)$ & \\multirow{2}{*}{$%+.2f$} & \\multirow{2}{*}{%s} \\\\",
            g("pam","RI","std"), sdv("pam","RI","std"), rg("pam","RI"),
            g("pam","FMI","std"), sdv("pam","FMI","std"), rg("pam","FMI"), fmtp(pv("pam","FMI"))),
    sprintf(" & Bagging & $%.3f\\,(%.3f)$ & & $%.3f\\,(%.3f)$ & & \\\\",
            g("pam","RI","bag"), sdv("pam","RI","bag"), g("pam","FMI","bag"), sdv("pam","FMI","bag")),
    sprintf("\\multirow{2}{*}{Shape-rankkmed} & Standard & $%.3f\\,(%.3f)$ & \\multirow{2}{*}{$%+.2f$} & $%.3f\\,(%.3f)$ & \\multirow{2}{*}{$%+.2f$} & \\multirow{2}{*}{%s} \\\\",
            g("rkm","RI","std"), sdv("rkm","RI","std"), rg("rkm","RI"),
            g("rkm","FMI","std"), sdv("rkm","FMI","std"), rg("rkm","FMI"), fmtp(pv("rkm","FMI"))),
    sprintf(" & Bagging & $%.3f\\,(%.3f)$ & & $%.3f\\,(%.3f)$ & & \\\\",
            g("rkm","RI","bag"), sdv("rkm","RI","bag"), g("rkm","FMI","bag"), sdv("rkm","FMI","bag")),
    "\\bottomrule", "\\end{tabular}", "\\end{table}", "")
}
writeLines(c(mk(r$sand,"Sand grains",50,49), mk(r$orang,"Orang-utan crania",8,54)),
           "paper/tab_real.tex")
cat("real tables written\n")
