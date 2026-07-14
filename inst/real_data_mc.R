suppressPackageStartupMessages({library(ShapeClust); library(shapes)})
B <- 100L; R_SEED <- 100L

analyse <- function(X, truth, name) {
  Z <- landmarks_to_preshape(X); D <- rho_matrix(Z)
  v <- riemannian_variances(X, truth)
  M <- t(vapply(seq_len(R_SEED), function(s) {
    set.seed(s)
    p_std <- shape_pam(D, K=2)$cluster
    p_bag <- bagclust1(Z, K=2, base="pam", B=B)$cluster
    r_std <- shape_rankkmed(D, K=2, m=5)$cluster
    r_bag <- bagclust1(Z, K=2, base="rankkmed", B=B, m=5)$cluster
    c(pam_RI_std=rand_index(truth,p_std),  pam_RI_bag=rand_index(truth,p_bag),
      pam_FMI_std=fmi_index(truth,p_std),  pam_FMI_bag=fmi_index(truth,p_bag),
      rkm_RI_std=rand_index(truth,r_std),  rkm_RI_bag=rand_index(truth,r_bag),
      rkm_FMI_std=fmi_index(truth,r_std),  rkm_FMI_bag=fmi_index(truth,r_bag))
  }, numeric(8)))
  cat(sprintf("\n===== %s =====\nVEG=%.4f  VIG=%.4f  ratio=%.2f\n", name, v$VEG, v$VIG, v$ratio))
  for (a in c("pam","rkm")) for (i in c("RI","FMI")) {
    s <- M[,paste0(a,"_",i,"_std")]; b <- M[,paste0(a,"_",i,"_bag")]
    cat(sprintf("%-4s %-3s  std %.3f (%.3f)  bag %.3f (%.3f)  RG=%+6.2f%%  p=%.2e\n",
        a, i, median(s), sd(s), median(b), sd(b),
        relative_gain(median(s), median(b)), wilcoxon_paired_gain(s,b)))
  }
  list(M=M, var=v)
}
data(sand); rs <- analyse(sand$x, as.integer(sand$group), "SAND GRAINS (k=50, n=49)")
data(pongof.dat); data(pongom.dat)
X <- array(0, dim=c(8,2,54)); X[,,1:24] <- pongof.dat; X[,,25:54] <- pongom.dat
ro <- analyse(X, c(rep(1L,24), rep(2L,30)), "ORANG-UTAN CRANIA (k=8, n=54)")
saveRDS(list(sand=rs, orang=ro), "paper/results_real_mc.rds")
