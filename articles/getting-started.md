# Getting started with ShapeClust

## Overview

**ShapeClust** provides medoid-based partitional clustering on Kendall’s
2D shape space, together with a bagging ensemble wrapper. This vignette
walks through the typical workflow.

## From landmarks to pre-shapes

A landmark configuration is a `k × 2` matrix of Cartesian coordinates. A
sample of `n` configurations is stored as an array of shape `k × 2 × n`.

``` r
library(ShapeClust)

set.seed(1)
k <- 4; n <- 30
X <- array(rnorm(k * 2 * n), dim = c(k, 2, n))

# Convert to complex pre-shape (columns lie on S^{2k-3}_C)
Z <- landmarks_to_preshape(X)
dim(Z)                    # (k - 1) x n
#> [1]  3 30
head(colSums(Mod(Z)^2))   # all 1
#> [1] 1 1 1 1 1 1
```

The Helmert sub-matrix `H` used internally satisfies
`H^T H = I - 1/k * 1 1^T`, i.e. it is centring in orthonormal
coordinates:

``` r
H <- helmert_submatrix(4)
round(t(H) %*% H, 12)
#>       [,1]  [,2]  [,3]  [,4]
#> [1,]  0.75 -0.25 -0.25 -0.25
#> [2,] -0.25  0.75 -0.25 -0.25
#> [3,] -0.25 -0.25  0.75 -0.25
#> [4,] -0.25 -0.25 -0.25  0.75
round(diag(4) - matrix(1/4, 4, 4), 12)
#>       [,1]  [,2]  [,3]  [,4]
#> [1,]  0.75 -0.25 -0.25 -0.25
#> [2,] -0.25  0.75 -0.25 -0.25
#> [3,] -0.25 -0.25  0.75 -0.25
#> [4,] -0.25 -0.25 -0.25  0.75
```

## Shape-space distances

Two intrinsic distances are provided: the Riemannian (geodesic) distance
`rho` and the full Procrustes distance `dF`. They satisfy
`dF = sin(rho)`:

``` r
D_rho <- rho_matrix(Z)
D_dF  <- dF_matrix(Z)
max(abs(D_dF - sin(D_rho)))
#> [1] 1.110223e-16
```

## Shape-PAM and Shape-rankkmed

The two base clusterers accept either a pre-shape matrix or a
precomputed distance matrix.

``` r
fit_pam <- shape_pam(Z, K = 3)
str(fit_pam)
#> List of 5
#>  $ cluster  : int [1:30] 2 3 1 1 1 1 1 3 2 2 ...
#>  $ medoids  : int [1:3] 6 11 18
#>  $ objective: num 16.9
#>  $ n_iter   : int 4
#>  $ distance : chr "rho"

fit_rkm <- shape_rankkmed(Z, K = 3, m = 5)
table(fit_pam$cluster, fit_rkm$cluster)
#>    
#>     1 2 3
#>   1 1 4 3
#>   2 3 5 0
#>   3 7 3 4
```

## Bagging ensemble

[`bagclust1()`](https://raydonal.github.io/ShapeClust/reference/bagclust1.md)
wraps any base clusterer in a bootstrap ensemble that stabilizes
partitions via majority voting after Hungarian label alignment.

``` r
fit_bag <- bagclust1(Z, K = 3, base = "pam", B = 50, seed = 1)
fit_bag
#> BagClust1 fit
#>   base clusterer: pam 
#>   B = 50 bootstrap replicates
#>   cluster sizes: 9 / 7 / 14
```

## Validation indices

``` r
# Simulate three well-separated clusters
mu_list <- list(
  matrix(c(0,0, 10,0, 10,10, 0,10), 4, 2, byrow = TRUE),
  matrix(c(0,0, 20,0, 20,10, 0,10), 4, 2, byrow = TRUE),
  matrix(c(0,0, 15,0, 10,10, 0,10), 4, 2, byrow = TRUE)
)
sim <- sim_landmark_gaussian(mu_list, c(30, 30, 30),
                             sigma = 0.5, regime = "isotropic",
                             seed = 1)
Z_sim <- landmarks_to_preshape(sim$X)
fit <- shape_pam(Z_sim, K = 3)
cat("RI  =", round(rand_index(sim$group, fit$cluster), 3), "\n")
#> RI  = 1
cat("FMI =", round(fmi_index(sim$group, fit$cluster), 3), "\n")
#> FMI = 1
```

## Between/within variances (VEG, VIG)

``` r
v <- riemannian_variances(sim$X, sim$group)
cat(sprintf("VEG = %.3f, VIG = %.3f, ratio = %.2f\n",
            v$VEG, v$VIG, v$ratio))
#> VEG = 0.040, VIG = 0.004, ratio = 10.32
```

## Model selection

``` r
res_K <- select_K(Z_sim, K_range = 2:5, method = "silhouette")
res_K$K_hat
#> [1] 3
res_K$scores
#>         2         3         4         5 
#> 0.6065260 0.6547472 0.5299256 0.5265576
```

## Reproducing the sand-grains analysis

``` r
library(shapes)
data(sand)
Z_sand <- landmarks_to_preshape(sand$x)
truth  <- as.integer(sand$group)

fit_std <- shape_pam(Z_sand, K = 2)
fit_bag <- bagclust1(Z_sand, K = 2, base = "pam", B = 100, seed = 1)

cat("Standard PAM  : RI =", round(rand_index(truth, fit_std$cluster), 3),
    " FMI =", round(fmi_index(truth, fit_std$cluster), 3), "\n")
#> Standard PAM  : RI = 0.495  FMI = 0.487
cat("Bagged PAM    : RI =", round(rand_index(truth, fit_bag$cluster), 3),
    " FMI =", round(fmi_index(truth, fit_bag$cluster), 3), "\n")
#> Bagged PAM    : RI = 0.551  FMI = 0.544
cat("Relative gain : RI =",
    round(relative_gain(rand_index(truth, fit_std$cluster),
                        rand_index(truth, fit_bag$cluster)), 2), "%\n")
#> Relative gain : RI = 11.34 %
```

## Session info

``` r
sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.5 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.12.0 
#> LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.12.0  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
#>  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
#>  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#>  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
#> 
#> time zone: America/Bahia
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] shapes_1.2.8     ShapeClust_0.4.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] fitdistrplus_1.2-6   cli_3.6.6            knitr_1.52          
#>  [4] rlang_1.3.0          xfun_0.61            otel_0.2.0          
#>  [7] textshaping_1.0.5    clue_0.3-68          jsonlite_2.0.0      
#> [10] htmltools_0.5.9      ragg_1.5.2           sass_0.4.10         
#> [13] scatterplot3d_0.3-45 rmarkdown_2.32       grid_4.6.1          
#> [16] evaluate_1.0.5       jquerylib_0.1.4      MASS_7.3-66         
#> [19] fastmap_1.2.0        base64enc_0.1-6      yaml_2.3.12         
#> [22] lifecycle_1.0.5      rgl_1.3.36           cluster_2.1.8.3     
#> [25] compiler_4.6.1       fs_2.1.0             Rcpp_1.1.2          
#> [28] htmlwidgets_1.6.4    systemfonts_1.3.2    lattice_0.23-1      
#> [31] digest_0.6.39        R6_2.6.1             splines_4.6.1       
#> [34] minpack.lm_1.2-4     magrittr_2.0.5       Matrix_1.7-6        
#> [37] bslib_0.12.0         tools_4.6.1          survival_3.8-12     
#> [40] pkgdown_2.2.1        cachem_1.1.0         desc_1.4.3
```
