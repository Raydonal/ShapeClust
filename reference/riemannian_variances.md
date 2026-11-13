# Riemannian between- and within-group variances (VEG, VIG)

Estimates `VEG = sum_i rho(mu_i, mu)^2 / (K-1)` and
`VIG = sum_i sum_j rho(X_ij, mu_i)^2 / (n-K)` using the Riemannian
distance rho on Kendall's shape space.

## Usage

``` r
riemannian_variances(X, cluster)
```

## Arguments

- X:

  Landmark array (k x 2 x n).

- cluster:

  Vector of group labels, length n.

## Value

A list with VEG, VIG, and the ratio VEG/VIG.

## Details

The group and overall means are full Procrustes means: the dominant
eigenvector of the complex sum of squares \\\sum_j z_j z_j^\*\\ of the
pre-shapes (Kent, 1994). This is the Frechet mean under the full
Procrustes distance.

## References

Kent, J. T. (1994). The complex Bingham distribution and shape analysis.
*Journal of the Royal Statistical Society B*, 56(2), 285–299.
[doi:10.1111/j.2517-6161.1994.tb01978.x](https://doi.org/10.1111/j.2517-6161.1994.tb01978.x)

## Examples

``` r
mu_list <- list(
  matrix(c(0,0, 10,0, 10,10, 0,10), 4, 2, byrow = TRUE),
  matrix(c(0,0, 20,0, 20,10, 0,10), 4, 2, byrow = TRUE)
)
out <- sim_landmark_gaussian(mu_list, c(15, 15), sigma = 1, seed = 1)
riemannian_variances(out$X, out$group)
#> $VEG
#> [1] 0.06036169
#> 
#> $VIG
#> [1] 0.01257175
#> 
#> $ratio
#> [1] 4.801375
#> 
```
