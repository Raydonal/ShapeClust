# Select K by silhouette width or the gap statistic

Select K by silhouette width or the gap statistic

## Usage

``` r
select_K(
  Z,
  K_range = 2:8,
  method = c("silhouette", "gap"),
  base = c("pam", "rankkmed"),
  B_gap = 20L,
  distance = c("rho", "dF")
)
```

## Arguments

- Z:

  Complex pre-shape matrix (k-1) x n, or an n x n distance matrix
  (silhouette only: the gap statistic needs the pre-shapes to draw its
  null reference).

- K_range:

  Integer vector of candidate values of K.

- method:

  One of "silhouette" (default) or "gap".

- base:

  Base clusterer to use for each candidate K. One of "pam" (default) or
  "rankkmed".

- B_gap:

  Number of null reference samples for the gap statistic.

- distance:

  One of "rho" (default), "dF".

## Value

A list with:

- K_hat:

  Selected K.

- scores:

  Named numeric vector of scores (silhouette width or gap statistic) at
  each candidate K.

- se:

  For `"gap"`, the simulation standard errors \\s_K\\; `NULL` otherwise.

- method:

  Selection method used.

## Details

The gap statistic compares \\\log W_K\\ with its mean under `B_gap`
samples of the same size drawn uniformly on the shape space (uniform
pre-shapes, i.e. a Complex Bingham with zero concentration). The
selected K is the smallest one with \\\mathrm{Gap}(K) \ge
\mathrm{Gap}(K') - s\_{K'}\\, where \\K'\\ is the next candidate
(Tibshirani, Walther and Hastie, 2001).

## References

Rousseeuw, P. J. (1987). Silhouettes: a graphical aid to the
interpretation and validation of cluster analysis. *Journal of
Computational and Applied Mathematics*, 20, 53–65.
[doi:10.1016/0377-0427(87)90125-7](https://doi.org/10.1016/0377-0427%2887%2990125-7)

Tibshirani, R., Walther, G. and Hastie, T. (2001). Estimating the number
of clusters in a data set via the gap statistic. *Journal of the Royal
Statistical Society B*, 63(2), 411–423.
[doi:10.1111/1467-9868.00293](https://doi.org/10.1111/1467-9868.00293)

## Examples

``` r
set.seed(1)
# Three clusters
mu_list <- list(
  matrix(c(0,0, 10,0, 10,10, 0,10), 4, 2, byrow = TRUE),
  matrix(c(0,0, 20,0, 20,10, 0,10), 4, 2, byrow = TRUE),
  matrix(c(0,0, 15,0, 10,10, 0,10), 4, 2, byrow = TRUE)
)
out <- sim_landmark_gaussian(mu_list, c(15, 15, 15), sigma = 0.5)
Z <- landmarks_to_preshape(out$X)
res <- select_K(Z, K_range = 2:5, method = "silhouette")
res$K_hat
#> [1] 3
```
