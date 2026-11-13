# Full Procrustes distance matrix

\$\$d_F(z_1, z_2) = \sqrt{1 - \|z_1^\* z_2\|^2}\$\$

## Usage

``` r
dF_matrix(Z)
```

## Arguments

- Z:

  Complex pre-shape matrix of size (k-1) x n.

## Value

A symmetric numeric matrix.

## Examples

``` r
set.seed(1)
Z <- landmarks_to_preshape(array(rnorm(4 * 2 * 5), dim = c(4, 2, 5)))
all.equal(dF_matrix(Z), sin(rho_matrix(Z)))
#> [1] TRUE
```
