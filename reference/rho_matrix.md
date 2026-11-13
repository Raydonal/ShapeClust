# Riemannian geodesic distance matrix on Kendall's 2D shape space

Returns the n x n symmetric matrix of intrinsic (Riemannian) distances
between the columns of a complex pre-shape matrix Z, in radians in
\\\[0, \pi/2\]\\.

## Usage

``` r
rho_matrix(Z)
```

## Arguments

- Z:

  Complex matrix of size (k-1) x n, unit-norm columns.

## Value

A symmetric numeric matrix.

## Details

\$\$\rho(z_1, z_2) = \arccos \|z_1^\* z_2\|\$\$

Cost: \\\mathcal{O}(n^2 (k-1))\\ via a single BLAS-backed Hermitian Gram
matrix.

## Examples

``` r
k <- 4; n <- 5
set.seed(1)
X <- array(rnorm(k * 2 * n), dim = c(k, 2, n))
Z <- landmarks_to_preshape(X)
D <- rho_matrix(Z)
isSymmetric(D)
#> [1] TRUE
```
