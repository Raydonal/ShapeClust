# Convert 2D landmark configurations to complex pre-shape vectors

Given a 3D array of shape (k, 2, n) of Cartesian coordinates, returns
the (k-1) x n complex matrix Z whose columns are the corresponding
unit-norm complex pre-shapes.

## Usage

``` r
landmarks_to_preshape(X)
```

## Arguments

- X:

  Numeric array of shape (k, 2, n). Translation is removed with the
  Helmert sub-matrix and scale by normalising to unit norm.

## Value

A complex matrix of size (k-1) x n.

## Examples

``` r
set.seed(1)
k <- 4; n <- 5
X <- array(rnorm(k * 2 * n), dim = c(k, 2, n))
Z <- landmarks_to_preshape(X)
# Columns have unit modulus:
round(colSums(Mod(Z)^2), 12)
#> [1] 1 1 1 1 1
```
