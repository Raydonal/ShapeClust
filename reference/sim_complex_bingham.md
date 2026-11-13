# Simulate a single Complex Bingham sample

Uses the Simplex-truncation method of Kent, Constable and Er (2004).

## Usage

``` r
sim_complex_bingham(n, lambda, seed = NULL)
```

## Arguments

- n:

  Sample size.

- lambda:

  Numeric vector of concentration parameters \\(\lambda_1, \dots,
  \lambda\_{k-1})\\; note the last eigenvalue is fixed at 0 by
  convention.

- seed:

  Optional integer seed; the caller's random-number stream is left as it
  was.

## Value

A complex matrix of size (k-1) x n whose columns are unit-norm Complex
Bingham draws.

## Examples

``` r
Z <- sim_complex_bingham(n = 5, lambda = rep(9, 3))
round(colSums(Mod(Z)^2), 12)  # unit norm
#> [1] 1 1 1 1 1
```
