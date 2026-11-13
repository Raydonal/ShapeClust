# Align cluster labels via optimal permutation

Given two label vectors on the same set of n objects, returns the
relabelling of `a` to `b` that maximises the number of matches.

## Usage

``` r
align_labels(a, b, K)
```

## Arguments

- a:

  Integer vector of labels in 1:K.

- b:

  Integer vector of labels in 1:K, same length as `a`.

- K:

  Number of clusters.

## Value

Integer vector of the relabelled `a`.

## Details

Uses the Hungarian algorithm via
[`clue::solve_LSAP`](https://rdrr.io/pkg/clue/man/solve_LSAP.html) if
available; otherwise falls back to exhaustive enumeration of \\S_K\\
(feasible for K \<= 8).

## Examples

``` r
a <- c(2, 2, 3, 3, 1, 1)
b <- c(1, 1, 2, 2, 3, 3)
align_labels(a, b, K = 3)  # equals b
#> [1] 1 1 2 2 3 3
```
