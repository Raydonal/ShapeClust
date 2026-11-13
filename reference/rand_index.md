# Rand Index

Rand (1971).

## Usage

``` r
rand_index(truth, pred)
```

## Arguments

- truth:

  Integer vector of ground-truth labels.

- pred:

  Integer vector of predicted labels, same length.

## Value

A numeric in \[0, 1\].

## Examples

``` r
rand_index(c(1,1,2,2,3,3), c(1,1,2,2,3,3))
#> [1] 1
rand_index(c(1,1,2,2,3,3), c(1,2,2,1,3,3))
#> [1] 0.7333333
```
