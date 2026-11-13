# Fowlkes-Mallows Index

Fowlkes and Mallows (1983).

## Usage

``` r
fmi_index(truth, pred)
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
fmi_index(c(1,1,2,2,3,3), c(1,1,2,2,3,3))
#> [1] 1
fmi_index(c(1,1,2,2,3,3), c(1,2,2,1,3,3))
#> [1] 0.3333333
```
