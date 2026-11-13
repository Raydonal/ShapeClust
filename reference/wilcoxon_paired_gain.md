# Paired one-sided Wilcoxon signed-rank test

Convenience wrapper for the paired one-sided test of
`median(V_bag - V_std) > 0`.

## Usage

``` r
wilcoxon_paired_gain(V_std, V_bag)
```

## Arguments

- V_std:

  Numeric vector of standard-variant values.

- V_bag:

  Numeric vector of bagged-variant values, same length.

## Value

A numeric p-value.

## Examples

``` r
set.seed(1)
v_std <- runif(30, 0.5, 0.8)
v_bag <- v_std + rnorm(30, 0.05, 0.02)
wilcoxon_paired_gain(v_std, v_bag)
#> [1] 9.126857e-07
```
