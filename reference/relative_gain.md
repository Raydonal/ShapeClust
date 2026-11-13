# Relative Gain

Percentage increase of a validation measure from a standard to a bagged
variant.

## Usage

``` r
relative_gain(V_std, V_bag)
```

## Arguments

- V_std:

  Numeric. Value of the validation index for the standard variant.

- V_bag:

  Numeric. Value for the bagged variant.

## Value

Numeric percentage.

## Examples

``` r
relative_gain(0.80, 0.88)  # +10 percent
#> [1] 10
```
