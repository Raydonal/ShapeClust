# Helmert sub-matrix

Returns the (k-1) x k Helmert sub-matrix used to remove translation from
a landmark configuration.

## Usage

``` r
helmert_submatrix(k)
```

## Arguments

- k:

  Integer. Number of landmarks (\>= 2).

## Value

A real matrix of size (k-1) x k.

## Details

The full Helmert matrix \\H^F\\ is the k x k orthogonal matrix with
first row \\(1/\sqrt{k}, \dots, 1/\sqrt{k})\\ and, for \\j \ge 2\\, jth
row \\(-1, \dots, -1, j-1, 0, \dots, 0)/\sqrt{j(j-1)}\\ where the entry
\\j-1\\ occurs at column \\j\\. Dropping the first row yields the
Helmert sub-matrix \\H\\.

Applying \\H\\ to a configuration removes translation, and \\H^\top H =
I_k - k^{-1} \mathbf{1}\_k \mathbf{1}\_k^\top\\ is the centring matrix –
so \\H\\ is the orthonormal representation of centring.

## Examples

``` r
H <- helmert_submatrix(4)
# Rows are orthonormal:
round(H %*% t(H), 12)
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    1    0
#> [3,]    0    0    1
```
