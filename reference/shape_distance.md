# Shape distances for real landmark arrays

Pairwise distances in Kendall's shape space for configurations of \\k\\
landmarks in \\m = 2\\ or \\3\\ dimensions, computed in C++.

## Usage

``` r
shape_distance(
  X,
  distance = c("rho", "dF", "dP"),
  reflection = FALSE,
  Y = NULL
)
```

## Arguments

- X:

  Numeric array of dimension \\k \times m \times n\\, with \\m\\ equal
  to 2 or 3 and \\k \ge m + 1\\.

- distance:

  One of `"rho"` (Riemannian, default), `"dF"` (full Procrustes) or
  `"dP"` (partial Procrustes).

- reflection:

  Logical. If `FALSE` (default) reflections are not removed, so a
  configuration and its mirror image are different shapes. If `TRUE`,
  shapes are compared up to reflection too.

- Y:

  Optional second array of dimension \\k \times m \times n_y\\. If
  given, the cross-distances between `X` and `Y` are returned instead.

## Value

A numeric matrix: the symmetric \\n \times n\\ matrix of distances
between the configurations of `X`, or the \\n \times n_y\\ matrix
between those of `X` and `Y`.

## Details

Each configuration is centered with the Helmert sub-matrix and scaled to
unit norm, giving a \\(k-1) \times m\\ pre-shape \\z\\. If \\s_1 \ge
\dots \ge s_m\\ are the singular values of \\z^\top w\\, the cosine of
the Riemannian distance is \\c = s_1 + \dots + s\_{m-1} \pm s_m\\, with
the sign of \\s_m\\ equal to the sign of \\\det(z^\top w)\\ unless
`reflection = TRUE`, in which case all signs are positive. Then \$\$\rho
= \arccos c, \quad d_F = \sqrt{1 - c^2}, \quad d_P = \sqrt{2(1 -
c)},\$\$ with \\c\\ clamped to \\\[0, 1\]\\. For \\m = 2\\ without
reflection \\c\\ equals \\\|z^\* w\|\\, so the result agrees with
[`rho_matrix`](https://raydonal.github.io/ShapeClust/reference/rho_matrix.md)
and
[`dF_matrix`](https://raydonal.github.io/ShapeClust/reference/dF_matrix.md).
Two configurations that are exactly equal get distance 0.

## References

Dryden, I. L. and Mardia, K. V. (2016). *Statistical Shape Analysis,
with Applications in R*, 2nd ed. Wiley, sections 4.1 and 4.2.

## See also

[`shape_mean`](https://raydonal.github.io/ShapeClust/reference/shape_mean.md),
[`shape_kmeans`](https://raydonal.github.io/ShapeClust/reference/shape_kmeans.md)

## Examples

``` r
set.seed(1)
X <- array(rnorm(8 * 3 * 6), dim = c(8, 3, 6))
D <- shape_distance(X)
round(D[1:3, 1:3], 3)
#>       [,1]  [,2]  [,3]
#> [1,] 0.000 1.241 0.951
#> [2,] 1.241 0.000 1.141
#> [3,] 0.951 1.141 0.000
all.equal(shape_distance(X, "dF"), sin(D))
#> [1] TRUE

# A mirror image is a different shape unless reflections are allowed
Xm <- X; Xm[, 1, ] <- -Xm[, 1, ]
shape_distance(X[, , 1, drop = FALSE], Y = Xm[, , 1, drop = FALSE])
#>          [,1]
#> [1,] 0.909222
shape_distance(X[, , 1, drop = FALSE], Y = Xm[, , 1, drop = FALSE],
               reflection = TRUE)
#>      [,1]
#> [1,]    0
```
