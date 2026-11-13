# Mean shape of a sample of 3D or 2D landmark configurations

Two notions of mean shape for real landmark arrays: the full Procrustes
mean, and the MDS (or \\\Phi\\) mean, an extrinsic mean obtained from
the Veronese-Whitney embedding.

## Usage

``` r
shape_mean(X, type = c("procrustes", "phi"), tol = 1e-12, max_iter = 1000L)
```

## Arguments

- X:

  Numeric array of dimension \\k \times m \times n\\, with \\m\\ equal
  to 2 or 3.

- type:

  `"procrustes"` (default) or `"phi"`.

- tol, max_iter:

  Convergence tolerance (relative change of the objective) and iteration
  cap of the Procrustes iteration. Used for `type = "procrustes"` only:
  the `"phi"` mean has an explicit formula and needs the Procrustes mean
  only to choose between a shape and its mirror image, for which a
  looser fixed tolerance is used.

## Value

A \\k \times m\\ matrix: a centered configuration of unit size
(Frobenius norm 1) representing the mean shape. A warning is given if
the Procrustes iteration does not reach `tol` within `max_iter` steps.

## Details

**Procrustes.** The full Procrustes mean minimizes \\\sum_i d_F^2(x_i,
\mu)\\. It is computed by the fixed-point iteration \\\mu \leftarrow
\mathrm{normalize}(\sum_i \cos\rho_i \\ z_i R_i)\\, where \\R_i \in
SO(m)\\ aligns the pre-shape \\z_i\\ to \\\mu\\, started at the first
configuration. Each step does not decrease \\\sum_i \cos^2\rho_i\\. For
\\m = 2\\ it agrees with the dominant eigenvector of \\\sum_i z_i
z_i^\*\\. It is the mean that
[`shapes::procGPA()`](https://rdrr.io/pkg/shapes/man/procGPA.html)
returns in `$mshape`, up to the looser convergence tolerance of that
function.

**phi.** The MDS mean takes the \\m\\ leading eigenvectors \\v_j\\
(eigenvalues \\\lambda_j\\) of \\n^{-1}\sum_i z_i z_i^\top\\, forms the
pre-shape \\\[\sqrt{\lambda_1} v_1, \dots, \sqrt{\lambda_m} v_m\] /
\sqrt{\lambda_1 + \dots + \lambda_m}\\ and lifts it with the transposed
Helmert sub-matrix. This is `shapes::MDSshape(x)` with its defaults
(`alpha = 1`, `projalpha = 1/2`). Because \\z z^\top\\ does not change
under a reflection of \\z\\, the last axis is flipped if that brings the
mean closer to the Procrustes mean, exactly as in
[`shapes::MDSshape`](https://rdrr.io/pkg/shapes/man/shapes-internal.html).

## References

Kent, J. T. (1994). The complex Bingham distribution and shape analysis.
*Journal of the Royal Statistical Society B*, 56, 285–299.
[doi:10.1111/j.2517-6161.1994.tb01978.x](https://doi.org/10.1111/j.2517-6161.1994.tb01978.x)

Dryden, I. L. and Mardia, K. V. (2016). *Statistical Shape Analysis,
with Applications in R*, 2nd ed. Wiley, chapters 6 and 7.

## See also

[`shape_kmeans`](https://raydonal.github.io/ShapeClust/reference/shape_kmeans.md),
which updates its centers with these means.

## Examples

``` r
set.seed(2)
base <- matrix(rnorm(8 * 3), 8, 3)
X <- array(rep(base, 20), dim = c(8, 3, 20)) +
  array(rnorm(8 * 3 * 20, sd = 0.05), dim = c(8, 3, 20))
mu <- shape_mean(X)
nu <- shape_mean(X, type = "phi")
# Both means are close to the template
shape_distance(array(c(mu, nu, base), dim = c(8, 3, 3)))
#>              [,1]         [,2]       [,3]
#> [1,] 0.0000000000 0.0001997356 0.01004940
#> [2,] 0.0001997356 0.0000000000 0.01003735
#> [3,] 0.0100493984 0.0100373527 0.00000000
```
