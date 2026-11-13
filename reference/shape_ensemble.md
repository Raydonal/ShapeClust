# Consensus clustering of shapes in tangent PCA coordinates

Clusters \\k \times m \times n\\ landmark configurations with a
bootstrap ensemble of k-means fits in the leading principal components
of the tangent coordinates, a co-association matrix, and a spectral
partition of that matrix. Unlike
[`bagclust1`](https://raydonal.github.io/ShapeClust/reference/bagclust1.md),
no label alignment is needed, and the result carries a real observed
shape as representative of each cluster.

## Usage

``` r
shape_ensemble(
  X,
  K,
  B = 30L,
  var_explained = 0.9,
  max_dim = 10L,
  seed = NULL,
  base = c("kmeans", "vii"),
  ...
)
```

## Arguments

- X:

  Numeric array of dimension \\k \times m \times n\\, with \\m\\ equal
  to 2 or 3.

- K:

  Number of clusters (\>= 2).

- B:

  Number of bootstrap replicates (default 30).

- var_explained:

  Proportion of tangent variance that fixes the number of principal
  components (default 0.90).

- max_dim:

  Largest number of components kept (default 10). The number used is at
  least `max(2, K - 1)` and at most `max_dim`.

- seed:

  Optional integer seed. The caller's random-number stream is restored
  on exit. With `NULL` the current stream is used.

- base:

  Base clusterer of each replicate: `"kmeans"` (default, the original
  behavior) or `"vii"`, the spherical-variance-per-cluster Gaussian
  mixture of
  [`shape_mixture`](https://raydonal.github.io/ShapeClust/reference/shape_mixture.md)
  (EM with 3 starts), which assigns all \\n\\ observations by maximum
  posterior. Use `"vii"` when clusters have unequal spread.

- ...:

  With `base = "kmeans"`, further arguments for
  [`kmeans`](https://rdrr.io/r/stats/kmeans.html) in each replicate
  (defaults `nstart = 1`, `iter.max = 30`). Ignored for `"vii"`.

## Value

An object of class `"shape_ensemble"`, a list with

- cluster:

  Integer vector of length n, values in 1:K.

- coassoc:

  The n x n co-association matrix.

- dim:

  Number of principal components used.

- representative:

  Integer vector of length K: index of the real observation closest to
  the Procrustes mean of each cluster.

- centers:

  Array \\k \times m \times K\\ with the cluster means.

- B, K, var_explained, base:

  The options used.

- n_fallback:

  With `base = "vii"`, the number of replicates whose mixture fit fell
  back to k-means.

## Details

The tangent coordinates are those of the Procrustes mean shape
([`shape_mean`](https://raydonal.github.io/ShapeClust/reference/shape_mean.md)):
each configuration is centered, scaled to unit size, rotated onto the
mean without reflection, and projected onto the tangent space at the
mean. A principal component analysis gives the scores; the smallest
number of components with cumulative variance of at least
`var_explained` is used, between `max(2, K - 1)` and `max_dim`: K groups
need at least K - 1 dimensions to be separated by k-means. With K = 3
this is the earlier floor of 2.

Each replicate draws a bootstrap sample of the scores, drops duplicate
rows, runs k-means with a single start, and assigns all \\n\\
observations to the nearest of the \\K\\ centers. The co-association
matrix holds the share of replicates in which two observations get the
same label. The partition is the k-means (20 starts) of the
row-normalized first \\K\\ eigenvectors of the symmetrically normalized
matrix.

The representative of a cluster is the member with the smallest
Riemannian distance \\\rho\\ to the full Procrustes mean of that
cluster.

## See also

[`bagclust1`](https://raydonal.github.io/ShapeClust/reference/bagclust1.md),
[`shape_kmeans`](https://raydonal.github.io/ShapeClust/reference/shape_kmeans.md)

## Examples

``` r
sim <- sim_shapes3d(n = 40, sigma = 1, k = 8, seed = 1)
fit <- shape_ensemble(sim$X, K = 2, B = 10, seed = 1)
fit$dim
#> [1] 2
fit$representative
#> [1] 15 55
table(sim$group, fit$cluster)
#>    
#>      1  2
#>   1 40  0
#>   2  0 40
```
