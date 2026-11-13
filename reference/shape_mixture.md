# Gaussian mixture clustering of shapes with a variance per cluster

Clusters \\k \times m \times n\\ landmark configurations with a Gaussian
mixture with spherical covariance \\\sigma_h^2 I\\ for each cluster
(model `"VII"` of mclust) in the leading principal components of the
tangent coordinates. The EM algorithm is implemented here, with
safeguards for duplicated rows, tiny groups and degenerate components.

## Usage

``` r
shape_mixture(
  X,
  K,
  var_explained = 0.9,
  max_dim = 10L,
  model = "VII",
  nstart = 10L,
  init = c("kmeans++", "kmeans", "random"),
  iter.max = 200L,
  tol = 1e-08,
  var_floor = 1e-06,
  seed = NULL,
  ...
)
```

## Arguments

- X:

  Numeric array of dimension \\k \times m \times n\\, with \\m\\ equal
  to 2 or 3.

- K:

  Number of clusters (\>= 2).

- var_explained:

  Proportion of tangent variance that fixes the number of principal
  components (default 0.90).

- max_dim:

  Largest number of components kept (default 10). The number used is at
  least `max(2, K - 1)`, as in
  [`shape_ensemble`](https://raydonal.github.io/ShapeClust/reference/shape_ensemble.md).

- model:

  Covariance model. Only `"VII"` (spherical, unequal volume) is
  implemented.

- nstart:

  Number of EM starts (default 10); the best log-likelihood among starts
  without a variance at the floor is kept.

- init:

  Initialization of each start: `"kmeans++"` (default, k-means++
  seeding, nearest-seed labels), `"kmeans"` (seeding followed by Lloyd
  iterations) or `"random"` (random labels).

- iter.max:

  Maximum EM iterations per start (default 200).

- tol:

  Relative log-likelihood tolerance for convergence (default 1e-8).

- var_floor:

  Variance floor, as a fraction of the mean per-dimension variance of
  the scores (default 1e-6).

- seed:

  Optional integer seed. The caller's random-number stream is restored
  on exit. With `NULL` the current stream is used.

- ...:

  Not used; reserved.

## Value

An object of class `"shape_mixture"`, a list with

- cluster:

  Integer vector of length n, maximum posterior label.

- posterior:

  n x K matrix of posterior probabilities.

- means:

  K x dim matrix of component means in score coordinates.

- variances:

  Length-K vector of spherical variances.

- proportions:

  Length-K vector of mixing proportions.

- loglik:

  Log-likelihood of the fit.

- dim:

  Number of principal components used.

- representative:

  Integer vector of length K: index of the real observation closest
  (Riemannian \\\rho\\) to the Procrustes mean of each cluster.

- centers:

  Array \\k \times m \times K\\ with the cluster means.

- converged:

  Whether EM met the tolerance.

- fallback:

  `TRUE` if the k-means fallback was used.

- K, model, var_explained:

  The options used.

## Details

The scores are those of `shape_ensemble`: tangent coordinates at the
Procrustes mean, then PCA, with the smallest number of components
reaching `var_explained`, between `max(2, K - 1)` and `max_dim`.

Why a variance per cluster: k-means (and so a k-means ensemble) assigns
to the nearest center, which is the Bayes rule only when all clusters
have the same spread. With unequal spread the Bayes boundary is a sphere
around the compact cluster and the bisector misassigns the tail of the
diffuse one. This model estimates that spread.

Safeguards. Variances are floored at `var_floor` times the mean
per-dimension variance of the scores, so duplicated rows (bootstrap
resamples) or a group of identical points never make the likelihood
singular. A component with less than one observation of total weight is
restarted at the observation worst explained by the current model (at
most `5 * K` restarts per start, then that start is abandoned). Starts
that end with a variance at the floor are treated as likelihood spikes
and are used only if no regular start exists. If no start succeeds, the
partition is the k-means (5 starts) partition with per-group spherical
variances, `fallback` is `TRUE` and `converged` is `FALSE`.

## See also

[`shape_ensemble`](https://raydonal.github.io/ShapeClust/reference/shape_ensemble.md),
[`shape_kmeans`](https://raydonal.github.io/ShapeClust/reference/shape_kmeans.md)

## Examples

``` r
sim <- sim_shapes3d(n = 40, sigma = 1, k = 8, seed = 1)
fit <- shape_mixture(sim$X, K = 2, seed = 1)
fit$variances
#> [1] 0.001442267 0.002163340
table(sim$group, fit$cluster)
#>    
#>      1  2
#>   1  0 40
#>   2 40  0
```
