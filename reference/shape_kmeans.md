# Shape k-means for real landmark arrays

Lloyd's algorithm for clustering configurations of \\k\\ landmarks in
two or three dimensions. Observations are assigned to the nearest center
in shape space and each center is then replaced by the mean shape of its
cluster, with several random restarts.

## Usage

``` r
shape_kmeans(
  X,
  K,
  mean = c("procrustes", "phi"),
  alg_steps = 10L,
  niter = 10L,
  tol = 1e-04,
  distance = c("rho", "dF", "dP")
)
```

## Arguments

- X:

  Numeric array of dimension \\k \times m \times n\\, with \\m\\ equal
  to 2 or 3.

- K:

  Number of clusters (\>= 2).

- mean:

  Mean shape used to update the centers: `"procrustes"` (full Procrustes
  mean, as
  [`shapes::procGPA()`](https://rdrr.io/pkg/shapes/man/procGPA.html);
  default) or `"phi"` (MDS mean, as
  [`shapes::MDSshape()`](https://rdrr.io/pkg/shapes/man/shapes-internal.html)).
  See
  [`shape_mean`](https://raydonal.github.io/ShapeClust/reference/shape_mean.md).

- alg_steps:

  Maximum number of Lloyd steps per restart (default 10).

- niter:

  Number of random restarts (default 10).

- tol:

  Stopping rule: a restart stops when the relative decrease of the
  objective between two steps is below `tol` (default 1e-4).

- distance:

  Distance in which the objective is measured, one of `"rho"` (default),
  `"dF"` or `"dP"`. All three are increasing functions of the same
  cosine, so assignments do not depend on this choice; the objective
  value and the stopping rule do.

## Value

A list with components:

- cluster:

  Integer vector of length n, values in 1:K.

- centers:

  Numeric array \\k \times m \times K\\ with the cluster means, each
  centered with unit norm.

- objective:

  Mean squared distance between each observation and the center it was
  assigned to, at the step that gave the returned clustering (the
  centers in force at the start of that step).

- n_iter:

  Index of the Lloyd step, within the restart that gave the returned
  clustering, at which the objective was lowest.

- initials:

  List with the K observations drawn as initial centers in each restart.

- mean, distance:

  The options used.

## Details

The algorithm follows `LloydShapes()` of the scripts used for the 3D
bagging study. Initial centers are `K` observations drawn with
[`sample.int()`](https://rdrr.io/r/base/sample.html) at the start of
each restart. In each step every observation goes to its nearest center
([`max.col()`](https://rdrr.io/r/base/maxCol.html) with random tie
breaking, as in the scripts), the centers are updated, and the objective
is evaluated with the distances of the assignment step, i.e. before the
update. The clustering and centers kept are those of the step with the
lowest objective over all steps of all restarts. A cluster with a single
member has that member as its center. If a cluster becomes empty its
center is left unchanged (the script stops with an error in that case).
Reflections are not removed (see
[`shape_distance`](https://raydonal.github.io/ShapeClust/reference/shape_distance.md)).

The caller's random-number stream is advanced; use
[`set.seed()`](https://rdrr.io/r/base/Random.html) for reproducibility.

## References

Dryden, I. L. and Mardia, K. V. (2016). *Statistical Shape Analysis,
with Applications in R*, 2nd ed. Wiley.

## See also

[`shape_mean`](https://raydonal.github.io/ShapeClust/reference/shape_mean.md),
[`bagclust1`](https://raydonal.github.io/ShapeClust/reference/bagclust1.md)
with `base = "kmeans"`.

## Examples

``` r
sim <- sim_shapes3d(n = 10, sigma = 0.5, seed = 1)
set.seed(1)
fit <- shape_kmeans(sim$X, K = 2)
table(fit$cluster, sim$group)
#>    
#>      1  2
#>   1  0 10
#>   2 10  0
fit$objective
#> [1] 0.004469681
```
