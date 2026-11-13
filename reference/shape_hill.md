# Hill-climbing medoid search on a shape distance matrix

Random-neighbor hill climbing for K medoids, ported from
`hill_climbing_shapes()` of the 3D bagging scripts.

## Usage

``` r
shape_hill(
  D,
  K,
  max_iter = 300L,
  from_best = FALSE,
  initial = NULL,
  distance = c("rho", "dF", "dP"),
  reflection = FALSE
)
```

## Arguments

- D:

  Numeric \\n \times n\\ distance matrix, or a real landmark array \\k
  \times m \times n\\, or a complex pre-shape matrix, from which
  distances are computed with `distance`.

- K:

  Number of clusters (\>= 2).

- max_iter:

  Number of neighbors tried (default 300, as in the scripts).

- from_best:

  Logical. If `FALSE` (default), faithful to the script, each neighbor
  is a modification of the *initial* medoid set. If `TRUE`, each
  neighbor modifies the best set found so far, which is hill climbing in
  the usual sense.

- initial:

  Optional vector of length K with initial medoid indices (1-based). If
  `NULL` they are drawn with
  [`sample.int()`](https://rdrr.io/r/base/sample.html).

- distance:

  Distance used if `D` is not a distance matrix: one of `"rho"`
  (default), `"dF"`, `"dP"`.

- reflection:

  Passed to
  [`shape_distance`](https://raydonal.github.io/ShapeClust/reference/shape_distance.md)
  when `D` is a landmark array.

## Value

A list with components:

- cluster:

  Integer vector of length n, values in 1:K: the index of the nearest
  medoid.

- medoids:

  Integer vector of length K, 1-based object indices.

- objective:

  Sum of the distances of the objects to their nearest medoid.

- n_iter:

  Number of neighbors that improved the objective.

- distance:

  Distance name, or `NA` if `D` was a distance matrix.

## Details

Starting from `initial`, each of the `max_iter` iterations replaces one
randomly chosen medoid by a randomly chosen object, redrawing
(cumulatively, on the same candidate set) while the set contains a
repeated object, and keeps the candidate set if its objective is
strictly lower. Objects are assigned with
[`max.col()`](https://rdrr.io/r/base/maxCol.html)-style random tie
breaking, as in the script, using R's generator, so for the same seed
and distances the result equals that of the script.

In the script every neighbor is built from the initial set, never from
the best set found, so it is a random search around the starting point
and not a hill climb in the strict sense. This is the default so that
results can be compared with the published scripts; `from_best = TRUE`
gives the usual hill climbing.

## See also

[`shape_pam`](https://raydonal.github.io/ShapeClust/reference/shape_pam.md),
[`shape_clarans`](https://raydonal.github.io/ShapeClust/reference/shape_clarans.md),
[`bagclust1`](https://raydonal.github.io/ShapeClust/reference/bagclust1.md)

## Examples

``` r
sim <- sim_shapes3d(n = 10, sigma = 0.5, seed = 1)
D <- shape_distance(sim$X)
set.seed(1)
fit <- shape_hill(D, K = 2, max_iter = 100)
table(fit$cluster, sim$group)
#>    
#>      1  2
#>   1 10  0
#>   2  0 10
```
