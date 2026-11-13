# BagClust1 for shape clustering

Each bootstrap replicate fits the base clusterer to a resample and then
assigns every one of the n original observations to its nearest
bootstrap medoid or center, so all observations vote in every replicate
(Dudoit and Fridlyand, 2003, BagClust1).

## Usage

``` r
bagclust1(
  Z,
  K,
  base = c("pam", "rankkmed", "kmeans", "clarans", "hill"),
  B = 100L,
  distance = c("rho", "dF", "dP"),
  seed = NULL,
  parallel = FALSE,
  n_cores = NULL,
  reflection = FALSE,
  ...
)
```

## Arguments

- Z:

  Either a complex pre-shape matrix of size (k-1) x n (2D shapes, from
  [`landmarks_to_preshape`](https://raydonal.github.io/ShapeClust/reference/landmarks_to_preshape.md)),
  or a real landmark array of dimension k x m x n with m = 2 or 3. Base
  clusterer `"kmeans"` needs the real array.

- K:

  Number of clusters.

- base:

  Base clusterer. One of "pam" (default), "rankkmed", "kmeans"
  ([`shape_kmeans`](https://raydonal.github.io/ShapeClust/reference/shape_kmeans.md)),
  "clarans"
  ([`shape_clarans`](https://raydonal.github.io/ShapeClust/reference/shape_clarans.md))
  or "hill"
  ([`shape_hill`](https://raydonal.github.io/ShapeClust/reference/shape_hill.md)).
  The first, second, fourth and fifth fit on the distance matrix of the
  resample and label observations by the nearest medoid; "kmeans" fits
  on the landmarks of the resample and labels observations by the
  nearest center (mean shape).

- B:

  Number of bootstrap replicates (default 100).

- distance:

  One of "rho" (default), "dF" or "dP" ("dP", the partial Procrustes
  distance, is new in the development version). It builds the distance
  matrix the medoid methods work on and the distance to the centers of
  "kmeans".

- seed:

  Optional integer seed. The result is then reproducible and identical
  with and without `parallel`; the caller's random-number stream is left
  as it was.

- parallel:

  Logical. If TRUE, use
  [`parallel::mclapply`](https://rdrr.io/r/parallel/mclapply.html)
  (Unix/macOS) or
  [`parallel::parLapply`](https://rdrr.io/r/parallel/clusterApply.html)
  (Windows).

- n_cores:

  Number of cores to use if `parallel = TRUE`. Default =
  `max(1, parallel::detectCores() - 1)`.

- reflection:

  Logical, for real landmark arrays only: if TRUE, reflections are
  removed too (see
  [`shape_distance`](https://raydonal.github.io/ShapeClust/reference/shape_distance.md)).
  Not available for `base = "kmeans"`.

- ...:

  Additional arguments passed to the base clusterer (e.g., `m = 7` for
  rankkmed, `mean = "phi"` for kmeans, `numlocal` and `maxneighbor` for
  clarans, `max_iter` for hill).

## Value

A list of class `"bagclust1"` with components:

- cluster:

  Integer vector of length n, the majority-vote consensus labels.

- votes:

  n x K integer matrix of vote counts (rows sum to B).

- reference:

  Integer vector, the reference (non-bagged) clustering used for label
  alignment.

- B:

  Number of bootstrap replicates.

- base:

  Base clusterer name.

## References

Dudoit, S. and Fridlyand, J. (2003). Bagging to improve the accuracy of
a clustering procedure. *Bioinformatics*, 19(9), 1090–1099.
[doi:10.1093/bioinformatics/btg038](https://doi.org/10.1093/bioinformatics/btg038)

## Examples

``` r
set.seed(1)
k <- 4; n <- 60
# 3 clusters via mean shift
X <- array(0, dim = c(k, 2, n))
for (i in seq_len(n)) {
  g <- ((i - 1) %% 3) + 1
  mu <- matrix(c(0,0, 1,0, 0,1, 1,1), 4, 2, byrow = TRUE)
  X[,,i] <- mu + matrix(rnorm(k*2, sd = 0.1 * g), k, 2)
}
Z <- landmarks_to_preshape(X)
fit <- bagclust1(Z, K = 3, base = "pam", B = 30, seed = 1)
table(fit$cluster)
#> 
#>  1  2  3 
#>  9 20 31 

# 3D shapes: a real k x 3 x n array, here with the k-means base clusterer
sim <- sim_shapes3d(n = 10, sigma = 0.5, seed = 1)
fit3 <- bagclust1(sim$X, K = 2, base = "kmeans", B = 10, seed = 1)
table(fit3$cluster, sim$group)
#>    
#>      1  2
#>   1  0 10
#>   2 10  0
```
