# Shape-rankkmed: Ranked K-medoids on Kendall's 2D shape space

Shape-rankkmed: Ranked K-medoids on Kendall's 2D shape space

## Usage

``` r
shape_rankkmed(
  Z,
  K,
  distance = c("rho", "dF"),
  m = 5L,
  max_iter = 30L,
  initial = NULL
)
```

## Arguments

- Z:

  Complex pre-shape matrix of size (k-1) x n, a numeric distance matrix
  of size n x n, or a real landmark array of dimension k x m x n with m
  = 2 or 3 (distances from
  [`shape_distance`](https://raydonal.github.io/ShapeClust/reference/shape_distance.md)).

- K:

  Number of clusters (\>= 2).

- distance:

  Character. One of "rho" (default), "dF".

- m:

  Neighbourhood size for cluster formation (default 5; Zadegan et
  al. (2013) recommend 5-15).

- max_iter:

  Integer. Maximum number of iterations (default 30). The loop stops
  earlier at a fixed point (medoids unchanged and no relocation), from
  which further iterations cannot move.

- initial:

  Optional vector of length K giving initial medoid indices (1-based).

## Value

A list with components:

- cluster:

  Integer vector of length n, values in 1:K.

- medoids:

  Integer vector of length K.

- objective:

  Sum of distances of each object to its medoid.

- n_iter:

  Number of iterations performed.

- distance:

  Distance name used.

## Details

With the same distance matrix, initial medoids and random-number stream,
the partition equals that of
[`kmed::rankkmed`](https://rdrr.io/pkg/kmed/man/rankkmed.html) run with
`iterate = max_iter`.

## References

Zadegan, S. M. R., Mirzaie, M. and Sadoughi, F. (2013). Ranked
k-medoids: a fast and accurate rank-based partitioning algorithm for
clustering large datasets. *Knowledge-Based Systems*, 39, 133–143.
[doi:10.1016/j.knosys.2012.10.012](https://doi.org/10.1016/j.knosys.2012.10.012)

## Examples

``` r
set.seed(1)
k <- 4; n <- 30
X <- array(rnorm(k * 2 * n), dim = c(k, 2, n))
Z <- landmarks_to_preshape(X)
fit <- shape_rankkmed(Z, K = 3, m = 5)
table(fit$cluster)
#> 
#>  1  2  3 
#> 10  8 12 
```
