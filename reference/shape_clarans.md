# CLARANS on a shape distance matrix

Clustering Large Applications based on RANdomized Search (Ng and Han,
2002), in the faster variant of Schubert and Rousseeuw (2021): random
restarts of a randomized medoid-swap search. Own implementation in C++
that uses R's random-number stream.

## Usage

``` r
shape_clarans(
  D,
  K,
  numlocal = 2L,
  maxneighbor = 0.0125,
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

- numlocal:

  Number of random restarts (default 2).

- maxneighbor:

  Search effort of each restart. A value below 1 is a rate and means
  `r = ceiling(maxneighbor * (n - K))`; a value of 1 or more is taken as
  `r` itself. A restart ends after `r - 1` consecutive unsuccessful swap
  trials, so `r = 1` tries no swap at all. Default 0.0125.

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

  Sum of the distances of the objects to their medoid.

- n_iter:

  Number of accepted swaps in the restart that gave the returned
  medoids.

- distance:

  Distance name, or `NA` if `D` was a distance matrix.

## Details

A restart draws K distinct random medoids. It then repeatedly draws a
random non-medoid \\h\\, evaluates the change of the total cost for
swapping \\h\\ with each of the K medoids at once, and performs the best
swap if it lowers the cost; the restart ends after \\r - 1\\ consecutive
trials without a swap, with \\r\\ defined at `maxneighbor`. The best
restart is returned. The defaults are those of the 3D bagging scripts,
which called `fastkmedoids::fastclarans()` with `numlocal = 2` and
`maxneighbor = 0.0125`. This function follows the same algorithm but
draws from R's generator, so for the same
[`set.seed()`](https://rdrr.io/r/base/Random.html) it does not return
the partition of that package, and nothing beyond base R is needed.

With the default `maxneighbor`, \\r = 1\\ whenever \\n - K \le 80\\, so
no swap is ever tried and the result is the better of `numlocal` random
medoid sets; for \\n = 100\\ and \\K = 2\\ (\\r = 2\\) one failed trial
ends a restart. Raise `maxneighbor` (for instance to 0.3, or to a count
such as 100) for a real search.

## References

Ng, R. T. and Han, J. (2002). CLARANS: a method for clustering objects
for spatial data mining. *IEEE Transactions on Knowledge and Data
Engineering*, 14(5), 1003–1016.
[doi:10.1109/TKDE.2002.1033770](https://doi.org/10.1109/TKDE.2002.1033770)

Schubert, E. and Rousseeuw, P. J. (2021). Fast and eager k-medoids
clustering: O(k) runtime improvement of the PAM, CLARA, and CLARANS
algorithms. *Information Systems*, 101, 101804.
[doi:10.1016/j.is.2021.101804](https://doi.org/10.1016/j.is.2021.101804)

## See also

[`shape_pam`](https://raydonal.github.io/ShapeClust/reference/shape_pam.md),
[`shape_hill`](https://raydonal.github.io/ShapeClust/reference/shape_hill.md),
[`bagclust1`](https://raydonal.github.io/ShapeClust/reference/bagclust1.md)

## Examples

``` r
sim <- sim_shapes3d(n = 10, sigma = 0.5, seed = 1)
D <- shape_distance(sim$X)
set.seed(1)
fit <- shape_clarans(D, K = 2)
table(fit$cluster, sim$group)
#>    
#>     1 2
#>   1 7 1
#>   2 3 9
```
