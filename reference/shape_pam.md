# Shape-PAM: PAM on Kendall's 2D shape space

Shape-PAM: PAM on Kendall's 2D shape space

## Usage

``` r
shape_pam(Z, K, distance = c("rho", "dF"), initial = NULL, max_iter = 100L)
```

## Arguments

- Z:

  Complex pre-shape matrix of size (k-1) x n produced by
  [`landmarks_to_preshape`](https://raydonal.github.io/ShapeClust/reference/landmarks_to_preshape.md),
  OR a numeric distance matrix of size n x n, OR a real landmark array
  of dimension k x m x n with m = 2 or 3 (distances from
  [`shape_distance`](https://raydonal.github.io/ShapeClust/reference/shape_distance.md)).
  If a distance matrix is supplied, it is used directly.

- K:

  Number of clusters (\>= 2).

- distance:

  Character. One of "rho" (Riemannian; default), "dF" (full Procrustes).
  Ignored if `Z` is already a distance matrix.

- initial:

  Optional vector of length K giving initial medoid indices (1-based).
  If NULL, medoids are sampled uniformly at random from 1:n.

- max_iter:

  Integer. Maximum number of swap iterations (default 100).

## Value

A list with components:

- cluster:

  Integer vector of length n, values in 1:K.

- medoids:

  Integer vector of length K.

- objective:

  Sum of distances of each object to its medoid.

- n_iter:

  Number of accepted swaps performed.

- distance:

  Distance name used.

## Examples

``` r
set.seed(1)
k <- 4; n <- 30
X <- array(rnorm(k * 2 * n), dim = c(k, 2, n))
Z <- landmarks_to_preshape(X)
fit <- shape_pam(Z, K = 3)
table(fit$cluster)
#> 
#>  1  2  3 
#>  8  8 14 
```
