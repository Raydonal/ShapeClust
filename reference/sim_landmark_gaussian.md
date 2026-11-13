# Simulate landmark configurations from a Gaussian perturbation

Perturbs a template configuration by a multivariate Gaussian in landmark
space, following Amaral et al. (2007) and Preston & Wood (2010).

## Usage

``` r
sim_landmark_gaussian(
  mu_list,
  n_per_group,
  sigma = 1,
  regime = c("isotropic", "anisotropic"),
  gamma = 4,
  seed = NULL
)
```

## Arguments

- mu_list:

  List of template configurations, each a k x 2 matrix. Length = number
  of groups G.

- n_per_group:

  Integer vector of length G, sample sizes per group.

- sigma:

  Standard deviation of the perturbation.

- regime:

  One of "isotropic" (default) or "anisotropic".

- gamma:

  Anisotropy parameter (default 4).

- seed:

  Optional integer seed; the caller's random-number stream is left as it
  was.

## Value

A list with:

- X:

  Numeric array of shape (k, 2, n_total).

- group:

  Integer vector of group labels in 1..G.

## Examples

``` r
mu_list <- list(
  matrix(c(0,0, 10,0, 10,10, 0,10), 4, 2, byrow = TRUE),
  matrix(c(0,0, 20,0, 20,10, 0,10), 4, 2, byrow = TRUE),
  matrix(c(0,0, 15,0, 10,10, 0,10), 4, 2, byrow = TRUE)
)
out <- sim_landmark_gaussian(mu_list, c(10, 10, 10), sigma = 1)
dim(out$X); table(out$group)
#> [1]  4  2 30
#> 
#>  1  2  3 
#> 10 10 10 
```
