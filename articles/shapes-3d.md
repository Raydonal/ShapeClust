# Clustering 3D shapes with ShapeClust

## Overview

Landmark configurations in two or three dimensions are stored as a real
`k × m × n` array: `k` landmarks, `m = 2` or `3` coordinates, `n`
configurations. The functions below work on such arrays directly. The
complex 2D interface described in the getting-started vignette is
unchanged.

## Simulating 3D shapes

[`sim_shapes3d()`](https://raydonal.github.io/ShapeClust/reference/sim_shapes3d.md)
draws two groups of configurations around a cube and a parallelepiped
template with Gaussian landmark noise, after a random rotation of each
configuration. Here `n` is the number of configurations per group.

``` r
library(ShapeClust)

sim <- sim_shapes3d(n = 20, sigma = 0.5, seed = 1)
dim(sim$X)
#> [1]  8  3 40
table(sim$group)
#> 
#>  1  2 
#> 20 20
```

## Distances

[`shape_distance()`](https://raydonal.github.io/ShapeClust/reference/shape_distance.md)
returns the Riemannian (`"rho"`), full Procrustes (`"dF"`) or partial
Procrustes (`"dP"`) distance matrix. Reflections are kept as distinct
shapes unless `reflection = TRUE`.

``` r
D <- shape_distance(sim$X, "rho")
dim(D)
#> [1] 40 40
round(range(D[upper.tri(D)]), 3)
#> [1] 0.044 1.107
```

## Clustering

[`shape_kmeans()`](https://raydonal.github.io/ShapeClust/reference/shape_kmeans.md)
runs Lloyd’s algorithm with a Procrustes mean shape.
[`shape_pam()`](https://raydonal.github.io/ShapeClust/reference/shape_pam.md)
and
[`shape_rankkmed()`](https://raydonal.github.io/ShapeClust/reference/shape_rankkmed.md)
accept the array as well (they compute
[`shape_distance()`](https://raydonal.github.io/ShapeClust/reference/shape_distance.md)
internally) or a distance matrix.

``` r
set.seed(1)
fit_km  <- shape_kmeans(sim$X, K = 2)
fit_pam <- shape_pam(sim$X, K = 2)

table(group = sim$group, kmeans = fit_km$cluster)
#>      kmeans
#> group  1  2
#>     1 20  0
#>     2  0 20
table(group = sim$group, pam = fit_pam$cluster)
#>      pam
#> group  1  2
#>     1 20  0
#>     2  0 20
cat("FMI, k-means:", round(fmi_index(sim$group, fit_km$cluster), 3), "\n")
#> FMI, k-means: 1
cat("FMI, PAM    :", round(fmi_index(sim$group, fit_pam$cluster), 3), "\n")
#> FMI, PAM    : 1
```

## Bagging

[`bagclust1()`](https://raydonal.github.io/ShapeClust/reference/bagclust1.md)
takes the same array. With `base = "kmeans"`, each bootstrap replicate
fits the centers on the resample and labels all `n` configurations by
the nearest center; the medoid-based bases (`"pam"`, `"rankkmed"`,
`"clarans"`, `"hill"`) work on the distance matrix.

``` r
fit_bag <- bagclust1(sim$X, K = 2, base = "pam", B = 30, seed = 1)
table(group = sim$group, bagged = fit_bag$cluster)
#>      bagged
#> group  1  2
#>     1 20  0
#>     2  0 20
cat("FMI, bagged PAM:", round(fmi_index(sim$group, fit_bag$cluster), 3), "\n")
#> FMI, bagged PAM: 1

fit_bagk <- bagclust1(sim$X, K = 2, base = "kmeans", B = 10, seed = 1)
cat("FMI, bagged k-means:",
    round(fmi_index(sim$group, fit_bagk$cluster), 3), "\n")
#> FMI, bagged k-means: 1
```

## Clustering tendency

[`hopkins_shape()`](https://raydonal.github.io/ShapeClust/reference/hopkins_shape.md)
tests whether a sample has cluster structure, with a Monte Carlo p-value
under a tangent-space Gaussian null.

``` r
hopkins_shape(sim$X, B = 49, seed = 1)
#> Hopkins test of clustering tendency (tangent null)
#> H = 0.7347, p = 0.02  (m = 4, power = 1, B = 49)
```

## A real data example: macaque skulls

The `macaques` data of the **shapes** package are 7 landmarks in 3D on
18 skulls, nine female and nine male.

``` r
data(macaques, package = "shapes")
X <- macaques$x
sex <- macaques$group
dim(X)
#> [1]  7  3 18
table(sex)
#> sex
#> f m 
#> 9 9

set.seed(1)
fit_pam <- shape_pam(X, K = 2)
fit_bag <- bagclust1(X, K = 2, base = "pam", B = 50, seed = 1)

table(sex = sex, pam = fit_pam$cluster)
#>    pam
#> sex 1 2
#>   f 1 8
#>   m 6 3
table(sex = sex, bagged = fit_bag$cluster)
#>    bagged
#> sex 1 2
#>   f 1 8
#>   m 6 3
cat("FMI, PAM       :", round(fmi_index(as.integer(sex), fit_pam$cluster), 3), "\n")
#> FMI, PAM       : 0.622
cat("FMI, bagged PAM:", round(fmi_index(as.integer(sex), fit_bag$cluster), 3), "\n")
#> FMI, bagged PAM: 0.622
```

## Consensus clustering with a representative shape

[`shape_ensemble()`](https://raydonal.github.io/ShapeClust/reference/shape_ensemble.md)
clusters in the leading principal components of the tangent coordinates
(enough components for 90% of the variance, between 2 and 10). It runs
k-means on bootstrap samples, builds the co-association matrix of the
replicates and cuts it with a spectral partition. Each cluster also gets
a representative, the observed shape closest to the cluster’s Procrustes
mean.

``` r
fit_ens <- shape_ensemble(sim$X, K = 2, B = 30, seed = 1)
fit_ens
#> Shape ensemble (tangent PCA, co-association, spectral)
#>   n = 40  K = 2  B = 30  PCA dimension = 2 
#>   cluster sizes: 20 20 
#>   representatives: 15 35
table(truth = sim$group, ensemble = fit_ens$cluster)
#>      ensemble
#> truth  1  2
#>     1 20  0
#>     2  0 20
```
