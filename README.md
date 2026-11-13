# ShapeClust

<!-- badges: start -->
[![License: GPL-3](https://img.shields.io/badge/License-GPL_3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

Medoid-based partitional clustering on **Kendall's two-dimensional shape
space**, together with the **BagClust1** bootstrap ensemble that stabilizes
partitions by majority voting over resampled clusterings.

Once translation, scale and rotation are removed, planar landmark
configurations lie on a curved manifold rather than in a flat Euclidean space,
so clustering algorithms that assume straight-line distances misjudge how close
two shapes really are. `ShapeClust` implements Partitioning Around Medoids and
Ranked *K*-medoids with the intrinsic Riemannian geodesic distance of that
manifold, and wraps either of them in a bagging ensemble. The computationally
intensive routines are written in C++ via `Rcpp` and `RcppArmadillo`.

This repository holds the R package and the scripts that reproduce the study in
the accompanying manuscript. The manuscript itself is not in this repository.

## Installation

```r
# install.packages("remotes")
remotes::install_github("raydonal/ShapeClust")
```

The package is not on CRAN.

## Quick start

```r
library(ShapeClust)

# Three shape groups -- a square, a rectangle and a trapezoid -- perturbed
mu <- list(
  matrix(c(0,0,  0,10, 10,10, 10,0),  4, 2, byrow = TRUE),
  matrix(c(0,0, 10,0,   0,20, 10,20), 4, 2, byrow = TRUE),
  matrix(c(0,0, 40,0,   0,30, 20,30), 4, 2, byrow = TRUE)
)
sim <- sim_landmark_gaussian(mu, n_per_group = c(50, 50, 50),
                             sigma = 9, regime = "isotropic", seed = 1)

# Landmark configurations -> complex pre-shapes on the manifold
Z <- landmarks_to_preshape(sim$X)

# The base clusterer, and the same clusterer inside the bagging ensemble
fit_base <- shape_pam(Z, K = 3)
fit_bag  <- bagclust1(Z, K = 3, base = "pam", B = 100, seed = 1)

fmi_index(sim$group, fit_base$cluster)   # base
fmi_index(sim$group, fit_bag$cluster)    # bagged
```

Not sure how many groups there are? `select_K()` computes the average silhouette
width and the gap statistic in the Riemannian metric. Note that bagging does not
estimate *K*: given a fixed *K*, it stabilizes the partition that the base
clusterer produces.

## 3D shapes

Configurations in two or three dimensions can also be given as a real
`k x m x n` array (m = 2 or 3). `shape_pam()`, `shape_rankkmed()` and
`bagclust1()` accept it directly, and `shape_kmeans()` clusters around
Procrustes mean shapes.

```r
# Two groups of 3D configurations (a cube and a parallelepiped), n per group
sim <- sim_shapes3d(n = 20, sigma = 0.5, seed = 1)

hopkins_shape(sim$X, B = 49, seed = 1)    # clustering tendency
D <- shape_distance(sim$X, "rho")         # Riemannian distance matrix

fit_km  <- shape_kmeans(sim$X, K = 2)
fit_pam <- shape_pam(sim$X, K = 2)
fit_bag <- bagclust1(sim$X, K = 2, base = "kmeans", B = 100, seed = 1)

fmi_index(sim$group, fit_bag$cluster)
```

The `shapes-3d` vignette walks through this and a real-data example (the
macaque skulls of the **shapes** package).

## What is in here

```
ShapeClust/
├── R/                  # R layer: geometry, algorithms, ensemble, validation,
│                       #   simulation, model selection
├── src/                # C++ kernels (Rcpp + RcppArmadillo)
├── man/                # documentation
├── tests/testthat/
├── vignettes/          # getting-started and 3D vignettes
└── inst/               # scripts that reproduce the whole study
```

### The package

| Function | Purpose |
|---|---|
| `landmarks_to_preshape()`, `helmert_submatrix()` | landmark configurations to complex pre-shapes |
| `rho_matrix()`, `dF_matrix()` | Riemannian and full Procrustes distance matrices |
| `shape_pam()`, `shape_rankkmed()` | the two base clusterers on the shape manifold |
| `bagclust1()`, `align_labels()` | the bagging ensemble and its label alignment |
| `rand_index()`, `fmi_index()`, `relative_gain()`, `riemannian_variances()`, `wilcoxon_paired_gain()` | validation |
| `sim_complex_bingham()`, `sim_landmark_gaussian()` | the data-generating regimes |
| `select_K()` | choosing the number of clusters |
| `shape_distance()`, `shape_mean()`, `shape_kmeans()`, `shape_clarans()`, `shape_hill()` | 3D shapes: distances, mean shapes and further clusterers on real `k x m x n` arrays |
| `shape_ensemble()` | consensus clustering of 3D shapes: bootstrap k-means in tangent PCA coordinates, co-association, spectral partition, real representative per cluster |
| `sim_shapes3d()`, `hopkins_shape()` | 3D simulator and Hopkins test of clustering tendency |

## Reproducing the study

The full design is 24 scenarios x 100 Monte Carlo replications x 100 bootstrap
replicates, for two base clusterers. On a single core that is roughly 70
minutes, so `inst/driver.R` runs it in resumable chunks that cache to disk and
skip whatever is already done. The scripts write to `paper/` under the working
directory, so create that directory first:

```bash
mkdir -p paper/figures
Rscript inst/driver.R        # repeat until it reports 96 cached chunks
Rscript inst/assemble.R      # 48 violin plots + the summary table
Rscript inst/finish_figs.R   # 24 biplots
Rscript inst/missing-figures.R
Rscript inst/real_data_mc.R  # real data, 100 random initializations
Rscript inst/make_tables.R   # regenerates the LaTeX tables
```

## Citation

```r
citation("ShapeClust")
```

> Ospina, R., Nascimento, I. R. A. and Amaral, G. J. A. (2026). Bagged
> medoid-based clustering on Kendall's two-dimensional shape space. Preprint.

## Authors

- **Raydonal Ospina** (maintainer) — Departamento de Estatística, LInCa,
  Universidade Federal da Bahia; Departamento de Estatística, CASTLab,
  Universidade Federal de Pernambuco · <raydonal@de.ufpe.br>
- **Inácio Robson Alves do Nascimento** — Departamento de Estatística, CASTLab,
  Universidade Federal de Pernambuco
- **Getúlio José Amorim do Amaral** — Departamento de Estatística, CASTLab,
  Universidade Federal de Pernambuco

## License

GPL-3. See the [license text](https://www.gnu.org/licenses/gpl-3.0).
