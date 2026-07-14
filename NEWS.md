# ShapeClust 0.2.0 (2026-07-12)

## Major changes

* Full pipeline rebuilt on top of Rcpp + RcppArmadillo. Hot paths
  (Riemannian distance matrix, PAM swap-cost, rank matrix, hostility
  update) are now O(k) per pair and BLAS-backed.
* New `bagclust1()` function wrapping any base clusterer in the
  BagClust1 bootstrap ensemble of Dudoit & Fridlyand (2003), with
  Hungarian label alignment via `clue::solve_LSAP` when available.
* New `select_K()` helper implementing silhouette width and gap
  statistic in shape space.
* New `sim_complex_bingham()` using the Simplex-truncation generator of
  Kent et al. (2004).
* Full unit-test suite (40 tests, testthat).
* Package now passes `R CMD check --as-cran` with zero warnings and
  zero notes.

## New functions

* `helmert_submatrix()`, `landmarks_to_preshape()`
* `rho_matrix()`, `dF_matrix()`
* `shape_pam()`, `shape_rankkmed()`
* `bagclust1()`, `align_labels()`
* `rand_index()`, `fmi_index()`, `relative_gain()`,
  `riemannian_variances()`, `wilcoxon_paired_gain()`
* `sim_complex_bingham()`, `sim_landmark_gaussian()`
* `select_K()`

# ShapeClust 0.1.0

* Initial pure-R prototype.
