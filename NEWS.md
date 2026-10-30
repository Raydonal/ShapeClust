# ShapeClust 0.3.0

## Portability and checks

* Added `src/Makevars` and `src/Makevars.win`: complex matrix products in the
  C++ layer need BLAS at link time.
* Added `tests/testthat.R`, so `R CMD check` runs the test suite.
* `bagclust1()` now handles `parallel::detectCores()` returning `NA` and stops
  with the original message when a parallel replicate fails. The comparison of
  serial and parallel runs is now tested on Windows too.

## Bug fixes (change numerical results)

* `shape_rankkmed()` now follows Zadegan et al. (2013) as implemented in
  `kmed::rankkmed()`: hostility is computed over the m nearest objects of
  each medoid, and colliding medoids are relocated at random to objects
  outside every group. Previously the hostility was summed over the whole
  assigned cluster, where it is constant whenever the cluster consists of
  mutual nearest neighbours, so the update picked the lowest-indexed member;
  and no relocation step existed. Partitions now equal those of
  `kmed::rankkmed()` for the same distances and random stream (tested).
* `bagclust1()` now lets every replicate label all n observations by their
  nearest bootstrap medoid, as in Dudoit and Fridlyand (2003). Previously
  only the in-bag observations (about 63%) voted in each replicate.
* `select_K(method = "gap")` drew its null reference by permuting the rows
  and columns of the distance matrix, which leaves the clustering problem
  unchanged and made the gap close to zero for every K. The reference is now
  a uniform sample on the shape space, and K is chosen by the one-standard-
  error rule of Tibshirani et al. (2001). The gap needs pre-shapes, not
  distances.
* `riemannian_variances()` now uses full Procrustes means (Kent 1994) instead
  of the sample medoid, and no longer claims to use 'shapes'.

## Other changes

* `bagclust1()` draws one seed per replicate, so results are identical with
  and without `parallel = TRUE`.
* Functions with a `seed` argument restore the caller's random-number stream.
* Removed the no-op arguments `center` (`landmarks_to_preshape()`) and
  `verbose` (`shape_pam()`).
* Examples for every exported function; regression tests for each fix.

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
