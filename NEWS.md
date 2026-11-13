# ShapeClust 0.4.0

* `shape_mixture()`: Gaussian mixture with a spherical variance per cluster
  (mclust model VII) in tangent PCA coordinates (same scores and K - 1 floor
  as `shape_ensemble()`), with its own EM (k-means++ multiple starts,
  variance floor, empty-component restarts, flagged k-means fallback), no new
  dependency. For clusters of unequal spread, where k-means and k-means
  ensembles misassign the tail of the diffuse cluster. Returns posteriors,
  means, variances, proportions, log-likelihood and real representatives.

* `shape_ensemble()` gains `base = c("kmeans", "vii")`. The default keeps all
  earlier results identical; `"vii"` fits the mixture on each bootstrap
  resample and assigns all observations by maximum posterior.

* `shape_ensemble()` and `shape_attention()` keep at least `K - 1` principal
  components (before: at least 2). On 260 primate crania with K = 4 the 90%
  variance rule kept 2 components and split the largest species. Results with
  K <= 3 are unchanged.

* `shape_attention()`: self-attention smoothing of shapes before clustering.
  Each configuration is replaced by a softmax-weighted average of the others
  in tangent (PCA) coordinates, with per-row temperatures set to a fixed
  perplexity, as a transformer layer with no trained projections; layers can
  be stacked (blurring mean shift) and a residual weight keeps part of each
  point. The smoothed points are clustered by k-means or PAM. Reports the
  share of attention inside each point's cluster as a diagnostic.
  Experimental: in simulations it helped only with many landmarks and moderate
  sample sizes.

* `shape_ensemble()`: consensus clustering of landmark arrays. Bootstrap
  k-means in tangent PCA coordinates (90% variance rule, 2 to 10 components),
  co-association matrix, spectral partition, and the real observation closest
  to each cluster's Procrustes mean as representative. Reproduces the
  reference implementation used in the authors' 3D simulation study for the
  same RNG state.

* `hopkins_shape()`: Hopkins test of clustering tendency for 2D/3D landmark
  arrays with a Monte Carlo p-value (tangent Gaussian or uniform null).
* `shape_pam()` and `shape_rankkmed()` accept a real `k x m x n` landmark
  array (distances from `shape_distance()`). 2D inputs give the same results
  as before.
* New vignette `shapes-3d`: simulation, clustering, bagging, `hopkins_shape()`
  and the macaque skulls of 'shapes'.

## 3D shapes

The package now clusters configurations of landmarks in two or three
dimensions given as real `k x m x n` arrays. The complex 2D interface
(`landmarks_to_preshape()`, `rho_matrix()`, `dF_matrix()`, `shape_pam()`,
`shape_rankkmed()`) is unchanged, and `bagclust1()` returns the same results
as 0.3.0 for the same 2D calls (tested).

* `shape_distance()`: Riemannian, full Procrustes and partial Procrustes
  distance matrices (`"rho"`, `"dF"`, `"dP"`) in C++, from the singular values
  of the cross-product of Helmert pre-shapes. Reflections are kept unless
  `reflection = TRUE`. Agrees with `shapes::riemdist()` and
  `shapes::procdist()` (tested), and with `rho_matrix()` in 2D.
* `shape_mean()`: full Procrustes mean (`shapes::procGPA()`) and MDS/phi mean
  (`shapes::MDSshape()`).
* `shape_kmeans()`: Lloyd's algorithm with a Procrustes or phi mean, a port of
  `LloydShapes()` and `LloydShapes.phi()` of the 3D bagging scripts.
* `shape_clarans()`: FastCLARANS with R's random-number stream, defaults
  `numlocal = 2`, `maxneighbor = 0.0125` as in the scripts. No dependency on
  'fastkmedoids'. With these defaults no swap is tried when
  `n - K <= 80` (see the help page).
* `shape_hill()`: the hill-climbing search of the scripts, with the same
  random draws (tested against a transcription). As in the script every
  neighbor is a modification of the initial medoids; `from_best = TRUE` gives
  hill climbing in the usual sense.
* `bagclust1()` accepts a real landmark array as well as complex pre-shapes,
  `base = "kmeans"`, `"clarans"` and `"hill"`, `distance = "dP"` and, for
  real arrays, `reflection`. For `"kmeans"` every replicate labels all
  observations by the nearest bootstrap center.
* `sim_shapes3d()`: noisy cubes and parallelepipeds (8 or 34 landmarks,
  isotropic or anisotropic), reproducing the generators of the scripts
  draw for draw, with no dependency on 'mvtnorm' or 'CircStats'.
* 'mvtnorm' added to Suggests (a test compares the simulator with
  `mvtnorm::rmvt()`).

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
