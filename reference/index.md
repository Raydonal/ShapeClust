# Package index

## Geometry

From landmark configurations to complex pre-shapes and shape-space
distances.

- [`helmert_submatrix()`](https://raydonal.github.io/ShapeClust/reference/helmert_submatrix.md)
  : Helmert sub-matrix
- [`landmarks_to_preshape()`](https://raydonal.github.io/ShapeClust/reference/landmarks_to_preshape.md)
  : Convert 2D landmark configurations to complex pre-shape vectors
- [`rho_matrix()`](https://raydonal.github.io/ShapeClust/reference/rho_matrix.md)
  : Riemannian geodesic distance matrix on Kendall's 2D shape space
- [`dF_matrix()`](https://raydonal.github.io/ShapeClust/reference/dF_matrix.md)
  : Full Procrustes distance matrix

## Geometry for 3D shapes

Distances and mean shapes for real k x m x n landmark arrays (m = 2 or
3).

- [`shape_distance()`](https://raydonal.github.io/ShapeClust/reference/shape_distance.md)
  : Shape distances for real landmark arrays
- [`shape_mean()`](https://raydonal.github.io/ShapeClust/reference/shape_mean.md)
  : Mean shape of a sample of 3D or 2D landmark configurations

## Clustering algorithms

PAM and Ranked K-medoids adapted to Kendall’s 2D shape space.

- [`shape_pam()`](https://raydonal.github.io/ShapeClust/reference/shape_pam.md)
  : Shape-PAM: PAM on Kendall's 2D shape space
- [`shape_rankkmed()`](https://raydonal.github.io/ShapeClust/reference/shape_rankkmed.md)
  : Shape-rankkmed: Ranked K-medoids on Kendall's 2D shape space

## Clustering algorithms for 3D shapes

k-means, CLARANS and hill climbing for real landmark arrays or distance
matrices.

- [`shape_kmeans()`](https://raydonal.github.io/ShapeClust/reference/shape_kmeans.md)
  : Shape k-means for real landmark arrays
- [`shape_clarans()`](https://raydonal.github.io/ShapeClust/reference/shape_clarans.md)
  : CLARANS on a shape distance matrix
- [`shape_hill()`](https://raydonal.github.io/ShapeClust/reference/shape_hill.md)
  : Hill-climbing medoid search on a shape distance matrix

## Ensembles and mixtures

Consensus ensemble in tangent PCA, Gaussian mixture with one variance
per cluster, and experimental self-attention smoothing.

- [`shape_ensemble()`](https://raydonal.github.io/ShapeClust/reference/shape_ensemble.md)
  : Consensus clustering of shapes in tangent PCA coordinates
- [`shape_mixture()`](https://raydonal.github.io/ShapeClust/reference/shape_mixture.md)
  : Gaussian mixture clustering of shapes with a variance per cluster
- [`shape_attention()`](https://raydonal.github.io/ShapeClust/reference/shape_attention.md)
  : Self-attention smoothing and clustering of shapes

## Bagging ensemble

BagClust1 wrapper and label alignment.

- [`bagclust1()`](https://raydonal.github.io/ShapeClust/reference/bagclust1.md)
  : BagClust1 for shape clustering
- [`align_labels()`](https://raydonal.github.io/ShapeClust/reference/align_labels.md)
  : Align cluster labels via optimal permutation

## Validation

External validation indices and paired Wilcoxon tests.

- [`rand_index()`](https://raydonal.github.io/ShapeClust/reference/rand_index.md)
  : Rand Index
- [`fmi_index()`](https://raydonal.github.io/ShapeClust/reference/fmi_index.md)
  : Fowlkes-Mallows Index
- [`relative_gain()`](https://raydonal.github.io/ShapeClust/reference/relative_gain.md)
  : Relative Gain
- [`riemannian_variances()`](https://raydonal.github.io/ShapeClust/reference/riemannian_variances.md)
  : Riemannian between- and within-group variances (VEG, VIG)
- [`wilcoxon_paired_gain()`](https://raydonal.github.io/ShapeClust/reference/wilcoxon_paired_gain.md)
  : Paired one-sided Wilcoxon signed-rank test

## Simulation

- [`sim_complex_bingham()`](https://raydonal.github.io/ShapeClust/reference/sim_complex_bingham.md)
  : Simulate a single Complex Bingham sample
- [`sim_landmark_gaussian()`](https://raydonal.github.io/ShapeClust/reference/sim_landmark_gaussian.md)
  : Simulate landmark configurations from a Gaussian perturbation
- [`sim_shapes3d()`](https://raydonal.github.io/ShapeClust/reference/sim_shapes3d.md)
  : Simulate noisy cubes and parallelepipeds in 3D

## Model selection and clustering tendency

- [`select_K()`](https://raydonal.github.io/ShapeClust/reference/select_K.md)
  : Select K by silhouette width or the gap statistic
- [`hopkins_shape()`](https://raydonal.github.io/ShapeClust/reference/hopkins_shape.md)
  : Hopkins test of clustering tendency for landmark shapes
