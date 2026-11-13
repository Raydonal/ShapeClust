# Self-attention smoothing and clustering of shapes

Replaces each configuration by an attention-weighted average of the
others in tangent coordinates, then clusters the smoothed points. The
attention is the scaled dot-product attention of transformers with no
trained projections: queries, keys and values are the tangent
coordinates themselves, and the weights are a softmax of minus the
squared distances.

## Usage

``` r
shape_attention(
  X,
  K,
  perplexity = NULL,
  layers = 1L,
  alpha = 0.5,
  reduce = c("pca", "none"),
  var_explained = 0.9,
  max_dim = 10L,
  base = c("kmeans", "pam"),
  nstart = 10L,
  seed = NULL
)
```

## Arguments

- X:

  Numeric array of dimension \\k \times m \times n\\, with \\m\\ equal
  to 2 or 3.

- K:

  Number of clusters (\>= 2).

- perplexity:

  Effective number of neighbors each observation attends to. Each row's
  temperature is set by bisection so that \\\exp(-\sum_j a\_{ij}\log
  a\_{ij})\\ equals this value. Default
  `min(30, max(5, floor(n / (2 K))))`.

- layers:

  Number of attention layers (default 1). Weights are recomputed from
  the smoothed points at each layer.

- alpha:

  Residual weight in \\(0, 1\]\\: each layer returns \\(1-\alpha) x_i +
  \alpha \sum_j a\_{ij} x_j\\. Default 0.5.

- reduce:

  `"pca"` (default) works in the leading principal components of the
  tangent coordinates, chosen as in
  [`shape_ensemble`](https://raydonal.github.io/ShapeClust/reference/shape_ensemble.md);
  `"none"` uses all tangent coordinates.

- var_explained, max_dim:

  PCA rule, as in
  [`shape_ensemble`](https://raydonal.github.io/ShapeClust/reference/shape_ensemble.md).

- base:

  Clusterer applied to the smoothed points: `"kmeans"` (default,
  [`kmeans`](https://rdrr.io/r/stats/kmeans.html) with `nstart` starts)
  or `"pam"`
  ([`shape_pam`](https://raydonal.github.io/ShapeClust/reference/shape_pam.md)
  on their Euclidean distances).

- nstart:

  Number of k-means starts (default 10).

- seed:

  Optional integer seed. The caller's random-number stream is restored
  on exit.

## Value

An object of class `"shape_attention"`, a list with

- cluster:

  Integer vector of length n, values in 1:K.

- smoothed:

  n x d matrix of smoothed coordinates.

- attention:

  n x n attention matrix of the first layer.

- beta:

  Inverse temperature of each row in the first layer.

- self_mass:

  Mean first-layer attention inside the own cluster.

- dim:

  Number of coordinates used.

- perplexity, layers, alpha, base, K:

  The options used.

## Details

The attention matrix has zero diagonal, so each point is replaced by an
average of other points; with `alpha < 1` part of the point is kept, as
in the residual connection of a transformer. Stacking layers with
recomputed weights is the blurring mean shift of Cheng (1995), and it is
also the self-attention dynamics studied by Geshkovski et al. (2023),
under which points merge into clusters.

The smoothing helps when the attention mass of a point falls mostly on
its own group. The returned `self_mass` is the average share of
attention that falls inside the point's final cluster, a label-free
proxy for that condition. Values near \\1/K\\ mean the attention did not
separate the groups.

Under a mixture of independent observations no function of the whole
sample classifies a point better than the Bayes rule applied to that
point alone, so the smoothing cannot remove overlap between groups; it
can only reduce the loss due to estimating centers and searching
partitions. In Gaussian simulations it did not beat k-means with 10
starts in any setting tried. The function is experimental.

## References

Vaswani, A. et al. (2017). Attention is all you need. NeurIPS 30.

Cheng, Y. (1995). Mean shift, mode seeking, and clustering. IEEE TPAMI
17(8), 790–799.

Geshkovski, B., Letrouit, C., Polyanskiy, Y. and Rigollet, P. (2023).
The emergence of clusters in self-attention dynamics. NeurIPS 36.

## See also

[`shape_ensemble`](https://raydonal.github.io/ShapeClust/reference/shape_ensemble.md),
[`shape_kmeans`](https://raydonal.github.io/ShapeClust/reference/shape_kmeans.md)

## Examples

``` r
sim <- sim_shapes3d(n = 60, sigma = 1, k = 8, seed = 1)
fit <- shape_attention(sim$X, K = 2, seed = 1)
fit
#> Shape attention clustering
#>   n = 120  K = 2  dim = 2  perplexity = 30  layers = 1  alpha = 0.5  base = kmeans 
#>   cluster sizes: 60 60 
#>   attention mass inside own cluster: 1.000 (1/K = 0.500)
table(sim$group, fit$cluster)
#>    
#>      1  2
#>   1  0 60
#>   2 60  0
```
