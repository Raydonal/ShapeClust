# Hopkins test of clustering tendency for landmark shapes

Hopkins statistic with a Monte Carlo reference distribution, for a
sample of shapes given as a real \\k \times d \times n\\ landmark array
(\\d = 2\\ or \\3\\).

## Usage

``` r
hopkins_shape(
  X,
  null = c("tangent", "uniform"),
  m = NULL,
  power = NULL,
  B = 199,
  distance = "rho",
  seed = NULL
)
```

## Arguments

- X:

  Numeric array \\k \times d \times n\\, \\d \in \\2, 3\\\\.

- null:

  `"tangent"` (default) or `"uniform"`; see Details.

- m:

  Number of artificial and of sampled points. Default `round(0.1 * n)`,
  at least 1.

- power:

  Exponent \\p\\; default 1.

- B:

  Number of Monte Carlo reference draws.

- distance:

  Shape distance (`"rho"`, `"dF"`, `"dP"`) for `null = "uniform"`. The
  tangent null uses Euclidean distance in the aligned coordinates.

- seed:

  Optional integer seed; the caller's random-number stream is restored
  on exit.

## Value

An object of class `"hopkins_shape"`: a list with `statistic`,
`p.value`, `null`, `power`, `m`, `B`, and `reference` (the `B` reference
values).

## Details

The statistic is \\H = \sum u^p / (\sum u^p + \sum w^p)\\, where \\u\\
are the nearest-neighbor distances from \\m\\ artificial points to the
data and \\w\\ those from \\m\\ sampled data points to the other data.
Large values indicate clustering.

With `null = "tangent"`, the artificial points come from a unimodal
Gaussian with the mean and covariance of the Procrustes-aligned data
coordinates (centered, unit-size configurations rotated onto the
Procrustes mean), and distances are Euclidean in those coordinates. This
tests unimodality against clustering. With `null = "uniform"`, the
artificial points are uniform on shape space (Gaussian landmarks) and
distances are `distance`. **The uniform null rejects for any
concentrated sample, clustered or not, so it says nothing about
clusters.**

The p-value does not assume \\H = 0.5\\. For each of `B` draws a full
data set of \\n\\ objects is simulated from the null (for the tangent
null, the Gaussian fitted to the data; for the uniform null, \\n\\
uniform shapes), \\H\\ is computed on it by the same procedure, and \\p
= (1 + \\\\H_b \ge H\\) / (B + 1)\\.

The exponent usually recommended is the dimension of the space. With a
Monte Carlo reference the choice affects power, not validity.

## See also

[`shape_distance`](https://raydonal.github.io/ShapeClust/reference/shape_distance.md),
[`shape_mean`](https://raydonal.github.io/ShapeClust/reference/shape_mean.md)

## Examples

``` r
sim <- sim_shapes3d(n = 30, sigma = 1, k = 8, seed = 1)
hopkins_shape(sim$X, B = 19, seed = 2)
#> Hopkins test of clustering tendency (tangent null)
#> H = 0.6938, p = 0.05  (m = 6, power = 1, B = 19)
```
