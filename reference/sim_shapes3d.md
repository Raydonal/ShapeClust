# Simulate noisy cubes and parallelepipeds in 3D

Simulator of the 3D clustering study: each group is a sample of noisy
copies of a template configuration (a cube, or a box with one edge twice
as long) with 8 or 34 landmarks, with isotropic or anisotropic
perturbation, rotated by a common random angle about the z axis.

## Usage

``` r
sim_shapes3d(
  n,
  sigma = 3,
  type = c("iso", "aniso"),
  k = 8L,
  groups = c("cube", "parallelepiped"),
  templates = NULL,
  gamma = 4,
  df = 99,
  rotate = TRUE,
  independent = FALSE,
  seed = NULL
)
```

## Arguments

- n:

  Number of configurations per group: a single number, or a vector of
  length `K`.

- sigma:

  Dispersion (default 3). The covariance scale is `sigma^2`.

- type:

  `"iso"` (independent Gaussian-like noise on every coordinate, default)
  or `"aniso"` (landmark and coordinate covariance of Kronecker form,
  see Details).

- k:

  Number of landmarks of the built-in templates, 8 (default) or 34.
  Ignored if `templates` is given.

- groups:

  Character vector naming the template of each group, from `"cube"` and
  `"parallelepiped"` (default: both, so K = 2). Names may be repeated or
  reordered.

- templates:

  Optional list of \\k \times 3\\ matrices used instead of the built-in
  ones, one per group.

- gamma:

  Anisotropy parameter of `type = "aniso"` (default 4, as in the
  scripts).

- df:

  Degrees of freedom of the multivariate t noise (default 99, as in the
  scripts).

- rotate:

  Logical. If `TRUE` (default) each group is rotated by one random angle
  about the z axis, as in the scripts. Shape distances do not depend on
  it.

- independent:

  Logical. The scripts reseed the generator with the same seed before
  each group, so that all groups share the same noise realization and
  the same rotation angle. This is reproduced when `independent = FALSE`
  (default). With `TRUE` the seed is set once and the groups get
  independent noise.

- seed:

  Optional integer seed; the caller's random-number stream is left as it
  was. In the scripts, replicate number `c` used `semente = c`. If
  `NULL` and `independent = FALSE` one seed is drawn from the current
  stream.

## Value

A list with components:

- X:

  Numeric array \\k \times 3 \times N\\ of the rotated configurations,
  groups stacked in order.

- group:

  Integer vector of length \\N\\: the true group (1 to K).

- X_norot:

  The configurations before rotation.

- K:

  Number of groups.

- type, k, sigma, groups:

  The options used.

## Details

For group \\g\\ with template \\M_g\\ (a \\k \times 3\\ matrix),
configurations are \\M_g + E\\, where \\\mathrm{vec}(E)\\ is drawn from
a multivariate t with `df` degrees of freedom and scale matrix
\\\Sigma\\: \\\Sigma = \sigma^2 I\_{3k}\\ for `"iso"`, and \\\Sigma =
\sigma^2 (J_k + (\gamma - 1) I_k) \otimes (J_3 + (\gamma - 1)
I_3)/\gamma^2\\ for `"aniso"`, with \\J\\ a matrix of ones, the
landmark-major ordering \\(x_1, y_1, z_1, x_2, \dots)\\. The random
draws follow the scripts: first a concentration (one of `seq(-pi, pi)`)
and a mean angle for the von Mises draw, then the noise, then the angle.
The noise reproduces
[`mvtnorm::rmvt()`](https://rdrr.io/pkg/mvtnorm/man/Mvt.html) (version
0.9-9994 or later, whose random numbers fill the matrix by rows) and the
angle reproduces `CircStats::rvm()`, without depending on those
packages. Note that some concentrations in `seq(-pi, pi)` are negative
and the mean angle is a large integer; this is kept as it is in the
scripts.

Built-in 8-landmark templates are the vertices of a cube \\\[0, 10\]^3\\
and of the box \\\[0, 10\] \times \[0, 20\] \times \[0, 10\]\\ with a
vertex labelling that is the mirror image of the cube's (reflections
matter, see
[`shape_distance`](https://raydonal.github.io/ShapeClust/reference/shape_distance.md));
the 34-landmark templates add edge and face points. They are the
coordinates of `cube8landm`, `parallelep8landm`, `cube34landm` and
`parallelep34landm` in the 'Anthropometry' package.

## See also

[`shape_kmeans`](https://raydonal.github.io/ShapeClust/reference/shape_kmeans.md),
[`shape_distance`](https://raydonal.github.io/ShapeClust/reference/shape_distance.md),
[`sim_landmark_gaussian`](https://raydonal.github.io/ShapeClust/reference/sim_landmark_gaussian.md)
for 2D data.

## Examples

``` r
sim <- sim_shapes3d(n = 5, sigma = 1, type = "aniso", k = 8, seed = 7)
dim(sim$X)
#> [1]  8  3 10
sim$group
#>  [1] 1 1 1 1 1 2 2 2 2 2
# Reproducible, and the caller's random-number stream is not touched
identical(sim$X, sim_shapes3d(n = 5, sigma = 1, type = "aniso", seed = 7)$X)
#> [1] TRUE
```
