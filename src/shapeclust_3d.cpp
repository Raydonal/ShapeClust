// -----------------------------------------------------------------------------
//  ShapeClust: real-array (m = 2 or 3) shape geometry and medoid searches.
//
//  Landmark configurations are k x m real matrices. They are mapped to
//  pre-shapes with the Helmert sub-matrix (translation) and unit Frobenius
//  norm (scale), giving (k-1) x m matrices. For two pre-shapes z and w the
//  singular values s_1 >= ... >= s_m of z'w give the full Procrustes
//  alignment. Reflections are excluded unless asked for: when det(z'w) < 0
//  the smallest singular value enters the sum with a negative sign
//  (Dryden and Mardia 2016, section 4.1.2).
// -----------------------------------------------------------------------------

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <R_ext/Random.h>
#include <cmath>
#include <algorithm>
#include <limits>
#include <vector>

using namespace Rcpp;

namespace {

// (k-1) x k Helmert sub-matrix, same entries as helmert_submatrix() in R.
arma::mat helmert_sub(arma::uword k) {
  arma::mat H(k - 1, k, arma::fill::zeros);
  for (arma::uword j = 2; j <= k; ++j) {
    const double denom = std::sqrt(static_cast<double>(j) * (j - 1));
    for (arma::uword c = 0; c < j - 1; ++c) H(j - 2, c) = -1.0 / denom;
    H(j - 2, j - 1) = (j - 1) / denom;
  }
  return H;
}

// Pre-shapes of all slices of X: (k-1) x m matrices with unit norm.
std::vector<arma::mat> preshapes(const arma::cube& X, const arma::mat& H) {
  std::vector<arma::mat> Z(X.n_slices);
  for (arma::uword i = 0; i < X.n_slices; ++i) {
    arma::mat z = H * X.slice(i);
    const double nz = arma::norm(z, "fro");
    if (nz > 0.0) z /= nz;
    Z[i] = z;
  }
  return Z;
}

// Sum of the singular values of the m x m matrix A, with the sign
// correction on the smallest one unless reflections are allowed.
inline double sv_sum(const arma::mat& A, bool reflection) {
  const arma::vec s = arma::svd(A);
  double tot = arma::accu(s);
  if (!reflection && arma::det(A) < 0.0) tot -= 2.0 * s(s.n_elem - 1);
  return tot;
}

inline double clamp01(double c) {
  if (c > 1.0) return 1.0;
  if (c < 0.0) return 0.0;
  return c;
}

// type: 0 = rho, 1 = dF, 2 = dP
inline double dist_from_cos(double c, int type) {
  c = clamp01(c);
  switch (type) {
    case 0: return std::acos(c);
    case 1: return std::sqrt(1.0 - c * c);
    default: return std::sqrt(2.0 * (1.0 - c));
  }
}

inline bool identical_config(const arma::cube& X, arma::uword i,
                             const arma::cube& Y, arma::uword j) {
  const double* a = X.slice_memptr(i);
  const double* b = Y.slice_memptr(j);
  const arma::uword len = X.n_rows * X.n_cols;
  for (arma::uword t = 0; t < len; ++t) if (a[t] != b[t]) return false;
  return true;
}

// Full Procrustes mean in Helmert coordinates ((k-1) x m, unit norm).
// Fixed-point iteration mu <- normalise(sum_i cos(rho_i) z_i R_i), where R_i
// is the rotation aligning z_i to mu. It is an ascent algorithm for
// sum_i cos^2(rho_i), i.e. a descent algorithm for sum_i d_F^2(z_i, mu).
arma::mat procrustes_mean_helmert(const std::vector<arma::mat>& Z,
                                  double tol, int max_iter, int* iters,
                                  bool* converged) {
  const arma::uword n = Z.size();
  const arma::uword m = Z[0].n_cols;
  arma::mat mu = Z[0];
  double Fold = -1.0;
  int it = 0;
  bool conv = false;
  arma::mat U, V, N, R;
  arma::vec s;
  for (it = 1; it <= max_iter; ++it) {
    arma::mat S(mu.n_rows, m, arma::fill::zeros);
    double F = 0.0;
    for (arma::uword i = 0; i < n; ++i) {
      N = mu.t() * Z[i];
      arma::svd(U, s, V, N);
      arma::vec d(m, arma::fill::ones);
      if (arma::det(U) * arma::det(V) < 0.0) d(m - 1) = -1.0;
      R = V * arma::diagmat(d) * U.t();
      const double c = arma::dot(s, d);
      S += c * (Z[i] * R);
      F += c * c;
    }
    const double nS = arma::norm(S, "fro");
    if (nS > 0.0) mu = S / nS;
    if (std::abs(F - Fold) <= tol * std::max(1.0, F)) { conv = true; break; }
    Fold = F;
  }
  if (iters) *iters = std::min(it, max_iter);
  if (converged) *converged = conv;
  return mu;
}

}  // namespace

// Pairwise shape distances between the slices of X and the slices of Y.
// type: 0 = rho, 1 = dF, 2 = dP. If `same`, Y is ignored and the symmetric
// n x n matrix of X is returned.
// [[Rcpp::export]]
arma::mat shape_dist_cpp(const arma::cube& X, const arma::cube& Y,
                         int type, bool reflection, bool same) {
  const arma::uword k = X.n_rows;
  const arma::mat H = helmert_sub(k);
  const std::vector<arma::mat> Zx = preshapes(X, H);
  if (same) {
    const arma::uword n = X.n_slices;
    arma::mat D(n, n, arma::fill::zeros);
    arma::mat A;
    for (arma::uword i = 0; i < n; ++i) {
      for (arma::uword j = i + 1; j < n; ++j) {
        double d = 0.0;
        if (!identical_config(X, i, X, j)) {
          A = Zx[i].t() * Zx[j];
          d = dist_from_cos(sv_sum(A, reflection), type);
        }
        D(i, j) = d; D(j, i) = d;
      }
    }
    return D;
  }
  const std::vector<arma::mat> Zy = preshapes(Y, H);
  arma::mat D(Zx.size(), Zy.size(), arma::fill::zeros);
  arma::mat A;
  for (arma::uword i = 0; i < Zx.size(); ++i) {
    for (arma::uword j = 0; j < Zy.size(); ++j) {
      double d = 0.0;
      if (!identical_config(X, i, Y, j)) {
        A = Zx[i].t() * Zy[j];
        d = dist_from_cos(sv_sum(A, reflection), type);
      }
      D(i, j) = d;
    }
  }
  return D;
}

// Full Procrustes mean shape of the slices of X, as a centred k x m
// configuration of unit size.
// [[Rcpp::export]]
Rcpp::List procrustes_mean_cpp(const arma::cube& X, double tol, int max_iter) {
  const arma::uword k = X.n_rows;
  const arma::mat H = helmert_sub(k);
  const std::vector<arma::mat> Z = preshapes(X, H);
  int iters = 0;
  bool conv = false;
  const arma::mat mu = procrustes_mean_helmert(Z, tol, max_iter, &iters, &conv);
  return Rcpp::List::create(
    Rcpp::Named("mean")      = arma::mat(H.t() * mu),
    Rcpp::Named("n_iter")    = iters,
    Rcpp::Named("converged") = conv
  );
}

// MDS / phi mean shape (extrinsic mean in the Veronese-Whitney embedding):
// leading m eigenvectors of the average of z z'. The embedding cannot tell
// a shape from its mirror image, so the orientation closest to the
// Procrustes mean `mu` (a k x m configuration) is returned, as in
// shapes::MDSshape.
// [[Rcpp::export]]
arma::mat phi_mean_cpp(const arma::cube& X, const arma::mat& mu) {
  const arma::uword k = X.n_rows;
  const arma::uword m = X.n_cols;
  const arma::mat H = helmert_sub(k);
  const std::vector<arma::mat> Z = preshapes(X, H);
  arma::mat G(k - 1, k - 1, arma::fill::zeros);
  for (arma::uword i = 0; i < Z.size(); ++i) G += Z[i] * Z[i].t();
  G /= static_cast<double>(Z.size());
  arma::vec lam;
  arma::mat vec;
  arma::eig_sym(lam, vec, G);          // ascending
  const arma::uword p = k - 1;
  arma::mat M(p, m);
  double tot = 0.0;
  for (arma::uword j = 0; j < m; ++j) tot += std::max(lam(p - 1 - j), 0.0);
  for (arma::uword j = 0; j < m; ++j) {
    M.col(j) = vec.col(p - 1 - j) * std::sqrt(std::max(lam(p - 1 - j), 0.0)) /
               std::sqrt(tot);
  }
  arma::mat zmu = H * mu;
  const double nz = arma::norm(zmu, "fro");
  if (nz > 0.0) zmu /= nz;
  const arma::mat A = M.t() * zmu;
  if (sv_sum(A, false) < sv_sum(A, true)) M.col(m - 1) *= -1.0;
  return H.t() * M;
}

// -----------------------------------------------------------------------------
//  Medoid searches on a distance matrix. Both use R's random-number stream
//  (R_unif_index draws what base::sample() would draw), so that results are
//  reproducible with set.seed().
// -----------------------------------------------------------------------------
namespace {

// max.col(m, ties.method = "random") for one row, as in R's R_max_col:
// values within RELTOL * max|value| of the running maximum count as ties
// and are resolved with unif_rand(). The caller owns GetRNGstate().
inline int max_col_random(const double* row, int nc) {
  const double RELTOL = 1e-5;
  double large = 0.0;
  for (int c = 0; c < nc; ++c) {
    const double a = row[c];
    if (!R_FINITE(a)) continue;
    large = std::max(large, std::fabs(a));
  }
  int mx = 0;
  double a = row[0];
  const double tol = RELTOL * large;
  int ntie = 1;
  for (int c = 1; c < nc; ++c) {
    const double b = row[c];
    if (b > a + tol) {
      a = b; mx = c; ntie = 1;
    } else if (b >= a - tol) {
      ++ntie;
      if (ntie * unif_rand() < 1.0) mx = c;
    }
  }
  return mx;
}

// Hill-climbing criterion for candidate medoids `cand` (0-based): sum over
// objects of the distance to the nearest candidate, and the assignment.
double hill_criterion(const arma::mat& D, const std::vector<int>& cand,
                      std::vector<int>& belong) {
  const arma::uword n = D.n_rows;
  const int K = static_cast<int>(cand.size());
  std::vector<double> negd(K);
  double est = 0.0;
  for (arma::uword i = 0; i < n; ++i) {
    double mn = std::numeric_limits<double>::infinity();
    for (int j = 0; j < K; ++j) {
      const double d = D(i, cand[j]);
      negd[j] = -d;
      if (d < mn) mn = d;
    }
    est += mn;
    belong[i] = max_col_random(negd.data(), K);
  }
  return est;
}

inline bool has_duplicates(const std::vector<int>& v) {
  for (std::size_t a = 0; a < v.size(); ++a)
    for (std::size_t b = a + 1; b < v.size(); ++b)
      if (v[a] == v[b]) return true;
  return false;
}

// Assignment of every object to its nearest medoid, with the nearest and
// second-nearest distances (after Schubert and Rousseeuw 2021). A medoid is
// always assigned to itself. Returns the total cost.
double assign_nearest(const arma::mat& D, const std::vector<int>& med,
                      std::vector<int>& asg, std::vector<double>& nearest,
                      std::vector<double>& second) {
  const int n = static_cast<int>(D.n_rows);
  const int K = static_cast<int>(med.size());
  double total = 0.0;
  for (int id = 0; id < n; ++id) {
    double d1 = std::numeric_limits<double>::max();
    double d2 = std::numeric_limits<double>::max();
    int i1 = -1;
    for (int i = 0; i < K; ++i) {
      const double d = D(id, med[i]);
      if (id == med[i] || d < d1) { d2 = d1; d1 = d; i1 = i; }
      else if (d < d2) { d2 = d; }
    }
    asg[id] = i1; nearest[id] = d1; second[id] = d2;
    total += d1;
  }
  return total;
}

}  // namespace

// Hill climbing as in hill_climbing_shapes() (Nascimento et al.): `init`
// holds the K initial medoids (0-based). Each iteration replaces one random
// medoid by a random object, redrawing while the candidate set has repeats,
// and keeps the best set seen. With from_best = FALSE the neighbor is always
// built from `init` (as the original script does); with TRUE it is built
// from the best set so far.
// [[Rcpp::export]]
Rcpp::List shape_hill_cpp(const arma::mat& D, const arma::uvec& init,
                          int max_iter, bool from_best) {
  const int n = static_cast<int>(D.n_rows);
  const int K = static_cast<int>(init.n_elem);
  std::vector<int> start(K), best(K), cand(K);
  for (int j = 0; j < K; ++j) start[j] = best[j] = static_cast<int>(init(j));
  std::vector<int> best_belong(n), belong(n);
  double best_est = hill_criterion(D, best, best_belong);
  int n_acc = 0;
  for (int it = 0; it < max_iter; ++it) {
    cand = from_best ? best : start;
    int ridx = static_cast<int>(R_unif_index(K));
    int neigh = static_cast<int>(R_unif_index(n));
    cand[ridx] = neigh;
    while (has_duplicates(cand)) {
      ridx = static_cast<int>(R_unif_index(K));
      neigh = static_cast<int>(R_unif_index(n));
      cand[ridx] = neigh;
    }
    const double est = hill_criterion(D, cand, belong);
    if (est < best_est) {
      best_est = est; best = cand; best_belong = belong; ++n_acc;
    }
  }
  IntegerVector med(K), cl(n);
  for (int j = 0; j < K; ++j) med[j] = best[j] + 1;
  for (int i = 0; i < n; ++i) cl[i] = best_belong[i] + 1;
  return Rcpp::List::create(
    Rcpp::Named("medoids")   = med,
    Rcpp::Named("cluster")   = cl,
    Rcpp::Named("objective") = best_est,
    Rcpp::Named("n_iter")    = n_acc
  );
}

// FastCLARANS (Schubert and Rousseeuw 2021, as in ELKI and fastkmedoids):
// numlocal random restarts; in each, a random non-medoid is tried against
// all K medoids at once, the best swap is kept if it lowers the cost, and
// the search stops once the failure counter j, which starts at 1 and is reset
// to 1 after each swap, reaches `retries`, i.e. after retries - 1 consecutive
// failures, with retries = ceil(maxneighbor * (n - K)) if maxneighbor < 1 and
// ceil(maxneighbor) otherwise.
// [[Rcpp::export]]
Rcpp::List shape_clarans_cpp(const arma::mat& D, int K, int numlocal,
                             double maxneighbor) {
  const int n = static_cast<int>(D.n_rows);
  const int retries = static_cast<int>(
    std::ceil(maxneighbor < 1.0 ? maxneighbor * (n - K) : maxneighbor));
  std::vector<int> subsampler(n);
  for (int i = 0; i < n; ++i) subsampler[i] = i;

  std::vector<int> med(K), asg(n), best_med(K), best_asg(n);
  std::vector<double> nearest(n), second(n), cost(K);
  double best_score = std::numeric_limits<double>::max();
  int best_swaps = 0;

  for (int loc = 0; loc < numlocal; ++loc) {
    // random initial medoids: K distinct draws
    std::vector<char> used(n, 0);
    int got = 0;
    while (got < K) {
      const int r = static_cast<int>(R_unif_index(n));
      if (!used[r]) { used[r] = 1; med[got++] = r; }
    }
    double total = assign_nearest(D, med, asg, nearest, second);
    int swaps = 0;
    int j = 1;
    while (j < retries) {
      // choose a random non-medoid (Fisher-Yates, no repeats within a draw)
      int cand = -1;
      bool skip = false;
      for (int r = 0; r < n; ++r) {
        const int rnd = static_cast<int>(R_unif_index(n - r)) + r;
        std::swap(subsampler[r], subsampler[rnd]);
        cand = subsampler[r];
        if (nearest[cand] > 0.0) break;
        if (second[cand] == 0.0) { ++j; skip = true; break; }
        bool is_med = false;
        for (int q = 0; q < K; ++q) if (med[q] == cand) { is_med = true; break; }
        if (!is_med) break;
        if (r >= 1000) Rcpp::stop("Failed to choose a non-medoid in 1000 attempts.");
      }
      if (skip) continue;
      // cost change of swapping `cand` with each medoid
      std::fill(cost.begin(), cost.end(), 0.0);
      for (int o = 0; o < n; ++o) {
        if (o == cand) continue;
        const double dcur = nearest[o];
        const double dh = D(cand, o);
        const int jcur = asg[o];
        cost[jcur] += std::min(dh, second[o]) - dcur;
        const double change = dh - dcur;
        if (change < 0.0) {
          for (int q = 0; q < K; ++q) if (q != jcur) cost[q] += change;
        }
      }
      int lastbest = 0;
      double mn = cost[0];
      for (int q = 1; q < K; ++q) if (cost[q] < mn) { mn = cost[q]; lastbest = q; }
      if (!(mn < -1e-12 * total)) { ++j; continue; }
      med[lastbest] = cand;
      total = assign_nearest(D, med, asg, nearest, second);
      ++swaps;
      j = 1;
    }
    if (total < best_score) {
      best_score = total; best_med = med; best_asg = asg; best_swaps = swaps;
    }
  }
  IntegerVector mv(K), cl(n);
  for (int q = 0; q < K; ++q) mv[q] = best_med[q] + 1;
  for (int i = 0; i < n; ++i) cl[i] = best_asg[i] + 1;
  return Rcpp::List::create(
    Rcpp::Named("medoids")   = mv,
    Rcpp::Named("cluster")   = cl,
    Rcpp::Named("objective") = best_score,
    Rcpp::Named("n_iter")    = best_swaps
  );
}
