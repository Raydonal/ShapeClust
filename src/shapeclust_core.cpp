// -----------------------------------------------------------------------------
//  ShapeClust: performance-critical routines in C++ via Rcpp + RcppArmadillo
// -----------------------------------------------------------------------------

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <algorithm>
#include <limits>
#include <vector>

using namespace Rcpp;

// Riemannian distance matrix
// [[Rcpp::export]]
arma::mat rho_matrix_cpp(const arma::cx_mat& Z) {
  const arma::uword n = Z.n_cols;
  arma::mat D(n, n, arma::fill::zeros);
  arma::cx_mat G = Z.t() * Z;
  for (arma::uword i = 0; i < n; ++i) {
    for (arma::uword j = i + 1; j < n; ++j) {
      double c = std::abs(G(i, j));
      if (c > 1.0) c = 1.0;
      const double d = std::acos(c);
      D(i, j) = d; D(j, i) = d;
    }
  }
  return D;
}

// Full Procrustes distance matrix
// [[Rcpp::export]]
arma::mat dF_matrix_cpp(const arma::cx_mat& Z) {
  const arma::uword n = Z.n_cols;
  arma::mat D(n, n, arma::fill::zeros);
  arma::cx_mat G = Z.t() * Z;
  for (arma::uword i = 0; i < n; ++i) {
    for (arma::uword j = i + 1; j < n; ++j) {
      double c = std::abs(G(i, j));
      if (c > 1.0) c = 1.0;
      const double d = std::sqrt(1.0 - c * c);
      D(i, j) = d; D(j, i) = d;
    }
  }
  return D;
}

// Nearest and second-nearest medoid
// [[Rcpp::export]]
arma::umat pam_nearest_two_cpp(const arma::mat& D,
                               const arma::uvec& medoids) {
  const arma::uword n = D.n_rows;
  const arma::uword K = medoids.n_elem;
  arma::umat NS(n, 2, arma::fill::zeros);
  for (arma::uword j = 0; j < n; ++j) {
    double d1 = std::numeric_limits<double>::infinity();
    double d2 = std::numeric_limits<double>::infinity();
    arma::uword i1 = 0, i2 = 0;
    for (arma::uword r = 0; r < K; ++r) {
      const double d = D(j, medoids(r));
      if (d < d1) { d2 = d1; i2 = i1; d1 = d; i1 = r; }
      else if (d < d2) { d2 = d; i2 = r; }
    }
    NS(j, 0) = i1; NS(j, 1) = i2;
  }
  return NS;
}

// PAM swap cost (four cases of Ng & Han 2002)
// [[Rcpp::export]]
double pam_swap_cost_cpp(const arma::mat& D,
                         const arma::uvec& medoids,
                         const arma::uvec& nearest,
                         const arma::uvec& second,
                         arma::uword s,
                         arma::uword h) {
  const arma::uword n = D.n_rows;
  const arma::uword os_idx = medoids(s);
  double TCsh = 0.0;
  for (arma::uword j = 0; j < n; ++j) {
    if (j == os_idx || j == h) continue;
    const arma::uword oi_pos = nearest(j);
    const double d_j_h  = D(j, h);
    const double d_j_os = D(j, os_idx);
    if (oi_pos == s) {
      const arma::uword oi2 = medoids(second(j));
      const double d_j_oi = D(j, oi2);
      TCsh += (d_j_h >= d_j_oi) ? (d_j_oi - d_j_os) : (d_j_h - d_j_os);
    } else {
      const arma::uword oi_med = medoids(oi_pos);
      const double d_j_oi = D(j, oi_med);
      if (d_j_h < d_j_oi) TCsh += (d_j_h - d_j_oi);
    }
  }
  return TCsh;
}

// Row-wise ranks with ties averaged, matching base::rank(ties.method = "average")
// [[Rcpp::export]]
arma::mat rank_matrix_cpp(const arma::mat& D) {
  const arma::uword n = D.n_rows;
  arma::mat R(n, n, arma::fill::zeros);
  for (arma::uword i = 0; i < n; ++i) {
    const arma::rowvec row = D.row(i);
    const arma::uvec ord = arma::stable_sort_index(row);
    arma::uword q = 0;
    while (q < n) {
      arma::uword r = q;
      while (r + 1 < n && row(ord(r + 1)) == row(ord(q))) ++r;
      const double avg = 0.5 * static_cast<double>(q + r) + 1.0;
      for (arma::uword t = q; t <= r; ++t) R(i, ord(t)) = avg;
      q = r + 1;
    }
  }
  return R;
}

// Row-wise ordering (1-based, ties kept in index order, as base::order)
// [[Rcpp::export]]
Rcpp::IntegerMatrix order_matrix_cpp(const arma::mat& D) {
  const arma::uword n = D.n_rows;
  Rcpp::IntegerMatrix SI(n, n);
  for (arma::uword i = 0; i < n; ++i) {
    const arma::uvec ord = arma::stable_sort_index(D.row(i));
    for (arma::uword q = 0; q < n; ++q) SI(i, q) = static_cast<int>(ord(q)) + 1;
  }
  return SI;
}

// -----------------------------------------------------------------------------
//  8. Full Shape-PAM loop in C++ (assignment + best-swap search + update).
//     Doing the whole sweep here avoids one R->C++ call per candidate pair,
//     which dominates the cost when K(n-K) is large.
//     `medoids` is 0-based on entry; the returned first row is the final
//     0-based medoid set, the second block is the 0-based cluster assignment.
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List shape_pam_cpp(const arma::mat& D,
                         arma::uvec medoids,
                         int max_iter) {
  const arma::uword n = D.n_rows;
  const arma::uword K = medoids.n_elem;
  int n_iter = 0;

  arma::uvec nearest(n), second(n);

  for (int it = 0; it < max_iter; ++it) {
    // assignment: nearest and second-nearest medoid (positional)
    for (arma::uword j = 0; j < n; ++j) {
      double d1 = std::numeric_limits<double>::infinity();
      double d2 = std::numeric_limits<double>::infinity();
      arma::uword i1 = 0, i2 = 0;
      for (arma::uword r = 0; r < K; ++r) {
        const double d = D(j, medoids(r));
        if (d < d1) { d2 = d1; i2 = i1; d1 = d; i1 = r; }
        else if (d < d2) { d2 = d; i2 = r; }
      }
      nearest(j) = i1; second(j) = i2;
    }

    // best swap over all (s, h)
    double best_TC = 0.0;
    arma::uword best_s = 0, best_h = 0;
    bool found = false;

    std::vector<bool> is_med(n, false);
    for (arma::uword r = 0; r < K; ++r) is_med[medoids(r)] = true;

    for (arma::uword s = 0; s < K; ++s) {
      const arma::uword os_idx = medoids(s);
      for (arma::uword h = 0; h < n; ++h) {
        if (is_med[h]) continue;
        double TCsh = 0.0;
        for (arma::uword j = 0; j < n; ++j) {
          if (j == os_idx || j == h) continue;
          const arma::uword oi_pos = nearest(j);
          const double d_j_h  = D(j, h);
          if (oi_pos == s) {
            const double d_j_os = D(j, os_idx);
            const double d_j_oi = D(j, medoids(second(j)));
            TCsh += (d_j_h >= d_j_oi) ? (d_j_oi - d_j_os) : (d_j_h - d_j_os);
          } else {
            const double d_j_oi = D(j, medoids(oi_pos));
            if (d_j_h < d_j_oi) TCsh += (d_j_h - d_j_oi);
          }
        }
        if (TCsh < best_TC) { best_TC = TCsh; best_s = s; best_h = h; found = true; }
      }
    }

    if (!found) break;
    medoids(best_s) = best_h;
    ++n_iter;
  }

  // final assignment + objective
  arma::uvec cluster(n);
  double obj = 0.0;
  for (arma::uword j = 0; j < n; ++j) {
    double d1 = std::numeric_limits<double>::infinity();
    arma::uword i1 = 0;
    for (arma::uword r = 0; r < K; ++r) {
      const double d = D(j, medoids(r));
      if (d < d1) { d1 = d; i1 = r; }
    }
    cluster(j) = i1;
    obj += d1;
  }

  return Rcpp::List::create(
    Rcpp::Named("medoids")   = medoids,
    Rcpp::Named("cluster")   = cluster,
    Rcpp::Named("objective") = obj,
    Rcpp::Named("n_iter")    = n_iter
  );
}
