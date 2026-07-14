// -----------------------------------------------------------------------------
//  ShapeClust: performance-critical routines in C++ via Rcpp + RcppArmadillo
// -----------------------------------------------------------------------------

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <algorithm>

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

// Row-wise ranking
// [[Rcpp::export]]
arma::umat rank_matrix_cpp(const arma::mat& D) {
  const arma::uword n = D.n_rows;
  arma::umat R(n, n, arma::fill::zeros);
  for (arma::uword i = 0; i < n; ++i) {
    arma::rowvec row = D.row(i);
    arma::uvec ord = arma::sort_index(row);
    for (arma::uword q = 0; q < n; ++q) R(i, ord(q)) = q + 1;
  }
  return R;
}

// Hostility-based medoid update
// [[Rcpp::export]]
arma::uvec hostility_new_medoids_cpp(const arma::umat& R,
                                     const arma::uvec& clusters,
                                     arma::uword K) {
  const arma::uword n = R.n_rows;
  arma::uvec new_medoids(K, arma::fill::zeros);
  for (arma::uword r = 0; r < K; ++r) {
    std::vector<arma::uword> members;
    members.reserve(n);
    for (arma::uword j = 0; j < n; ++j) if (clusters(j) == r + 1) members.push_back(j);
    const arma::uword m = members.size();
    if (m == 0) { new_medoids(r) = 0; continue; }
    double best_hv = -1.0;
    arma::uword best_i = members[0];
    for (arma::uword a = 0; a < m; ++a) {
      const arma::uword i = members[a];
      double hv = 0.0;
      for (arma::uword b = 0; b < m; ++b) hv += static_cast<double>(R(i, members[b]));
      if (hv > best_hv) { best_hv = hv; best_i = i; }
    }
    new_medoids(r) = best_i;
  }
  return new_medoids;
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
