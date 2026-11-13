#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// For every new sample: min over dominant balls of (d / radius)^expo.
// [[Rcpp::export]]
NumericVector min_ratio_cpp(NumericMatrix x_dom,
                            NumericMatrix x_new,
                            NumericVector radii,
                            NumericVector expo) {
  const int nd = x_dom.nrow();
  const int n = x_new.nrow();
  const int p = x_dom.ncol();
  NumericVector out(n, R_PosInf);

  for (int i = 0; i < n; i++) {
    double best = R_PosInf;
    for (int j = 0; j < nd; j++) {
      double s = 0.0;
      for (int k = 0; k < p; k++) {
        double t = x_new(i, k) - x_dom(j, k);
        s += t * t;
      }
      double v = std::sqrt(s) / radii[j];
      if (expo[j] != 1.0) v = std::pow(v, expo[j]);
      if (v < best) best = v;
    }
    out[i] = best;
  }
  return out;
}

// Euclidean distance matrix between rows of a (n_a x p) and rows of b (n_b x p).
// [[Rcpp::export]]
NumericMatrix cross_dist_cpp(NumericMatrix a, NumericMatrix b) {
  const int na = a.nrow(), nb = b.nrow(), p = a.ncol();
  NumericMatrix out(na, nb);
  for (int j = 0; j < nb; j++) {
    for (int i = 0; i < na; i++) {
      double s = 0.0;
      for (int k = 0; k < p; k++) {
        double t = a(i, k) - b(j, k);
        s += t * t;
      }
      out(i, j) = std::sqrt(s);
    }
  }
  return out;
}
