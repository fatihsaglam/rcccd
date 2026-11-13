#include <Rcpp.h>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

// Random walk cover. For every target sample the distances to all samples
// (target first, then non-target) are sorted ONCE; the order does not change
// between iterations, only the weights do (covered samples get weight 0).
// [[Rcpp::export]]
List cover_rwcccd_cpp(NumericMatrix x_main,
                      NumericMatrix x_other,
                      NumericVector w,
                      double proportion,
                      double m,
                      double p) {
  const int n_main = x_main.nrow();
  const int n_other = x_other.nrow();
  const int n = n_main + n_other;
  const int dim = x_main.ncol();

  std::vector<double> xt((size_t) n * dim);
  for (int k = 0; k < dim; k++) {
    for (int i = 0; i < n_main; i++) xt[(size_t) i * dim + k] = x_main(i, k);
    for (int i = 0; i < n_other; i++) xt[(size_t) (n_main + i) * dim + k] = x_other(i, k);
  }

  std::vector<int> ord((size_t) n_main * n);
  std::vector<double> dsort((size_t) n_main * n);
  std::vector<double> pen((size_t) n_main * n);
  std::vector<double> dmax_main(n_main, 0.0);

  {
    std::vector<double> d(n);
    std::vector<int> idx(n);
    for (int j = 0; j < n_main; j++) {
      const double* a = &xt[(size_t) j * dim];
      double mx = 0.0;
      for (int i = 0; i < n; i++) {
        const double* b = &xt[(size_t) i * dim];
        double s = 0.0;
        for (int k = 0; k < dim; k++) {
          double t = a[k] - b[k];
          s += t * t;
        }
        d[i] = std::sqrt(s);
        if (i < n_main && d[i] > mx) mx = d[i];
      }
      dmax_main[j] = mx;
      std::iota(idx.begin(), idx.end(), 0);
      std::sort(idx.begin(), idx.end(), [&](int u, int v) {
        return d[u] < d[v] || (d[u] == d[v] && u < v);
      });
      const size_t base = (size_t) j * n;
      for (int k = 0; k < n; k++) {
        ord[base + k] = idx[k];
        dsort[base + k] = d[idx[k]];
        pen[base + k] = m * std::pow(d[idx[k]], p);
      }
    }
  }

  std::vector<double> ww(w.begin(), w.end());
  std::vector<char> cover_main(n_main, 0), code_dom(n_main, 0);
  NumericVector radii(n_main), T_score(n_main);
  int n_uncovered = n_main;
  int n_dom = 0;
  double cover_proportion = 0.0;

  std::vector<double> T_scores(n_main), rs(n_main);

  while (true) {
    int sel = -1;
    double best = R_NegInf;
    int first_cand = -1;
    for (int j = 0; j < n_main; j++) {
      if (code_dom[j] || cover_main[j]) continue;
      if (first_cand < 0) first_cand = j;
      const size_t base = (size_t) j * n;
      double R = 0.0, best_val = R_NegInf, R_at = 0.0;
      int k_at = 0;
      for (int k = 0; k < n; k++) {
        R += ww[ord[base + k]];
        double v = R - pen[base + k];
        if (v > best_val) {
          best_val = v;
          R_at = R;
          k_at = k;
        }
      }
      const double r = dsort[base + k_at];
      const double ratio = dmax_main[j] > 0 ? r / dmax_main[j] : 0.0;
      const double T = R_at - ratio * (n_uncovered / 2.0);
      rs[j] = r;
      T_scores[j] = T;
      if (T > best) {
        best = T;
        sel = j;
      }
    }
    if (first_cand < 0) break;     // no candidate left
    if (sel < 0) sel = first_cand; // all scores NaN

    const double r_sel = rs[sel];
    const size_t base = (size_t) sel * n;
    for (int k = 0; k < n; k++) {
      if (!(dsort[base + k] < r_sel)) break;
      const int idx = ord[base + k];
      if (idx < n_main && !cover_main[idx]) {
        cover_main[idx] = 1;
        n_uncovered--;
      }
      ww[idx] = 0.0;
    }

    code_dom[sel] = 1;
    n_dom++;
    radii[sel] = r_sel;
    T_score[sel] = T_scores[sel];
    cover_proportion = (double) (n_main - n_uncovered) / n_main;

    if (cover_proportion >= proportion) break;
    if (n_dom == n_main) break;
    if (n_uncovered == 0) break;
  }

  IntegerVector code_dominant(n_main);
  for (int j = 0; j < n_main; j++) code_dominant[j] = code_dom[j];

  return List::create(Named("code_dominant") = code_dominant,
                      Named("radii") = radii,
                      Named("T_score") = T_score,
                      Named("cover_proportion") = cover_proportion);
}
