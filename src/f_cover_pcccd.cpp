// Legacy greedy cover of rcccd 0.3.2. Kept only for backward compatibility with
// packages that call rcccd:::f_cover_pcccd() (not exported, not used by rcccd).
#include <Rcpp.h>
#include <vector>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
List f_cover_pcccd(IntegerVector cover,
                  double thresh,
                  NumericMatrix M,
                  NumericMatrix dist_main2main,
                  NumericVector dist_main2other) {
  IntegerVector i_dominant;
  int n_main = M.nrow();
  int counter = 0;
  int sum_cover_old = sum(cover);
  int sum_cover = sum(cover);
  while (sum_cover < thresh) {
    sum_cover_old = sum_cover;
    IntegerVector k_covered = rep(0, n_main);
    for (int i = 0; i < n_main; i++) {
      for (int j = 0; j < n_main; j++) {
        if (cover[j] == 0) {
          k_covered[i] = k_covered[i] + M(i,j);
        }
      }
    }
    i_dominant.push_back(which_max(k_covered));
    LogicalVector i_covered = M(i_dominant[counter], _) > 0;
    for (int i = 0; i < n_main; i++) {
      if (i_covered[i] == 1) {
        cover[i] = 1;
      }
      cover[i_dominant[counter]] = 1;
    }
    sum_cover = sum(cover);
    if (sum_cover == sum_cover_old) {
      break;
    }
    counter ++;
  }
  return List::create(Named("i_dominant") = i_dominant + 1,
                      Named("cover_proportion") = sum_cover/thresh);
}


// ---------------------------------------------------------------------------

// Pure and proper cover. Ball i has center x_main[i, ] and radius radii[i] and
// covers sample j if d(i, j) < radii[i]. Distances are computed on the fly, so
// no n x n distance matrix is stored. Greedy rule: the ball covering the most
// samples is selected first (first one in case of ties); covered samples are
// not selected again.
// [[Rcpp::export]]
List cover_pcccd_cpp(double thresh, NumericMatrix x_main, NumericVector radii) {
  const int n = x_main.nrow();
  const int p = x_main.ncol();

  // transposed copy: samples are contiguous
  std::vector<double> xt((size_t) n * p);
  for (int k = 0; k < p; k++)
    for (int i = 0; i < n; i++)
      xt[(size_t) i * p + k] = x_main(i, k);

  auto dist = [&](int i, int j) {
    const double* a = &xt[(size_t) i * p];
    const double* b = &xt[(size_t) j * p];
    double s = 0.0;
    for (int k = 0; k < p; k++) {
      double d = a[k] - b[k];
      s += d * d;
    }
    return std::sqrt(s);
  };

  std::vector<int> n_covers(n, 0);
  for (int i = 0; i < n; i++) {
    const double r = radii[i];
    int cnt = 0;
    for (int j = 0; j < n; j++) {
      if (dist(i, j) < r) cnt++;
    }
    n_covers[i] = cnt;
  }

  std::vector<char> covered(n, 0);
  std::vector<int> i_dominant;
  int n_covered = 0;

  while (n_covered < thresh && (int) i_dominant.size() < n) {
    int sel = 0;
    int best = -2;
    for (int i = 0; i < n; i++) {
      if (covered[i]) continue;
      if (n_covers[i] > best) {
        best = n_covers[i];
        sel = i;
      }
    }
    if (best == -2) break;
    i_dominant.push_back(sel + 1);

    const double r = radii[sel];
    for (int j = 0; j < n; j++) {
      if (!covered[j] && dist(sel, j) < r) {
        covered[j] = 1;
        n_covered++;
      }
    }
    if (!covered[sel]) {
      covered[sel] = 1;
      n_covered++;
    }
  }

  return List::create(Named("i_dominant") = wrap(i_dominant),
                      Named("cover_proportion") = (double) n_covered / n);
}
