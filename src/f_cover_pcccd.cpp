#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List f_cover_pcccd(double thresh, IntegerMatrix M) {
  int n_main = M.nrow();
  IntegerVector cover(n_main, 0);
  IntegerVector i_dominant;

  IntegerVector n_covers(n_main);
  for (int i = 0; i < n_main; i++) {
    int sum_row = 0;
    for (int j = 0; j < n_main; j++) {
      sum_row += M(i, j);
    }
    n_covers[i] = sum_row;
  }

  while (sum(cover) < thresh) {
    // Zero out cover counts for already covered
    for (int i = 0; i < n_main; i++) {
      if (cover[i] == 1) n_covers[i] = 0;
    }

    int i_selected = which_max(n_covers);
    i_dominant.push_back(i_selected + 1);  // R is 1-based

    // Mark all it covers (and itself)
    for (int j = 0; j < n_main; j++) {
      if (M(i_selected, j) == 1) {
        cover[j] = 1;
      }
    }
    cover[i_selected] = 1;
  }

  double proportion = (double) sum(cover) / n_main;

  return List::create(
    Named("i_dominant") = i_dominant,
    Named("cover_proportion") = proportion
  );
}
