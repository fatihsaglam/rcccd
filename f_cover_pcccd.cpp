#include <Rcpp.h>
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

