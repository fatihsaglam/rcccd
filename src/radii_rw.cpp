#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec radii_rw(const arma::mat& dist_main2all,
                   const arma::mat& dist_main2main,
                   const arma::vec& w,
                   const double m,
                   const int p) {
  int n_main = dist_main2main.n_rows;
  arma::vec radii(n_main);

  for (int j = 0; j < n_main; ++j) {
    arma::uvec i_order = sort_index(dist_main2all.row(j));
    arma::vec w_ordered = w.elem(i_order);
    arma::vec R_x_all = cumsum(w_ordered);
    arma::vec r = sort(dist_main2main.row(j)).t();
    arma::vec k_x = R_x_all.subvec(0, n_main - 1) - m * arma::pow(r, p);
    int max_idx = k_x.index_max();

    radii[j] = r[max_idx];
  }

  return radii;
}
