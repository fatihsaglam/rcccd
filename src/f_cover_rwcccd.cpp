#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List f_cover_rwcccd(
    const arma::mat& dist_main2all,
    const arma::mat& dist_main2main,
    const arma::vec& w,
    double proportion,
    double m,
    double p) {

  int n_main = dist_main2all.n_rows;
  int n = dist_main2all.n_cols;

  vec cover_main(n_main, fill::zeros);
  vec cover_all(n, fill::zeros);
  vec code_dominant(n_main, fill::zeros);
  vec radii(n_main);
  vec T_score(n_main);
  double cover_proportion;

  vec dd;
  vec T_scores(n_main);
  vec rs(n_main);
  int i_selected;
  double T_score_selected;
  double r_selected;
  vec d_m2m_i_selected;
  vec d_m2a_i_selected;
  vec ww = w;

  while (sum((code_dominant + cover_main) == 0) >= 1) {
    T_scores = rep(-1e200, n_main);
    rs = rep(-1e200, n_main);

    uvec i_potentialdom = find((code_dominant == 0) && (cover_main == 0));

    uvec ii_sort;
    vec dd_sorted;
    vec w_sorted;
    vec R_all;
    vec P_all;
    vec R_penalized_all;
    int i_R_penalized_max;
    int n_uncovered = sum(cover_main == 0);

    for (int j : i_potentialdom) {
      dd = dist_main2all.row(j).t();
      ii_sort = sort_index(dd);
      dd_sorted = dd(ii_sort);
      w_sorted = ww(ii_sort);
      R_all = cumsum(w_sorted);
      P_all = m*pow(dd_sorted, p);
      R_penalized_all = R_all - P_all;
      i_R_penalized_max = R_penalized_all.index_max();
      rs(j) = dd_sorted(i_R_penalized_max);
      T_scores(j) = R_all(i_R_penalized_max) - rs(j)/max(dist_main2main.row(j)) * (n_uncovered/2);
    }

    i_selected = T_scores.index_max();
    T_score_selected = T_scores(i_selected);
    r_selected = rs(i_selected);

    d_m2m_i_selected = dist_main2main.row(i_selected).t();
    cover_main(find(d_m2m_i_selected < r_selected)) = ones<vec>(sum(d_m2m_i_selected < r_selected));

    d_m2a_i_selected = dist_main2all.row(i_selected).t();
    cover_all(find(d_m2a_i_selected < r_selected)) = ones<vec>(sum(d_m2a_i_selected < r_selected));
    ww(find(d_m2a_i_selected < r_selected)) = zeros<vec>(sum(d_m2a_i_selected < r_selected));

    code_dominant(i_selected) = 1;
    radii(i_selected) = r_selected;
    T_score(i_selected) = T_score_selected;
    cover_proportion = sum(cover_main)/n_main;

    if (cover_proportion >= proportion) {
      break;
    }
    if (sum(code_dominant) == n_main) {
      break;
    }
    if (sum(cover_main) == n_main) {
      break;
    }
  }

  return List::create(Named("code_dominant") = code_dominant,
                      Named("radii") = radii,
                      Named("T_score") = T_score,
                      Named("cover_proportion") = cover_proportion);
}

