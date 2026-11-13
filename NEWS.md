# rcccd 0.4.0

## Speed and memory

* PCCCD cover: the greedy dominant-ball search is rewritten in C++ and computes
  distances on the fly, so the n x n distance and indicator matrices are no
  longer built. Memory use is now linear in the number of samples.
* RWCCCD cover: the sorting of distances is done once per sample instead of in
  every iteration of the greedy search; only the weights change between
  iterations. This removes most of the computation time on larger data.
* `predict()` methods compute the nearest-ball ratios in C++ instead of building
  a full dominant-by-new-sample distance matrix with 'proxy'. The 'proxy'
  dependency is dropped, and the 'RcppArmadillo' dependency is no longer needed.

## Bug fixes

* `rwcccd_classifier()`: the T score used an integer division (`n_uncovered/2`)
  in the C++ code, so it differed from `rwcccd_classifier_2()` whenever the
  number of uncovered samples was odd. Both functions now give the same result.
  Dominant balls (and therefore predictions with the default `e = 0`) are
  unchanged in most cases, but can differ in a few data sets because of this fix.
* A predictor matrix with a single column or a class with a single observation
  no longer fails (`x[i, ]` dropped the dimensions).
* Ties in `predict()` are now broken deterministically (first class) instead of
  randomly.
* Classes without any observation in the training data are allowed (e.g. in a
  bootstrap sample of the ensemble); they get probability zero.
* `pcccd_ensemble_classifier()` resamples when a sample contains a single class,
  and its `predict()` is vectorized.

## Other changes

* Input checks for all classifiers: `x` must be numeric without missing values,
  `y` must be a factor of matching length without missing values, and at least
  two classes must be observed.
* `rwcccd_classifier_2()` is kept as a slow pure R reference implementation.
* Added tests (testthat), including checks of the C++ covers against plain R
  implementations.
* `f_cover_pcccd()` (internal, not exported) keeps its old signature so that
  packages calling it, such as 'imbalanceDatRel', keep working.
