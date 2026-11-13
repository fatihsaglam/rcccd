## Release summary

rcccd 0.4.0. Faster C++ cover algorithms (lower memory use), 'proxy' and
'RcppArmadillo' dependencies dropped, bug fixes, input checks and tests. See
NEWS.md.

## Test environments

* local Windows 11, R (release)
* GitHub Actions / R-hub (add results here if run)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

The only reverse dependency was 'imbalanceDatRel', which calls the internal
function `rcccd:::f_cover_pcccd()`. That function is kept with its previous
signature and behaviour, and 'imbalanceDatRel' is superseded by 'imbalanceKit'.
