# circacompare 0.1.1.9001

## Improvements

* allow per-sample weights in `circacompare()` and `circa_single()`. This allows to downweight individual samples, for example in the case of outliers rather than hard-filtered them which reduces power. Such per-sample weights can reproducibly be estimated using approaches such as `arrayWeights()` in the [limma](https://bioconductor.org/packages/release/bioc/html/limma.html) package for datasets with many genes/observations, such as RNA-seq or microarrays.

# circacompare 0.1.1.9000

## Improvements

* better error messages - also makes the shiny app more verbose for users

* allow the user to use `suppress_all` argument to show/hide messages during model fitting across all circacompare functions.

* allow user to specify use of linear model (non-mixed) for fitting models to each group when using circacompare_mixed(). (Faster computation!)

# circacompare 0.1.1

## Improvements

* User-definable parameters for all functions. Notably: period can be estimated rather than always user-defined, rhythmic characteristics can be selectively shared between groups, all parameters can have additional decay terms and those decay terms can have group-effects to determine differences in decay between groups.

* More informative summary tables from all models.

## Bug fixes

* More sensible default values for random effects.

* `nlme` is used for estimating presence of rhythmicity of each group when `circacompare_mixed()` is used, rather than `nls`.
