# circacompare 0.1.1.9001

## Improvements

* allow per-sample weights in `circacompare()`. This is experimental as it has only been added to the
`nls()` function that runs the actual differential analysis, not to the ones in the `model_each_group` function which models both groups separately to assess per-group rhythmicity.

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
