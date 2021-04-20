# circacompare 0.1.1

## Improvements

* User-definable parameters for all functions. Notably: period can be estimated rather than always user-defined, rhythmic characteristics can be selectively shared between groups, all parameters can have additional decay terms and those decay terms can have group-effects to determine differences in decay between groups.

* More informative summary tables from all models.

## Bug fixes

* More sensible default values for random effects.

* `nlme` is used for estimating presence of rhythmicity of each group when `circacompare_mixed()` is used, rather than `nls`.
