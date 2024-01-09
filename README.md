<!-- badges: start -->
[![DOI](https://img.shields.io/badge/doi-10.1093/bioinformatics/btz730-green.svg)](https://doi.org/10.1093/bioinformatics/btz730)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/last-month/circacompare)](https://www.r-pkg.org/pkg/circacompare)
[![R-CMD-check](https://github.com/RWParsons/circacompare/workflows/R-CMD-check/badge.svg)](https://github.com/RWParsons/circacompare/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/license/mit/)
<!-- badges: end -->

# circacompare
`circacompare` is an R package that allows for the statistical analyses and comparison of two circadian rhythms.
This work is published [here](https://academic.oup.com/bioinformatics/article-abstract/doi/10.1093/bioinformatics/btz730/5582266) and can be cited as: 


Rex Parsons, Richard Parsons, Nicholas Garner, Henrik Oster, Oliver Rawashdeh, CircaCompare: A method to estimate and statistically support differences in mesor, amplitude, and phase, between circadian rhythms, Bioinformatics, https://doi.org/10.1093/bioinformatics/btz730

There have been several improvements to the package since initial release. In addition to what was available in 1.0.0 and described in the publication, the package offers approaches to:
* Perform analysis on a single rhythmic dataset to estimate its mesor, amplitude and phase.
* Choose to use a known period (user-determined) or to let the model estimate the period from the data.
* Add parameters to estimate the exponential decay in any of the rhythmic characteristics.
* Use a mixed-model instead of a fixed effects model to take into account within-subject correlation regarding any set of rhythmic parameters.
* Perform a comparison between groups all or a subset of rhythmic characteristics.

Please see the vignette for full details examples of these features.

```
browseVignettes(package="circacompare")
```

# Installation

### Installing from CRAN

```
install.packages("circacompare")
```



### Installing development version from GitHub

If you have not done so already, install `devtools` using the following code:

```
install.packages("devtools")
```

Then install circacompare directly using the following code:
```
devtools::install_github("RWParsons/circacompare")
```
# Help

Once loaded into R, load the documentation using `?circacompare`.  As per the example, you can use the other function, `make_data`, to generate example data appropriate for use in the `circacompare` function.

If you're having further troubles or suggestions for improvement, please create an issue or email me (rex.parsons94@gmail.com).

# ShinyR application

An implementation of the `circacompare` program is available as a Shiny app here: https://rwparsons.shinyapps.io/circacompare/
The data uploaded should be in csv format. The file which you upload ought to have columns for:
1. a time variable (which should be numeric and in hours)
2. a grouping variable (which can be of any format but must have only two possible values)
3. an outcome variable (which should be numeric)

Upload your csv file and select the respective columns from the dropdown menu.  Click 'run' to conduct the comparison.

# Python implementation

An implementation of this package in Python is available [here](https://github.com/RWParsons/circacompare_py). This package comes with the added functionality of specification for the loss function.
