# circacompare
`circacompare` is an R package that allows for the statistical analyses and comparison of two circadian rhythms.
This work is published [here](https://academic.oup.com/bioinformatics/article-abstract/doi/10.1093/bioinformatics/btz730/5582266) and can be cited as: 


Rex Parsons, Richard Parsons, Nicholas Garner, Henrik Oster, Oliver Rawashdeh, CircaCompare: A method to estimate and statistically support differences in mesor, amplitude, and phase, between circadian rhythms, Bioinformatics, https://doi.org/10.1093/bioinformatics/btz730


# Installation

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
The data uploaded should be in csv format. The columns in the file are required:
1. a time variable (which should be numeric and in hours)
2. a grouping variable (which can be of any format but must have only two possible values)
3. a outcome variable (which should be numeric)

Upload your csv file and select the respective columns from the dropdown menu.  Click 'run' to conduct the comparison.
