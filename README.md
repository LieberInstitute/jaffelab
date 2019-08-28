
<!-- README.md is generated from README.Rmd. Please edit that file -->

# jaffelab

[![Travis-CI build
status](https://travis-ci.org/LieberInstitute/jaffelab.svg?branch=master)](https://travis-ci.org/LieberInstitute/jaffelab)
[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![Codecov test
coverage](https://codecov.io/gh/LieberInstitute/jaffelab/branch/master/graphs/badge.svg)](https://codecov.io/gh/LieberInstitute/jaffelab?branch=master)
[![DOI](https://zenodo.org/badge/70074284.svg)](https://zenodo.org/badge/latestdoi/70074284)

This package contains custom functions that are frequently used by the
Jaffe lab. Please check the package [documentation
website](http://lieberinstitute.github.io/jaffelab) for more
information.

# Installation instructions

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `jaffelab` using the
following code:

``` r
## If needed:
if (!requireNamespace("remotes", quietly = TRUE))
   install.packages("remotes")

## Install with:
remotes::install_github('LieberInstitute/jaffelab')
```

# Citation

Below is the citation output from using `citation('jaffelab')` in R.
Please run this yourself to check for any updates on how to cite
**jaffelab**.

``` r
citation('jaffelab')
#> 
#> To cite package 'jaffelab' in publications use:
#> 
#>   Leonardo Collado-Torres, Andrew E. Jaffe and Emily E. Burke
#>   (2019). jaffelab: Commonly used functions by the Jaffe lab. R
#>   package version 0.99.27.
#>   https://github.com/LieberInstitute/jaffelab
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {jaffelab: Commonly used functions by the Jaffe lab},
#>     author = {Leonardo Collado-Torres and Andrew E. Jaffe and Emily E. Burke},
#>     year = {2019},
#>     note = {R package version 0.99.27},
#>     url = {https://github.com/LieberInstitute/jaffelab},
#>   }
```

# Testing

Testing on Bioc-devel is feasible thanks to [R
Travis](http://docs.travis-ci.com/user/languages/r/).
