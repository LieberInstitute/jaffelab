
<!-- README.md is generated from README.Rmd. Please edit that file -->

# jaffelab

<!-- badges: start -->

[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![R-CMD-check-bioc](https://github.com/LieberInstitute/jaffelab/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/LieberInstitute/jaffelab/actions)
[![Codecov test
coverage](https://codecov.io/gh/LieberInstitute/jaffelab/branch/master/graph/badge.svg)](https://codecov.io/gh/LieberInstitute/jaffelab?branch=master)
<!-- badges: end -->

`jaffelab` is a package initially developed by Andrew E Jaffe’s team at
the Lieber Institute for Brain Development. It contains several helper
functions used in their analyses.

## Installation instructions

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `jaffelab` using the
following code:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("LieberInstitute/jaffelab")
```

## Citation

Below is the citation output from using `citation('jaffelab')` in R.
Please run this yourself to check for any updates on how to cite
**jaffelab**.

``` r
print(citation('jaffelab'), bibtex = TRUE)
#> 
#> To cite package 'jaffelab' in publications use:
#> 
#>   Leonardo Collado-Torres, Andrew E. Jaffe and Emily E. Burke (2021).
#>   jaffelab: Commonly used functions by the Jaffe lab. R package version
#>   0.99.31. https://github.com/LieberInstitute/jaffelab
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {jaffelab: Commonly used functions by the Jaffe lab},
#>     author = {Leonardo Collado-Torres and Andrew E. Jaffe and Emily E. Burke},
#>     year = {2021},
#>     note = {R package version 0.99.31},
#>     url = {https://github.com/LieberInstitute/jaffelab},
#>   }
```

Please note that the `jaffelab` was only made possible thanks to many
other R and bioinformatics software authors, which are cited either in
the vignettes and/or the paper(s) describing this package.

## Code of Conduct

Please note that the `jaffelab` project is released with a [Contributor
Code of Conduct](http://bioconductor.org/about/code-of-conduct/). By
contributing to this project, you agree to abide by its terms.

## Development tools

-   Continuous code testing is possible thanks to [GitHub
    actions](https://www.tidyverse.org/blog/2020/04/usethis-1-6-0/)
    through *[usethis](https://CRAN.R-project.org/package=usethis)*,
    *[remotes](https://CRAN.R-project.org/package=remotes)*, and
    *[rcmdcheck](https://CRAN.R-project.org/package=rcmdcheck)*
    customized to use [Bioconductor’s docker
    containers](https://www.bioconductor.org/help/docker/) and
    *[BiocCheck](https://bioconductor.org/packages/3.13/BiocCheck)*.
-   Code coverage assessment is possible thanks to
    [codecov](https://codecov.io/gh) and
    *[covr](https://CRAN.R-project.org/package=covr)*.
-   The [documentation
    website](http://LieberInstitute.github.io/jaffelab) is automatically
    updated thanks to
    *[pkgdown](https://CRAN.R-project.org/package=pkgdown)*.
-   The code is styled automatically thanks to
    *[styler](https://CRAN.R-project.org/package=styler)*.
-   The documentation is formatted thanks to
    *[devtools](https://CRAN.R-project.org/package=devtools)* and
    *[roxygen2](https://CRAN.R-project.org/package=roxygen2)*.

For more details, check the `dev` directory.

This package was developed using
*[biocthis](https://bioconductor.org/packages/3.13/biocthis)*.
