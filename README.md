Status: Travis CI [![Build Status](https://travis-ci.org/LieberInstitute/jaffelab.svg?branch=master)](https://travis-ci.org/LieberInstitute/jaffelab), Codecov [![codecov.io](https://codecov.io/github/LieberInstitute/jaffelab/coverage.svg?branch=master)](https://codecov.io/github/LieberInstitute/jaffelab?branch=master)

jaffelab
========
[![DOI](https://zenodo.org/badge/70074284.svg)](https://zenodo.org/badge/latestdoi/70074284)

This package contains custom functions that are frequently used by the Jaffe lab. This is a work in progress.

# Installation instructions

Get R 3.4.x from [CRAN](http://cran.r-project.org/).

```R
## If needed:
# install.packages('devtools')

## Install with:
library('devtools')
install_github('LieberInstitute/jaffelab')
```


# Citation

Below is the citation output from using `citation('jaffelab')` in R. Please 
run this yourself to check for any updates on how to cite __jaffelab__.

To cite the __jaffelab__ package in publications use:

Leonardo Collado-Torres and Andrew E. Jaffe (2017). jaffelab: Commonly used functions by the Jaffe lab. R package version 0.99.16. https://github.com/LieberInstitute/jaffelab

@Manual{,
    title = {jaffelab: Commonly used functions by the Jaffe lab},
    author = {Leonardo Collado-Torres and Andrew E. Jaffe},
    year = {2017},
    note = {R package version 0.99.16},
    url = {https://github.com/LieberInstitute/jaffelab},
}

# Testing

Testing on Bioc-devel is feasible thanks to [R Travis](http://docs.travis-ci.com/user/languages/r/).
