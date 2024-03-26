# jaffelab 0.99.33

* `merge_rse_metrics()` can now also handle `ERCCsumLogErr`. Note that its an
approximation given properties of logarithms. To compute the actual correct
values, you would need access to the original output ERCC quantification output
files.

# jaffelab 0.99.32

* Added the `skipLines` argument to `junctionCount` to resolve
https://github.com/LieberInstitute/jaffelab/issues/8.

# jaffelab 0.99.30

* `cleaningY()` now supports `NA`s in the input `y`.

# jaffelab 0.99.29

* Added the `googledrive_csv()` function for uploading CSV files to Google
Drive. It used to be internal code of `lab_presenters()` but could be more
useful outside of it.

# jaffelab 0.99.28

* Added the `lab_presenters()` function for randomizing lab presenters. Could
also be used for randomizing the lab snacks rotation.

# jaffelab 0.99.27

* Added a `NEWS.md` file to track changes to the package, thus replacing the
existing `NEWS` file.
* Created a `pkgdown::build_site()` version for the documentation of the
package.


# jaffelab 0.99.26

NEW FEATURES

* Added examples to `agePlotter()`, `cleaningY()`, `getR2()`, `edge.pvalue()`
and noted that `junctionCount()` is missing a full example.


# jaffelab 0.99.25

NEW FEATURES

* Added the functions `n_max()` and `n_min()` which return the top n maximum
or minimum values of a vector.


# jaffelab 0.99.23

NEW FEATURES

* Added the function `merge_rse_metrics()` that merges some of the
sample metadata metrics for samples that were sequenced in more than
one lane. Useful for the RangedSummarizedExperiment objects we
typically produce for RNA-seq experiments at the Jaffe lab.


# jaffelab 0.99.22

NEW FEATURES

* Added the `granges_to_ucsc()` function, which is the reverse of
`ucsc_to_granges()`.


# jaffelab 0.99.21

NEW FEATURES

* Added the function `corner()` which takes head of rows and columns of objects


# jaffelab 0.99.20

SIGNIFICANT USER-VISIBLE CHANGES

* `agePlotter()`: changed how the X-axis ticks and range are made for the
prenatal samples.


# jaffelab 0.99.19

NEW FEATURES

* Introduced the `ucsc_to_granges()` function that converts a character vector
of coordinates from the UCSC genome browser to a GRanges object.


# jaffelab 0.99.18

SIGNIFICANT USER-VISIBLE CHANGES

* `agePlotter()` is now more flexible and can handle up a vector of
`pointColor` and `lineColor`. Also, the breaks
are now flexible and support more scenarios (you need at least 3 groups).


# jaffelab 0.99.16

SIGNIFICANT USER-VISIBLE CHANGES

* GitHub repo is now public.


# jaffelab 0.99.15

NEW FEATURES

* Added the function `expression_cutoff()` which finds suggested expression
cutoffs for an RPKM matrix.


# jaffelab 0.99.14

BUG FIXES

* Corrected the documentation of` cleaningY()`. This is after Stephen
Semick and Emily Burke reported issues with the function.



# jaffelab 0.99.13

SIGNIFICANT USER-VISIBLE CHANGES

* Moved `coverage_bwtool()` to the
[recount.bwtool](https://github.com/LieberInstitute/recount.bwtool) package.


# jaffelab 0.99.12

NEW FEATURES

* Added the function `coverage_bwtool()` for computing coverage matrices
using bwtool for a user-specified list of bigWig files and a
user-specified phenotype table. It's heavily based on recount.bwtool


# jaffelab 0.99.11

SIGNIFICANT USER-VISIBLE CHANGES

* `junctionCount()` uses `check.names = FALSE` when calling `DataFrame()`.
Related to a bug Badoi Phan reported.


# jaffelab 0.99.8

BUG FIXES

* Added code by Emily Burke that makes `junctionCount()` work for single-end
libraries. Still need to trace why ? are being added as the strand.


# jaffelab 0.99.10

NEW FEATURES

* Added the `getT()` function.


# jaffelab 0.99.7

NEW FEATURES

* Added `agePlotter()` -- missing examples and tests


# jaffelab 0.99.6

NEW FEATURES

* Added `junctionCount()` -- missing examples and tests


# jaffelab 0.99.5

NEW FEATURES

* Added `getR2()` -- missing examples and tests


# jaffelab 0.99.4

NEW FEATURES

* Added `cleaningY()` -- missing examples and tests


# jaffelab 0.99.3

NEW FEATURES

* Added `getOR()`


# jaffelab 0.99.2

SIGNIFICANT USER-VISIBLE CHANGES

* Dropped `splitit()` since it's available in rafalib.
* Now depending on rafalib making those functions available to the user.


# jaffelab 0.99.0

NEW FEATURES

* Added `ss()` as a wrapper for string `split()` and `sapply()`.
* Added `splitit()` and `split0()` that splits a vector into a list and return
the integer indexes of the elements.
* Added `getPcaVars()` which computes the percent of variance explained
for the principal components of prcomp object.
* Added `getF()`.
