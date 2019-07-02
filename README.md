
<!-- README.md is generated from README.Rmd. Please edit that file -->
schex
=====

<!-- badges: start -->
[![Travis build status](https://travis-ci.org/SaskiaFreytag/schex.svg?branch=master)](https://travis-ci.org/SaskiaFreytag/schex) [![Codecov test coverage](https://codecov.io/gh/SaskiaFreytag/schex/branch/master/graph/badge.svg)](https://codecov.io/gh/SaskiaFreytag/schex?branch=master) <!-- badges: end -->

The goal of schex is to provide easy plotting of hexagon cell representations of single cell data stored in `SingleCellExperiment` or `Seurat` objects.

<img src="misc/schex_hex.png" style="width:25.0%" />

Installation
------------

You can install the development version of schex with:

``` r
# install.packages("devtools")
devtools::install_github("SaskiaFreytag/schex")
```

Why you need schex?
-------------------

Did you know that order in which points are plotted depends on their location in the data frame? For example when plotting the expression of CD19, a B-cell maker, you may get the following three plots depending on how you order your observations.

<img src="misc/example_schex_files/figure-html/ggplot-decreasing-1.png" alt="Observations in decreasing order with regrads to their CD19 expression" style="width:49.0%" /> <img src="misc/example_schex_files/figure-html/ggplot-increasing-1.png" alt="Observations in increasing order with regrads to their CD19 expression" style="width:49.0%" /> <img src="misc/example_schex_files/figure-html/ggplot-random-1.png" alt="Observation in random order" style="width:49.0%" />

Using the first plot you would not decide to call the central cluster a B-cell population. Using the second plot you would probably decide to call the same cluster a B-cell population. Using the last plot you might be undecided.

The solution
------------

Instead of plotting points on top of each other, schex summarizes points into hexagon cells. Hence avoiding confusion due to observation order.

<img src="misc/example_schex_files/figure-html/schex-1.png" alt="schex plotting" style="width:49.0%" />
