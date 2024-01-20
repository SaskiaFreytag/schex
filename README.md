
<!-- README.md is generated from README.Rmd. Please edit that file -->

# schex

<!-- badges: start -->

[![](https://img.shields.io/badge/download-129/total-blue.svg)](https://bioconductor.org/packages/stats/bioc/schex)
[![](https://img.shields.io/badge/download-81/total-yellow.svg)](https://bioconductor.org/packages/stats/bioc/schex)
[![Codecov test
coverage](https://codecov.io/gh/SaskiaFreytag/schex/branch/master/graph/badge.svg)](https://codecov.io/gh/SaskiaFreytag/schex?branch=master)
<!-- badges: end -->
<img src='man/figures/logo.png' align="right" height="139" />

The goal of schex is to provide easy plotting of hexagon cell
representations of single cell data stored in `SingleCellExperiment` objects.

## Installation

You can install schex using the [Bioconductor
project](https://bioconductor.org/):

``` r
# install.packages("BiocManager")
BiocManager::install("schex")
```

You can install the development version of schex with:

``` r
# install.packages("devtools")
devtools::install_github("SaskiaFreytag/schex")
```

## Why you need schex?

Did you know that the order in which points are plotted depends on their
location in the data frame? For example when plotting the expression of
CD19, a B-cell maker, you may get the following three plots depending on
how you order your observations.

![Observations in decreasing order with regrads to their CD19
expression](man/figures/figure-html/ggplot-decreasing-1.png)
![Observations in increasing order with regrads to their CD19
expression](man/figures/figure-html/ggplot-increasing-1.png)
![Observation in random
order](man/figures/figure-html/ggplot-random-1.png)

Using the first plot you would not decide to call the central cluster a
B-cell population. Using the second plot you would probably decide to
call the same cluster a B-cell population. Using the last plot you might
be undecided.

## The solution

Instead of plotting points on top of each other, schex summarizes points
into hexagon cells. Hence avoiding confusion due to observation order.

![schex plotting](man/figures/figure-html/schex-1.png)

Check out the vignettes to learn about how to get started. Or for a
quick start, use the following code.

    library(TENxPBMCData)
    library(scater)
    tenx_pbmc3k <- TENxPBMCData(dataset = "pbmc3k")
    rm_ind <- calculateAverage(tenx_pbmc3k) < 0.1
    tenx_pbmc3k <- tenx_pbmc3k[!rm_ind, ]
    tenx_pbmc3k <- logNormCounts(tenx_pbmc3k)
    tenx_pbmc3k <- runPCA(tenx_pbmc3k)
    tenx_pbmc3k <- make_hexbin(tenx_pbmc3k, 10, dimension_reduction = "PCA")
    plot_hexbin_density(tenx_pbmc3k)
