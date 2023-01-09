
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Crossdome <a href=''><img src="man/figures/logo.png" align="right" height="139"/></a>

<!-- badges: start -->

[![R-CMD-check](https://github.com/oandrefonseca/crossdome/workflows/R-CMD-check/badge.svg)](https://github.com/oandrefonseca/crossdome/actions)

<!-- badges: end -->

<br>

## Abstract

**Crossdome: An R package to measure cross-reactivity risk on the
sequence-space**

Currently, several clinical protocols are leveraging on distinct immune
mechanisms, such as adoptive T-cell therapy and peptide-based vaccines.
However, multiple factors can impact the accuracy of these immune-based
applications, such as expression heterogeneity, immunogenicity, and
cross-reactivity (CR) risk. Crossdome was created to measure
cross-reactivity potential based on biochemical properties. Our approach
aims to rank potential agnostic peptides and measure cross-reactivity
risk using TCR binding, and MHC presentation probability. Additionally,
we provide the expression profile related to each CR candidate.

## Features

- `cross_universe()` Peptide database spanning eluted candidates
  (experimentally validated) and custom (user-defined).
- `cross_pair_summary()` Calculates relatedness score between peptides
- `cross_compose()` Predicts relatedness among peptides in a given
  database. Low values are associated with cross-reactive candidates.
- `cross_browser()` Opens an interactive shiny application

## Installation

``` r

devtools::install_github("oandrefonseca/crossdome")
```

## Basic Usage

``` r

library(crossdome)

database <- cross_universe(off_targets = 'ESDPIVAQY', allele = "HLA-A*01:01")
result <- cross_compose(query = 'EVDPIGHLY', background = database)
#> ##------ Sun Jan  8 18:48:20 2023 ------##
str(result)
#> Formal class 'xrResult' [package "crossdome"] with 7 slots
#>   ..@ query          : chr "EVDPIGHLY"
#>   ..@ result         :'data.frame':  37656 obs. of  10 variables:
#>   .. ..$ index            : int [1:37656] 1 2 3 4 5 6 7 8 9 10 ...
#>   .. ..$ query            : chr [1:37656] "EVDPIGHLY" "EVDPIGHLY" "EVDPIGHLY" "EVDPIGHLY" ...
#>   .. ..$ subject          : chr [1:37656] "EVDPIGHLY" "EVDPIGHVY" "EADPTGHSY" "EVDPTSHSY" ...
#>   .. ..$ n_positive       : int [1:37656] 9 9 6 6 6 8 6 7 6 5 ...
#>   .. ..$ n_mismatch       : int [1:37656] 0 1 3 3 4 4 4 5 4 4 ...
#>   .. ..$ relatedness_score: num [1:37656] 0 1.35 8.72 10.43 10.95 ...
#>   .. ..$ zscore           : num [1:37656] -5.48 -5.25 -3.99 -3.69 -3.6 ...
#>   .. ..$ pvalue           : num [1:37656] 2.12e-08 7.66e-08 3.36e-05 1.11e-04 1.57e-04 ...
#>   .. ..$ hla_allele       : chr [1:37656] "HLA-A*01:01" "HLA-A*01:01" "HLA-A*01:01" "HLA-A*01:01" ...
#>   .. ..$ percentile_rank  : num [1:37656] 0 0.00266 0.00531 0.00797 0.01062 ...
#>   ..@ allele         : chr "HLA-A*01:01"
#>   ..@ expression     : list()
#>   ..@ analysis       : list()
#>   ..@ position_weight: num [1:9] 1 1 1 1 1 1 1 1 1
#>   ..@ timestamp      : chr "##------ Sun Jan  8 18:48:20 2023 ------##"
```

## Included datasets

``` r
data("hla_database")
data("hpa_database")
data("peptide_annotation")
data("mage_off_targets")
```
