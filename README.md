
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Crossdome <a href='https://XXX.XXXXXX.org'><img src='man/figures/logo.png' align="right" height="139" /></a>

<!-- badges: start -->

[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/oandrefonseca/crossdome?branch=main&svg=true)](https://ci.appveyor.com/project/oandrefonseca/crossdome)
[![R-CMD-check](https://github.com/oandrefonseca/crossdome/workflows/R-CMD-check/badge.svg)](https://github.com/oandrefonseca/crossdome/actions)
[![Codecov test
coverage](https://codecov.io/gh/oandrefonseca/crossdome/branch/main/graph/badge.svg)](https://codecov.io/gh/oandrefonseca/crossdome?branch=main)
[![CRAN
status](https://www.r-pkg.org/badges/version/crossdome)](https://CRAN.R-project.org/package=crossdome)
<!-- badges: end --> <br>

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

## Description

-   `cross_universe()` Description
-   `cross_bp_summary()` Description
-   `cross_compose()` Description
-   `cross_browser()` Description

## Installation

``` r
devtools::install_github("oandrefonseca/crossdome")
```

## Usage

``` r
library(crossdome)
#> Warning: replacing previous import 'shiny::dataTableOutput' by
#> 'DT::dataTableOutput' when loading 'crossdome'
#> Warning: replacing previous import 'shiny::renderDataTable' by
#> 'DT::renderDataTable' when loading 'crossdome'
```

## Getting help
