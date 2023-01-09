
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Crossdome <a href='https://XXX.XXXXXX.org'><img src='man/figures/logo.png' align="right" height="139" /></a>

<!-- badges: start -->

[![R-CMD-check](https://github.com/oandrefonseca/crossdome/workflows/R-CMD-check/badge.svg)](https://github.com/oandrefonseca/crossdome/actions)
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

## Usage

``` r

library(crossdome)
```

## Getting help
