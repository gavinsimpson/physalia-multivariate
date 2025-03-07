# Multivariate data analysis with R and vegan

### Physalia-Courses 

https://www.physalia-courses.org/

### Gavin Simpson

#### 10th &ndash; 13th February, 2025

## Overview

The R statistical language has enjoyed wide and rapid adoption by many
ecologists, and is used across many ecological subdisciplines for statistical
analyses and the production of publication-quality figures. For community
ecologists using R, one of the most-used, and useful, add-on packages is vegan,
which provides a wide range of functionality covering inter alia ordination,
diversity analysis, and ecological simulation. This workshop will offer
participants a practical introduction to some of the most useful functions
available within vegan. We will focus on the use of ordination methods and on 
the use of restricted permutations to test a range of experimental designs.

In particular, we discuss when and how to use multivariate methods
including unconstrained and constrained ordination (CCA, RDA, Constrained PCoA),
as well as between-group tests such as PERMANOVA. We will cover concepts such
as design- and model-based permutations and the exchangeability of samples in
tests. We will also discuss the use of vegan to go beyond simply fitting a
constrained ordination model, to diagnostics, plotting, etc.

## Slides

* [Monday](https://gavinsimpson.github.io/physalia-multivariate/01-monday/slides.html)
    * [PDF](https://gavinsimpson.github.io/physalia-multivariate/01-monday/slides.pdf)
* [Tuesday](https://gavinsimpson.github.io/physalia-multivariate/02-tuesday/slides.html)
    * [PDF](https://gavinsimpson.github.io/physalia-multivariate/02-tuesday/slides.pdf)
* [Wednesday](https://gavinsimpson.github.io/physalia-multivariate/03-wednesday/slides.html)
    * [PDF](https://gavinsimpson.github.io/physalia-multivariate/03-wednesday/slides.pdf)
    * [Ponds example](https://gavinsimpson.github.io/physalia-multivariate/03-wednesday/constrained-ordination.html)
    * [Spring meadows example](https://gavinsimpson.github.io/physalia-multivariate/03-wednesday/spring-meadows.html)
    * [PERMANOVA](https://gavinsimpson.github.io/physalia-multivariate/05-friday/permanova.html)
* [Thursday](https://gavinsimpson.github.io/physalia-multivariate/04-thursday/slides.html)
    * [PDF](https://gavinsimpson.github.io/physalia-multivariate/04-thursday/slides.pdf)
    * [Randomised Block](https://gavinsimpson.github.io/physalia-multivariate/04-thursday/randomised-complete-block.html)
    * [Permutation test examples](https://gavinsimpson.github.io/physalia-multivariate/04-thursday/permutation-tests-solutions.html)
<!--
* [Friday](https://gavinsimpson.github.io/physalia-multivariate/05-friday/slides.html)
    * [PDF](https://gavinsimpson.github.io/physalia-multivariate/05-friday/slides.pdf)
    * [PERMANOVA](https://gavinsimpson.github.io/physalia-multivariate/05-friday/permanova.html)
    * [dbRDA](https://gavinsimpson.github.io/physalia-multivariate/05-friday/dbrda.html)
    * [PRC](https://gavinsimpson.github.io/physalia-multivariate/05-friday/prc.html)
    * [Co-CA](https://gavinsimpson.github.io/physalia-multivariate/05-friday/cocorrespondence-analysis.html)
    * [Variation partitioning](https://gavinsimpson.github.io/physalia-multivariate/05-friday/variation-partitioning.html)
    * [Linear discriminant analysis](https://gavinsimpson.github.io/physalia-multivariate/05-friday/linear-discriminants.html)
-->

## Target audience and assumed background

This course is suitable for PhD students (including senior thesis-based masters
students) and researchers working with multivariate data sets in biology
(inter alia ecology, animal science agriculture, microbial
ecology/microbiology), with limited statistical knowledge but a willingness to
learn more.

Participants should be familiar with RStudio and have some fluency in
programming R code, including being able to import, manipulate (e.g. modify
variables) and visualise data. There will be a mix of lectures, and hands-on
practical exercises throughout the course.

## Learning outcomes

1. Have a good introductory understanding of the main approaches used in the
   analysis of multivariate data sets,
2. Be able to choose an appropriate method to use to analyse a data set
3. Understand how to use restricted permutation tests with constrained
   ordination methods to test the effects of predictor variables or
   experimental treatments,
4. Be able to use the R statistical software to analyse multivariate data

## Pre-course preparation

### R

Please insure you have at version 4.4.x of R installed. Note that R and RStudio are two different things:
it is not sufficient to just update RStudio, you also need to update R by
installing new versions as they are released.

To download R go to the [CRAN Download](https://cran.r-project.org/) page and 
follow the links to download R for your operating system:

* [Windows](https://cran.r-project.org/bin/windows/)
* [MacOS X](https://cran.r-project.org/bin/macosx/)
* [Linux](https://cran.r-project.org/bin/linux/)

To check what version of R you have installed, you can run

```r
version
```

in R and look at the `version.string` entry (or the `major` and `minor`
entries).

### RStudio

We will use RStudio as the environment for interacting with R. You are free to
use your preferred IDE or R interface if you wish. To install RStudio, go to
<https://posit.co/download/rstudio-desktop/> and follow the instructions.

### Required R Packages

We will make use of several R packages that you'll need to have installed.
Prior to the start of the course, please run the following code to update
your installed packages and then install the required packages:

```r
# how many CPU cores are available?
cores <- parallel::detectCores()
# update any installed R packages
update.packages(ask = FALSE, checkBuilt = TRUE, Ncpus = cores - 1)

# packages to install
pkgs <- c("vegan", "tidyverse", "cocorresp", "permute", "here", "tibble",
          "readxl", "janitor", "stringr", "ggrepel", "e1071")

# install those packages
install.packages(pkgs, Ncpus = cores - 1)
```

Please also install my in-development *ggvegan* package as it will be useful
to know how to draw ordination diagrams using the *ggplot2* package. As
*ggvegan* is not yet on CRAN, we need to install from my R Universe:

```r
# Enable repository from gavinsimpson
options(repos = c(
  gavinsimpson = 'https://gavinsimpson.r-universe.dev',
  CRAN = 'https://cloud.r-project.org'))
# Download and install ggvegan in R
install.packages("ggvegan")
```

Setting the options like this only affects the current R session, which is all
we need to install *ggvegan*.

Ideally, we will also make use of the *DESeq2* package, which has to be
installed from the *Bioconductor* project:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
```

But don't worry if you can't get this installed, we won't be using it for any
of the longer computer examples, just to illustrate a couple of workflows for
handling high-throughput data.
