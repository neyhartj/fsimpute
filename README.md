# fsimpute
An R package for imputing sparse genotyping data in bi-parental populations
of partially inbred individuals

## Author
[Jeff Neyhart](neyhartje@gmail.com)

## Introduction

Imputation software abounds for low-coverage genotype data generated from 
platforms such as genotyping-by-sequencing, particularly in the field
of plant breeding and genetics. However, there are not many easy-to-use
imputation software packages that can handle populations that have undergone
finite generations of inbreeding.

`fsimpute` (**f**ull **s**ibling / **f**inite **s**elfing impute) is a simple package 
for imputation of low-coverage genotype data in bi-parental populations of 
partially inbred individuals. The software can handle bi-allelic SNP loci.
The methodology is based partly on functions implemented in the 
[mpMap](https://github.com/behuang/mpMap) package and their imputation procedure 
for multi-parent populations (Huang et al. 2014).

## Simulations

To test the efficacy of `fsimpute`, I ran a few simulations using fake and real
data. The fake data was generated through `R/qtl` and the real data consisted
of genotype data on barley breeding lines. This information was downloaded from
the [T3](https://triticeaetoolbox.org/barley/) database.

I tested 3 family sizes (30, 100, and 200 individuals), 3 levels of inbreeding
(2 selfing generations = $F_3$, 4 selfing generations = $F_5$, and 6 selfing
generations = $F_7$), and 4 per-SNP missing marker cutoffs (20%, 40%, 60%, 
and 80%). 50 replications of each scenario were performed.

### Fake Data Simulations



## Procedure

You will need three pieces of information to use `fsimpute`. 1) genotype data on 
parents, 2) genotype data on the progeny, and 3) a genetic map of the markers.


## Installation
The package can be installed via the `devtools` package:
```
devtools::install_github("neyhartj/fsimpute")
```


