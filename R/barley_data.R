#' Barley Genotype Data
#'
#' @description
#' A \code{matrix} of genotype data for barley breeding lines from
#' the University of Minnesota and North Dakota State University.
#'
#' @format A matrix with 764 rows and 1590 columns, where each row is a breeding
#' line and each column is a bi-allelic SNP marker. Markers are encoded as
#' z {0, 1, 2} where z is the number of "A" alleles.
#'
#' @source \url{https://triticeaetoolbox.org/barley/}
"barley_SNPs"


#' Genetic Map of Barley SNP Markers
#'
#' @description
#' A \code{list} of genetic map positions of 1590 bi-allelic SNP markers on
#' the 7 barley chromosomes. Can be used as a genetic map in the package
#' \code{qtl}.
#'
#' @format A list with 7 elements:
#' \describe{
#'   \item{Chromosome 1}{198 markers}
#'   \item{Chromosome 2}{259 markers}
#'   \item{Chromosome 3}{212 markers}
#'   \item{Chromosome 4}{219 markers}
#'   \item{Chromosome 5}{262 markers}
#'   \item{Chromosome 6}{216 markers}
#'   \item{Chromosome 7}{224 markers}
#' }
#'
#' @source \url{https://triticeaetoolbox.org/barley/}
"barley_map"
