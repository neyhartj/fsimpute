#' Create a cross object
#'
#' @description
#' Takes raw founder and final genotypes and assembles a \code{cross}
#' object for downstream analysis. The genotypic data is first filtered
#' for unambiguous genotypes (see \code{Details}), then the genotypes are
#' recoded based on the unambiguous genotypes. Finally a \code{cross} object
#' is assembled and returned.
#'
#' @param map The genetic map, formatted as a \code{list} of chromosomes, where
#' each chromosome is a named vector of marker positions (in cM) and the names
#' are the marker names.
#' @param founders A list of founder genotypes, where each element in the list
#' is a \code{matrix} of genotypes on a single chromosome. Genotypes should be
#' coded as \code{z {0, 1, 2, NA}} where \code{z} is the number of reference alleles.
#' Column names should be marker names. This is the observed genotype data.
#' @param finals A list of progeny genotypes. The encoding and formatting should
#' be identical to the argument \code{founders}.
#' @param selfing.gen The number of selfing generations that the \code{finals}
#' have undergone. If the generation of the finals is \emph{F_t}, then the
#' argument \code{selfing.gen} would be \emph{t} - 1.
#'
#' @details
#' To force genotype data from a bi-parental family into a \code{cross} object
#' in \code{\link[qtl]{read.cross}}, the genotypes must be recoded into parental
#' states. Say two inbred parents have the observed gentypes \code{parent1 = 0} and
#' \code{parent2 = 2}, the parental states would be recoded as
#' \code{parent1 = 1} and \code{parent2 = 3}. Parent 1 is also given a parental
#' state of \code{1} and parent 2 is always given a parental state of \code{3}.
#' Among the progeny, any genotype call that is identical to the that of parent 1
#' would received a recoded genotype of \code{1}, and any gentype cal that is
#' identical to that of parent 2 would receive a recoded genotype of \code{3}.
#' Heterozygous genotype calls are recoded as \code{2}.
#'
#' Of course, in observed genotype data, the parental states are inherently
#' unknown (otherwise imputation would be easy). To determine parental states at
#' a marker, the marker must be unambiguous. That is, the parental states must
#' be easily inferred. To do this, we must only look at markers for which the
#' parents are both observed (i.e. no NA) and polymorphic between (i.e.
#' \code{0} and \code{2}). Any ambiguous markers are set to missing. This is ok,
#' because the imputation step is easily able to impute these genotypes.
#'
#' @return
#' An object of class \code{cross} with the normal elements \code{geno} and
#' \code{pheno}. The \code{geno} element is a list of chromosomes, each with
#' the elements:
#' \describe{
#' \item{data}{The recoded progeny genotypes}
#' \item{map}{The genetic map on that chromosome}
#' \item{founders}{The original founder genotypes}
#' \item{finals}{The original final genotypes}
#' \item{founders_unamb}{The unambiguous founder genotypes}
#' \item{finals_unamb}{The unambiguous final genotypes}
#' }
#'
#'
#' @import qtl
#' @import purrr
#'
#' @export
#'
prep_cross <- function(map, founders, finals, selfing.gen) {

  ## Error handling

  # Make sure the founders and finals are coded propery
  if (!all(unique(unlist(founders)) %in% c(0, 1, 2, NA)))
    stop("The gentypes in the argument 'founders' must be coded as
         z {0, 1, 2, NA}.")

  if (!all(unique(unlist(finals)) %in% c(0, 1, 2, NA)))
    stop("The gentypes in the argument 'founders' must be coded as
         z {0, 1, 2, NA}.")

  # Make sure the number of individuals in the founders and finals lists
  # are the same within those lists
  if (length(unique(sapply(founders, nrow))) != 1)
    stop("The number of individuals in the 'founders' arguments is not
         the same across elements of the list.")

  if (length(unique(sapply(finals, nrow))) != 1)
    stop("The number of individuals in the 'founders' arguments is not
         the same across elements of the list.")

  # First find the number of markers in the map, finals, and founders
  nMarkers.map <- sapply(map, length)

  nMarkers.founders <- sapply(founders, ncol)

  nMarkers.finals <- sapply(finals, ncol)

  # Are the number of chromosomes identical
  if (length(unique(c(length(nMarkers.map), length(nMarkers.founders), length(nMarkers.finals)))) != 1)
    stop("The number of chromosomes in the arguments 'map,' 'founders,' and 'finals'
         is not equal.")

  # Are the number of markers on each chromosome the same?
  same.markers <- mapply(nMarkers.map, nMarkers.founders, nMarkers.finals,
                         FUN = function(map.mar, fou.mar, fin.mar)
                           length(unique(c(map.mar, fou.mar, fin.mar))) == 1)

  if (!any(same.markers))
    stop("The number of markers in the arguments 'map,' 'founders,' and 'finals'
         is not equal.")


  # Find the unambiguous markers and set ambiguous markers to NA
  # Unambiguous markers are defined as those that are homozygous within a
  # parent but heterozygous between parents
  founders.unamb <- founders %>%
    map(function(chr) {

      # Apply over markers
      apply(X = chr, MARGIN = 2, FUN = function(snp) {

        # If one or both parents are missing, or both are the same, it is ambiguous
        if (any(is.na(snp)) | length(unique(snp)) == 1) {
          return(c(NA, NA))

          # If both parents are non-missing and different, it's good
        } else {
          return(snp)
        }}) })

  # Find the index of ambigous markers on each chromosome
  amb.markers <- founders.unamb %>%
    map(function(chr) {

      # Apply a function over a logical matrix
      apply(X = is.na(chr), MARGIN = 2, FUN = all) })

  # Using this index, set the ambigous markers in the finals to NA
  finals.unamb <- list(finals, amb.markers) %>%
    pmap(function(chr, mar.ind) {
      # Set NA
      chr[,mar.ind] <- NA
      # Return
      return(chr) })

  # Recode the final and founder genotypes as 1,2,3 format
  finals.recode <- list(founders.unamb, finals.unamb) %>%
    pmap(function(fou.chr, fi.chr) {
      # Combine the fou.chr and fi.chr for manipulation
      chr <- rbind(fou.chr, fi.chr)

      # Iterate over columns
      apply(X = chr, MARGIN = 2, FUN = function(snp) {
        # Positions 1 and 2 are parents
        fous <- snp[1:2]
        # The remainder are the finals
        fins <- snp[-1:-2]

        # Recode
        ifelse(test = fins == fous[1], yes = 1,
               no = ifelse(test = fins == 1, yes = 2,
                           no = ifelse(test = fins == fous[2], yes = 3, no = NA) ))
      }) })

  # Create the geno list for use in the cross object
  geno <- list(finals.recode, founders, finals, founders.unamb, finals.unamb, map) %>%
    pmap(function(geno.chr, fou.chr, fin.chr, fou.unamb.chr, fin.unamb.chr, map.chr) {

      # Return a list
      chr <- list(
        data = geno.chr,
        map = map.chr,
        founders = fou.chr,
        finals = fin.chr,
        founders_unamb = fou.unamb.chr,
        finals_unamb = fin.unamb.chr
      )

      # Change the list class
      class(chr) <- "A"
      # return
      return(chr)
      })

  # Create a nonsense cross object
  cross <- sim.cross(map = map, type = "bcsft", cross.scheme = c(0, selfing.gen + 1))

  # Assign the elements
  cross$geno <- geno
  cross$pheno <- data.frame(phenotype = rep(NA, nrow(finals$`1`)))

  # Return the cross
  return(cross)

} # Close the fuction






