#' Impute genotypes in bi-parental populations of finite selfing
#'
#' @description
#' Imputes missing genotype data using conditional genotype probabilities
#' derived from a Hidden Markov Model. Founder genotypes are then imputed
#' from the conditional genotype probabilities in the progeny. Then progeny
#' genotypes are imputed using the imputed founder genotypes.
#'
#'
#' @param cross An object of class \code{cross} from the package
#' \code{\link[qtl]{read.cross}}. This object should be generated using the
#' function \code{\link[fsimpute]{prep_cross}}.
#' @param prob.threshold The minimum conditional genotype probability to declare
#' a site as having been inherited from a certain parent.
#' @param type The genotype output type. Can be \code{discrete} to code each
#' imputed genotype as the most likely genotype, or \code{continuous} to code
#' each imputed genotype as a weighted average based on the genotype probabilities.
#' @param step See \code{\link[qtl]{calc.genoprob}}.
#' @param off.end See \code{\link[qtl]{calc.genoprob}}.
#' @param error.prob See \code{\link[qtl]{calc.genoprob}}.
#' @param map.functions See \code{\link[qtl]{calc.genoprob}}.
#' @param stepwidth See \code{\link[qtl]{calc.genoprob}}.
#'
#' @return
#' The \code{cross} object is returned with the additional item \code{imputed},
#' which contains the following information:
#' \describe{
#'   \item{$stats}{Some measures on the imputed genotypes, including the proportion of missing data}
#'   \item{$founders}{The imputed founder genotypes}
#'   \item{$finals}{The imputed progeny genotypes}
#'   }
#'
#'
#' @import qtl
#' @import purrr
#'
#' @export
#'
fsimpute <- function(cross, prob.threshold = 0.7, type = "discrete", step = 0, off.end = 0,
                     error.prob = 1e-4, map.function = c("haldane","kosambi","c-f","morgan"),
                     stepwidth=c("fixed", "variable", "max")) {

  ## Error handling

  # Check if cross is the correct format
  iscross <- inherits(cross, "cross")
  if (!iscross)
    stop("The argument 'cross' must be an object of class 'cross.'")

  # Check to make sure the prob.threshold is between 0 and 1
  if (prob.threshold < 0 | prob.threshold > 1)
    stop("The argument 'prob.threshold' must be between 0 and 1.")

  # Check to see if the type is a correct input
  if (!type %in% c("discrete", "continuous"))
    stop("The argument 'type' must be 'discrete' or 'continuous.'")


  ## First step to is to calculate conditional genotype probabilities

  cross.probs <- calc.genoprob(cross = cross, step = step, off.end = off.end,
                               error.prob = error.prob, map.function = map.function,
                               stepwidth = stepwidth)

  # Extract the probabilities
  probs <- lapply(X = cross.probs$geno, FUN = function(chr) chr$prob)

  # Find the most likely parental genotype, given the threshold
  finals.likely <- lapply(X = probs, FUN = function(chr)
    # Apply over rows and columns of the array (i.e. look down the 3rd dimension)
    apply(X = chr, MARGIN = c(1,2), FUN = function(geno)
      ifelse(max(geno) >= prob.threshold, which.max(geno), NA) ))


  ## Once likely final genotypes are determined, the next step is to impute
  ## the founder genotypes For each founder, assess the progeny to find the the
  ## genotype that most likely originated from the founder
  founders.imputed <- list(cross$geno, finals.likely) %>%
    pmap(function(cross.chr, fl.chr) {

      # Extract the original final genotypes
      fi.chr <- cross.chr$finals

      # Apply a function over the markers
      founders.likely <- sapply(X = seq_len(ncol(fi.chr)), FUN = function(i)

        # Identify the index of progeny inheriting each parental allele
        sapply(X = c(1,3), FUN = function(fou) {
          # Identify the index of progeny inheriting that founder allele
          finals.ind <- fl.chr[,i] == fou
          # Count the alleles inherited from that parent
          tabl <- fi.chr[finals.ind,i] %>%
            table()
          # Assign the most frequent genotype
          freq.genotype <- which.max(tabl) %>%
            names() %>%
            as.numeric()
          # The length may be zero, if so assign NA
          if (length(freq.genotype) == 0) {
            return(NA)

          } else {
            return(freq.genotype)
          } }) )

      # Add column names
      colnames(founders.likely) <- colnames(fi.chr)
      return(founders.likely) })

  ## Once founders are imputed to a higher density, we can use that information
  ## to impute the finals
  finals.imputed <- list(founders.imputed, finals.likely) %>%
    pmap(function(fou.imp.chr, fi.chr) {

      # Combine matricies
      chr <- rbind(fou.imp.chr, fi.chr)

      # Apply a function over the markers
      apply(X = chr, MARGIN = 2, FUN = function(snp) {

        # The parents are the first two
        fou <- snp[1:2]
        # The finals are the remaining
        fin <- snp[-1:-2]

        # If the founders are identical, the het is homozygous for the founder
        if (length(unique(fou)) == 1 ) {
          het <- fou[1]
        } else {
          het <- 1
        }

        # Impute using founder genotypes - make sure NAs are not included
        ifelse(test = fin == 1, yes = fou[1],
                      no = ifelse(test = fin == 2, yes = het,
                                  no = ifelse(test = fin == 3, yes = fou[2], no = NA))) }) })


  ## Run some stats
  # Calculate missingness proportion
  finals.missing <- sapply(X = finals.imputed, FUN = function(chr) mean(is.na(chr)))
  finals.missing["mean"] = mean(finals.missing)

  founders.missing <- sapply(X = founders.imputed, FUN = function(chr) mean(is.na(chr)))
  founders.missing["mean"] = mean(founders.missing)

  # Add the imputed data back to the cross object
  cross$imputed <- list(
    stats = list(missingness = rbind(founders.missing, finals.missing)),
    founders = founders.imputed,
    finals = finals.imputed
  )

  return(cross)

} # Close the function
