#' Impute genotypes in bi-parental populations of finite selfing
#'
#' @description
#' Imputes missing genotype data using conditional genotype probabilities
#' derived from a Hidden Markov Model. Founder genotypes are then imputed
#' from the conditional genotype probabilities in the progeny. Then progeny
#' genotypes are imputed using the imputed founder genotypes.
#'
#'
#' @param cross An object of class \code{cross} from the package \code{qtl}. Genotype
#' data should only include those of the progeny in the family.
#' @param founders A list of founder genotypes, where each element in the list
#' is a \code{matrix} of genotypes on a single chromosome. Genotypes should be
#' coded as \code{z {0, 1, 2, NA}} where \code{z} is the number of reference alleles.
#' Column names should be marker names. This is the observed genotype data.
#' @param finals A list of progeny genotypes. The encoding and formatting should
#' be identical to the argument \code{founders}.
#' @param prob.threshold The minimum conditional genotype probability to declare
#' a site as having been inherited from a certain parent.
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
#' @import dplyr
#' @import stringr
#'
#' @export
#'
fsimpute <- function(cross, founders, finals, prob.threshold = 0.7, step = 0,
                     off.end = 0, error.prob = 1e-4,
                     map.function = c("haldane","kosambi","c-f","morgan"),
                     stepwidth=c("fixed", "variable", "max")) {

  ## Error handling

  # Check if cross is the correct format
  iscross <- inherits(cross, "cross")
  if (!iscross)
    stop("The argument 'cross' must be an object of class 'cross.'")

  # Check to make sure the finals and founders have the same number of markers
  founders.dim.check <- sapply(founders, dim)
  finals.dim.check <- sapply(finals, dim)
  cross.dim.check <- sapply(cross$geno, function(chr) dim(chr$data))


  if(!all.equal(founders.dim.check[2,], finals.dim.check[2,]))
    stop("The number of markers in each chromosome is not the same between
         the 'founders' and 'finals' arguments.")

  # Make sure founders has only two rows
  if (unique(founders.dim.check[1,]) != 2)
    stop("The number of entries in the 'founders' argument must be exactly 2.")


  # Now make sure that the number of markers is the same between finals,
  # founders, and genotypes in the cross object
  if (!all.equal(founders.dim.check[2,], cross.dim.check[2,]))
    stop("The number of markers in the 'founders' and 'finals' arguments is
         not the same as in the 'cross' object genotype data.")

  # Make sure the finals and cross geno data have the same dimension
  if (!identical(finals.dim.check, cross.dim.check))
    stop("The 'finals' argument does not have the same dimensions as the genotype
         data in the 'cross' argument.")


  # Check the founders and finals matrices for correct formatting
  founder.elements <- unlist(founders) %>%
    unique()
  final.elements <- unlist(finals) %>%
    unique()

  if (!founder.elements %in% c(0, 1, 2, NA))
    stop("The elements of the 'founders' argument must be in the set z {0, 1, 2, NA}.")

  if (!final.elements %in% c(0, 1, 2, NA))
    stop("The elements of the 'finals' argument must be in the set z {0, 1, 2, NA}.")


  # Check to make sure the prob.threshold is between 0 and 1
  if (prob.threshold < 0 | prob.threshold > 1)
    stop("The argument 'prob.threshold' must be between 0 and 1.")





  ## First step to is to calculate conditional genotype probabilities

  cross.probs <- calc.genoprob(cross = cross, step = step, off.end = off.end,
                               error.prob = error.prob, map.function = map.function,
                               stepwidth = stepwidth)

  # Extract the probabilities
  probs <- lapply(X = cross.probs$geno, FUN = function(chr) chr$prob)

  # Find the most likely genotype, given the threshold
  finals.likely <- lapply(X = probs, FUN = function(chr)
    # Apply over rows and columns of the array (i.e. look down the 3rd dimension)
    apply(X = chr, MARGIN = c(1,2), FUN = function(geno)
      ifelse(max(geno) >= prob.threshold, which.max(geno), NA) ))


  ## Once likely final genotypes are determined, the next step is to impute
  ## the founder genotypes For each founder, assess the progeny to find the the
  ## genotype that most likely originated from the founder

  founders.imputed <-
    mapply(founders, finals.likely, finals, FUN = function(fo.chr, fl.chr, fi.chr) {

      # Iterate over the number of markers
      for (i in seq_len(ncol(fo.chr))) {

        # Identify which founder is missing
        missing.founder <- which(is.na(fo.chr[,i]))

        # Are there missing founder? if not, skip
        if (length(missing.founder) < 1)
          next

        # Get the recoded founder identifier (i.e. 1 or 3)
        missing.founder.id <- c(1,3)[missing.founder]

        # Apply a function over the vector of missing founders
        founders.likely <- sapply(X = missing.founder.id, FUN = function(fou) {

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
          } })

        fo.chr[missing.founder,i] <- founders.likely

      }

      return(list(fo.chr)) })

  ## Once founders are imputed to a higher density, we can use that information
  ## to impute the finals
  finals.imputed <- mapply(founders.imputed, finals.likely, FUN = function(fi.chr, fl.chr) {

    # Iterate over the number of markers
    for (i in seq_len(ncol(fi.chr))) {

      # Extract parent genos
      pars <- fi.chr[,i]

      # Extract final genos
      fin <- fl.chr[,i]

      # Impute using founder genotypes - make sure NAs are not included
      fin[fin == 1 & !is.na(fin)] <- pars[1]
      fin[fin == 2 & !is.na(fin)] <- 1
      fin[fin == 3 & !is.na(fin)] <- pars[2]

      # Add back to matrix
      fl.chr[,i] <- fin

    }

    return(list(fl.chr)) })


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
