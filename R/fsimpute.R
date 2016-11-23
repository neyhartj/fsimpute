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
#' @import tidyr
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

  # Grab the founder and final names
  founder.names <- cross$geno$`1`$founders %>%
    row.names()
  final.names <- cross$geno$`1`$finals %>%
    row.names()


  ## First step to is to calculate conditional genotype probabilities

  cross.probs <- calc.genoprob(cross = cross, step = step, off.end = off.end,
                               error.prob = error.prob, map.function = map.function,
                               stepwidth = stepwidth)

  # Extract the probabilities
  probs <- lapply(X = cross.probs$geno, FUN = function(chr) chr$prob)
  
  finals.likely <- probs %>%
    map(function(chr) {
      
      # Create a list from the array - separates into the number of individuals
      array_tree(chr, margin = c(1,2)) %>%
        map(function(indiv) {
          # Determine the most likely parental state
          map_dbl(indiv, function(prob) ifelse(max(prob) >= prob.threshold, which.max(prob), NA))
        }) %>%
        do.call("rbind", .) })


  ## Once likely final genotypes are determined, the next step is to impute
  ## the founder genotypes.
  ## For each founder, assess the progeny to find the the
  ## genotype that most likely originated from the founder
  founders.imputed <- list(cross$geno, finals.likely) %>%
    pmap(function(cross.chr, fl.chr) {

      # Extract the original final genotypes
      fi.chr <- cross.chr$finals
      # Extract the original founder genotypes
      fou.chr <- cross.chr$founders
      # Combine the original genotypes
      orig.chr <- rbind(fou.chr, fi.chr)
      
      # Total number of individuals
      nInd <- nrow(orig.chr)
      
      # Combine a matrix of parent states and the likely parent states in the finals
      state.chr <- rbind(
        matrix(data = c(1,3), nrow = 2, ncol = ncol(orig.chr)),
        fl.chr
      )
      
      # Combine the original and the state matrices
      founders.likely <- rbind(state.chr, orig.chr) %>%
        # Apply a function over the markers
        apply(MARGIN = 2, FUN = function(snp) {
          
          # Separate the orig and state
          state <- snp[seq_len(nInd)]
          orig <- snp[-seq_len(nInd)]
          
          # Audit the progeny for those inheriting each parental state
          fou.likely <- list(state[1:2]) %>%
            pmap_dbl(function(par) {
              
              # Find the index of progeny inheriting that parental state
              finals.ind <- which(state[-1:-2] == par)
              # Which observed genotype most frequently corresponds to that state?
              tabl <- table(orig[-1:-2][finals.ind])
              most.freq <- which.max(tabl) %>%
                names() %>% as.numeric()
              
              ifelse(length(most.freq) == 0, NA, most.freq) }) %>%
            
            # Convert to matrix
            as.matrix()
          
          # Is one of the likely founders NA and are all of the likely finals not NA?
          ifelse( sum(is.na(fou.likely)) == 1 & sum(is.na(state[-1:-2])) == 0,
                yes = return(rep(na.omit(fou.likely), 2)), no = return(fou.likely) ) })


      # Add column and row names
      dimnames(founders.likely) <- dimnames(fou.chr)
      return(founders.likely) })

  ## Once founders are imputed to a higher density, we can use that information
  ## to impute the finals
  finals.imputed <- list(founders.imputed, finals.likely) %>%
    pmap(function(fou.imp.chr, fi.chr) {

      # Combine matricies
      chr <- rbind(fou.imp.chr, fi.chr)

      # Apply a function over the markers
      chr.new <- apply(X = chr, MARGIN = 2, FUN = function(snp) {

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
                                  no = ifelse(test = fin == 3, yes = fou[2], no = NA))) })

      # Add row names back
      row.names(chr.new) <- final.names
      return(chr.new) })
  
  # Add the imputed data back to the cross object
  cross$geno <- mapply(founders.imputed, finals.imputed, cross$geno, FUN = function(fou.chr, fi.chr, chr) {
    chr$founders_imputed <- fou.chr
    chr$finals_imputed <- fi.chr
    return(list(chr)) })
    
  ## Run some stats
  # Calculate missingness beforehand
  missing <- cross$geno %>%
    map_df(function(chr) {
      # Calculate overall missingness proportion for finals and founders simultaneously
      chr[c("founders", "finals", "founders_imputed", "finals_imputed")] %>%
        map_dbl(function(genos) mean(is.na(genos)) ) }) %>%
    mutate(group = rep(c("founders", "finals"), 2),
           type = rep(c("pre_imputation", "post_imputation"), each = 2 )) %>%
    gather(chrom, missingness, -group, -type) %>%
    spread(type, missingness) %>%
    select(group, chrom, pre_imputation, post_imputation) %>%
    arrange(chrom)
    
  # Add it back in to the cross object
  cross$missingness_stats <- missing

  return(cross)

} # Close the function
