## fsimpute simulation experiments
## For running on MSI

# Set the working directory
setwd("/home/smithkp/neyhartj/GBS/fsimpute/simulations")

# Load libraries
library(fsimpute)
library(qtl)
library(stringr)
library(tidyverse)
library(parallel)

# Detect cores
n.cores <- detectCores()

## De-novo simulations

# Size of each family
fam.size.list <- c(30, 100, 200)
# Maximum Proportion of missing marker data
missing.prob.list <- c(0.2, 0.4, 0.6, 0.8)
# Different generations of selfing
selfing.gen.list <- c(2, 4, 6)


## Non-varying parameters
# Number of chromosomes
n.chr <- 2
# Length of chromosomes
chr.len <- seq(150, by = 10, length.out = n.chr )
# Number of markers per chromosome
n.markers <- 1000
# Number of simulation replications
n.rep <- 50
# Heterozygosity in the parents
founder.het <- 0

# Genetic map
map <- sim.map(len = chr.len, n.mar = n.markers, include.x = F, eq.spacing = F)

# Build a data.frame for the simulation design
design <- expand.grid(
  rep = 1:n.rep,
  max_miss = missing.prob.list,
  fam_size = fam.size.list,
  self_gen = selfing.gen.list
) %>%
  select(self_gen, fam_size, max_miss, rep)

results <- design %>%
  # Add empty results columns
  mutate(founder_accuracy = NA,
         final_accuracy = NA,
         founder_missing = NA,
         final_missing = NA)

# Break up the design into n.core even chunks
results.split <- results %>%
  split(seq(n.cores))

## Run the simulations

# Parallelize the simulations
sim.results <- mclapply(X = results.split, FUN = function(df) {

  # Iterate over rows
  for (i in seq(nrow(df))) {

    # Simulate the family
    fam.cross <- sim.cross(map = map, n.ind = df[i,]$fam_size, type = "bcsft",
                           cross.scheme = c(0, (df[i,]$self_gen + 1)))

    # Generate founder genotypes
    # Adding residual heterozygosity
    founders.recode <- lapply(X = map, FUN = function(chr)
      sapply(X = chr, FUN = function(site)
        sample(c(0,1,2), size = 2, replace = T,
               prob = c(((1 - founder.het) / 2), founder.het, ((1 - founder.het) / 2)))) )

    # Extract the progeny genotypes from the cross
    finals <- lapply(X = fam.cross$geno, FUN = function(chr) chr$data)

    # Recode the finals to z {0, 1, 2} format
    finals.recode <- mapply(founders.recode, finals, FUN = function(p.chr, f.chr) {

      # Iterate over markers
      new.genos <- sapply(X = seq_len(ncol(p.chr)), FUN = function(i) {
        # Pull out parental genotypes
        pars <- p.chr[,i]
        names(pars) <- c(1,3)

        # If the parents are the same, the het is the same as the homs
        if (all(pars[1] == pars)) {
          het <- pars[1]
        } else {
          het <- 1
        }

        # Convert the genotypes
        ifelse(f.chr[,i] == 1, pars["1"],
               ifelse(f.chr[,i] == 2, het,
                      ifelse(f.chr[,i] == 3, pars["3"], NA))) })

      colnames(new.genos) <- colnames(p.chr)

      return(list(new.genos)) })

    # Combine matricies and simulate missing data
    missing.recode <- list(founders.recode, finals.recode) %>%
      pmap(function(fo.chr, fi.chr) {

        # Combine matrices
        chr <- rbind(fo.chr, fi.chr)

        # First generate the missing data proportions from a uniform distribution
        miss.prob <- runif(n = ncol(chr), min = 0, max = df[i,]$max_miss)

        miss.mat <- sapply(X = miss.prob, FUN = function(prob)
          sample(x = c(T, F), size = nrow(chr), replace = T, prob = c(prob, 1 - prob)))

        chr[miss.mat] <- NA

        fo.chr <- chr[1:2,]
        fi.chr <- chr[-1:-2,]

        list(founders = fo.chr, finals = fi.chr) })

    # Grab the founders and finals matrices
    founders.missing <- lapply(X = missing.recode, FUN = function(chr) chr$founders)
    finals.missing <- lapply(X = missing.recode, FUN = function(chr) chr$finals)

    # Create the cross object
    new.cross <- prep_cross(map = map, founders = founders.missing,
                            finals = finals.missing, selfing.gen = df[i,]$self_gen)


    # Impute
    impute.cross <- fsimpute(cross = new.cross, prob.threshold = 0.7)


    # Founder accuracy
    fi.combined <- impute.cross$geno %>%
      map(function(chr) chr$founders_imputed) %>%
      do.call("cbind", .)
    f.combined <- founders.recode %>%
      do.call("cbind", .)

    tab <- table(fi.combined, f.combined, dnn = list("imputed", "original"))

    fou_acc <- sum(diag(tab)) / sum(tab)

    # Finals accuracy
    fi.combined <- impute.cross$geno %>%
      map(function(chr) chr$finals_imputed) %>%
      do.call("cbind", .)
    f.combined <- do.call("cbind", finals.recode)

    tab <- table(fi.combined, f.combined, dnn = list("imputed", "original"))

    fin_acc <- sum(diag(tab)) / sum(tab)

    # Missingness
    fou_miss <- impute.cross$missingness_stats %>%
      filter(group == "founders") %>%
      summarize(mean = mean(post_imputation)) %>%
      unlist()
    fin_miss <- impute.cross$missingness_stats %>%
      filter(group == "finals") %>%
      summarize(mean = mean(post_imputation)) %>%
      unlist()

    df[i,]$founder_accuracy <- fou_acc
    df[i,]$final_accuracy <- fin_acc
    df[i,]$founder_missing <- fou_miss
    df[i,]$final_missing <- fin_miss

  }

  # Return the filled df
  return(df)

}, mc.cores = n.cores)

# Save the results
save.file <- "de_novo_simulation_results.RData"
save("sim.results", file = save.file)





