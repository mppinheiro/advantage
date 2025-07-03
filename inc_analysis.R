#############################################################
# Estimating true producer advantage in a hunter population #
#############################################################


# Initialize and define functions -----------------------------------------

# Load packages and set theme
if (!require("pacman")) install.packages("pacman")
pacman::p_load(lhs, tidyr, dplyr, ggplot2, moments, rstatix, boot, ggpubr, 
               forcats, overlapping)
theme_set(theme_bw())
# Set up random seed for reproducibility
set.seed(-1556080447)

# Function to generate a list of random matrices; where weights are household weights, 
# diag_values is the fixed diagonal and num_samples the number of matrices to be sampled
generate_random_matrices <- function(weights, diag_values, num_samples) {
  n <- length(diag_values) # matrix size
  sampled_matrices <- vector("list", num_samples)  # store matrices
  for (i in seq_len(num_samples)) {
    mat <- matrix(0, n, n) # create empty matrix
    for (j in 1:n) {
      # Obtain total sum of off-diagonal elements from row j
      remaining_sum <- 1 - diag_values[j]
      # Generate random sample for n-1 variables
      random_sample <- runif(n - 1)
      # Introduce weights excluding the diagonal element
      biased_sample <- random_sample * weights[-j]
      # Adjust sampled values to sum to the remaining_sum
      scaled_values <- biased_sample / sum(biased_sample) * remaining_sum
      # Fill row j with scaled values in off-diagonal positions
      col_idx <- setdiff(1:n, j)
      mat[j, col_idx] <- scaled_values
      mat[j, j] <- diag_values[j]
    }
    sampled_matrices[[i]] <- mat
  }
  return(sampled_matrices)
}

# Function to calculate producer advantage index; where x is an income distribution, 
# v the acquisition vector and k the acquisition rate of the focal producer
pad_fun <- function(x, v, k) {
  set.seed(-1596028447)
  prod <- mean(x[which(v==k)])
  other <- mean(x[which(v!=k)])
  res <- (prod - other) / min(prod, other) # positive advantage, negative disadvantage
  return(round(res, 4))
}

# Statistic functions for calculating CIs
f_med <- function(data, indices) {  # median
  return(median(data[indices]))
}
f_p05 <- function(data, indices) {  # 5th percentile
  return(quantile(data[indices], 0.05))
}
f_p95 <- function(data, indices) {  # 95th percentile
  return(quantile(data[indices], 0.95))
}

# Define production and self-sharing sets ---------------------------------

# Define acquisition rates for each producer type
at = 16818  # acquisition rate of top married hunters
ar = 2735  # acquisition rate of regular married hunters
ap = 1367  # acquisition rate of poor married hunters

# Generate the production vector

# Create a set of n-sized vectors with all possible combinations of 
# acquisition rates; where n is the desired population size
v <- c(at, ar, ap)
pv <- expand.grid(v, v, v, v, v, v, v) # for a population of 7 married hunters
# Sort in descending order and eliminate duplicates 
pv <- t(apply(pv, 1, FUN = function(x) sort(x, decreasing = TRUE)))
pv <- pv[!duplicated(pv),]
# Apply conditions for right skewed distribution of hunter quality
pv <-
  pv[sapply(seq_len(nrow(pv)), function(x) length(which(pv[x,] == at)) > 0) &
       sapply(seq_len(nrow(pv)), function(x) length(which(pv[x,] == ar)) > 0) &
       sapply(seq_len(nrow(pv)), function(x) length(which(pv[x,] == ap)) > 0) &
       sapply(seq_len(nrow(pv)), function(x) length(which(pv[x,] == at)) <= 
                length(which(pv[x,] == ar))) &
       sapply(seq_len(nrow(pv)), function(x) length(which(pv[x,] == at)) <= 
                length(which(pv[x,] == ap)))
     ,]

# Select the vector that most closely approaches the
# average acquisition rate of the hunter population
# a = 4447 # average adult male acquisition rate
a = 5187 # average married male acquisition rate
pv <- pv[which(abs(rowMeans(pv) - a) == min(abs(rowMeans(pv) - a))),]

# Define self-sharing rates for each producer type
st = round(0.2, 2)  # self-sharing rate of top married hunters
sr = round(0.44, 2)  # self-sharing rate of regular married hunters
sp = round(0.59, 2) # self-sharing rate of poor married hunters

# Generate the self-sharing vector from the production vector
sv <- pv
sv[sv == at] <- st
sv[sv == ar] <- sr
sv[sv == ap] <- sp


# Generate sharing set ----------------------------------------------------

# Inputs for the matrix function
weights <- c(rep(1,length(pv))) # household weights
diag_values <- sv  # self-sharing vector as the diagonal
num_samples <- 10000  # number of matrices to sample

# Generate the sharing set
share <- generate_random_matrices(weights, diag_values, num_samples)

# Output
 cat("Generated", length(share), 
     "random matrices in the sharing set.\nThis is the first sharing matrix:\n")
 print(share[[1]])  # Show the first generated matrix


# Derive income set -------------------------------------------------------

 # Calculate income set
 income <- lapply(share, FUN = function(x) as.vector(t(x) %*% pv))
 
 # Calculate production advantage indexes
 pad <- sapply(income, FUN = function(x) pad_fun(x,pv,at)) # based on the income set
 psa <- pad_fun(sv*pv,pv,at)  # based on post-sharing acquisition rates


# Calculate CIs -----------------------------------------------------------

 # Perform bootstrapping (10,000 resamples)
 boot_med <- boot(pad, statistic = f_med, R = 10000)
 boot_p05 <- boot(pad, statistic = f_p05, R = 10000)
 boot_p95 <- boot(pad, statistic = f_p95, R = 10000)
 
 # collect the selected statistic
 ci_med <- round(boot.ci(boot_med, type = "perc")$percent[4:5], 4)
 ci_p05 <- round(boot.ci(boot_p05, type = "perc")$percent[4:5], 4)
 ci_p95 <- round(boot.ci(boot_p95, type = "perc")$percent[4:5], 4)
 

# Save output -------------------------------------------------------------

# choose the appropriate path
save(pad, psa, ci_med, ci_p05, ci_p95,
     file = "/Users/admin/Desktop/Textos/2024/Vantagem do produtor/data.Rdata")

 
# Expand to full household population -------------------------------------

 # define acquisition and self-sharing rates for additional producers
 aw <- 3000 # acquisition rate of adult women
 ay <- 2000 # acquisition rate of young women
 as <- 1854 # acquisition rate of unmarried hunters
 sw <- 0.71 # self-sharing rate of women
 ss <- 0.54 # self-sharing rate of unmarried hunters

 # adjust acquisition and self-sharing rates of married couple households
 # based on women's acquisition and self-sharing rates
 married <- seq_along(pv) # select married hunters
 pv_full <- pv + aw
 sv_full <- sapply(married, function(x)
   round((1 - aw / pv_full[x]) * sv[x] + (aw / pv_full[x]) * sw, 2))

 # Add custom households to the set of vectors
 # (2 unmarried men, 1 elderly woman, 2 unmarried women)
 pv_full <- c(pv_full, 2*as, 1*aw, 2*ay)
 sv_full <- c(sv_full, ss, sw, sw)


 # Generate sharing set ----------------------------------------------------

 # Inputs for the matrix function
 weights <- c(rep(3,length(pv)), 2, 1, 2) # household weights
 diag_values <- sv_full  # self-sharing vector as the diagonal
 num_samples <- 10000  # number of matrices to sample

 # Generate sharing set
 share_full <- generate_random_matrices(weights, diag_values, num_samples)

 # Output
 cat("Generated", length(share_full),
     "random matrices in the sharing set.\nThis is the first sharing matrix:\n")
 print(share_full[[1]])  # Show the first generated matrix


 # Derive income set -------------------------------------------------------

 # Calculate income set
 income_full <- lapply(share_full, FUN = function(x) as.vector(t(x) %*% pv_full))

 # Calculate production advantage indexes
 # based on the income set
 pad_full <- sapply(income_full, FUN = function(x) pad_fun(x,pv_full[married],at+aw))
 # based on post-sharing acquisition rates
 psa_full <- pad_fun(sv_full*pv_full,pv_full[married], at+aw)


 # Calculate CIs -----------------------------------------------------------

 # Perform bootstrapping (10,000 resamples)
 boot_med_full <- boot(pad_full, statistic = f_med, R = 10000)
 boot_p05_full <- boot(pad_full, statistic = f_p05, R = 10000)
 boot_p95_full <- boot(pad_full, statistic = f_p95, R = 10000)

 # collect the selected statistic
 ci_med_full <- round(boot.ci(boot_med_full, type = "perc")$percent[4:5], 4)
 ci_p05_full <- round(boot.ci(boot_p05_full, type = "perc")$percent[4:5], 4)
 ci_p95_full <- round(boot.ci(boot_p95_full, type = "perc")$percent[4:5], 4)


 # Save output -------------------------------------------------------------

 # choose the appropriate path
 save(pad_full, psa_full, ci_med_full, ci_p05_full, ci_p95_full,
      file = "/Users/admin/Desktop/Textos/2024/Vantagem do produtor/data_full.Rdata")
 

# Additional statistics ---------------------------------------------------

 # t <- t.test(pad_full, pad) # t-test
 # boot_over <- boot.overlap(list(pad,pad_full), B = 1000 ) # bootstrap overlap 
 # over <- round(boot_over$OVboot_stats$estOV, 4) # overlap estimate
 # ci_over <- c(round(quantile(boot_over$OVboot_dist, 0.05)[[1]],4), # overlap CI
 #              round(quantile(boot_over$OVboot_dist, 0.95)[[1]],4))
 # 
 # save(t, over, ci_over,
 #      file = "/Users/admin/Desktop/Textos/2024/Vantagem do produtor/stats.Rdata")
 