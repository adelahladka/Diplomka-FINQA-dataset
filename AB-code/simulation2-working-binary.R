#Simulation from ground up
rm(list = ls())
library(ltm)
library(difR)

set.seed(123)

load("simulation-abe.RData")

###################################################################################
#Hyperparameters

#N = 10 #number of simulations
#I = 5 #number of items (should change to I)
#n_tot = 1000 # number of respondents (700, 1000, 1200, 1500) 
#rat_n = c(1:1) #(c(5,2), c(1,1), c(5,1), c(2,1))
#------------------------------------------------------------------------------

# ============================
# Function: generate_group_sizes
# Purpose: Calculate group sizes for focal and reference groups given total N and ratio
# ============================
generate_group_sizes <- function(n_total, ratio = c(1, 1)) { # Reference: Focal
  if (length(ratio) != 2) stop("Ratio must be of length 2, like c(1, 3) or 1:3.")
  
  total_parts <- sum(ratio)
  n2 <- round(n_total * ratio[1] / total_parts)  # Size of focal group (first in ratio)
  n1 <- n_total - n2                             # Size of reference group (remaining subjects)
  
  return(list(n = n_total, n1 = n1, n2 = n2))    # Return total and both group sizes
}

# ============================
# Function: generate_param
# Purpose: Generate item parameters for 3PL model: b (difficulty), a (discrimination), c (guessing)
# ============================
generate_param <- function(I) {
  b <- rnorm(I, mean = 0, sd = 1)                  # Difficulty parameters
  z <- rnorm(I, mean = 0, sd = sqrt(0.1225))       # Latent variable for log-discrimination
  a <- exp(z)                                      # Discrimination (log-normal)
  c <- rep(0.2, I)                                 # Guessing parameter fixed at 0.2
  
  param <- cbind(b, a, c)
  rownames(param) <- paste0("Item", 1:I)
  colnames(param) <- c("b", "a", "c")
  return(param)
}

params <- generate_param(50)

# ============================
# Function: simulation_abe2
# Purpose: Run a DIF detection simulation with multiple methods
# ============================
simulation_abe2 <- function(N, n1, n2, skup, param1, param2, theta1, theta2, 
                            statistics = list(MantelB = TRUE, MantelNUB = TRUE, BresB = TRUE, LogB = TRUE, SIBB = TRUE, cSIBB = TRUE)) {
  
  # Initialize result storage
  Mantel <- MantelLow <- MantelHigh <- Bres <- Bres2 <- Log <- SIB <- cSIB <- rep(NA, N)
  MantelStat <- MantelLowStat <- MantelHighStat <- BresStat <- Bres2Stat <- LogStat <- SIBStat <- cSIBStat <- rep(NA, N)
  
  for (i in 1:N) {
    # Simulate response data for each group using IRT model
    data1 <- rmvlogis(n1, param1, IRT = TRUE, link = "logit", z.vals = matrix(theta1)) # Reference
    data2 <- rmvlogis(n2, param2, IRT = TRUE, link = "logit", z.vals = matrix(theta2)) # Focal
    data  <- rbind(data1, data2)
    
    # Compute total test score and split data by score
    total_scores <- rowSums(data)
    mean_score <- mean(total_scores)
    data_low  <- data[total_scores <= mean_score, ]
    skup_low  <- skup[total_scores <= mean_score]
    data_high <- data[total_scores > mean_score, ]
    skup_high <- skup[total_scores > mean_score]
    
    # Mantel-Haenszel test (full sample)
    if (statistics$MantelB | statistics$MantelNUB) {
      help <- try(difMH(data, skup, focal.name = 1, correct = TRUE, MHstat = "MHChisq", exact = FALSE)[[1]][1], silent = TRUE)
      MantelStat[i] <- help
      Mantel[i] <- if (is.numeric(help)) 1 - pchisq(help, 1) else NA
    }
    
    # Mantel-Haenszel test (low ability group)
    if (statistics$MantelNUB) {
      help <- try(difMH(data_low, skup_low, focal.name = 1, correct = TRUE, MHstat = "MHChisq", exact = FALSE)[[1]][1], silent = TRUE)
      MantelLowStat[i] <- help
      MantelLow[i] <- if (is.numeric(help)) 1 - pchisq(help, 1) else NA
    }
    
    # Mantel-Haenszel test (high ability group)
    if (statistics$MantelNUB) {
      help <- try(difMH(data_high, skup_high, focal.name = 1, correct = TRUE, MHstat = "MHChisq", exact = FALSE)[[1]][1], silent = TRUE)
      MantelHighStat[i] <- help
      MantelHigh[i] <- if (is.numeric(help)) 1 - pchisq(help, 1) else NA
    }
    
    # Breslow-Day test for homogeneity of odds ratios
    if (statistics$BresB) {
      help <- try(difBD(data, skup, focal.name = 1)[[1]][1, 3], silent = TRUE)
      BresStat[i] <- help
      Bres[i] <- if (is.numeric(help)) help else NA
    }
    
    # Breslow-Day trend test
    if (statistics$BresB) {
      help <- try(difBD(data, skup, focal.name = 1, BDstat = "trend")[[1]][1, 3], silent = TRUE)
      Bres2Stat[i] <- help
      Bres2[i] <- if (is.numeric(help)) help else NA
    }
    
    # Logistic regression DIF test
    if (statistics$LogB) {
      help <- try(difLogistic(data, skup, focal.name = 1)[[1]][1], silent = TRUE)
      LogStat[i] <- help
      Log[i] <- if (is.numeric(help)) 1 - pchisq(help, 2) else NA
    }
    
    # SIBTEST for uniform DIF
    if (statistics$SIBB) {
      help <- try(difSIBTEST(data, skup, focal.name = 1, type = 'udif', purify = FALSE)[[3]][1], silent = TRUE)
      SIBStat[i] <- help
      SIB[i]    <- if (is.numeric(help)) 1 - pchisq(help, 1) else NA
    }
    
    # SIBTEST for non-uniform DIF
    if (statistics$cSIBB) {
      help <- try(difSIBTEST(data, skup, focal.name = 1, type = 'nudif', purify = FALSE)[[3]][1], silent = TRUE)
      cSIBStat[i] <- help
      cSIB[i]    <- if (is.numeric(help)) 1 - pchisq(help, 1) else NA
    }
    
    # Print progress
    message(sprintf("Completed iteration %d at %s", i, date()))
  }
  
  # Combine results
  mat     <- rbind(Mantel, MantelLow, MantelHigh, Mantel, Bres, Bres2, Log, SIB, cSIB)
  matStat <- rbind(MantelStat, MantelLowStat, MantelHighStat, MantelStat, BresStat, Bres2Stat, LogStat, SIBStat, cSIBStat)
  
  rownames(mat)     <- c("MantelNormal", "MantelLow", "MantelHigh", "Mantel", "Bres", "Bres2", "Log", "SIB", "cSIB")
  rownames(matStat) <- c("MantelStatNormal", "MantelLowStat", "MantelHighStat", "MantelStat", "BresStat", "Bres2Stat", "LogStat", "SIBStat", "cSIBStat")
  
  return(list(pvalues = mat, statistics = matStat)) # Return p-values and statistics
}


# ============================
# Function: simul_total2
# Purpose: Wrapper to run multiple simulation scenarios (Type I, uniform DIF, non-uniform DIF)
# ============================
simul_total2 <- function(N, n_total, rat_n, I, mu_R, mu_F, type = "alpha", diffs = NULL,
                         statistics = list(MantelB = TRUE, MantelNUB = TRUE, BresB = TRUE, LogB = TRUE, SIBB = TRUE, cSIBB = TRUE),
                         timeMeasure = FALSE) {
  
  start_time <- Sys.time() # For timing
  
  # Generate group sizes and ability distributions
  list_ns <- generate_group_sizes(n_total, rat_n)
  n1 <- list_ns$n1
  n2 <- list_ns$n2
  
  theta1_ab <- rnorm(n1, mean = mu_R, sd = 1)
  theta2_ab <- rnorm(n2, mean = mu_F, sd = 1)
  
  # Generate item parameters and group labels
  param1_ab <- generate_param(I)
  param2_ab <- param1_ab
  skup <- c(rep(0, n1), rep(1, n2))
  
  results <- list()
  
  if (type == "alpha") {
    # Type I error — identical parameters
    results[["alpha"]] <- simulation_abe2(N,n1,n2,skup, param1_ab, param2_ab, theta1_ab, theta2_ab, statistics)
    
  } else if (type == "unif") {
    # Uniform DIF — increase difficulty of first item for focal group
    for (i in seq_along(diffs)) {
      d <- diffs[i]
      param2_ab_mod <- param1_ab
      param2_ab_mod[1, 1] <- param1_ab[1, 1] + d # increase b
      
      label <- paste0("Unif_DIF_delta", d)
      results[[label]] <- simulation_abe2(N,n1,n2,skup, param1_ab, param2_ab_mod, theta1_ab, theta2_ab, statistics)
    }
    
  } else if (type == "nonunif") {
    # Non-uniform DIF — modify discrimination of first item
    for (i in seq_along(diffs)) {
      delta <- diffs[i]
      param2_ab_mod <- param1_ab
      a_orig <- param1_ab[1, 2]
      c_orig <- param1_ab[1, 3]
      
      # Formula to adjust 'a' for DIF
      param2_ab_mod[1, 2] <- 2 * a_orig / (2 + delta * a_orig / ((1 - c_orig) * log(2)))
      
      label <- paste0("NonUnif_DIF_delta", delta)
      results[[label]] <- simulation_abe2(N,n1,n2,skup, param1_ab, param2_ab_mod, theta1_ab, theta2_ab, statistics)
    }
  }
  
  end_time <- Sys.time()
  time_taken <- end_time - start_time
  #print(end_time, start_time)
  message(sprintf("Total time for simulation: %s", time_taken))
  
  if (timeMeasure) return(time_taken)
  return(results)
}
#################################################################################
#Evaluating------------------------------------------------------------

# Type I error
res_alpha <- simul_total2(N = 1000, n_total = 500, rat_n = c(1,2), I = 40, mu_R = 0, mu_F = -1, type = "alpha",statistics = list(SIBB = TRUE, MantelB = TRUE, MantelNUB = TRUE, BresB = TRUE, LogB = TRUE, cSIBB = TRUE)
)

# Uniform DIF
res_unif <- simul_total2(N = 1000, n_total = 500, rat_n = c(1,2), I = 5, mu_R = 0, mu_F = -1, 
                        type = "unif", diffs = c(0.5,1,2,4))

# Non-uniform DIF
res_nonunif <- simul_total2(N = 1000, n_total = 500, rat_n = c(1,2), I = 5, mu_R = 0, mu_F = -1, 
                           type = "nonunif", diffs = c(0.4,0.6,0.8,1))


###################################################################################
# Time measuring function: Runs each DIF detection method separately to measure runtime
run_all_statistics <- function() {
  # Define the names of the DIF methods to test
  stat_names <- c("SIBB", "MantelB", "MantelNUB", "BresB", "LogB", "cSIBB")
  
  # Initialize a numeric vector to store time results
  results <- numeric(length(stat_names))
  
  # Loop through each statistical method
  for (i in seq_along(stat_names)) {
    # Create a list of all methods set to FALSE
    stats_list <- setNames(as.list(rep(FALSE, length(stat_names))), stat_names)
    
    # Activate only the current method
    stats_list[[stat_names[i]]] <- TRUE
    
    # Run the simulation with the selected method and record time
    result <- simul_total2(
      N = 100,             # Number of replications
      n_total = 500,        # Total sample size
      rat_n = c(1, 2),      # Group size ratio (focal:reference)
      I = 40,               # Number of items
      mu_R = 0,             # Reference group ability mean
      mu_F = -1,            # Focal group ability mean
      type = "alpha",       # Type I error condition (no DIF)
      statistics = stats_list, # Use only the current method
      timeMeasure = TRUE    # Return time taken
    )
    
    # Store the result (time taken)
    results[i] <- result
  }
  
  # Combine method names and times into a data frame
  df_results <- data.frame(
    Method = stat_names,
    Value = results
  )
  
  return(df_results)
}
#Evaluation ---------------------------------------------------------------------
df_results <- run_all_statistics()


#################################################################################
#################################################################################
#Estimate of Type I error

calculate_alpha_estimates <- function(res_alpha) {
  # Calculate est_alpha (proportion of p-values <= 0.01 for rows 1:3)
  est_alpha <- apply(res_alpha$alpha$pvalues[1:3,], 1, function(x) {
    sum(x <= 0.01, na.rm = TRUE) / sum(!is.na(x))
  })
  
  # Calculate est_alpha2 (proportion of p-values <= 0.05 for rows 4:9)
  est_alpha2 <- apply(res_alpha$alpha$pvalues[4:9,], 1, function(x) {
    sum(x <= 0.05, na.rm = TRUE) / sum(!is.na(x))
  })
  
  # Combine both into a single vector (est_alpha_total)
  est_alpha_total <- c(est_alpha, est_alpha2)
  
  # Calculate the number of NAs in each row
  alphaNA <- apply(res_alpha$alpha$pvalues, 1, function(x) { sum(is.na(x)) })
  
  # Calculate confidence intervals for each estimate
  confintAlpha <- rep(NA, length(est_alpha_total))
  for (i in 1:length(est_alpha_total)) {
    confintAlpha[i] <- paste(
      "(", 
      round(binom.test(round((1000 - alphaNA[i]) * est_alpha_total[i]), n = (1000 - alphaNA[i]))$conf.int[1], 3), ", ",
      round(binom.test(round((1000 - alphaNA[i]) * est_alpha_total[i]), n = (1000 - alphaNA[i]))$conf.int[2], 3), 
      ")", sep = ""
    )
  }
  # Return the results as a list
  return(list(est_alpha_total = est_alpha_total,
              alphaNA = alphaNA, 
              confintAlpha = confintAlpha))
}


est_alpha_result<- calculate_alpha_estimates(res_alpha)
est_alpha_frame <- data.frame(
  est_alpha_total = est_alpha_result$est_alpha_total, 
  alphaNA = est_alpha_result$alphaNA, 
  confintAlpha = est_alpha_result$confintAlpha
)

# Set the row names to match the names in alphaNA
rownames(est_alpha_frame) <- names(est_alpha_result$alphaNA)
print(est_alpha_frame)


################################################################################
#Estimation of power of the test

calculate_power_and_NA <- function(data, diffs = c(0.4, 0.6, 0.8, 1), type = c("nonunif", "unif")) {
  
  # Ensure 'type' is either 'nonunif' or 'unif'
  if (!type %in% c("nonunif", "unif")) {
    stop("Invalid type. Please specify either 'nonunif' or 'unif'.")
  }
  
  # Get the prefix based on 'type' argument
  prefix <- ifelse(type == "nonunif", "NonUnif_DIF_delta", "Unif_DIF_delta")
  
  # Get the number of rows from the first pvalues matrix (assuming all have the same structure)
  num_rows <- nrow(data[[paste0(prefix, diffs[1])]]$pvalues)
  
  # Initialize the result matrices for both sets of rows
  est_power <- matrix(NA, nrow = num_rows, ncol = length(diffs)) # Dynamic number of rows
  NAValues <- matrix(NA, nrow = num_rows, ncol = length(diffs)) # Dynamic number of rows
  
  # Row labels (methods names)
  method_labels <- rownames(data[[paste0(prefix, diffs[1])]]$pvalues)
  
  # Loop through each delta in diffs
  for (i in 1:length(diffs)) {
    delta <- diffs[i]
    
    # Access the pvalues for the current delta
    pvalues <- data[[paste0(prefix, delta)]]$pvalues
    
    # Calculate est_power (proportion of p-values <= 0.01 for rows 1:3, <= 0.05 for rows 4:9)
    est_power[1:3, i] <- apply(pvalues[1:3, ], 1, function(x) { mean(x <= 0.01, na.rm = TRUE) })
    est_power[4:9, i] <- apply(pvalues[4:9, ], 1, function(x) { mean(x <= 0.05, na.rm = TRUE) })
    
    # Calculate NA400 (count of NA values)
    NAValues[1:3, i] <- apply(pvalues[1:3, ], 1, function(x) { sum(is.na(x)) })
    NAValues[4:9, i] <- apply(pvalues[4:9, ], 1, function(x) { sum(is.na(x)) })
  }
  
  # Add row labels
  rownames(est_power) <- method_labels
  rownames(NAValues) <- method_labels
  
  # Add column labels (diffs values)
  colnames(est_power) <- paste0("Delta_", diffs)
  colnames(NAValues) <- paste0("Delta_", diffs)
  
  # Return the two matrices
  return(list(est_power = est_power, NAValues = NAValues))
}


diffsNU <- c(0.4, 0.6, 0.8, 1)  # Example diffs
diffsU <- c(0.5,1,2,4)
est_power_nonufif <- calculate_power_and_NA(res_nonunif, diffsNU, type = "nonunif")
est_power_unif <- calculate_power_and_NA(res_unif, diffsU, type = "unif")

################################################################################
################################################################################
#Parameter testing


# Define the I values
I_values <- c(5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80)
n_total_value <- c(100,200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1400)
# Prepare lists to store the results
table12_list <- list()
table22_list <- list()

# Loop through I values
for (n_val in n_total_value) {
  print(n_val)
  res_alpha2 <- simul_total2(
    N = 1000,
    n_total = n_val,
    rat_n = c(1, 2),
    I = 30,
    mu_R = 0,
    mu_F = -1,
    type = "alpha",
    statistics = list(
      SIBB = TRUE,
      MantelB = TRUE,
      MantelNUB = TRUE,
      BresB = TRUE,
      LogB = TRUE,
      cSIBB = TRUE
    )
  )
est_alpha_result2 <- calculate_alpha_estimates(res_alpha2)
  
  # Add current I to the tables
table12 <- est_alpha_result2$est_alpha_total
table22 <- est_alpha_result2$alphaNA
  
table12$n <- n_val
table22$n <- n_val
  
  table12_list[[as.character(n_val)]] <- table12
  table22_list[[as.character(n_val)]] <- table22
}

# Combine all into final tables
final_table12 <- do.call(rbind, table12_list)
final_table22 <- do.call(rbind, table22_list)

# View results
print(final_table1) # I values variable, n = 500 (I think)
print(final_table2)

#Total time for simulation: 1.04047406858868

print(final_table12) # n values variable, I = 30
print(final_table22)
#Total time for simulations 45.715


#################################################################################
save.image(file = "simulation-abe.RData")



############################################################################
#Plotting

library(dplyr)
library(ggplot2)
library(tidyr)

#Type I error

est_alpha_frame2 <- est_alpha_frame %>%
  tibble::rownames_to_column("Method") %>%
  separate(confintAlpha, into = c("Lower", "Upper"), sep = ", ", remove = FALSE) %>%
  mutate(
    Lower = as.numeric(gsub("[()]", "", Lower)),
    Upper = as.numeric(gsub("[()]", "", Upper))
  )

ggplot(est_alpha_frame2, aes(x = Method, y = est_alpha_total)) +
  geom_bar(stat = "identity", fill = "steelblue", width = 0.6) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, color = "black") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", size = 0.8) +
  labs(title = "Estimated Type I Error with 95% Confidence Intervals",
       x = "Method",
       y = "Estimated Type I Error Rate") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#################################################################################
est_powerUnif <- as.data.frame(est_power_unif$est_power)

est_power_long_unif <- est_powerUnif %>%
  tibble::rownames_to_column("Method") %>%      # move rownames to a column
  pivot_longer(cols = starts_with("Delta"),     # reshape to long format
               names_to = "Delta", 
               values_to = "Power")


ggplot(est_power_long_unif, aes(x = Delta, y = Power, color = Method, group = Method)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(title = "Power by Method Across Different Effect Sizes",
       x = "Effect Size", y = "Estimated Power") +
  theme_minimal(base_size = 14) +
  scale_color_brewer(palette = "Set2") +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "grey30")



est_powerNONUnif <- as.data.frame(est_power_nonufif$est_power)



est_power_long_nonunif <- est_powerNONUnif %>%
  tibble::rownames_to_column("Method") %>%  # Add method as a column
  pivot_longer(cols = starts_with("Delta"), 
               names_to = "Delta", 
               values_to = "Power") %>%
  filter(!is.na(Power))  # Remove missing values

# Plot the data
ggplot(est_power_long_nonunif, aes(x = Delta, y = Power, color = Method, group = Method)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(title = "Power by Method Across Different Effect Sizes",
       x = "Effect Size (Delta)", y = "Estimated Power") +
  theme_minimal(base_size = 14) +
  scale_color_brewer(palette = "Set3") +  # Use a palette that supports more than 8 colors
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "grey30")  # Optional: Add a threshold line for power = 0.8


##############################################################
final_table1$I <- as.numeric(final_table1$I)
final_table1F <- as.data.frame(final_table1)
# Reshape the data from wide format to long format
final_table1_long <- final_table1F %>%
  pivot_longer(cols = -I, names_to = "Method", values_to = "Type_I_Error")  # Reshape the data

# Plot the data
ggplot(final_table1_long, aes(x = I, y = Type_I_Error, color = Method, group = Method)) +
  geom_line(size = 1) +  # Line for each method
  geom_point(size = 2) +  # Points for each method
  labs(title = "Type I Error by Number of Items (I)",
       x = "Number of Items (I)", y = "Type I Error") +
  theme_minimal(base_size = 14) +
  scale_color_brewer(palette = "Set3") +  # Use a color palette that supports many categories
  theme(legend.position = "bottom")  # Position the legend at the bottom

