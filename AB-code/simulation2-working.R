#Simulation from ground up
library(ltm)
library(difR)
rm(list = ls())
###################################################################################
#Hyperparameters

#N = 10 #number of simulations
#I = 5 #number of items (should change to I)
#n_tot = 1000 # number of respondents (700, 1000, 1200, 1500) 
#rat_n = c(1:1) #(c(5,2), c(1,1), c(5,1), c(2,1))
set.seed(123)
#------------------------------------------------------------------------------
#Functions by piece
generate_group_sizes <- function(n_total, ratio = c(1, 1)) {#Reference: Focal
  if (length(ratio) != 2) stop("Ratio must be of length 2, like c(1, 3) or 1:3.")
  
  total_parts <- sum(ratio)
  n2 <- round(n_total * ratio[1] / total_parts)  # focal group (first in ratio)
  n1 <- n_total - n2                             # reference group (second in ratio)
  
  return(list(n = n_total, n1 = n1, n2 = n2))
}

generate_param <- function(I) {
  b <- rnorm(I, mean = 0, sd = 1)
  z <- rnorm(I, mean = 0, sd = sqrt(0.1225))
  a <- exp(z)
  c <- rep(0.2, I)
  
  param <- cbind(b, a, c)
  rownames(param) <- paste0("Item", 1:I)
  colnames(param) <- c("b", "a", "c")
  return(param)
}



simulation_abe2 <- function(N, n1, n2, skup, param1, param2, theta1, theta2, statistics = list(MantelB = TRUE, MantelNUB = TRUE, BresB = TRUE, LogB = TRUE, SIBB = TRUE, cSIBB = TRUE)) {
  
  # Initialize result vectors
  Mantel     <- MantelLow <- MantelHigh <-  Bres     <- Bres2     <- Log     <- SIB     <- cSIB <-  rep(NA, N)
  MantelStat  <- MantelLowStat <- MantelHighStat <- BresStat <- Bres2Stat <- LogStat <- SIBStat <- cSIBStat <- rep(NA, N)
  
  for (i in 1:N) {
    # Generate data for each group
    data1 <- rmvlogis(n1, param1, IRT = TRUE, link = "logit", z.vals = matrix(theta1))
    data2 <- rmvlogis(n2, param2, IRT = TRUE, link = "logit", z.vals = matrix(theta2))
    data  <- rbind(data1, data2)
    
    # Compute total score
    total_scores <- rowSums(data)
    mean_score <- mean(total_scores)
    
    # Split data and grouping vector
    data_low  <- data[total_scores <= mean_score, ]
    skup_low  <- skup[total_scores <= mean_score]
    data_high <- data[total_scores >  mean_score, ]
    skup_high <- skup[total_scores >  mean_score]
    
    # Mantel-Haenszel Test
    if (statistics$MantelB | statistics$MantelNUB) {
      help <- try(difMH(data, skup, focal.name = 1, correct = TRUE, MHstat = "MHChisq", exact = FALSE)[[1]][1], silent = TRUE)
      MantelStat[i] <- help
      Mantel[i] <- if (is.numeric(help)) 1 - pchisq(help, 1) else NA
    }
    
    # Mantel-Haenszel Low
    if (statistics$MantelNUB) {
      help <- try(difMH(data_low, skup_low, focal.name = 1, correct = TRUE, MHstat = "MHChisq", exact = FALSE)[[1]][1], silent = TRUE)
      MantelLowStat[i] <- help
      MantelLow[i] <- if (is.numeric(help)) 1 - pchisq(help, 1) else NA
    }
    
    # Mantel-Haenszel High
    if (statistics$MantelNUB) {
      help <- try(difMH(data_high, skup_high, focal.name = 1, correct = TRUE, MHstat = "MHChisq", exact = FALSE)[[1]][1], silent = TRUE)
      MantelHighStat[i] <- help
      MantelHigh[i] <- if (is.numeric(help)) 1 - pchisq(help, 1) else NA
    }
    
    # Breslow-Day Test
    if (statistics$BresB) {
      help <- try(difBD(data, skup, focal.name = 1)[[1]][1, 3], silent = TRUE)
      BresStat[i] <- help
      Bres[i] <- if (is.numeric(help)) help else NA
    }
    
    # Breslow-Day (trend) Test
    if (statistics$BresB) {
      help <- try(difBD(data, skup, focal.name = 1, BDstat = "trend")[[1]][1, 3], silent = TRUE)
      Bres2Stat[i] <- help
      Bres2[i] <- if (is.numeric(help)) help else NA
    }
    
    # Logistic Regression DIF Test
    if (statistics$LogB) {
      help <- try(difLogistic(data, skup, focal.name = 1)[[1]][1], silent = TRUE)
      LogStat[i] <- help
      Log[i] <- if (is.numeric(help)) 1 - pchisq(help, 2) else NA
    }
    
    # SIBTEST (assumes first item is tested; adapt if needed)
    #print(statistics$SIBB)
    if (statistics$SIBB) {
      help <- try(difSIBTEST(data, skup, focal.name = 1, type = 'udif', purify = FALSE)[[3]][1], silent = TRUE)
      SIBStat[i] <- help
      SIB[i]    <- if (is.numeric(help)) 1 - pchisq(help, 1) else NA
    }
    if (statistics$cSIBB) {
      help <- try(difSIBTEST(data, skup, focal.name = 1, type = 'nudif', purify = FALSE)[[3]][1], silent = TRUE)
      cSIBStat[i] <- help
      cSIB[i]    <- if (is.numeric(help)) 1 - pchisq(help, 1) else NA
    }
    
    # Progress
    message(sprintf("Completed iteration %d at %s", i, date()))
  }
  
  # Combine results
  mat     <- rbind(Mantel, MantelLow, MantelHigh, Mantel, Bres, Bres2, Log, SIB, cSIB)
  matStat <- rbind(MantelStat, MantelLowStat, MantelHighStat,MantelStat, BresStat, Bres2Stat, LogStat, SIBStat, cSIBStat)
  
  rownames(mat)     <- c("MantelNormal", "MantelLow", "MantelHigh", "Mantel", "Bres", "Bres2", "Log", "SIB", "cSIB")
  rownames(matStat) <- c("MantelStatNormal", "MantelLowStat", "MantelHighStat","MantelStat", "BresStat", "Bres2Stat", "LogStat", "SIBStat", "cSIBStat")
  
  return(list(pvalues = mat, statistics = matStat))
}


simul_total2 <- function(N, n_total, rat_n, I, mu_R, mu_F, type = "alpha", diffs = NULL, statistics = list(MantelB = TRUE, MantelNUB = TRUE, BresB = TRUE, LogB = TRUE, SIBB = TRUE, cSIBB = TRUE), timeMeasure = FALSE) {
  # Start time for the entire function
  start_time <- Sys.time()
  # Group sizes
  list_ns <- generate_group_sizes(n_total, rat_n)
  n <- list_ns$n
  n1 <- list_ns$n1
  n2 <- list_ns$n2
  
  # Abilities
  theta1_ab <- rnorm(n1, mean = mu_R, sd = 1)
  theta2_ab <- rnorm(n2, mean = mu_F, sd = 1)
  
  # Item parameters (same initially)
  param1_ab <- generate_param(I)
  param2_ab <- param1_ab
  skup <- c(rep(0, n1), rep(1, n2))
  
  # Init result list
  results <- list()
  
  if (type == "alpha") {
    # Type I error — no DIF introduced
    result_alpha <- simulation_abe2(N,n1,n2,skup, param1_ab, param2_ab, theta1_ab, theta2_ab, statistics)
    results[["alpha"]] <- result_alpha
    
  } else if (type == "unif") {
    # Uniform DIF — manipulate difficulty (b)
    for (i in seq_along(diffs)) {
      d <- diffs[i]
      param2_ab_mod <- param1_ab
      param2_ab_mod[1, 1] <- param1_ab[1, 1] + d # shift difficulty
      
      label <- paste0("Unif_DIF_delta", d)
      results[[label]] <- simulation_abe2(N,n1,n2,skup, param1_ab, param2_ab_mod, theta1_ab, theta2_ab, statistics)
    }
    
  } else if (type == "nonunif") {
    # Non-uniform DIF — manipulate discrimination (a)
    for (i in seq_along(diffs)) {
      delta <- diffs[i]
      param2_ab_mod <- param1_ab
      a_orig <- param1_ab[1, 2]
      c_orig <- param1_ab[1, 3]
      param2_ab_mod[1, 2] <- 2 * a_orig / (2 + delta * a_orig / ((1 - c_orig) * log(2)))
      
      label <- paste0("NonUnif_DIF_delta", delta)
      results[[label]] <- simulation_abe2(N,n1,n2,skup, param1_ab, param2_ab_mod, theta1_ab, theta2_ab, statistics)
    }
  }
  # End time for the entire function
  end_time <- Sys.time()
  time_taken <- end_time - start_time
  message(sprintf("Total time for simulation: %s", time_taken))
  if (timeMeasure) {
    return(time_taken)
  }
  return(results)
}
#Evaluating------------------------------------------------------------

# Type I error
res_alpha <- simul_total2(N = 10, n_total = 500, rat_n = c(1,2), I = 40, mu_R = 0, mu_F = -1, type = "alpha",statistics = list(SIBB = TRUE, MantelB = TRUE, MantelNUB = TRUE, BresB = TRUE, LogB = TRUE, cSIBB = TRUE)
)
res_alpha$alpha$pvalues
# Uniform DIF
res_unif <- simul_total2(N = 10, n_total = 1000, rat_n = c(1,2), I = 5, mu_R = 0, mu_F = -1, 
                        type = "unif", diffs = c(0.5,1,2,4))

# Non-uniform DIF
res_nonunif <- simul_total2(N = 10, n_total = 1000, rat_n = c(1,2), I = 5, mu_R = 0, mu_F = -1, 
                           type = "nonunif", diffs = c(0.4,0.6,0.8,1))



# Time measuring (Need to change arguments from results to time_taken)

run_all_statistics <- function() {
  stat_names <- c("SIBB", "MantelB", "MantelNUB", "BresB", "LogB", "cSIBB")
  results <- numeric(length(stat_names))
  
  for (i in seq_along(stat_names)) {
    # Create statistics list with one TRUE
    stats_list <- setNames(as.list(rep(FALSE, length(stat_names))), stat_names)
    stats_list[[stat_names[i]]] <- TRUE
    
    # Run simulation
    result <- simul_total2(
      N = 10, n_total = 500, rat_n = c(1, 2), I = 40,
      mu_R = 0, mu_F = -1, type = "alpha",
      statistics = stats_list, timeMeasure = TRUE
    )
    
    # Extract numeric result
    results[i] <- result
  }
  
  # Create data frame
  df_results <- data.frame(
    Method = stat_names,
    Value = results
  )
  
  return(df_results)
}

# Run it and store the result

df_results <- run_all_statistics()
print(df_results)


# to do --------------------------------------------------------------
# multiple items ???
# make a time and sample graph

# look into specific libraries 

#####################################################################
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

# Print the dataframe
print(est_alpha_frame)
############################################################################
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
