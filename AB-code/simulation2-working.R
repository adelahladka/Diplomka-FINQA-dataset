#Simulation from ground up
library(ltm)
library(difR)
rm(list = ls())
###################################################################################
#Hyperparameters

N = 10 #number of simulations
m = 5 #number of items (should change to I)
n_tot = 1000 # number of respondents (700, 1000, 1200, 1500) 
rat_n = c(1:1) #(c(5,2), c(1,1), c(5,1), c(2,1))
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



simulation_abe2 <- function(N, n1, n2, skup, param1, param2, theta1, theta2) {
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
    
    # Mantel-Haenszel
    help <- try(difMH(data, skup, focal.name = 1, correct = TRUE, MHstat = "MHChisq", exact = FALSE)[[1]][1], silent = TRUE)
    MantelStat[i] <- help
    Mantel[i] <- if (is.numeric(help)) 1 - pchisq(help, 1) else NA
    
    # Mantel-Haenszel Low
    help <- try(difMH(data_low, skup_low, focal.name = 1, correct = TRUE, MHstat = "MHChisq", exact = FALSE)[[1]][1], silent = TRUE)
    MantelLowStat[i] <- help
    MantelLow[i] <- if (is.numeric(help)) 1 - pchisq(help, 1) else NA
    
    # Mantel-Haenszel High
    help <- try(difMH(data_high, skup_high, focal.name = 1, correct = TRUE, MHstat = "MHChisq", exact = FALSE)[[1]][1], silent = TRUE)
    MantelHighStat[i] <- help
    MantelHigh[i] <- if (is.numeric(help)) 1 - pchisq(help, 1) else NA
    
    # Breslow-Day Test
    help <- try(difBD(data, skup, focal.name = 1)[[1]][1, 3], silent = TRUE)
    BresStat[i] <- help
    Bres[i] <- if (is.numeric(help)) help else NA
    
    # Breslow-Day (trend) Test
    help <- try(difBD(data, skup, focal.name = 1, BDstat = "trend")[[1]][1, 3], silent = TRUE)
    Bres2Stat[i] <- help
    Bres2[i] <- if (is.numeric(help)) help else NA
    
    # Logistic Regression DIF Test
    help <- try(difLogistic(data, skup, focal.name = 1)[[1]][1], silent = TRUE)
    LogStat[i] <- help
    Log[i] <- if (is.numeric(help)) 1 - pchisq(help, 2) else NA
    
    # SIBTEST (assumes first item is tested; adapt if needed)
    help <- try(difSIBTEST(data, skup, focal.name = 1, type = 'udif', purify = FALSE)[[3]][1], silent = TRUE)
    SIBStat[i] <- help
    SIB[i]    <- if (is.numeric(help)) 1 - pchisq(help, 1) else NA
    #SIB[i]     <- safe_stat(quote(difSIBTEST(data, skup, focal.name = 1, purify = FALSE)[[4]][1]))
    help <- try(difSIBTEST(data, skup, focal.name = 1, type = 'nudif', purify = FALSE)[[3]][1], silent = TRUE)
    cSIBStat[i] <- help
    cSIB[i]    <- if (is.numeric(help)) 1 - pchisq(help, 1) else NA
    # Progress
    message(sprintf("Completed iteration %d at %s", i, date()))
  }
  
  # Combine results
  mat     <- rbind(Mantel,MantelLow, MantelHigh, Bres, Bres2, Log, SIB, cSIB)
  matStat <- rbind(MantelStat,MantelLowStat, MantelHighStat, BresStat, Bres2Stat, LogStat, SIBStat, cSIBStat)
  
  rownames(mat)     <- c("Mantel","MantelLow", "MantelHigh", "Bres", "Bres2", "Log", "SIB", "cSIB")
  rownames(matStat) <- c("MantelStat","MantelLowStat", "MantelHighStat", "BresStat", "Bres2Stat", "LogStat", "SIBStat", "cSIBTSTat")
  
  return(list(pvalues = mat, statistics = matStat))
}


simul_total2 <- function(N, n_total, rat_n, I, mu_R, mu_F, type = "alpha", diffs = NULL) {
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
    result_alpha <- simulation_abe2(N,n1,n2,skup, param1_ab, param2_ab, theta1_ab, theta2_ab)
    results[["alpha"]] <- result_alpha
    
  } else if (type == "unif") {
    # Uniform DIF — manipulate difficulty (b)
    for (i in seq_along(diffs)) {
      d <- diffs[i]
      param2_ab_mod <- param1_ab
      param2_ab_mod[1, 1] <- param1_ab[1, 1] + d # shift difficulty
      
      label <- paste0("Unif_DIF_d", d)
      results[[label]] <- simulation_abe2(N,n1,n2,skup, param1_ab, param2_ab_mod, theta1_ab, theta2_ab)
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
      results[[label]] <- simulation_abe2(N,n1,n2,skup, param1_ab, param2_ab_mod, theta1_ab, theta2_ab)
    }
  }
  
  return(results)
}
#Evaluating------------------------------------------------------------

# Type I error
res_alpha <- simul_total2(N = 10, n_total = 1000, rat_n = c(1,1), I = 5, mu_R = 0, mu_F = -1, type = "alpha")

# Uniform DIF
res_unif <- simul_total2(N = 10, n_total = 1000, rat_n = c(5,2), I = 5, mu_R = 0, mu_F = -1, 
                        type = "unif", diffs = c(0.5,1,2,4))

# Non-uniform DIF
res_nonunif <- simul_total2(N = 10, n_total = 1000, rat_n = c(2,1), I = 5, mu_R = 0, mu_F = -1, 
                           type = "nonunif", diffs = c(0.4,0.6,0.8,1))


# to do --------------------------------------------------------------
# multiple items ???
# add NU variants

# look into specific libraries 

?difMH
?difSIBTEST
?
difSIBTEST
?mantelHaenszel
#####################################################################
