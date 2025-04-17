#Simulation from ground up
  rm(list = ls())
library(ltm)
library(difR)

set.seed(123)

#load("simulation2-abe.RData")

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


# ============================
# Function: simulation_abe2
# Purpose: Run a DIF detection simulation with multiple methods
simulation_abe3 <- function(N, n1, n2, skup, param1, param2, theta1, theta2,
                            statistics = list(MantelB = TRUE, MantelNUB = TRUE, BresB = TRUE, LogB = TRUE, SIBB = TRUE, cSIBB = TRUE), I) {
  
  new_matrix <- function() matrix(NA_real_, nrow = I, ncol = N)
  
  results <- list(
    Mantel = new_matrix(), MantelStat = new_matrix(),
    MantelNormal = new_matrix(), MantelNormalStat = new_matrix(),
    MantelLow = new_matrix(), MantelLowStat = new_matrix(),
    MantelHigh = new_matrix(), MantelHighStat = new_matrix(),
    Bres = new_matrix(), BresStat = new_matrix(),
    Bres2 = new_matrix(), Bres2Stat = new_matrix(),
    Log = new_matrix(), LogStat = new_matrix(),
    SIB = new_matrix(), SIBStat = new_matrix(),
    cSIB = new_matrix(), cSIBStat = new_matrix()
  )
  
  for (i in 1:N) {
    data1 <- rmvlogis(n1, param1, IRT = TRUE, link = "logit", z.vals = matrix(theta1))
    data2 <- rmvlogis(n2, param2, IRT = TRUE, link = "logit", z.vals = matrix(theta2))
    data <- rbind(data1, data2)
    
    total_scores <- rowSums(data)
    mean_score <- mean(total_scores)
    data_low <- data[total_scores <= mean_score, ]
    skup_low <- skup[total_scores <= mean_score]
    data_high <- data[total_scores > mean_score, ]
    skup_high <- skup[total_scores > mean_score]
    
    if (statistics$MantelB | statistics$MantelNUB) {
      resstat <- try(difMH(data, skup, focal.name = 1, purify = FALSE), silent = TRUE)
      if (inherits(resstat, "try-error")) {
        results$MantelStat[, i] <- rep(NA_real_, I)
        results$Mantel[, i] <- rep(NA_real_, I)
        
        results$MantelNormalStat[, i] <- rep(NA_real_, I)
        results$MantelNormal[, i] <- rep(NA_real_, I)
      } else {
        results$MantelStat[, i] <- resstat[[1]]
        results$Mantel[, i] <- resstat$p.value
        
        
        results$MantelNormalStat[, i] <- resstat[[1]]
        results$MantelNormal[, i] <- resstat$p.value
      }
    }
    
    if (statistics$MantelNUB) {
      stat_low <- try(difMH(data_low, skup_low, focal.name = 1), silent = TRUE)
      if (inherits(stat_low, "try-error")) {
        results$MantelLowStat[, i] <- rep(NA_real_, I)
        results$MantelLow[, i] <- rep(NA_real_, I)
      } else {
        results$MantelLowStat[, i] <- stat_low[[1]]
        results$MantelLow[, i] <- stat_low$p.value
      }
      
      stat_high <- try(difMH(data_high, skup_high, focal.name = 1), silent = TRUE)
      if (inherits(stat_high, "try-error")) {
        results$MantelHighStat[, i] <- rep(NA_real_, I)
        results$MantelHigh[, i] <- rep(NA_real_, I)
      } else {
        results$MantelHighStat[, i] <- stat_high[[1]]
        results$MantelHigh[, i] <- stat_high$p.value
      }
    }
    
    if (statistics$BresB) {
      stat1 <- try(difBD(data, skup, focal.name = 1), silent = TRUE)
      stat2 <- try(difBD(data, skup, focal.name = 1, BDstat = "trend"), silent = TRUE)
      if (inherits(stat1, "try-error")) {
        results$BresStat[, i] <- rep(NA_real_, I)
        results$Bres[, i] <- rep(NA_real_, I)
      } else {
        results$BresStat[, i] <- stat1[[1]][,1]
        results$Bres[, i] <- stat1$p.value
      }
      if (inherits(stat2, "try-error")) {
        results$Bres2Stat[, i] <- rep(NA_real_, I)
        results$Bres2[, i] <- rep(NA_real_, I)
      } else {
        results$Bres2Stat[, i] <- stat2[[1]][,1]
        results$Bres2[, i] <- stat2$p.value
      }
    }
    
    if (statistics$LogB) {
      stat <- try(difLogistic(data, skup, focal.name = 1), silent = TRUE)
      if (inherits(stat, "try-error")) {
        results$LogStat[, i] <- rep(NA_real_, I)
        results$Log[, i] <- rep(NA_real_, I)
      } else {
        results$LogStat[, i] <- stat[[1]]
        results$Log[, i] <- stat$p.value
      }
    }
    
    if (statistics$SIBB) {
      stat <- try(difSIBTEST(data, skup, focal.name = 1, type = 'udif', purify = FALSE), silent = TRUE)
      if (inherits(stat, "try-error")) {
        results$SIBStat[, i] <- rep(NA_real_, I)
        results$SIB[, i] <- rep(NA_real_, I)
      } else {
        results$SIBStat[, i] <- stat[[1]]
        results$SIB[, i] <- stat$p.value
      }
    }
    
    if (statistics$cSIBB) {
      stat <- try(difSIBTEST(data, skup, focal.name = 1, type = 'nudif', purify = FALSE), silent = TRUE)
      if (inherits(stat, "try-error")) {
        results$cSIBStat[, i] <- rep(NA_real_, I)
        results$cSIB[, i] <- rep(NA_real_, I)
      } else {
        results$cSIBStat[, i] <- stat[[1]]
        results$cSIB[, i] <- stat$p.value
      }
    }
    
    message(sprintf("Completed iteration %d at %s", i, date()))
  }
  
  results_df <- lapply(results, function(x) {
    df <- as.data.frame(x)
    rownames(df) <- paste0("Item_", 1:I)
    colnames(df) <- paste0("Sim_", 1:N)
    df
  })
  
  return(results_df)
}





simul_total3 <- function(N, n_total, rat_n, I, mu_R, mu_F, 
                         type = c(0, 0),  # proportions of unif and nonunif DIF
                         diffs_unif = 0.5, diffs_nonunif = 1,
                         statistics = list(MantelB = TRUE, MantelNUB = TRUE, BresB = TRUE, LogB = TRUE, SIBB = TRUE, cSIBB = TRUE))
                         {
  
  start_time <- Sys.time()
  
  # Setup
  list_ns <- generate_group_sizes(n_total, rat_n)
  n1 <- list_ns$n1
  n2 <- list_ns$n2
  theta1_ab <- rnorm(n1, mean = mu_R, sd = 1)
  theta2_ab <- rnorm(n2, mean = mu_F, sd = 1)
  param1_ab <- generate_param(I)
  param2_ab <- param1_ab
  skup <- c(rep(0, n1), rep(1, n2))
  
  # Determine how many items get DIF
  n_unif <- if (type[1] > 0) round(type[1] * I) else 0
  n_nonunif <- if (type[2] > 0) round(type[2] * I) else 0
  
  # Select DIF items
  dif_items <- if ((n_unif + n_nonunif) > 0) sample(1:I, n_unif + n_nonunif, replace = FALSE) else integer(0)
  unif_items <- if (n_unif > 0) dif_items[1:n_unif] else integer(0)
  nonunif_items <- if (n_nonunif > 0) dif_items[(n_unif + 1):(n_unif + n_nonunif)] else integer(0)
  
  # Apply Uniform DIF (adjust difficulty b → param[,1])
  if (n_unif > 0) {
    diffs_u <- rep(diffs_unif, length.out = n_unif)
    param2_ab[unif_items, 1] <- param1_ab[unif_items, 1] + diffs_u
  }
  
  # Apply Non-Uniform DIF (adjust discrimination a → param[,2])
  if (n_nonunif > 0) {
    diffs_n <- rep(diffs_nonunif, length.out = n_nonunif)#TODO: can be added as string, or randomly in some range
    for (j in seq_along(nonunif_items)) {
      idx <- nonunif_items[j]
      a_orig <- param1_ab[idx, 2]  # original discrimination
      c_orig <- param1_ab[idx, 3]  # guessing parameter (usually 0 for 2PL)
      delta <- diffs_n[j]
      
      # Adjusted discrimination for focal group (2PL logic: more/less sensitive to ability)
      new_a <- 2 * a_orig / (2 + delta * a_orig / ((1 - c_orig) * log(2)))
      
      param2_ab[idx, 2] <- new_a
    }
  }
  #print(I)
  res <- simulation_abe3(N, n1, n2, skup, param1_ab, param2_ab, theta1_ab, theta2_ab, statistics,I)
  
  end_time <- Sys.time()
  
  # Metadata for output
  metadata <- list(
    N = N,
    n_total = n_total,
    n1 = n1,
    n2 = n2,
    rat_n = rat_n,
    I = I,
    mu_R = mu_R,
    mu_F = mu_F,
    type = type,
    n_unif = n_unif,
    n_nonunif = n_nonunif,
    unif_items = unif_items,
    nonunif_items = nonunif_items,
    diffs_unif = diffs_unif,
    diffs_nonunif = diffs_nonunif,
    time_taken = end_time - start_time,
    start_time = start_time,
    end_time = end_time,
    statistics = statistics
  )
  
  return(list(results = res, metadata = metadata))
}

#################################################################################
#Evaluating------------------------------------------------------------
res_alpha2 <- simul_total3(N = 1000, n_total = 500, rat_n = c(1,2), I = 50, mu_R = 0, mu_F = 0, 
                         type = c(0,0), diffs_unif = 0)

res_bothdif <- simul_total3(N = 1000, n_total = 500, rat_n = c(1,2), I = 50, mu_R = 0, mu_F = 0, 
                           type = c(0.1,0.1), diffs_unif = 0.6, diffs_nonunif = 0.8)

res_bothdif2 <- simul_total3(N = 5000, n_total = 500, rat_n = c(1,2), I = 50, mu_R = 0, mu_F = 0, 
                            type = c(0.1,0.1), diffs_unif = 0.6, diffs_nonunif = 0.8)


res_udif <- simul_total3(N = 1000, n_total = 500, rat_n = c(1,2), I = 50, mu_R = 0, mu_F = 0, 
                            type = c(0.1,0), diffs_unif = 0.6)

res_nondif <- simul_total3(N = 1000, n_total = 500, rat_n = c(1,2), I = 50, mu_R = 0, mu_F = 0, 
                            type = c(0,0.1), diffs_nonunif = 0.8)


run_all_simulations <- function() {
  print('fun1')
  res_alpha2 <- simul_total3(
    N = 1000, n_total = 500, rat_n = c(1, 2), I = 50,
    mu_R = 0, mu_F = 0, type = c(0, 0), diffs_unif = 0
  )
  print('fun2')
  res_bothdif <- simul_total3(
    N = 1000, n_total = 500, rat_n = c(1, 2), I = 50,
    mu_R = 0, mu_F = 0, type = c(0.1, 0.1), diffs_unif = 0.6, diffs_nonunif = 0.8
  )
  print('fun3')
  res_udif <- simul_total3(
    N = 1000, n_total = 500, rat_n = c(1, 2), I = 50,
    mu_R = 0, mu_F = 0, type = c(0.1, 0), diffs_unif = 0.6
  )
  print('fun4')
  res_nondif <- simul_total3(
    N = 1000, n_total = 500, rat_n = c(1, 2), I = 50,
    mu_R = 0, mu_F = 0, type = c(0, 0.1), diffs_nonunif = 0.8
  )
  
  return(list(
    alpha2 = res_alpha2,
    both_dif = res_bothdif,
    uniform_dif = res_udif,
    nonuniform_dif = res_nondif
  ))
}

results <- run_all_simulations()




######################################################################################
#----------------------------------------------------------------
#Debugging
group_size_test <- generate_group_sizes(500, c(1,2))
params_item <- generate_param(50)
n1<- group_size_test$n1
n2<- group_size_test$n2
theta1_ab <- rnorm(n1, mean = 0, sd = 1)
theta2_ab <- rnorm(n2, mean = 0, sd = 1)
param1_ab <- generate_param(50)
param2_ab <- param1_ab
data1 <- rmvlogis(n1, param1_ab, IRT = TRUE, link = "logit", z.vals = matrix(theta1_ab)) # reference
data2 <- rmvlogis(n2, param2_ab, IRT = TRUE, link = "logit", z.vals = matrix(theta2_ab)) # focal
data <- rbind(data1, data2)
skup <- c(rep(0, n1), rep(1, n2))
difSIBTEST(data, skup, focal.name = 1, type = 'udif', purify = FALSE)[[1]]
difSIBTEST(data, skup, focal.name = 1, type = 'udif', purify = FALSE)$p.value

?difSIBTEST
length(ress$X2)
ress <- try(difSIBTEST(data, skup, focal.name = 1, type = 'udif', purify = FALSE), silent = TRUE)

ress2 <- try(difBD(data, skup, focal.name = 1), silent = TRUE)[[1]]
ress2[[1]][,1]
ress2$BDstat
ress3 <- try(difLogistic(data, skup, focal.name = 1), silent  =TRUE)
ress3[[1]]
ress4 <- try(difSIBTEST(data, skup, focal.name = 1, type = 'udif', purify = FALSE), silent = TRUE)[[1]]
#---------------------------------------------------------------------------------------------------------
###################################################################################################
calculate_power_rate <- function(input) {
  methods <- c()
  if (input$metadata$statistics$MantelB) methods <- c(methods, "Mantel")
  if (input$metadata$statistics$MantelNUB) methods <- c(methods, "MantelNormal", "MantelLow", "MantelHigh")
  if (input$metadata$statistics$BresB) methods <- c(methods, "Bres", "Bres2")
  if (input$metadata$statistics$LogB) methods <- c(methods, "Log")
  if (input$metadata$statistics$SIBB) methods <- c(methods, "SIB")
  if (input$metadata$statistics$cSIBB) methods <- c(methods, "cSIB")
  
  dif_items <- sort(unique(c(input$metadata$unif_items, input$metadata$nonunif_items)))
  
  power_results <- list()
  
  for (method in methods) {
    #print(method)
    #print(dif_items)
    mat <- input$results[[method]]
    mat <- mat[dif_items,]
    alpha <- if (method %in% c("MantelNormal", "MantelLow", "MantelHigh")) 0.01 else 0.05
    
    detected <- mat <= alpha
    num_detected <- sum(detected, na.rm = TRUE)
    denom <- length(dif_items) * input$metadata$N
    #print(detected)
    #print(denom)
    power_results[[method]] <- num_detected / denom
  }
  
  return(power_results)
}


calculate_rejection_rate <- function(input) {
  methods <- c()
  if (input$metadata$statistics$MantelB) methods <- c(methods, "Mantel")
  if (input$metadata$statistics$MantelNUB) methods <- c(methods, "MantelNormal", "MantelLow", "MantelHigh")
  if (input$metadata$statistics$BresB) methods <- c(methods, "Bres", "Bres2")
  if (input$metadata$statistics$LogB) methods <- c(methods, "Log")
  if (input$metadata$statistics$SIBB) methods <- c(methods, "SIB")
  if (input$metadata$statistics$cSIBB) methods <- c(methods, "cSIB")
  
  #dif_items <- sort(unique(c(input$metadata$unif_items, input$metadata$nonunif_items)))
  #all_items <- seq_len(nrow(input$results[[methods[1]]]))
  #non_dif_items <- setdiff(all_items, dif_items)
  
  
  rejection_results <- list()
  
  for (method in methods) {
    mat <- input$results[[method]]
    alpha <- if (method %in% c("MantelNormal", "MantelLow", "MantelHigh")) 0.01 else 0.05
    
    detected <- mat <= alpha
    num_false_detected <- sum(detected, na.rm = TRUE)
    denom <- input$metadata$I * input$metadata$N
    
    rejection_results[[method]] <- num_false_detected / denom
  }
  
  return(rejection_results)
}

results$alpha2$metadata
results$both_dif$metadata
results$uniform_dif$metadata
results$nonuniform_dif$metadata
calculate_rejection_rate(
  results$alpha2)
calculate_power_rate(results$both_dif)

calculate_power_rate(results$nonuniform_dif)

calculate_power_rate(results$uniform_dif)

calculate_power_rate(res_bothdif2)

#Testing --------------------------------------------------------------------
methods <- c()
if (input$metadata$statistics$MantelB) methods <- c(methods, "Mantel")
if (input$metadata$statistics$MantelNUB) methods <- c(methods, "MantelNormal", "MantelLow", "MantelHigh")
if (input$metadata$statistics$BresB) methods <- c(methods, "Bres", "Bres2")
if (input$metadata$statistics$LogB) methods <- c(methods, "Log")
if (input$metadata$statistics$SIBB) methods <- c(methods, "SIB")
if (input$metadata$statistics$cSIBB) methods <- c(methods, "cSIB")

dif_items <- sort(unique(c(input$metadata$unif_items, input$metadata$nonunif_items)))

power_results <- list()

for (method in methods) {
  mat <- input$results$method[dif_items,]
  alpha <- if (method %in% c("MantelNormal", "MantelLow", "MantelHigh")) 0.01 else 0.05
  
  detected <- mat <= alpha
  num_detected <- sum(detected, na.rm = TRUE)
  denom <- length(dif_items) * input$metadata$N
  print(detected)
  print(denom)
  power_results <- num_detected / denom
}

#################################################################################

















################################################################################

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



######################################################################
# Define the I values
I_values <- c(5, 10, 15)#, 20, 25, 30, 40, 50, 60, 70, 80)
n_total_value <- c(100,200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1400)
# Prepare lists to store the results

# Loop through I values
for (n_val in n_total_value) {
  print(n_val)
  res_alpha2 <- simul_total3(
    N = 10, n_total = n_val, rat_n = c(1, 2), I = 50,
    mu_R = 0, mu_F = 0, type = c(0, 0), diffs_unif = 0
  )
  est_alpha_result2 <- calculate_rejection_rate(res_alpha2)
}#not working, to do


save.image(file = "simulation2-abe.RData")
