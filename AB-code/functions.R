# Functions for simulations

# ============================
# Packages
# ============================
library(ltm)
library(difR)

# ============================
# Function: generate_group_sizes
# Purpose: Calculate group sizes for focal and reference groups given total N and ratio
# ============================
generate_group_sizes <- function(n_total, ratio = c(1, 1)) { # Reference: Focal
  if (length(ratio) != 2) stop("Ratio must be of length 2, like c(1, 3) or 1:3.")

  total_parts <- sum(ratio)
  n2 <- round(n_total * ratio[1] / total_parts) # Size of focal group (first in ratio)
  n1 <- n_total - n2 # Size of reference group (remaining subjects)

  return(list(n = n_total, n1 = n1, n2 = n2)) # Return total and both group sizes
}

# ============================
# Function: generate_param
# Purpose: Generate item parameters for 3PL model: b (difficulty), a (discrimination), c (guessing)
# ============================
generate_param <- function(I) {
  b <- rnorm(I, mean = 0, sd = 1) # Difficulty parameters
  z <- rnorm(I, mean = 0, sd = sqrt(0.1225)) # Latent variable for log-discrimination
  a <- exp(z) # Discrimination (log-normal)
  c <- rep(0.2, I) # Guessing parameter fixed at 0.2

  param <- cbind(b, a, c)
  rownames(param) <- paste0("Item", 1:I)
  colnames(param) <- c("b", "a", "c")
  return(param)
}

# ============================
# Function: simulation_abe2
# Purpose: Run a DIF detection simulation with multiple methods
# ============================
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
        results$BresStat[, i] <- stat1[[1]][, 1]
        results$Bres[, i] <- stat1$p.value
      }
      if (inherits(stat2, "try-error")) {
        results$Bres2Stat[, i] <- rep(NA_real_, I)
        results$Bres2[, i] <- rep(NA_real_, I)
      } else {
        results$Bres2Stat[, i] <- stat2[[1]][, 1]
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
      stat <- try(difSIBTEST(data, skup, focal.name = 1, type = "udif", purify = FALSE), silent = TRUE)
      if (inherits(stat, "try-error")) {
        results$SIBStat[, i] <- rep(NA_real_, I)
        results$SIB[, i] <- rep(NA_real_, I)
      } else {
        results$SIBStat[, i] <- stat[[1]]
        results$SIB[, i] <- stat$p.value
      }
    }

    if (statistics$cSIBB) {
      stat <- try(difSIBTEST(data, skup, focal.name = 1, type = "nudif", purify = FALSE), silent = TRUE)
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


# ============================
# Function: simul_total3
# Purpose:
# ============================
simul_total3 <- function(N, n_total, rat_n, I, mu_R, mu_F,
                         type = c(0, 0), # proportions of unif and nonunif DIF
                         diffs_unif = 0.5, diffs_nonunif = 1,
                         statistics = list(MantelB = TRUE, MantelNUB = TRUE, BresB = TRUE, LogB = TRUE, SIBB = TRUE, cSIBB = TRUE)) {
  start_time <- Sys.time()
  #-------------------
  # generating data
  #-------------------
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
    diffs_n <- rep(diffs_nonunif, length.out = n_nonunif) # TODO: can be added as string, or randomly in some range
    for (j in seq_along(nonunif_items)) {
      idx <- nonunif_items[j]
      a_orig <- param1_ab[idx, 2] # original discrimination
      c_orig <- param1_ab[idx, 3] # guessing parameter (usually 0 for 2PL)
      delta <- diffs_n[j]

      # Adjusted discrimination for focal group (2PL logic: more/less sensitive to ability)
      new_a <- 2 * a_orig / (2 + delta * a_orig / ((1 - c_orig) * log(2)))

      param2_ab[idx, 2] <- new_a
    }
  }
  # print(I)
  res <- simulation_abe3(N, n1, n2, skup, param1_ab, param2_ab, theta1_ab, theta2_ab, statistics, I)

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

# ============================
# Function: calculate_power_rate
# Purpose: calculates power rate
# ============================
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
    mat <- input$results[[method]]
    mat <- mat[dif_items, ]
    alpha <- if (method %in% c("MantelNormal", "MantelLow", "MantelHigh")) 0.01 else 0.05

    detected <- mat <= alpha
    num_detected <- sum(detected, na.rm = TRUE)
    denom <- length(dif_items) * input$metadata$N

    power_results[[method]] <- num_detected / denom
  }

  return(power_results)
}

# ============================
# Function: calculate_rejection_rate
# Purpose: calculates rejection rate
# ============================
calculate_rejection_rate <- function(input) {
  methods <- c()
  if (input$metadata$statistics$MantelB) methods <- c(methods, "Mantel")
  if (input$metadata$statistics$MantelNUB) methods <- c(methods, "MantelNormal", "MantelLow", "MantelHigh")
  if (input$metadata$statistics$BresB) methods <- c(methods, "Bres", "Bres2")
  if (input$metadata$statistics$LogB) methods <- c(methods, "Log")
  if (input$metadata$statistics$SIBB) methods <- c(methods, "SIB")
  if (input$metadata$statistics$cSIBB) methods <- c(methods, "cSIB")

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
