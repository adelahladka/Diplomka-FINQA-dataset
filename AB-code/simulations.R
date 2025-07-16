# Evaluation

rm(list = ls())
source("functions.R")
set.seed(2025)
################################################################################
# ===========================
# Paramterers
# N - number of simulation (M in thesis)
# n_total - Sample size (N in thesis)
# rat_n - Reference: Focal group (rho in thesis)
# I - number of items 
# type = c(a, b) - Proportions of uniform (a) and non-uniform (b) DIF
# diffs_unif  - DIF effect size for uniform DIF (delta in thesis)
# diffs_nonunif - DIF effect size for non-uniform DIF (Delta in thesis)
# diff_random - Randomize the location of DIF items, in case of FALSE DIF items 
#will begin indexing from 1
#============================

################################################################################
#===========================
# Scenario 1: Sample size N
# Scenario 3: DIF effect size
#============================

tab <- c()
res_total <- list() 

samples_n <- c(10,20,30, 40, 50, 60, 70, 80, 100, 200, 400, 500, 1000,1200,2000)
difSize <- c(0.4, 0.6, 0.8, 1.0, 1.2, 1.5)

for(sample.size in samples_n) {
  res <- simul_total3(
    N = 1000,
    n_total = sample.size,
    rat_n = c(1, 1), 
    I = 40,
    mu_R = 0, mu_F = 0,
    type = c(0.025,0),  # Needs to be changed depending on type of DIF
    diffs_unif = 1, # Needs to be changed depending on type of DIF
    seed_arg = 2025, diff_random = FALSE,
   statistics =  list(MantelB = FALSE, MantelNUB = FALSE, BresB = FALSE, LogB = FALSE, SIBB = TRUE, cSIBB = FALSE))
  tab <- rbind(
    tab,
    cbind(
      N = 100,
      n = sample.size,
      I = 40 ,
      DIF_type = "uniform",
      DIF_proportion = 0.025,
      DIF_size = 1,
      power = unlist(calculate_power_rate(res)),
      rejection = calculate_rejection_rate(res)
    )
  )
  res_total[[as.character(sample.size)]] <- res
  message(sprintf("Completed iteration %.2f at %s", sample.size, date()))
}

###############################################################################
#=======================
# Scenario 2: Varying group ratio rho
#=======================

ratios_R <- cbind(c(1,1), c(2,1), c(3,1), c(4,1), c(1,2), c(1,3), c(1,4))
for (i in 1:ncol(ratios_R)) {
  ratio.size <- ratios_R[, i]  
  res <- simul_total3(
    N = 1000,
    n_total = 400,
    rat_n = ratio.size,
    I = 40,
    mu_R = 0,
    mu_F = 0,
    type = c(0, 0.025), # Needs to be changed depending on type of DIF
    diffs_nonunif = 1, # Needs to be changed depending on type of DIF
    seed_arg = 2025,
    diff_random = FALSE,
    statistics = list(MantelB = FALSE, MantelNUB = TRUE, BresB = FALSE, LogB = TRUE, SIBB = FALSE, cSIBB = TRUE)
  )
  
  tab <- rbind(
    tab,
    cbind(
      N = 1000,
      n = 400,
      rat = paste0(ratio.size[1], ":", ratio.size[2]),
      I = 40,
      DIF_type = "nonuniform",
      DIF_proportion = 0.025,
      DIF_size = 1,
      power = unlist(calculate_power_rate(res)),
      rejection = calculate_rejection_rate(res)
    )
  )
  
  res_total[[paste0(ratio.size[1], ":", ratio.size[2])]] <- res
  message(sprintf("Completed iteration %s at %s", paste(ratio.size, collapse = ":"), date()))
}



###############################################################################
#==========================
#Scenario 4
#=========================

samples_I <- c(5, 10, 20, 40, 50, 60, 70, 80, 100, 200)
set.seed(2025) 
max_I <- max(samples_I)
master_params <- generate_param(max_I, seed_arg = 2025)

# One DIF item in each case: 1/I proportion of DIF items
ratios <- list(
  c(1/5,0),
  c(1/10,0),
  c(1/20, 0),
  c(1/40, 0),
  c(1/50, 0),
  c(1/60, 0),
  c(1/70, 0),
  c(1/80, 0),
  c(1/100, 0),
  c(1/200, 0)
)

ratios2 <- list(
  c(0, 1/5),
  c(0, 1/10),
  c(0, 1/20),
  c(0, 1/40),
  c(0, 1/50),
  c(0, 1/60),
  c(0, 1/70),
  c(0, 1/80),
  c(0, 1/100),
  c(0, 1/200)
)


tab2 <- data.frame()
res_total2 <- list()

for (i in seq_along(samples_I)) {
  item.size <- samples_I[i]
  type.size <- ratios2[[i]] # Needs to be changed depending on type of DIF
  
  res <- simul_total3(
    N = 1000,
    n_total = 1000,
    rat_n = c(1, 1),
    I = item.size,
    mu_R = 0, mu_F = 0,
    type = type.size,
    diffs_nonunif = 1,
    seed_arg = 2025,
    diff_random = FALSE,
    statistics = list(
      MantelB = FALSE,
      MantelNUB = TRUE,
      BresB = FALSE,
      LogB = TRUE,
      SIBB = FALSE,
      cSIBB = TRUE
    ), master_params = master_params[1:item.size, ]
  )
  
  tab2 <- rbind(
    tab2,
    cbind(
      N_total = 1000,
      ratio_ = "1:1",
      n_sample = 1000,
      I = item.size,
      DIF_type = "nonuniform",
      DIF_proportion = 1 / item.size,
      DIF_size = 1,
      power = unlist(calculate_power_rate(res)),
      rejection = calculate_rejection_rate(res)
    )
  )
  
  res_total2[[item.size]] <- res
  message(sprintf("Completed iteration for I = %d at %s", item.size, date()))
}
##############################################################################
#==================
# Scenario 5: Penfield ratios
#==================


n_tot <- c(700, 1000, 1200, 1500)
ratios <- list(c(5, 2), c(1, 1), c(10, 2), c(2, 1))

tab2 <- data.frame()
res_total2 <- list()


for (i in seq_along(n_tot)) {
  sample.size <- n_tot[i]
  ratio <- ratios[[i]]
  
  res <- simul_total3(
    N = 1000,
    n_total = sample.size,
    rat_n = ratio,
    I = 40,
    mu_R = 0, mu_F = 0,
    type = c(0, 0.025), # Needs to be changed depending on type of DIF
    #diffs_unif = 1,
    diffs_nonunif = 1, # Needs to be changed depending on type of DIF
    seed_arg = 2025,
    statistics = list(
      MantelB = TRUE,
      MantelNUB = TRUE,
      BresB = FALSE,
      LogB = TRUE,
      SIBB = FALSE,
      cSIBB = TRUE
    )
  )
  
  tab2 <- rbind(
    tab2,
    cbind(
      N_total = sample.size,
      ratio_ref_foc = paste(ratio, collapse = ":"),
      I = I,
      DIF_type = "non-uniform",
      DIF_proportion = 1,
      DIF_size = 1,
      power = unlist(calculate_power_rate(res)),
      rejection = calculate_rejection_rate(res)
    )
  )
  
  res_total2[[paste0("N", sample.size, "_R", ratio[1], "_F", ratio[2])]] <- res
}















#=================================
# Saving results
#================================

#save(res_total, file = "results/variable-n-res-U-logreg.RData")
#save(tab, file = "results/variable-n-U-logreg.RData")

# For Scenario 4 and 5
#save(tab2, file = "results/variable-I-rat-NONU-5-10.RData")
#save(res_total2, file = "results/variable-I-rat-res-NONU-5-10.RData")



