rm(list = ls())

source("functions.R")

set.seed(2025)




tab <- c()
res_total <- list() 
#Useful for making some pretty graphs
samples_n <- c(10,20,30, 40, 50, 60, 70, 80, 100, 200, 400, 500, 1000,1200,2000)
samples_I <- c(20, 40, 50, 60, 70, 80, 100, 200)
difSize <- c(0.4, 0.6, 0.8, 1.0, 1.2, 1.5)


for (dif.size in difSize) {
  # základní nastavení, pro různé sample sizes zkoušíme jednoduchý setting
  # poměr reference:focal 1:1, 20 položek, N(0, 1) rozdělení traitu,
  # uniformní DIF o velikosti 1, 5% DIFových položek (tj. 1)
  
  res <- simul_total3(
    N =1,
    n_total =400, rat_n = c(1, 1), I = 40,
    mu_R = 0, mu_F = 0,
    type = c(0, 0.025), diffs_nonunif = dif.size, seed_arg = 2025, diff_random = FALSE,
   statistics =  list(MantelB = FALSE, MantelNUB = TRUE, BresB = FALSE, LogB = TRUE, SIBB = FALSE, cSIBB = TRUE))
  tab <- rbind(
    tab,
    cbind(
      N = 1000,
      n = 400,
      I = 40 ,
      DIF_type = "nonuniform",
      DIF_proportion = 0.025,
      DIF_size = dif.size,
      power = unlist(calculate_power_rate(res)),
      rejection = calculate_rejection_rate(res)
    )
  )
  res_total[[as.character(dif.size)]] <- res
  message(sprintf("Completed iteration %.2f at %s", dif.size, date()))
}
tab
res_total$`1.2`$metadata
calculate_itemwise_detection2(res_total$`10`, include_all_items = TRUE)
###############################################################################
ratios_R <- cbind(c(1,1), c(2,1), c(3,1), c(4,1), c(1,2), c(1,3), c(1,4))
for (i in 1:ncol(ratios_R)) {
  ratio.size <- ratios_R[, i]  # this gives you c(a, b)
  
  #print(ratio.size)
  
  res <- simul_total3(
    N = 1000,
    n_total = 400,
    rat_n = ratio.size,
    I = 40,
    mu_R = 0,
    mu_F = 0,
    type = c(0, 0.025),
    diffs_nonunif = 1,
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
#Ratio testing as per Penfiled THIS

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
    type = c(0, 0.025),
    #diffs_unif = 1,
    diffs_nonunif = 1,
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
###############################################################################
#One DIF item always
samples_I <- c(20, 40, 50, 60, 70, 80, 100, 200)
# Generate master item parameters once for the largest test
set.seed(2025)  # Use same seed as in your generate_param function
max_I <- max(samples_I)  # This will be 200
master_params <- generate_param(max_I, seed_arg = 2025)

# One DIF item always - using consistent parameters



# One DIF item in each case: 1/I proportion of DIF items
ratios <- list(
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
  type.size <- ratios[[i]]
  
  res <- simul_total3(
    N = 1000,
    n_total = 1000,
    rat_n = c(1, 1),
    I = item.size,
    mu_R = 0, mu_F = 0,
    type = type.size,
    diffs_unif = 1,
    seed_arg = 2025,
    diff_random = FALSE,
    statistics = list(
      MantelB = TRUE,
      MantelNUB = FALSE,
      BresB = FALSE,
      LogB = TRUE,
      SIBB = TRUE,
      cSIBB = FALSE
    ), master_params = master_params[1:item.size, ]  # Use first item.size items from master
  )
  
  tab2 <- rbind(
    tab2,
    cbind(
      N_total = 1000,
      ratio_ = "1:1",
      n_sample = 1000,
      I = item.size,
      DIF_type = "uniform",
      DIF_proportion = 1 / item.size,
      DIF_size = 1,
      power = unlist(calculate_power_rate(res)),
      rejection = calculate_rejection_rate(res)
    )
  )
  
  res_total2[[item.size]] <- res
  message(sprintf("Completed iteration for I = %d at %s", item.size, date()))
}

####
# Zkoumání parametru a
res_total2N$N700_R5_F2$metadata$paramR
sorted_matrix <- res_total2N$N700_R5_F2$metadata$paramF[
  order(res_total2N$N700_R5_F2$metadata$paramF[, "a"]), 
]



a_values <- sorted_matrix[, "a"]
row_labels <- rownames(sorted_matrix)

# Main plot
plot(a_values, type = "b", pch = 19, col = "blue",
     xaxt = "n", xlab = "Row Name", ylab = "a",
     main = "Plot of 'a' by Row Name")

# Add row names as x-axis labels
axis(1, at = 1:length(row_labels), labels = row_labels, las = 2, cex.axis = 0.7)

# Highlight Item7 in red
item7_index <- which(row_labels == "Item7")
points(item7_index, res_total2$N700_R5_F2$metadata$paramR['Item7','a'], col = "red", pch = 19, cex = 1.2)


res_total2$N700_R5_F2$metadata$paramF['Item7','b']
##############################################################################





##############################################################################
##Time measuring per method
n_tot <- c(700, 1000, 1200, 1500) #should be trigerred together
ratios <- c(c(5,2),c(1,1),c(10,2),c(2,1))
# Define the six methods and their names
method_names <- c("MantelB","MantelNUB", "LogB","BresB", "SIBB", "cSIBB")
# Preallocate a list to hold each run’s output
isolated_runs <- vector("list", length(method_names))
names(isolated_runs) <- method_names

for (m in method_names) {
  # build a statistics list with only m = TRUE
  stats_list <- setNames(as.list(rep(FALSE, length(method_names))), method_names)
  stats_list[[m]] <- TRUE
  
  # run simul_total3 with only that method turned on
  isolated_runs[[m]] <- simul_total3(
    N           = 1000,
    n_total     = 1000,
    rat_n       = c(1, 1),
    I           = 40,
    mu_R        = 0,
    mu_F        = 0,
    type        = c(0, 0.025),    # 5% uniform, 0% non-uniform
    diffs_unif  = 0,
    diffs_nonunif = 1,
    seed_arg    = 2025,
    statistics  = stats_list
  )
}




# Initialize or load your tracking data frame
if (!exists("timing_results")) {
  timing_results <- data.frame()
}

# Helper function to safely extract time
get_time <- function(run, method) {
  if (!is.null(run[[method]])) {
    return(run[[method]]$metadata$time_taken)
  } else {
    return(NA)
  }
}
calculate_itemwise_detection(isolated_runs$MantelNUB)
calculate_itemwise_detection(isolated_runs$LogB)
calculate_itemwise_detection(isolated_runs$cSIBB)


# Collect shared metadata
n_total <- isolated_runs$LogB
n_ratio <- paste(isolated_runs$cSIBB$metadata$rat_n, collapse = ":")
I <- isolated_runs$cSIBB$metadata$I

# Create a new row of results
new_row <- data.frame(
  n_total = n_total,
  n_ratio = n_ratio,
  I = I,
  Mantel = get_time(isolated_runs, "MantelB"),
  MantelNormal = get_time(isolated_runs, "MantelNUB"),
  Bres = get_time(isolated_runs, "BresB"),
  Log = get_time(isolated_runs, "LogB"),
  SIB = get_time(isolated_runs, "SIBB"),
  cSIB = get_time(isolated_runs, "cSIBB")
)

# Add to your results table
timing_results <- rbind(timing_results, new_row)

timing_results[c(1,8,9,10,11),c(1,5,7,9)]
tab

tab
tab2
res_total2$N700_R5_F2$metadata
# save results
save(tab2, file = "results/variable-I-rat-U.RData")

save(res_total2, file = "results/variable-I-rat-res-U.RData")
#save(res_total2, file = "results/penfieldRatio-resNU.RData")

save(res_total, file = "results/variable-ratio-res-NONU.RData")
save(tab, file = "results/variable-ratio-NONU.RData")
#I should do for non uniform
res_total2$


#Vote ratio --------------------------------------------------------------------
#Scenario 1: N varied
rm(list = ls())
source("functions.R")
samples_n <- c(10,20,30, 40, 50, 60, 70, 80, 100, 200, 400, 500, 1000,1200,2000)
load("results/variable-n-res-U.RData")
combined_results_n_U <- combine_vote_ratios(res_total, samples_n)
#-----------------------------------------------------------------------------
#rm(list = ls())
#source("functions.R")
#samples_n <- c(10,20,30, 40, 50, 60, 70, 80, 100, 200, 400, 500, 1000,1200,2000)
load("results/variable-n-res-NONU.RData")
combined_results_n_NONU <-combine_vote_ratios(res_total, samples_n)

################################################################################
#Scenario 2: I varied
load("results/variable-I-res-U.RData")
samples_I <- c(20, 40, 50, 60, 70, 80, 100, 200)
combined_results_I_U <- combine_vote_ratios(res_total, samples_I)
#-------------------------------------------------------------------------------7
load("results/variable-I-res-NONU.RData")
combined_results_I_NONU <- combine_vote_ratios(res_total, samples_I)
combined_results_I_NONU$item
################################################################################
#Scenario 3: I varied gamma varried
load("results/variable-I-rat-res-U.RData")
samples_Irat <- c('I20', 'I40', 'I50', 'I60', 'I70', 'I80', 'I100', 'I200')
combined_results_I_rat_U <- combine_vote_ratios(res_total2, samples_Irat)
#-------------------------------------------------------------------------------7
load("results/variable-I-rat-res-NONU.RData")
combined_results_I_rat_NONU <- combine_vote_ratios(res_total2, samples_Irat)

################################################################################
#Scenario 4: rho varied
ratios_R <- c("1:1", "2:1", "3:1", "4:1", "1:2", "1:3", "1:4")
load("results/variable-ratio-res-U.RData")
combined_results_rho_U <- combine_vote_ratios(res_total, ratios_R)
#-------------------------------------------------------------------------------7
load("results/variable-ratio-res-NONU.RData")
combined_results_rho_NONU <- combine_vote_ratios(res_total, ratios_R)


################################################################################
#Scenario 2: dfiSIze varied
difSize <- c(0.4, 0.6, 0.8, 1.0, 1.2, 1.5)
load("results/variable-difSize-res-U.RData")
combined_results_difSize_U <- combine_vote_ratios(res_total, difSize)
#-------------------------------------------------------------------------------7
load("results/variable-difSize-res-NONU.RData")
combined_results_difSize_NONU <- combine_vote_ratios(res_total, difSize)





item_1 <- combined_results_I_U[combined_results_I_U$item == 1, ]
item_1_N <- combined_results_I_NONU[combined_results_I_NONU$item == 1, ]
item_2 <- combined_results_I_U[combined_results_I_U$item == 2, ]
  item_2_N <- combined_results_I_NONU[combined_results_I_NONU$item == 2, ]
item_3 <- combined_results_I_U[combined_results_I_U$item == 3, ]
item_3_N <- combined_results_I_NONU[combined_results_I_NONU$item == 3, ]
item_4 <- combined_results_I_U[combined_results_I_U$item == 4, ]
item_4_N <- combined_results_I_NONU[combined_results_I_NONU$item == 4, ]
item_5 <- combined_results_I_U[combined_results_I_U$item == 5, ]
item_5_N <- combined_results_I_NONU[combined_results_I_NONU$item == 5, ]
item_6 <- combined_results_I_U[combined_results_I_U$item == 6, ]
item_6_N <- combined_results_I_NONU[combined_results_I_NONU$item == 6, ]
item_7 <- combined_results_I_U[combined_results_I_U$item == 7, ]
item_7_N <- combined_results_I_NONU[combined_results_I_NONU$item == 7, ]
item_8 <- combined_results_I_U[combined_results_I_U$item == 8, ]
item_8_N <- combined_results_I_NONU[combined_results_I_NONU$item == 8, ]
item_9 <- combined_results_I_U[combined_results_I_U$item == 9, ]
item_9_N <- combined_results_I_NONU[combined_results_I_NONU$item == 9, ]
item_10 <- combined_results_I_U[combined_results_I_U$item == 10, ]
item_10_N <- combined_results_I_NONU[combined_results_I_NONU$item == 10, ]

