rm(list = ls())

source("functions.R")

set.seed(2025)
#tab <- c()
#res_total <- list() 
#Useful for making some pretty graphs
samples_n <- c(10,20,30, 40, 50, 60, 70, 80, 100, 200, 400, 500, 1000,1200,2000)
samples_I <- c(10,20,30, 40, 50, 60, 70, 80, 100, 200)

#From article S&N 1990 
n_tot <- c(700, 1000, 1200, 1500) #should be trigerred together
ratios <- c(c(5,2),c(1,1),c(10,2),c(2,1))


for (sample.size in samples_I) {
  # základní nastavení, pro různé sample sizes zkoušíme jednoduchý setting
  # poměr reference:focal 1:1, 20 položek, N(0, 1) rozdělení traitu,
  # uniformní DIF o velikosti 1, 5% DIFových položek (tj. 1)
  res <- simul_total3(
    N = 1000,
    n_total =1000, rat_n = c(1, 1), I = sample.size,
    mu_R = 0, mu_F = 0,
    type = c(0.05, 0), diffs_unif = 1, seed_arg = 2025,
   statistics =  list(MantelB = TRUE, MantelNUB = FALSE, BresB = FALSE, LogB = TRUE, SIBB = TRUE, cSIBB = FALSE))
  tab <- rbind(
    tab,
    cbind(
      N = 1000,
      n = 1000,
      I = sample.size,
      DIF_type = "uniform",
      DIF_proportion = 0.05,
      DIF_size = 1,
      power = unlist(calculate_power_rate(res)),
      rejection = calculate_rejection_rate(res)
    )
  )
  res_total[[as.character(sample.size)]] <- res
}

calculate_itemwise_detection2(res_total$`10`, include_all_items = TRUE)
###############################################################################
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
    type = c(0.025, 0),
    diffs_unif = 1,
    seed_arg = 2025,
    statistics = list(
      MantelB = TRUE,
      MantelNUB = FALSE,
      BresB = FALSE,
      LogB = TRUE,
      SIBB = TRUE,
      cSIBB = FALSE
    )
  )
  
  tab2 <- rbind(
    tab2,
    cbind(
      N_total = sample.size,
      ratio_ref_foc = paste(ratio, collapse = ":"),
      I = I,
      DIF_type = "uniform",
      DIF_proportion = diffs_unif,
      DIF_size = 1,
      power = unlist(calculate_power_rate(res)),
      rejection = calculate_rejection_rate(res)
    )
  )
  
  res_total2[[paste0("N", sample.size, "_R", ratio[1], "_F", ratio[2])]] <- res
}


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





# save results
save(tab, file = "results/timing_results.RData")






#I should do for non uniform
