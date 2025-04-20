rm(list = ls())

source("functions.R")

set.seed(2025)
samples_n <- c(10,20,30, 40, 50, 60, 70, 80, 100, 200, 400, 500, 1000,1200,2000)
samples_I <- c(10,20,30, 40, 50, 60, 70, 80, 100, 200)
tab <- c()
res_total <- list() 

for (sample.size in samples_I) {
  # základní nastavení, pro různé sample sizes zkoušíme jednoduchý setting
  # poměr reference:focal 1:1, 20 položek, N(0, 1) rozdělení traitu,
  # uniformní DIF o velikosti 1, 5% DIFových položek (tj. 1)
  res <- simul_total3(
    N = 1000,
    n_total =200, rat_n = c(1, 1), I = sample.size,
    mu_R = 0, mu_F = 0,
    type = c(0.05, 0), diffs_unif = 1, seed_arg = 2025,
   statistics =  list(MantelB = TRUE, MantelNUB = FALSE, BresB = FALSE, LogB = TRUE, SIBB = TRUE, cSIBB = FALSE))
  tab <- rbind(
    tab,
    cbind(
      N = sample.size,
      I = 20,
      DIF_type = "uniform",
      DIF_proportion = 0.05,
      DIF_size = 1,
      power = unlist(calculate_power_rate(res)),
      rejection = calculate_rejection_rate(res)
    )
  )
  res_total[[as.character(sample.size)]] <- res
}


# AH:
# 1. co se stane, když nějaká metoda nepůjde fitovat?
# 2. projděte poznámky v práci - jsou tam nějaké tipy na to, na co se v simulacích soustředit
# např. distribuce latent traits, nízký sample size, ...
# 3. zvolte různé sample sizes - např. od 50 do 5000 nějak v rozumném množství
# 4. jsou parametry pořád stejné? Toto je důležité - pokud je generujete uvnitř funkce bez seedu,
# pravděpodobně máte pro každý scénář/sample size (nebo dokonce simulaci?) jinou množinu parametrů,
# měly by být stejné - vygenerujte jednou a pak to použijte pro všechny stejné scénáře (např. stejný počet položek
# stejný typ DIFu, stejná velikost DIFu)
# 5. asi nemá smysl používat neuniformní metody na uniformní DIF - je možno je vypnout? Pokud si vyberu jen některé,
# funkce mi nefunguje.
# 6. funkce calculate_rejection_rate() špatně počítá rejection rate 
# - započítává tam i správně detekovanou položku (tj. bere všechny detekované položky)
# - mělo by se dělit počtem všech neDIFových položek (nyní se dělí počtem všech položek)
# 7. jak se určí, která položka bude DIFová?


#AB
#4. Solved? I have put another argument seed_arg to both generate_parameters and theta1, theta2. I am not truly sure,
# if the same seed should be put to theta1 and simultaneously theta2 (never mind, I noticed it is not correct)
#6. Solved (was confused for Type I error)
#5 Solved
#1 I think try(, silent = TRUE) + the if condition should prevent any cancellation of the simulation
#??? U rejection/power rate, neměl by denominator vyřazovat NA hodnoty?


#2
# Use purification and p correction
# NU particularly when the sample includes distinct high- and low-performing group


set.seed(2025)

# Expanded values
samples_n  <- c(10, 20, 30, 40, 50, 60, 70, 80, 100, 200, 400, 500, 1000, 1200, 2000)

items      <- c(10, 20, 30, 40, 50, 60, 70, 80, 100, 200)
ratios     <- list(c(1, 1), c(2, 1), c(3, 1))  # still three reference:focal ratios

# Build parameter grid
param_grid <- expand.grid(
  n_total = samples_n,
  I       = items,
  ratio   = seq_along(ratios)
)

# Initialize output
tab <- data.frame()
res_total <- list()

# Loop through all parameter combinations
for (row in seq_len(nrow(param_grid))) {
  n_total <- param_grid$n_total[row]
  I       <- param_grid$I[row]
  rat_n   <- ratios[[param_grid$ratio[row]]]
  
  # Run the simulation
  res <- simul_total3(
    N = 1000,
    n_total = n_total,
    rat_n   = rat_n,
    I       = I,
    mu_R = 0,
    mu_F = 0,
    type = c(0.05, 0),
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
  
  # Store result in list with descriptive key
  key <- paste0("n", n_total, "_I", I, "_rat", paste(rat_n, collapse = "_"))
  res_total[[key]] <- res
  
  # Add summary info to the table
  tab <- rbind(
    tab,
    data.frame(
      N_total = n_total,
      I       = I,
      rat_n1  = rat_n[1],
      rat_n2  = rat_n[2],
      DIF_type = "uniform",
      DIF_proportion = 0.05,
      DIF_size = 1,
      power = unlist(calculate_power_rate(res)),
      rejection = calculate_rejection_rate(res),
      row.names = NULL
    )
  )
}






# save results
save(tab, file = "results/res_I20_uniform_size1_5perc_normaltrait.RData")
