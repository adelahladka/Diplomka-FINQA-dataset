rm(list = ls())

source("functions.R")

set.seed(123)

tab <- c()
for (sample.size in c(500, 1000)) {
  # základní nastavení, pro různé sample sizes zkoušíme jednoduchý setting
  # poměr reference:focal 1:1, 20 položek, N(0, 1) rozdělení traitu,
  # uniformní DIF o velikosti 1, 5% DIFových položek (tj. 1)
  res <- simul_total3(
    N = 100,
    n_total = sample.size, rat_n = c(1, 1), I = 20,
    mu_R = 0, mu_F = 0,
    type = c(0.05, 0), diffs_unif = 1
  ) # ,
  # statistics = list(MantelB = TRUE, LogB = TRUE, SIBB = TRUE))
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
}

tab 

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


# save results
save(tab, file = "results/res_I20_uniform_size1_5perc_normaltrait.RData")
