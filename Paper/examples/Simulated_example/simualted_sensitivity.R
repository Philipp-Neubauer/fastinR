FAP.predators <- "even_good_FA_preds.csv"
FAP.preys <- "even_good_FA_preys.csv"
Conv.Coeffs.mean <- "even_good_FA_cc_means.csv"
Conv.Coeffs.var <- "even_good_FA_cc_var.csv"
fat.conts <- "even_good_fat_cont.csv"
even_good <- add_FA(FA.predators = FAP.predators,
FA.preys = FAP.preys,
fat.conts = fat.conts,
Conv.Coeffs.mean = Conv.Coeffs.mean,
Conv.Coeffs.var = Conv.Coeffs.var)
FAP.predators <- "uneven_good_FA_preds.csv"
FAP.preys <- "uneven_good_FA_preys.csv"
Conv.Coeffs.mean <- "uneven_good_FA_cc_means.csv"
Conv.Coeffs.var <- "uneven_good_FA_cc_var.csv"
fat.conts <- "uneven_good_fat_cont.csv"
uneven_good <- add_FA(FA.predators = FAP.predators,
FA.preys = FAP.preys,
fat.conts = fat.conts,
Conv.Coeffs.mean = Conv.Coeffs.mean,
Conv.Coeffs.var = Conv.Coeffs.var)
FAP.predators <- "even_bad_FA_preds.csv"
FAP.preys <- "even_bad_FA_preys.csv"
Conv.Coeffs.mean <- "even_bad_FA_cc_means.csv"
Conv.Coeffs.var <- "even_bad_FA_cc_var.csv"
fat.conts <- "even_bad_fat_cont.csv"
even_bad <- add_FA(FA.predators = FAP.predators,
FA.preys = FAP.preys,
fat.conts = fat.conts,
Conv.Coeffs.mean = Conv.Coeffs.mean,
Conv.Coeffs.var = Conv.Coeffs.var)

bad_prop <- read.csv('even_bad_props.csv',h=F)[,2:5]
good_prop <- read.csv('even_good_props.csv',h=F)[,2:5]
uneven_prop <- read.csv('uneven_good_props.csv',h=F)[,2:5]

save(list=c("uneven_prop","good_prop","bad_prop","uneven_good","even_good","even_bad"),file='sensitivity.RData')

