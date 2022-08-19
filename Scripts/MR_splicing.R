library(TwoSampleMR)
library(tidyverse)
library(ggforestplot)

## Splicing Proteins on Breast Cancer
exposure_dat <- extract_instruments(c("prot-a-2452", "prot-a-814", "prot-a-2838",
                                      "prot-a-2839", "prot-a-2706", "prot-a-2707",
                                      "prot-a-3130"))


outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes = "ieu-a-1126")

dat <- harmonise_data(exposure_dat, outcome_dat)

res <- mr(dat) %>% generate_odds_ratios()

forestplot(df = res, name = exposure,
           estimate = b, se=se, logodds=T, xlab = "OR")

### Breast Cancer on splicing related proteins
exposure_dat <- extract_instruments("ieu-a-1126")
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes = "prot-a-814")
dat <- harmonise_data(exposure_dat, outcome_dat)
res <- mr(dat)
mr_scatter_plot(res, dat)
res_single <- mr_singlesnp(dat)
mr_forest_plot(res_single)

