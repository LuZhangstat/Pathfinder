rm(list = ls())
## install the beta release version of R package posterior
# install.packages("posterior", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# check the dataset in posteriordb #

library(posteriordb)
library(posterior)
pd <- pdb_default() # Posterior database connection
pn <- posterior_names(pd)
head(pn)

# pick model: gp_pois_regr-gp_regr
po <- posterior("gp_pois_regr-gp_regr", pdb = pd)
po

# access data and model
sc <- stan_code(po)
sc

# Check data 
dat <- get_data(po)
dat
# The same dataset as "gp_regr.data.R" in stat_comp_benchmarks

# Access gold standard posterior draws and information on how those were computed as follows.
rp_info <- reference_posterior_info(po, type = "draws")
rp_info

gsd <- reference_posterior_draws(po)
gsd

# Error
# draws_df <- posterior::as_draws_df(gsd$draws)
# head(draws_df)
