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

posterior::summarize_draws(gsd)

### obtain the histogram of lp__ ##
library(rstan)
source("./sim.R")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
INV <- lp_Int_q_posteriordb(po, alpha = 0.01)
INV
#    0.5%        99.5% 
#   -149854.7079 -227.5184
