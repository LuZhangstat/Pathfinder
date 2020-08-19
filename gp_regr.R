rm(list = ls())
## install the beta release version of R package posterior
# install.packages("posterior", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# check the model: gp_pois_regr-gp_regr with new adaptation algorithm #

library(posteriordb)
pd <- pdb_default() # Posterior database connection
po <- posterior("gp_pois_regr-gp_regr", pdb = pd)
po

# access data and model
sc <- stan_code(po)
sc
dat <- get_data(po)

# Access gold standard posterior draws and information on how those were computed as follows.
rp_info <- reference_posterior_info(po, type = "draws")
rp_info
gsd <- reference_posterior_draws(po)
draws_df <- posterior::as_draws_df(gsd)
head(draws_df)

### test Bob's code ###
#library(cmdstanr)
library(rstan)
#library(posterior)
#library(bayesplot) # draw dis
set.seed(123)
m <- stan_model(model_code = sc)

opath <- opt_path_stan(m, data = list(N = dat$N, x = dat$x, y = dat$y), 
                       N = 75, init_bound = 5)
find_typical(m, data = list(N = dat$N, x = dat$x, y = dat$y), opath)



