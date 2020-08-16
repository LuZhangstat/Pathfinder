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
dat <- get_data(po)

# Access gold standard posterior draws and information on how those were computed as follows.
rp_info <- reference_posterior_info(po, type = "draws")
rp_info
gsd <- reference_posterior_draws(po)
draws_df <- posterior::as_draws_df(gsd)
head(draws_df)

### test with adaptive algorithm from bob ###
library(cmdstanr)
#library(rstan)
library(posterior)
#library(bayesplot) # draw dis
set.seed(123)
# m <- stan_model(model_code = sc)
# f <- optimizing(m, data = list(N = dat$N, x = dat$x, y = dat$y),
#                 verbose = TRUE, iter = 2000, save_iterations = TRUE,
#                 refresh = 1)
file <- file.path("./stat_comp_benchmarks/benchmarks/gp_regr/gp_regr.stan")
mod <-  cmdstan_model(file)
# check the optimization algorithm:
fit_mle <- mod$optimize(data = list(N = dat$N, x = dat$x, y = dat$y), 
                        seed = 123,
                        save_latent_dynamics = TRUE,
                        output_dir = "./")
# the exact value of parameters are rho = 5.5, alpha = 3, sigma = 2
fit_mle$summary()
fit_mle$mle("rho")
# mcmc_hist(fit$draws("rho")) +
#   vline_at(fit_mle$mle(), size = 1.5)

## test with campfire from ben ##
devtools::install_github("bbbales2/campfire")
library(campfire)
library(rstan)

