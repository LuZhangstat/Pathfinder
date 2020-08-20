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

m <- stan_model(model_code = sc)
INV <- lp_Int_q_posteriordb(po, alpha = 0.01)
INV
#    0.5%        99.5% 
#   -33.84823    -27.63794

###  run Stan with a long Phase I warmup time  ###
set.seed(123)
L = 2000
phI_sample <- sampling(m, data = get_data(po), 
                       seed = 1234,
                       iter = L + 1, 
                       warmup = L,
                       chains = 4, 
                       cores = 4,
                       algorithm ="NUTS",
                       control = list(adapt_init_buffer = L,
                                      adapt_term_buffer = 0,
                                      adapt_window = 0),
                       save_warmup = TRUE)

###  record the number of iterations required to find an lp__ value  ###
list_of_draws <- extract(phI_sample)
str(list_of_draws)
summary(list_of_draws$lp__)
summary(list_of_draws$rho)
summary(list_of_draws$alpha)
summary(list_of_draws$sigma)

