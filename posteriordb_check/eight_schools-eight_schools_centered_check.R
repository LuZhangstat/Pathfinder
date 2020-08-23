rm(list = ls())
## install the beta release version of R package posterior
# install.packages("posterior", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# check the dataset in posteriordb #

library(posteriordb)
library(posterior)
pd <- pdb_default() # Posterior database connection
pn <- posterior_names(pd)
head(pn)

# pick model: eight_schools-eight_schools_centered
po <- posterior("eight_schools-eight_schools_centered", pdb = pd)
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

### obtain the posterior interval of lp__ ##
library(rstan)
source("./utils/sim.R")
source("./utils/lp_utils.R")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

model <- stan_model(model_code = sc)
# lps <- lp_recover(model, data = dat, pos_draws = gsd)
alpha = 0.01
INV <- lp_Int_q_posteriordb(po, alpha)
INV
# 0.5%       99.5% 
# -29.43469  20.07446

###  run Stan with a long Phase I warmup time  ###
set.seed(123)
L = 2000
M = 12
mc.cores = 2
phiI_sample <- sampling(model, data = get_data(po), 
                       seed = 1234,
                       iter = L + 1, 
                       warmup = L,
                       chains = M, 
                       cores = mc.cores,
                       algorithm ="NUTS",
                       control = list(adapt_init_buffer = L,
                                      adapt_term_buffer = 0,
                                      adapt_window = 0),
                       save_warmup = TRUE)

###  record the number of iterations required to reach INV ###
lp_phI <- ls_lp_phI(phiI_sample, L)
summary(lp_phI)

# check the trace plot of lp__ #
L_p = 20
p_lp_trace = data.frame(iter = rep(1:L_p, M), 
                        chain = rep(paste(1:M), each = L_p),
                  lp__ = c(lp_phI[1:L_p, ]))
library(ggplot2)
p_lp <- ggplot(data = p_lp_trace, 
               aes(x=iter, y=lp__, group=chain, color=chain)) + geom_line() +
  geom_hline(yintercept = INV)
p_lp

# Get the number of iterations and  leapfrogs #
lp_explore_sum <- lp_explore(phiI_sample, INV, L, M)
lp_explore_sum$n_iters;
lp_explore_sum$n_sum_leapfrog
mean(lp_explore_sum$n_sum_leapfrog); max(lp_explore_sum$n_sum_leapfrog)
# [1] 24.08333
# [1] 63


