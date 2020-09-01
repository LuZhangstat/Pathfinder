rm(list = ls())
## install the beta release version of R package posterior
# install.packages("posterior", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# check the dataset in posteriordb #

library(posteriordb)
library(posterior)
pd <- pdb_default() # Posterior database connection
pn <- posterior_names(pd)

# pick model 41: mcycle_gp-accel_gp
po <- posterior("mcycle_gp-accel_gp", pdb = pd)

# access data and model
sc <- stan_code(po)
sc

# Check data 
dat <- get_data(po)
# The same dataset as "gp_regr.data.R" in stat_comp_benchmarks

# Access gold standard posterior draws and information on how those were computed as follows.
rp_info <- reference_posterior_info(po, type = "draws")
rp_info

gsd <- reference_posterior_draws(po)

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
# 0.5%     99.5% 
#   -668.4108 -633.6211


###  run Stan with a long Phase I warmup time  ###
set.seed(123)
L = 2000
M = 20
mc.cores = parallel::detectCores() - 2
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

# check the trace plot of lp__, take chain 6 off #
L_p1 = 16
L_p2 = 1000
p_lp_trace = data.frame(iter = rep(L_p1:L_p2, M-1), 
                        chain = rep(paste((1:M)[-6]), each = L_p2 - L_p1 + 1),
                  lp__ = c(lp_phI[L_p1:L_p2, -6]))
library(ggplot2)
p_lp <- ggplot(data = p_lp_trace, 
               aes(x=iter, y=lp__, group=chain, color=chain)) + geom_line() +
  geom_hline(yintercept = INV)
p_lp
width = 860; height = 740
jpeg(filename = paste0("./pics/No41-mcycle_gp-accel_gp_2.jpeg"),
     width = width, height = height, units = "px", pointsize = 12)
print(p_lp)
dev.off()

# Get the number of iterations and  leapfrogs #
lp_explore_sum <- lp_explore(phiI_sample, INV, L, M)
lp_explore_sum$n_iters
# [1]   21   79   23   22  118 2000   34   57  402   33   18   18   13  137   19   23
# [17]   21   27   21   27
lp_explore_sum$n_sum_leapfrog
# [1]     303    1711     640     193    2846 2037979     922    1317    8257     922
# [11]    1370     468     229    3183     236     311     723     708     493     859
mean(lp_explore_sum$n_sum_leapfrog); max(lp_explore_sum$n_sum_leapfrog)
# [1] 458.5
# [1] 1483


