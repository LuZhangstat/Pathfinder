rm(list = ls())

## test with campfire from ben ##
# devtools::install_github("bbbales2/campfire")
# devtools::load_all()
# devtools::build()
library(cmdstanr)
library(posterior)
library(campfire)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
model_path = "../campfire_test/accel_splines.stan"

data_env = new.env(parent = baseenv())
source("../campfire_test/accel_splines.data.R", local = data_env)
data = as.list.environment(data_env)

stan_fit = stan(model_path, data = data, iter = 1)
out = warmup(model_path, stan_fit, data = data, num_chains = 4,
             parallel_chains = 4, print_stdout = TRUE)

# And you can take the output of a warmup call and use it to sample the actual
# model. It’ll use the metric computed in campfire::warmup, but it’ll still do
# 50 steps of timestep adaptation at the beginning:

model = cmdstanr::cmdstan_model(model_path)
fit = do.call(model$sample, out$args)
