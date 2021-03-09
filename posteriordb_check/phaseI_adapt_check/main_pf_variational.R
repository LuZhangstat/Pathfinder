rm(list = ls())
setwd("./posteriordb") # set working dir to cloned package
library(rstan)
library(parallel)
library(foreach)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
## install the beta release version of R package posterior
# install.packages("posterior", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(posteriordb)
library(posterior)
source("../utils/sim_variational_pf_mean.R")
source("../utils/lp_utils.R")


pd <- pdb_local() # Posterior database connection
pn <- posterior_names(pd)
L_pn = length(pn)

# parameters settings 22#
alpha = 0.01
mc.cores = parallel::detectCores() - 2
MC = 20
init_bound = 2
width = 800; height = 700 # the size of the plot
seed_list = 1:MC 

N1 = 200    # Maximum iters in optimization
N_mode_max = 1
N_sam = 3   # ELBO
N_mass = 5
factr_tol = 1e2
MC = 20


load(file = "../results/lp_posteriordb_LBFGS.RData")
set.seed(123)

# preallocate results #
lp_opath <- c()

which(model_record == 48)
t_0 <- proc.time()
for(i in 1:49){ #length(model_record)
  modelname <- pn[model_record[i]]
  printf("model %d: %s", model_record[i], modelname)
  
  # pick model
  po <- posterior(modelname, pdb = pd)
  # get reference posterior samples
  gsd <- reference_posterior_draws(po)
  
  # compile the model
  sc <- stan_code(po)
  model <- stan_model(model_code = sc)
  
  ###  get the data  ###
  data <- get_data(po)
  
  # set up lmm #
  posterior <- to_posterior(model, data)
  D <- get_num_upars(posterior)
  lmm = min(max(D, 5), N1)
  cat("No. pars:", D," lmm in L-BFGS: ", lmm, "\n")
  
  # t <- proc.time()
  # opath <- opt_path_stan_parallel(seed_list, seed_list, mc.cores, model, data,
  #                                 N1, N_mode_max, N_sam, N_mass,
  #                                 init_bound = 15.0, factr_tol, lmm) # plot for 8school init_bound = 15
  # proc.time() - t

  t <- proc.time()
  opath <- opt_path_stan_init_parallel(
    initial_ls[[i]], mc.cores, model, data, N1, N_mode_max, N_sam, N_mass,
    init_bound, factr_tol, lmm, seed_list)
  print(proc.time() - t)
  # opath <- lp_opath[[i]]$opath
  
  # Check the scaled log prob mass and keep the index with relatively high prob mass
  pick_mode1 <- filter_mode(opath, fit_info)
  pick_mode2 <- filter_mode(opath, fit_info_DIV)
  pick_mode <- intersect(pick_mode1, pick_mode2)
  pick_samples <- filter_samples(opath[pick_mode])
  
  lp_opath[[i]] <- list(opath = opath, pick_mode = pick_mode,
                        pick_samples = pick_samples)
  
}
proc.time() - t_0

save(file = "../results/lp_posteriordb_phI_adapt_set22.RData",
     list = c("lp_opath"))

# load("../results/lp_posteriordb_phI_adapt_set22_center.RData")



 
  
  