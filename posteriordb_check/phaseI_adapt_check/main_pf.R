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
source("../utils/sim_pf.R")
source("../utils/lp_utils.R")

pd <- pdb_local() # Posterior database connection
pn <- posterior_names(pd)
L_pn = length(pn)

mc.cores = parallel::detectCores() - 2 # number of parallel cores
width = 800; height = 700 # the size of the plot

## tuning parameters (default) ##
init_bound = 2.0 # parameter for initial distribution 
N1 = 1000    # maximum iters in optimization
factr_tol = 1e2 # relative tolerance = 1-4 is not enough, should use at least 1e7
N_sam_DIV = 5   # samples for ELBO evaluation
N_sam = 100
lmm = 6 # maximum histogram size
# single pathfinder repeat for 100 times #
MC = 100
seed_list = 1:MC 

## sensitivity test ##
# shorter optimization path #
sen_test_L <- FALSE
if(sen_test_L){
  N1 = 200    # maximum iters in optimization
  factr_tol = 1e7 # relative tolerance = 1-4 is not enough, should use at least 1e7
}

# larger sample size for ELBO evaluation #
sen_test_K <- FALSE
if(sen_test_K){
  N_sam_DIV = 30  # samples for ELBO evaluation
}

# larger history #
sen_test_J <- FALSE
if(sen_test_J){
  lmm = 60
}

# load precalculated results #
load(file = "../results/lp_posteriordb_LBFGS_h6.RData")
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
  cat("No. pars:", D," lmm in L-BFGS: ", lmm, "\n")
  
  t <- proc.time()
  opath <- opt_path_stan_parallel(seed_list, seed_list, mc.cores, model, data,
                                  init_bound = init_bound, N1, N_sam_DIV, N_sam, 
                                  factr_tol, lmm) # plot for 8school init_bound = 15
  print(proc.time() - t)
  # opath <- lp_opath[[i]]$opath
  
  pick_samples <- random_sample_Each(opath, seed = 1) # randomly pick one approximating draws for each run of Pathfinder
  # pick_samples <- Imp_Resam_Each(opath, seed = 1)   # use PSIS-IR to pick one approximating draws for each run of Pathfinder (no longer used in the paper)
  # pick_samples <- Imp_Resam_WOR(opath, n_inits = 20, seed = 1) # use PSIS-IR to pick n_inits distinct approximating draws for multi-path Pathfinder

  lp_opath[[i]] <- list(opath = opath, pick_samples = pick_samples)
  
}
proc.time() - t_0

save(file = "../results/lp_posteriordb_phI_adapt_long_hist.RData",
     list = c("lp_opath"))
# _default
# _short_L
# _large_K
# _long_hist

lapply(lp_opath, f <- function(x){ncol(x$pick_samples)})


## multi-path Pathfinder ##
# preallocate results #
lp_multi_opath <- c()

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
  cat("No. pars:", D," lmm in L-BFGS: ", lmm, "\n")
  
  lp_multi_opath[[i]] <- list()
  t <- proc.time()
  for(j in 1:100){
    cat(j, "\t")
    seed_list = (j-1)*20 + 1:20
    opath <- opt_path_stan_parallel(seed_list, seed_list, mc.cores, model, data,
                                    init_bound = init_bound, N1, N_sam_DIV, N_sam, 
                                    factr_tol, lmm)
    pick_samples <- Imp_Resam_WOR(opath, n_inits = N_sam, seed = 1)
    lp_multi_opath[[i]][[j]] <- pick_samples   # record the 100 samples 
  }
  print(proc.time() - t)
  
}
proc.time() - t_0

save(file = "../results/multi_pf_samples_default.RData",
     list = c("lp_multi_opath"))


############### old code #######################
# random sample 4 pathfinder + SIR WOR #
load(file = "../results/lp_posteriordb_LBFGS_h6.RData")
load("../results/lp_posteriordb_phI_adapt_large_K.RData")
set.seed(123)

pf_SIR_WOR_4_draws <- list()
pf_4_each <- list()

which(model_record == 48)
t_0 <- proc.time()
for(i in 1:49){ #length(model_record)
  modelname <- pn[model_record[i]]
  printf("model %d: %s", model_record[i], modelname)
  
  # pick model
  
  pf_SIR_WOR_4_draws[[i]] <- list()
  pf_4_each[[i]] <- list()
  for(j in 1:25){
    pf_4_each[[i]][[j]] <- lp_opath[[i]]$
      pick_samples[, (((j - 1) * 4 + 1):(j * 4))]
    set.seed(j)
    pf_SIR_WOR_4_draws[[i]][[j]] <- 
      Imp_Resam_WOR(lp_opath[[i]]$opath[((j - 1) * 4 + 1):(j * 4)],
                    n_inits = 4, seed = j)
    
  }
}
proc.time() - t_0

save(file = "../results/pf_SIR_WOR_4_draws_large_K.RData",
     list = c("pf_SIR_WOR_4_draws", "pf_4_each"))
  
  