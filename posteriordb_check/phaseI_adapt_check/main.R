rm(list = ls())
setwd("./posteriordb") # set working dir to cloned package
library(rstan)
library(parallel)
library(foreach)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
## install the beta release version of R package posterior
# install.packages("posterior", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# check the dataset in posteriordb #
library(posteriordb)
library(posterior)
library(ggplot2)
source("../utils/sim.R")
source("../utils/lp_utils.R")

set.seed(123)
pd <- pdb_local() # Posterior database connection
pn <- posterior_names(pd)
L_pn = length(pn)

# parameters settings #
alpha = 0.01
N = 75    # Maximum iters in optimization
M = 6    # 20 iterations
init_bound = 2
width = 860; height = 740 # the size of the plot
mc.cores = parallel::detectCores() - 2
seed_list = 1:M

# preallocate results #
#lp_explore_n_iters <- array(data = NA, dim = c(M, L_pn))
#lp_explore_n_leapfrog <- array(data = NA, dim = c(M, L_pn))
lp_INV <- array(data = NA, dim = c(2, L_pn))
lp_opath <- c()

for(l in 1:L_pn){
  modelname <- pn[l]
  printf("model %d: %s", l, modelname)
  # pick model
  po <- posterior(modelname, pdb = pd)
  # get reference posterior samples
  skip_to_next <- FALSE
  tryCatch(gsd <- reference_posterior_draws(po), 
           error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next) { 
    print("Error in obtaining reference posterior for this posterior.")
    next }  
  
  # compile the model
  sc <- stan_code(po)
  model <- stan_model(model_code = sc)
  # obtain posterior interval of lp__
  INV <- lp_Int_q_posteriordb(po, alpha)
  lp_INV[, l] = INV
  ###  run Bob's Phase I  ###
  data <- get_data(po)
  opath <- opt_path_stan_parallel(seed_list, mc.cores, 
                                  model, data, N, init_bound)
  pick_ind <- mclapply(opath, find_typical, model = model, 
                       data = data, mc.cores = mc.cores) 
  lp_opath[[l]] <- list(opath = opath, pick_ind = pick_ind)
  
  # check the trace plot of lp__ #
  L_p = nrow(opath)
  label = rep("grey", L_p); label[pick_ind] = "green";
  p_lp_trace = data.frame(iter = 1:L_p, lp__ = opath[, ncol(opath)],
                          label = label)
  p_lp <- ggplot(data = p_lp_trace, 
                 aes(x=iter, y=lp__)) + geom_line() + 
    geom_point(colour = label) + 
    geom_hline(yintercept = INV) 
  jpeg(filename = paste0("../pics/phI_adapt/No",l,"-", modelname, ".jpeg"),
       width = width, height = height, units = "px", pointsize = 12)
  print(p_lp)
  dev.off()
}

# save(file = "../results/lp_posteriordb_explore.RData", 
#      list = c("lp_explore_n_iters", "lp_explore_n_leapfrog",
#               "lp_INV"))
# 
# load("../results/lp_posteriordb_explore.RData")
## check reference posterior ##

N_models = 0
model_record = c()
for(l in 1:L_pn){
  modelname <- pn[l]
  # printf("model %d: %s", l, modelname)

  # pick model
  po <- posterior(modelname, pdb = pd)
  # get reference posterior samples
  skip_to_next <- FALSE
  tryCatch(gsd <- reference_posterior_draws(po),
           error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next) {
    # print("Error in obtaining reference posterior for this posterior.")
    next }
  N_models = N_models + 1
  model_record = c(model_record, l)
}
N_models
# only 49 out of 97 models have reference posterior samples

## check the transformed parameters block ##
# for(id in model_record){
id = 24
  cat("id:", id)
  po <- posterior(pn[id], pdb = pd)
  sc <- stan_code(po)
  print(sc)
  gsd <- reference_posterior_draws(po)
  tt <- sapply(gsd[[1]], unlist)
  tt[1, ]
  
  model <- stan_model(model_code = sc)
  data <- get_data(po)
  posterior <- to_posterior(model, data)
  get_inits(posterior)[[1]] 
  # readline(prompt="Press [enter] to continue")
# }
# 21, 24 does not match
takeoff <- c(21, 24)

