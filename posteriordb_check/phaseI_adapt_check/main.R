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
mc.cores = parallel::detectCores() - 2
M = mc.cores    # 20 iterations
init_bound = 2
width = 860; height = 740 # the size of the plot
seed_list = 1:M

# preallocate results #
lp_INV <- array(data = NA, dim = c(2, L_pn))
lp_opath <- c()
takeoff <- c(21, 24)

N_models = 0
model_record = c()
for(l in 1:L_pn){
  if(any(l == takeoff)){next}
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


for(l in 1:L_pn){
  if(all(l != model_record)){next}
  modelname <- pn[l]
  printf("model %d: %s", l, modelname)
  
  # pick model
  po <- posterior(modelname, pdb = pd)
  # get reference posterior samples
  gsd <- reference_posterior_draws(po)
  
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
}

# save(file = "../results/lp_posteriordb_phI_adapt.RData",
#      list = c("lp_opath", "lp_INV", "model_record"))
#

# load("../results/lp_posteriordb_phI_adapt.RData")

## check the plots ##
for(l in model_record){
  modelname <- pn[l]
  printf("model %d: %s", l, modelname)
  Lp_list = unlist(sapply(lp_opath[[l]]$opath, 
                          f <- function(x){
                            if(is.vector(x)){
                              return(1)
                            }else if(is.null(dim(x))){
                              return(NA)
                            }else{return(nrow(x))}}))
  chain_id = (1:M)[!is.na(Lp_list)]
  p_lp_trace = 
    data.frame(iter = c(unlist(sapply(chain_id, f <- function(x){1:Lp_list[x]}))), 
               lp__ = c(unlist(sapply(lp_opath[[l]]$opath[chain_id], 
                                      f <- function(x){if(is.null(dim(x))){
                                        return(x[length(x)])
                                      }else{ return(x[, ncol(x)])}}))),
               chain = c(unlist(sapply(chain_id, f <- function(x){rep(x, Lp_list[x])}))),
               label = c(unlist(sapply(chain_id, f <- function(x){
                 if(Lp_list[x] == 1){
                   out = c("not preferred init")
                 } else{
                   out = rep("not preferred init", Lp_list[x])
                   out[lp_opath[[l]]$pick_ind[[x]]] = "preferred init"
                   if(length(out) > Lp_list[x]){
                     out = rep("error in finding init", Lp_list[x])
                   }} 
                 return(out)}))))
  
  p_lp <- ggplot(data = p_lp_trace, 
                 aes(x=iter, y=lp__, group=chain, color=label)) +
    geom_line(colour = "grey") + geom_point(size = 1) + 
    geom_hline(yintercept = lp_INV[, l], colour = "black") + 
    ylim(min(lp_INV[, l] - (lp_INV[2, l] - lp_INV[1, l]), 
             quantile(p_lp_trace$lp__, 0.1)),
         max(lp_INV[2, l] + (lp_INV[2, l] - lp_INV[1, l]), 
             max(p_lp_trace$lp__)))
  
  jpeg(filename = paste0("../pics/phI_adapt/No",l,"-", modelname, ".jpeg"),
       width = width, height = height, units = "px", pointsize = 12)
  print(p_lp)
  dev.off()
}


