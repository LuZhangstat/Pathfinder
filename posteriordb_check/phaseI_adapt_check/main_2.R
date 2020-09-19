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
library(gridExtra)
#source("../utils/sim.R")
source("../utils/sim_H3.R")
source("../utils/lp_utils.R")

set.seed(123)
pd <- pdb_local() # Posterior database connection
pn <- posterior_names(pd)
L_pn = length(pn)

# parameters settings #
alpha = 0.01
N = 40    # Maximum iters in optimization
mc.cores = 1 #parallel::detectCores() - 2
MC = mc.cores    # 10 iterations
init_bound = 2
width = 600; height = 500 # the size of the plot
seed_list = 1:MC
#takeoff <- c(21, 24)

## setting 1 ##
# iter <- 2
# max_treedepth <- 3
# stepsize <- 0.005
# M = 20

## setting 2 larger treedepth##
# iter <- 2
# max_treedepth <- 5
# stepsize <- 0.005
# M = 40
# N = 100

## setting 3 adapt stepsize##
# iter <- 3
# max_treedepth <- 4
# stepsize <- get_sampler_params(fit_0)[[1]][1, "stepsize__"] / 2^2
# M = 60
# N = 100

## setting 4 adapt stepsize##
# iter <- 1
# max_treedepth <- 6
# stepsize <- get_sampler_params(fit_0)[[1]][1, "stepsize__"] / 3*2
# M = 60
# N = 60
# 50% Center interval

## setting 5 adapt stepsize##
# iter <- 1
# max_treedepth <- 6
# stepsize <- get_sampler_params(fit_0)[[1]][1, "stepsize__"] / 3*2
# M = 60
# N = 60
# 80% Center interval


##setting 6 ##
# iter <- 1
# max_treedepth <- 10
# stepsize <- get_sampler_params(fit_0)[[1]][1, "stepsize__"] / 2

##setting 7 ##
# hamiltonian dynamic
# M = 4
# int_time = 6
# stepsize <- get_sampler_params(fit_0)[[1]][1, "stepsize__"] / 2
# lb = qbinom(0.1, (M * int_time), 0.5) / (M * int_time)
# ub = qbinom(0.9, (M * int_time), 0.5) / (M * int_time)

##setting 8##
# sim_H2
# int_time = 30

# N_models = 0
# model_record = c()
# for(l in 1:L_pn){
#   if(any(l == takeoff)){next}
#   modelname <- pn[l]
#   
#   # pick model
#   po <- posterior(modelname, pdb = pd)
#   # get reference posterior samples
#   skip_to_next <- FALSE
#   tryCatch(gsd <- reference_posterior_draws(po),
#            error = function(e) { skip_to_next <<- TRUE})
#   if(skip_to_next) {
#     # print("Error in obtaining reference posterior for this posterior.")
#     next }
#   N_models = N_models + 1
#   model_record = c(model_record, l)
#   printf("model %d: %s", l, modelname)
# }
# N_models

model_record <- c(1, 2, 3, 4, 6, 8, 11, 13, 15, 20, 23, 25, 26, 27, 29, 31, 33,
                  34, 36, 37, 40, 41, 43, 47, 51, 55, 61, 94, 95)
N_models <- length(model_record)


# preallocate results #
lp_INV <- array(data = NA, dim = c(2, length(model_record)))
lp_mean <- c()
lp_opath <- c()

for(i in 1:length(model_record)){
  modelname <- pn[model_record[i]]
  printf("model %d: %s", model_record[i], modelname)
  
  # pick model
  po <- posterior(modelname, pdb = pd)
  # get reference posterior samples
  gsd <- reference_posterior_draws(po)
  
  # compile the model
  sc <- stan_code(po)
  model <- stan_model(model_code = sc)
  # obtain posterior interval of lp__
  INV <- lp_Int_q_posteriordb(po, alpha)
  lp_INV[, i] = INV[1:2]
  lp_mean[i] = INV[3]
  ###  run Bob's Phase I  ###
  data <- get_data(po)
  opath <- opt_path_stan_parallel(seed_list, mc.cores, 
                                  model, data, N, init_bound)
  pick_records <- mclapply(opath, find_typical, model = model, 
                       data = data, mc.cores = mc.cores) 
  
  lp_opath[[i]] <- list(opath = opath, pick_records = pick_records)
}

test_ind <- find_typical(opath[[4]], model, data)
opath[[9]][test_ind, ncol(opath[[9]])]
init_param_unc <- opath[[1]][22, -ncol(opath[[1]])]
opath[[3]][39, ncol(opath[[1]])]
# tt <- get_sampler_params(fit_0)
# tt2 <- unlist(sapply(opath2, f <- function(x){ x[ , ncol(x)]}))
# hist(tt2[tt2>-2000])

save(file = "../results/lp_posteriordb_phI_adapt_set9.RData",
     list = c("lp_opath", "lp_INV", "model_record"))


# load("../results/lp_posteriordb_phI_adapt_set2.RData")

## check the plots ##
for(i in 1:length(model_record)){
  modelname <- pn[model_record[i]]
  printf("model %d: %s", model_record[i], modelname)
  #i = model_record[i]
  Lp_list = unlist(sapply(lp_opath[[i]]$opath, 
                          f <- function(x){
                            if(is.vector(x)){
                              return(1)
                            }else if(is.null(dim(x))){
                              return(NA)
                            }else{return(nrow(x))}}))
  chain_id = (1:MC)[!is.na(Lp_list)]
  p_lp_trace = 
    data.frame(iter = c(unlist(sapply(chain_id, f <- function(x){1:Lp_list[x]}))), 
               lp__ = c(unlist(sapply(lp_opath[[i]]$opath[chain_id], 
                                      f <- function(x){if(is.null(dim(x))){
                                        return(x[length(x)])
                                      }else{ return(x[, ncol(x)])}}))),
               chain = c(unlist(sapply(chain_id, f <- function(x){rep(x, Lp_list[x])}))),
               label = c(unlist(sapply(chain_id, f <- function(x){
                 if(Lp_list[x] == 1){
                   out = c("point on optim path")
                 } else{
                   out = rep("point on optim path", Lp_list[x])
                   out[lp_opath[[i]]$pick_records[[x]]$typical_index] = "marked point"
                   if(length(out) > Lp_list[x]){
                     out = rep("point on optim path", Lp_list[x]) # error in HMC sampling
                   }} 
                 return(out)}))))
  
  p_lp1 <- ggplot(data = p_lp_trace, 
                 aes(x=iter, y=lp__, group=chain, color=label)) +
    geom_line(colour = "grey") + geom_point(size = 2) + 
    geom_hline(yintercept = lp_INV[, i], colour = "black") + 
    geom_hline(yintercept = lp_mean[i], colour = "blue", linetype = 2)+
    ylim(min(p_lp_trace$lp__)
      #lp_INV[1, i] - 1.5*(lp_INV[2, i] - lp_INV[1, i])
       # min(lp_INV[1, i] - (lp_INV[2, i] - lp_INV[1, i]),
       #        quantile(p_lp_trace$lp__, 0.2))
      ,
         max(lp_INV[2, i] + (lp_INV[2, i] - lp_INV[1, i]), 
             max(p_lp_trace$lp__))) + ggtitle(paste("model:", modelname))+
    theme_bw() 
  
  p_lp2 <- ggplot(data = p_lp_trace, 
                 aes(x=iter, y=lp__, group=chain, color=label)) +
    geom_line(colour = "grey") + geom_point(size = 2) + 
    geom_hline(yintercept = lp_INV[, i], colour = "black") + 
    geom_hline(yintercept = lp_mean[i], colour = "blue", linetype = 2)+
    ylim(#min(p_lp_trace$lp__)
         lp_INV[1, i] - 1.5*(lp_INV[2, i] - lp_INV[1, i])
         # min(lp_INV[1, i] - (lp_INV[2, i] - lp_INV[1, i]),
         #        quantile(p_lp_trace$lp__, 0.2))
         ,
         max(lp_INV[2, i] + (lp_INV[2, i] - lp_INV[1, i]), 
             max(p_lp_trace$lp__))) + ggtitle(paste("model:", modelname))+
    theme_bw() 
  
  p_lp3 <- ggplot(data = p_lp_trace, 
                 aes(x=iter, y=lp__, group=chain, color=label)) +
    geom_line(colour = "grey") + geom_point(size = 2) + 
    geom_hline(yintercept = lp_INV[, i], colour = "black") + 
    geom_hline(yintercept = lp_mean[i], colour = "blue", linetype = 2)+
    ylim(#min(p_lp_trace$lp__)
         #lp_INV[1, i] - 1.5*(lp_INV[2, i] - lp_INV[1, i])
         min(lp_INV[1, i] - (lp_INV[2, i] - lp_INV[1, i]),
                quantile(p_lp_trace$lp__, 0.2))
         ,
         max(lp_INV[2, i] + (lp_INV[2, i] - lp_INV[1, i]), 
             max(p_lp_trace$lp__))) + ggtitle(paste("model:", modelname))+
    theme_bw() 
  
  jpeg(filename = paste0("../pics/phI_adapt/No", model_record[i], "-", 
                         modelname, ".jpeg"), #model_record[i]
       width = width*3, height = height, units = "px", pointsize = 12)
  grid.arrange(p_lp1, p_lp2, p_lp3, nrow = 1)
  dev.off()
  
}


## check the model with largest leapfrogs iterations
##No: 9: earnings-earn_height          ## similar to the fitting, not reach the INV
##No: 37: kilpisjarvi_mod-kilpisjarvi  ## zero
##No: 14: earnings-logearn_interaction ## almost zero
##No: 11: earnings-logearn_height_male ## almost zero
##No: 10: earnings-log10earn_height    ## zero
##No: 12: earnings-logearn_height      ## zero
##No: 48: mesquite-mesquite            ## multimodel? A lot
##No: 6: diamonds-diamonds             ## check 6 ??
##No: 3: bball_drive_event_0-hmm_drive_0 # multimodel. Also need checking
##No: 41: mcycle_gp-accel_gp           ## Not very good

i = which(model_record == 3)
p_lp <- ggplot(data = p_lp_trace, 
               aes(x=iter, y=lp__, group=chain, color=label)) +
  geom_line(colour = "grey") + geom_point(size = 1) + 
  geom_hline(yintercept = lp_INV[, i], colour = "black") + 
  ylim(lp_INV[1, i] - 10*(lp_INV[2, i] - lp_INV[1, i]),
       max(lp_INV[2, i] + (lp_INV[2, i] - lp_INV[1, i]), 
           max(p_lp_trace$lp__)))
p_lp


## errors in finding optimization path
## divergents in the MCMC sampling, stepsize
## slow

