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
source("../utils/sim_RHR.R")
source("../utils/lp_utils.R")
source("../utils/L_RHR.R")

set.seed(123)
pd <- pdb_local() # Posterior database connection
pn <- posterior_names(pd)
L_pn = length(pn)

# parameters settings #
alpha = 0.01
N1 = 60    # Maximum iters in optimization
N_mode_max = 20 
N_sam = 100
mc.cores = parallel::detectCores() - 2
MC = mc.cores    # 10 iterations
MC = 3
init_bound = 2
width = 600; height = 500 # the size of the plot
seed_list = 1:MC 


##setting 11##

model_record <- c(1, 2, 3, 4, 6, 8, 11, 13, 15, 20, 23, 25, 26, 27, 29, 31, 33,
                  34, 36, 37, 40, 41, 43, 47, 51, 55, 61, 94, 95)
N_models <- length(model_record)


# preallocate results #
lp_INV <- array(data = NA, dim = c(2, length(model_record)))
lp_mean <- c()
lp_opath <- c()

for(i in 1:length(model_record)){ #length(model_record)
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
  
  #opt_path_stan(model, data, N1 = 30, N_rep = 6, init_bound = 2)
  opath <- opt_path_stan_parallel(seed_list, mc.cores, model, data, 
                                  N1, N_mode_max, N_sam, init_bound)
  pick_ind <- mclapply(opath, find_indx, mc.cores = mc.cores)
  lp_opath[[i]] <- list(opath = opath, pick_ind = pick_ind)
}

param_path <- opath[[1]]
test_ind <- find_typical(opath[[1]], model, data)
opath[[9]][test_ind, ncol(opath[[9]])]
init_param_unc <- opath[[1]][22, -ncol(opath[[1]])]
opath[[3]][39, ncol(opath[[1]])]
# tt <- get_sampler_params(fit_0)
# tt2 <- unlist(sapply(opath2, f <- function(x){ x[ , ncol(x)]}))
# hist(tt2[tt2>-2000])

# save(file = "../results/lp_posteriordb_phI_adapt_set9.RData",
#      list = c("lp_opath", "lp_INV", "model_record"))


# load("../results/lp_posteriordb_phI_adapt_set9.RData")

## check the plots ##
for(i in 1:length(model_record)){ #length(model_record)
  modelname <- pn[model_record[i]]
  printf("model %d: %s", model_record[i], modelname)
  #i = model_record[i]
  Lp_list = unlist(sapply(lp_opath[[i]]$opath, 
                          f <- function(x){
                            if(length(x) == 1){
                              return(NA)
                            }else if(is.vector(x$y)){
                              return(1)
                            }else{return(nrow(x$y))}}))
  chain_id = (1:MC)[!is.na(Lp_list)]
  p_lp_trace = 
    data.frame(iter = c(unlist(sapply(chain_id, f <- function(x){1:Lp_list[x]}))), 
               lp__ = c(unlist(sapply(lp_opath[[i]]$opath[chain_id], 
                                      f <- function(x){if(is.null(dim(x$y))){
                                        return(x$y[length(x$y)])
                                      }else{ return(x$y[, ncol(x$y)])}}))),
               chain = c(unlist(sapply(chain_id, f <- function(x){rep(x, Lp_list[x])}))),
               label = c(unlist(sapply(chain_id, f <- function(x){
                 if(Lp_list[x] == 1){
                   out = c("point on optim path")
                 } else{
                   out = rep("point on optim path", Lp_list[x])
                   out[lp_opath[[i]]$pick_ind[[x]]] = "marked point"
                   if(length(out) > Lp_list[x]){
                     out = rep("point on optim path", Lp_list[x]) # error in HMC sampling
                   }} 
                 return(out)}))))
  
  p_lp <- ggplot(data = p_lp_trace, 
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
             max(p_lp_trace$lp__))) + 
    ggtitle(paste("model:", modelname,"\n", "estimated E(lp):", 
                  paste(sapply(lp_opath[[i]]$opath[chain_id], 
                         f <- function(x){ 
                           if(is.na(x$E_lp[length(x$E_lp)])){"NA"}else{
                             round(x$E_lp[length(x$E_lp)], digits = 1)
                           }}), 
                        collapse = " "), 
                  "VS E(lp):", round(lp_mean[i], digits = 1),
                  "\n", "counts in line search:", 
                  paste(sapply(lp_opath[[i]]$opath[chain_id], 
                               f <- function(x){ 
                                 if(is.na(x$E_lp[length(x$E_lp)])){"NULL"}else{
                                   sum(x$step_count)
                                 }}), 
                        collapse = " "), "No. of sample",
                  paste(sapply(lp_opath[[i]]$opath[chain_id], 
                               f <- function(x){ 
                                 if(is.na(x$E_lp[1])){"0"}else{
                                   N_sam * length(x$E_lp)
                                 }}), 
                        collapse = " "))) +
    theme_bw() 
  
  jpeg(filename = paste0("../pics/phI_adapt_setting11/No", model_record[i], "-", 
                         modelname, "_L.jpeg"), #model_record[i]
       width = width, height = height, units = "px", pointsize = 12)
  print(p_lp)
  dev.off()
  
  p_lp <- ggplot(data = p_lp_trace, 
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
             max(p_lp_trace$lp__))) + 
    ggtitle(paste("model:", modelname,"\n", "estimated E(lp):", 
                  paste(sapply(lp_opath[[i]]$opath[chain_id], 
                               f <- function(x){ 
                                 if(is.na(x$E_lp[length(x$E_lp)])){"NA"}else{
                                   round(x$E_lp[length(x$E_lp)], digits = 1)
                                 }}), 
                        collapse = " "), 
                  "VS E(lp):", round(lp_mean[i], digits = 1),
                  "\n", "counts in line search:", 
                  paste(sapply(lp_opath[[i]]$opath[chain_id], 
                               f <- function(x){ 
                                 if(is.na(x$E_lp[length(x$E_lp)])){"NULL"}else{
                                   sum(x$step_count)
                                 }}), 
                        collapse = " "), "No. of sample",
                  paste(sapply(lp_opath[[i]]$opath[chain_id], 
                               f <- function(x){ 
                                 if(is.na(x$E_lp[1])){"0"}else{
                                   N_sam * length(x$E_lp)
                                 }}), 
                        collapse = " "))) +
    theme_bw() 
  
  jpeg(filename = paste0("../pics/phI_adapt_setting11/No", model_record[i], "-", 
                         modelname, "_s.jpeg"), #model_record[i]
       width = width, height = height, units = "px", pointsize = 12)
  print(p_lp)
  dev.off()
  
  p_lp <- ggplot(data = p_lp_trace, 
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
             max(p_lp_trace$lp__))) + 
    ggtitle(paste("model:", modelname,"\n", "estimated E(lp):", 
                  paste(sapply(lp_opath[[i]]$opath[chain_id], 
                               f <- function(x){ 
                                 if(is.na(x$E_lp[length(x$E_lp)])){"NA"}else{
                                   round(x$E_lp[length(x$E_lp)], digits = 1)
                                 }}), 
                        collapse = " "), 
                  "VS E(lp):", round(lp_mean[i], digits = 1),
                  "\n", "counts in line search:", 
                  paste(sapply(lp_opath[[i]]$opath[chain_id], 
                               f <- function(x){ 
                                 if(is.na(x$E_lp[length(x$E_lp)])){"NULL"}else{
                                   sum(x$step_count)
                                 }}), 
                        collapse = " "), "No. of sample",
                  paste(sapply(lp_opath[[i]]$opath[chain_id], 
                               f <- function(x){ 
                                 if(is.na(x$E_lp[1])){"0"}else{
                                   N_sam * length(x$E_lp)
                                 }}), 
                        collapse = " "))) +
    theme_bw() 
  
  jpeg(filename = paste0("../pics/phI_adapt_setting11/No", model_record[i], "-", 
                         modelname, "_trun.jpeg"), #model_record[i]
       width = width, height = height, units = "px", pointsize = 12)
  print(p_lp)
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

i = which(model_record == 36)
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

# > lp_INV[2, ] - lp_INV[1, ]
# [1]  9.586457  7.194073  9.718831  9.121589 18.548415  6.288808  7.219317  8.229141  7.301282
# [10] 49.509147  7.993822  6.210289  7.550884 11.020195  7.933567  8.152397  8.178039  6.330099
# [19]  6.191472  6.535123  8.541336 34.789726  8.384959 11.455657 11.638807 11.501653  8.177425
# [28]  9.262253  9.736387

# i = 10, p = 60, 50
# i = 22, p = 66, 35
# i = 5,  p = 27, 19
# i = 24, p = 9, 11.4
# i = 1, p = 7, 9.6
# i = 19, p = 7, 6.2

plot(c(7, 9, 27, 60, 66), c(9.6, 11.4, 19, 50, 35), type = "l")

