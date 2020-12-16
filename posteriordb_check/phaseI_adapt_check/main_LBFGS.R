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
library(dplyr)
#source("../utils/sim_LBFGS.R")
source("../utils/sim_LBFGS_diag_up.R")
#source("../utils/sim_CG_diag_up.R")
source("../utils/lp_utils.R")
#source("../utils/L_RHR.R")

set.seed(123)
pd <- pdb_local() # Posterior database connection
pn <- posterior_names(pd)
L_pn = length(pn)

# parameters settings #
alpha = 0.01
N1 = 300    # Maximum iters in optimization
N_mode_max = 20 
N_sam = 200
mc.cores = parallel::detectCores() - 2
MC = mc.cores    # 10 iterations
MC = 5
init_bound = 2
width = 800; height = 700 # the size of the plot
seed_list = 1:MC 
factr_tol = 1e9

# parameters settings 14 #
N1 = 300    # Maximum iters in optimization
N_mode_max = 20 
N_sam = 200
factr_tol = 1e9

# parameters settings 15 #
N1 = 1000    # Maximum iters in optimization
N_mode_max = 20 
N_sam = 200
factr_tol = 1e2

# parameters settings 16 #
# all models
#source("../utils/sim_LBFGS_diag_up.R")
N1 = 300    # Maximum iters in optimization
N_mode_max = 20 
N_sam = 100
factr_tol = 1e2
lmm = 5
MC = 20
seed_list = 1:20


#lmm = min(2*as.integer(log(D)) + 5, D)

# parameters settings 17 #
# all models
# N1 = 300    # Maximum iters in optimization
# N_mode_max = 20 
# N_sam = 100
# factr_tol = 1e2
# lmm = 5
# MC = 20

# parameters setting 19 #
# source("../utils/sim_CG_diag_up.R")

# model_record <- c(1, 2, 3, 4, 6, 8, 9, 11, 13, 15, 20, 21, 23, 25, 26, 27, 29, 31, 33,
#                   34, 36, 37, 40, 41, 43, 47, 51, 55, 61, 94, 95)
# N_models <- length(model_record)

load(file = "../results/lp_posteriordb_LBFGS.RData")



#model_record <- c(16, 17, 18, 19)

# preallocate results #
# lp_INV <- array(data = NA, dim = c(2, length(model_record)))
# lp_mean <- c()
lp_opath <- c()

which(model_record == 41)

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
  # if(modelname == "eight_schools-eight_schools_noncentered"){
  #   INV <- lp_Int_q_posteriordb_8school_noncen(po, alpha)
  # } else if (modelname == "gp_pois_regr-gp_pois_regr") {
  #   INV <- lp_Int_q_posteriordb_gp_pois_regr(po, alpha)
  # } else {INV <- lp_Int_q_posteriordb(po, alpha)}
  # lp_INV[, i] = INV[1:2]
  # lp_mean[i] = INV[3]
  
  ###  run Bob's Phase I  ###
  data <- get_data(po)
  
  # set up lmm #
  posterior <- to_posterior(model, data)
  D <- get_num_upars(posterior)
  #lmm = min(as.integer(log(D))*2 + 5, D)
  #lmm = max(lmm, 5)
  lmm = 5
  cat("No. pars:", D," lmm in L-BFGS: ", lmm, "\n")
  
  #opt_path_stan(model, data, N1 = 30, N_rep = 6, init_bound = 2)
  # opath <- opt_path_stan_parallel(seed_list, mc.cores, model, data,
  #                                 N1, N_mode_max, N_sam,
  #                                 init_bound = 2.0, factr_tol, lmm)

  t <- proc.time()
  set.seed(123)
  opath <- opt_path_stan_init_parallel(
    initial_ls[[i]], mc.cores, model, data, N1, N_mode_max, N_sam,
    init_bound, factr_tol, lmm)
  proc.time() - t
  
  pick_ind <- mclapply(opath, find_indx, mc.cores = mc.cores)
  # lp_opath <- list(opath = opath, pick_ind = pick_ind)
  lp_opath[[i]] <- list(opath = opath, pick_ind = pick_ind)
  
}


param_path <- opath[[1]]
test_ind <- find_typical(opath[[1]], model, data)
opath[[9]][test_ind, ncol(opath[[9]])]
init_param_unc <- opath[[1]][22, -ncol(opath[[1]])]
opath[[3]][39, ncol(opath[[1]])]

# save(file = "../results/lp_posteriordb_phI_adapt_set17.RData",
#      list = c("lp_opath", "lp_INV", "model_record"))


# load("../results/lp_posteriordb_phI_adapt_set17.RData")
which(model_record ==94) # 9 21 37 61
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
                                 if(is.null(x$E_lp[which.max(x$log_MASS_c)])){"NULL"}else{
                                   round(x$E_lp[which.max(x$log_MASS_c)], digits = 1)
                                 }}), 
                        collapse = " "), 
                  "VS E(lp):", round(lp_mean[i], digits = 1),
                  # "\n", "counts in line search:", 
                  # paste(sapply(lp_opath[[i]]$opath[chain_id], 
                  #              f <- function(x){ 
                  #                if(is.null(x$E_lp[which.max(x$log_MASS_c)])){"NULL"}else{
                  #                  sum(x$step_count)
                  #                }}), 
                  #       collapse = " "), "No. of sample",
                  # paste(sapply(lp_opath[[i]]$opath[chain_id], 
                  #              f <- function(x){ 
                  #                if(is.null(x$E_lp[which.max(x$log_MASS_c)])){"0"}else{
                  #                  N_sam * length(x$E_lp)
                  #                }}), 
                  #       collapse = " "),
                  "\n", "estimated log scaled Prob Mass:", 
                  paste(sapply(lp_opath[[i]]$opath[chain_id], 
                               f <- function(x){ 
                                 if(is.null(x$E_lp[which.max(x$log_MASS_c)])){"NULL"}else{
                                   round(x$log_MASS_c[which.max(x$log_MASS_c)], 
                                         digits = 1)}}), 
                        collapse = " "),
                  # "\n", "-0.5logdet(Hk) around picked mode:",
                  # paste(sapply(lp_opath[[i]]$opath[chain_id], 
                  #              f <- function(x){ 
                  #                if(is.null(x$E_lp[which.max(x$log_MASS_c)])){"NULL"}else{
                  #                  round(x$lVol[which.max(x$log_MASS_c)], digits = 1)}}), 
                  #       collapse = " "),
                  "\n", "Pareto K:",
                  paste(sapply(lp_opath[[i]]$opath[chain_id], 
                               f <- function(x){ 
                                 if(is.null(x$E_lp[which.max(x$log_MASS_c)])){"NULL"}else{
                                   round(x$pareto_k[which.max(x$log_MASS_c)], digits = 2)}}), 
                        collapse = " ")
                  # ,"\n", "condition number:",
                  # paste(sapply(lp_opath[[i]]$opath[chain_id], 
                  #              f <- function(x){ 
                  #                if(is.null(x$E_lp[which.max(x$log_MASS_c)])){"NULL"}else{
                  #                  round(x$cond_num[which.max(x$log_MASS_c)], digits = 2)}}), 
                  #       collapse = " ")
                  )) +
    theme_bw() +
    theme(title =element_text(size=16, face='bold'))

  
  jpeg(filename = paste0("../pics/phI_adapt_setting17/No", model_record[i], "-", 
                         modelname, "_L.jpeg"), #model_record[i]
       width = width, height = height, units = "px", pointsize = 16)
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
             max(p_lp_trace$lp__)))  + 
    ggtitle(paste("model:", modelname,"\n", "estimated E(lp):", 
                  paste(sapply(lp_opath[[i]]$opath[chain_id], 
                               f <- function(x){ 
                                 if(is.null(x$E_lp[which.max(x$log_MASS_c)])){"NULL"}else{
                                   round(x$E_lp[which.max(x$log_MASS_c)], digits = 1)
                                 }}), 
                        collapse = " "), 
                  "VS E(lp):", round(lp_mean[i], digits = 1),
                  # "\n", "counts in line search:", 
                  # paste(sapply(lp_opath[[i]]$opath[chain_id], 
                  #              f <- function(x){ 
                  #                if(is.null(x$E_lp[which.max(x$log_MASS_c)])){"NULL"}else{
                  #                  sum(x$step_count)
                  #                }}), 
                  #       collapse = " "), "No. of sample",
                  # paste(sapply(lp_opath[[i]]$opath[chain_id], 
                  #              f <- function(x){ 
                  #                if(is.null(x$E_lp[which.max(x$log_MASS_c)])){"0"}else{
                  #                  N_sam * length(x$E_lp)
                  #                }}), 
                  #       collapse = " "),
                  "\n", "estimated log scaled Prob Mass:", 
                  paste(sapply(lp_opath[[i]]$opath[chain_id], 
                               f <- function(x){ 
                                 if(is.null(x$E_lp[which.max(x$log_MASS_c)])){"NULL"}else{
                                   round(x$log_MASS_c[which.max(x$log_MASS_c)], 
                                         digits = 1)}}), 
                        collapse = " "),
                  # "\n", "-0.5logdet(Hk) around picked mode:",
                  # paste(sapply(lp_opath[[i]]$opath[chain_id], 
                  #              f <- function(x){ 
                  #                if(is.null(x$E_lp[which.max(x$log_MASS_c)])){"NULL"}else{
                  #                  round(x$lVol[which.max(x$log_MASS_c)], digits = 1)}}), 
                  #       collapse = " "),
                  "\n", "Pareto K:",
                  paste(sapply(lp_opath[[i]]$opath[chain_id], 
                               f <- function(x){ 
                                 if(is.null(x$E_lp[which.max(x$log_MASS_c)])){"NULL"}else{
                                   round(x$pareto_k[which.max(x$log_MASS_c)], digits = 2)}}), 
                        collapse = " ")
                  # ,"\n", "condition number:",
                  # paste(sapply(lp_opath[[i]]$opath[chain_id], 
                  #              f <- function(x){ 
                  #                if(is.null(x$E_lp[which.max(x$log_MASS_c)])){"NULL"}else{
                  #                  round(x$cond_num[which.max(x$log_MASS_c)], digits = 2)}}), 
                  #       collapse = " ")
    )) + theme_bw() +
    theme(title =element_text(size=16, face='bold'))
  
  jpeg(filename = paste0("../pics/phI_adapt_setting19/No", model_record[i], "-", 
                         modelname, "_s.jpeg"), #model_record[i]
       width = width, height = height, units = "px", pointsize = 16)
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
             max(p_lp_trace$lp__)))  + 
    ggtitle(paste("model:", modelname,"\n", "estimated E(lp):", 
                  paste(sapply(lp_opath[[i]]$opath[chain_id], 
                               f <- function(x){ 
                                 if(is.null(x$E_lp[which.max(x$log_MASS_c)])){"NULL"}else{
                                   round(x$E_lp[which.max(x$log_MASS_c)], digits = 1)
                                 }}), 
                        collapse = " "), 
                  "VS E(lp):", round(lp_mean[i], digits = 1),
                  # "\n", "counts in line search:", 
                  # paste(sapply(lp_opath[[i]]$opath[chain_id], 
                  #              f <- function(x){ 
                  #                if(is.null(x$E_lp[which.max(x$log_MASS_c)])){"NULL"}else{
                  #                  sum(x$step_count)
                  #                }}), 
                  #       collapse = " "), "No. of sample",
                  # paste(sapply(lp_opath[[i]]$opath[chain_id], 
                  #              f <- function(x){ 
                  #                if(is.null(x$E_lp[which.max(x$log_MASS_c)])){"0"}else{
                  #                  N_sam * length(x$E_lp)
                  #                }}), 
                  #       collapse = " "),
                  "\n", "estimated log scaled Prob Mass:", 
                  paste(sapply(lp_opath[[i]]$opath[chain_id], 
                               f <- function(x){ 
                                 if(is.null(x$E_lp[which.max(x$log_MASS_c)])){"NULL"}else{
                                   round(x$log_MASS_c[which.max(x$log_MASS_c)], 
                                         digits = 1)}}), 
                        collapse = " "),
                  # "\n", "-0.5logdet(Hk) around picked mode:",
                  # paste(sapply(lp_opath[[i]]$opath[chain_id], 
                  #              f <- function(x){ 
                  #                if(is.null(x$E_lp[which.max(x$log_MASS_c)])){"NULL"}else{
                  #                  round(x$lVol[which.max(x$log_MASS_c)], digits = 1)}}), 
                  #       collapse = " "),
                  "\n", "Pareto K:",
                  paste(sapply(lp_opath[[i]]$opath[chain_id], 
                               f <- function(x){ 
                                 if(is.null(x$E_lp[which.max(x$log_MASS_c)])){"NULL"}else{
                                   round(x$pareto_k[which.max(x$log_MASS_c)], digits = 2)}}), 
                        collapse = " ")
                  # ,"\n", "condition number:",
                  # paste(sapply(lp_opath[[i]]$opath[chain_id], 
                  #              f <- function(x){ 
                  #                if(is.null(x$E_lp[which.max(x$log_MASS_c)])){"NULL"}else{
                  #                  round(x$cond_num[which.max(x$log_MASS_c)], digits = 2)}}), 
                  #       collapse = " ")
    ))  + theme_bw() +
  theme(title =element_text(size=16, face='bold'))
  
  jpeg(filename = paste0("../pics/phI_adapt_setting19/No", model_record[i], "-", 
                         modelname, "_trun.jpeg"), #model_record[i]
       width = width, height = height, units = "px", pointsize = 16)
  print(p_lp)
  dev.off()
  
  lgnorm_trace = 
    data.frame(iter = c(unlist(sapply(chain_id, f <- function(x){1:Lp_list[x]}))), 
               lgnorm = c(unlist(sapply(lp_opath[[i]]$opath[chain_id], 
                                      f <- function(x){if(is.null(dim(x$y))){
                                        return(x$lgnorms[1])
                                      }else{ return(x$lgnorms)}}))),
               chain = c(unlist(sapply(chain_id, 
                                       f <- function(x){
                                         rep(as.character(x), Lp_list[x])}))))
  
  p_lg <- ggplot(data = lgnorm_trace, 
                 aes(x=iter, y=lgnorm, group=chain, color=chain)) +
    geom_line(colour = "grey") + geom_point(size = 2) + 
    theme_bw() 
  
  jpeg(filename = paste0("../pics/phI_adapt_setting19/No", model_record[i], "-", 
                         modelname, "_zlgnorm.jpeg"), #model_record[i]
       width = width, height = height, units = "px", pointsize = 16)
  print(p_lg)
  dev.off()
  
  # lsk_list = unlist(sapply(lp_opath[[i]]$opath, 
  #                         f <- function(x){
  #                           if(length(x) == 1){
  #                             return(NA)
  #                           }else{return(length(x$sknorm_ls))}}))
  # 
  # sknorm_trace = 
  #   data.frame(iter = c(unlist(sapply(chain_id, 
  #                                     f <- function(x){1:lsk_list[x]}))), 
  #              sknorm = c(unlist(sapply(lp_opath[[i]]$opath[chain_id], 
  #                                       f <- function(x){ return(x$sknorm_ls)}))),
  #              thetak = c(unlist(sapply(lp_opath[[i]]$opath[chain_id], 
  #                                    f <- function(x){ return(x$thetak_ls)}))),
  #              chain = c(unlist(sapply(chain_id, 
  #                                      f <- function(x){
  #                                        rep(as.character(x), lsk_list[x])}))))
  # 
  # p_sk <- ggplot(data = sknorm_trace, 
  #                aes(x=iter, y=log(sknorm), group=chain, color=chain)) +
  #   geom_line(colour = "grey") + geom_point(size = 2) + 
  #   theme_bw() 
  # 
  # jpeg(filename = paste0("../pics/phI_adapt_setting17/No", model_record[i], "-", 
  #                        modelname, "_zlsknorm.jpeg"), #model_record[i]
  #      width = width, height = height, units = "px", pointsize = 16)
  # print(p_sk)
  # dev.off()
  
  # p_thetak <- ggplot(data = sknorm_trace, 
  #                aes(x=iter, y=log(thetak), group=chain, color=chain)) +
  #   geom_line(colour = "grey") + geom_point(size = 2) + 
  #   theme_bw() 
  # 
  # jpeg(filename = paste0("../pics/phI_adapt_setting17/No", model_record[i], "-", 
  #                        modelname, "_zlthetaknorm.jpeg"), #model_record[i]
  #      width = width, height = height, units = "px", pointsize = 16)
  # print(p_thetak)
  # dev.off()
  
}

pathfinder_fn_call <- 
  sapply(lp_opath, f <- function(x){ 
    sapply(x$opath,  g <- function(z){sum(z$fn_call)}) } )
pathfinder_gr_call <- 
  sapply(lp_opath, f <- function(x){
    sapply(x$opath, g <- function(z){sum(z$gr_call)}) } )

# boxplot of counts
df <- data.frame(fn_call = c(pathfinder_fn_call),
                 gr_call = c(pathfinder_gr_call),
                 #n_leapfrogs = c(lp_explore_n_leapfrog),
                 model = rep(pn[model_record], each = 20))

jpeg(filename = paste0("../pics/phI_adapt_setting17/box_pathfinder_fn_counts_log.jpeg"),
     width = 860*1.3, height = 740*2, units = "px", pointsize = 12)
p_box_iter <- df %>%
  ggplot(aes(y = reorder(model, fn_call, FUN = median), 
              x = fn_call)) + 
  geom_boxplot() + 
  scale_x_log10() + ylab("") + xlab("No. calls to fn") + 
  theme_grey(base_size = 26)  
print(p_box_iter)
dev.off()

jpeg(filename = paste0("../pics/phI_adapt_setting17/box_pathfinder_gr_counts_log.jpeg"),
     width = 860*1.3, height = 740*2, units = "px", pointsize = 12)
p_box_iter <- df %>%
  ggplot(aes(y = reorder(model, gr_call, FUN = median), 
             x = gr_call)) + 
  geom_boxplot() + 
  scale_x_log10() + ylab("") + xlab("No. calls to fn") + 
  theme_grey(base_size = 26)  
print(p_box_iter)
dev.off()

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

i = which(model_record == 26)
p_lp <- ggplot(data = p_lp_trace, 
               aes(x=iter, y=lp__, group=chain, color=label)) +
  geom_line(colour = "grey") + geom_point(size = 1) + 
  geom_hline(yintercept = lp_INV[, i], colour = "black") + 
  ylim(lp_INV[1, i] - 10*(lp_INV[2, i] - lp_INV[1, i]),
       max(lp_INV[2, i] + (lp_INV[2, i] - lp_INV[1, i]), 
           max(p_lp_trace$lp__)))
p_lp




# check covid model  No. 17
# check ovarian-logistic_regression_rhs No. 62
which(pn == "ovarian-logistic_regression_rhs")
  modelname <- pn[62]
  
  printf("model %s", modelname)
  #i = model_record[i]
  Lp_list = unlist(sapply(lp_opath$opath, 
                          f <- function(x){
                            if(length(x) == 1){
                              return(NA)
                            }else if(is.vector(x$y)){
                              return(1)
                            }else{return(nrow(x$y))}}))
  chain_id = (1:MC)[!is.na(Lp_list)]
  p_lp_trace = 
    data.frame(iter = c(unlist(sapply(chain_id, f <- function(x){1:Lp_list[x]}))), 
               lp__ = c(unlist(sapply(lp_opath$opath[chain_id], 
                                      f <- function(x){if(is.null(dim(x$y))){
                                        return(-x$y[length(x$y) - 1])
                                      }else{ return(-x$y[, ncol(x$y) - 1])}}))),
               chain = c(unlist(sapply(chain_id, f <- function(x){rep(x, Lp_list[x])}))),
               label = c(unlist(sapply(chain_id, f <- function(x){
                 if(Lp_list[x] == 1){
                   out = c("point on optim path")
                 } else{
                   out = rep("point on optim path", Lp_list[x])
                   out[lp_opath$pick_ind[[x]]] = "marked point"
                   if(length(out) > Lp_list[x]){
                     out = rep("point on optim path", Lp_list[x]) # error in HMC sampling
                   }} 
                 return(out)}))))
  
  p_lp <- ggplot(data = p_lp_trace, 
                 aes(x=iter, y=lp__, group=chain, color=label)) +
    geom_line(colour = "grey") + geom_point(size = 2) + 
    #geom_hline(yintercept = lp_INV[, i], colour = "black") + 
    #geom_hline(yintercept = lp_mean[i], colour = "blue", linetype = 2)+
    ylim(-3300,
         # quantile(p_lp_trace$lp__, 0.05)
         #lp_INV[1, i] - 1.5*(lp_INV[2, i] - lp_INV[1, i])
         # min(lp_INV[1, i] - (lp_INV[2, i] - lp_INV[1, i]),
         #        quantile(p_lp_trace$lp__, 0.2))
         #,
         #max(p_lp_trace$lp__)
         -2500
         ) + 
    ggtitle(paste("model:", modelname,"\n", "estimated E(lp):", 
                  paste(sapply(lp_opath$opath[chain_id], 
                               f <- function(x){ 
                                 if(is.null(x$E_lp[which.max(x$log_MASS_c)])){"NULL"}else{
                                   round(x$E_lp[which.max(x$log_MASS_c)], digits = 1)
                                 }}), 
                        collapse = " "),
                  # "\n", "counts in line search:", 
                  # paste(sapply(lp_opath$opath[chain_id], 
                  #              f <- function(x){ 
                  #                if(is.null(x$E_lp[which.max(x$log_MASS_c)])){"NULL"}else{
                  #                  sum(x$step_count)
                  #                }}), 
                  #       collapse = " "), "No. of sample",
                  # paste(sapply(lp_opath$opath[chain_id], 
                  #              f <- function(x){ 
                  #                if(is.null(x$E_lp[which.max(x$log_MASS_c)])){"0"}else{
                  #                  N_sam * length(x$E_lp)
                  #                }}), 
                  #       collapse = " "),
                  "\n", "estimated log scaled Prob Mass:", 
                  paste(sapply(lp_opath$opath[chain_id], 
                               f <- function(x){ 
                                 if(is.null(x$E_lp[which.max(x$log_MASS_c)])){"NULL"}else{
                                   round(x$log_MASS_c[which.max(x$log_MASS_c)], 
                                         digits = 1)}}), 
                        collapse = " "),
                  # "\n", "-0.5logdet(Hk) around picked mode:",
                  # paste(sapply(lp_opath$opath[chain_id], 
                  #              f <- function(x){ 
                  #                if(is.null(x$E_lp[which.max(x$log_MASS_c)])){"NULL"}else{
                  #                  round(x$lVol[which.max(x$log_MASS_c)], digits = 1)}}), 
                  #       collapse = " "),
                  "\n", "Pareto K:",
                  paste(sapply(lp_opath$opath[chain_id], 
                               f <- function(x){ 
                                 if(is.null(x$E_lp[which.max(x$log_MASS_c)])){"NULL"}else{
                                   round(x$pareto_k[which.max(x$log_MASS_c)], digits = 2)}}), 
                        collapse = " ")
                  # , "\n", "condition number:",
                  # paste(sapply(lp_opath$opath[chain_id], 
                  #              f <- function(x){ 
                  #                if(is.null(x$E_lp[which.max(x$log_MASS_c)])){"NULL"}else{
                  #                  round(x$cond_num[which.max(x$log_MASS_c)], digits = 2)}}), 
                  #       collapse = " ")
                  )) + theme_bw() +
    theme(title =element_text(size=16, face='bold'))
  
  jpeg(filename = paste0("../pics/phI_adapt_setting19/No_",
                         modelname, "_trun2.jpeg"), #model_record[i]
       width = width, height = height, units = "px", pointsize = 12)
  print(p_lp)
  dev.off()
  
  lgnorm_trace = 
    data.frame(iter = c(unlist(sapply(chain_id, f <- function(x){1:Lp_list[x]}))), 
               lgnorm = c(unlist(sapply(lp_opath$opath[chain_id], 
                                        f <- function(x){if(is.null(dim(x$y))){
                                          return(x$lgnorms[1])
                                        }else{ return(x$lgnorms)}}))),
               chain = c(unlist(sapply(chain_id, 
                                       f <- function(x){
                                         rep(as.character(x), Lp_list[x])}))))
  
  p_lg <- ggplot(data = lgnorm_trace, 
                 aes(x=iter, y=lgnorm, group=chain, color=chain)) +
    geom_line(colour = "grey") + geom_point(size = 2) + 
    theme_bw() 
  
  jpeg(filename = paste0("../pics/phI_adapt_setting17/No_",  
                         modelname, "_zlgnorm.jpeg"), #model_record[i]
       width = width, height = height, units = "px", pointsize = 16)
  print(p_lg)
  dev.off()

  posterior <- to_posterior(model, data)
  good_init <- constrain_pars(posterior, opath[[1]]$y[pick_ind[[1]][1], 
                                         1:get_num_upars(posterior)])
  str(good_init)
  good_init$mu
  good_init$phi
  good_init$tau
  
  
  constrain_pars(posterior, opath[[4]]$y[4, 
                                         1:get_num_upars(posterior)])

  
 
  
  
  