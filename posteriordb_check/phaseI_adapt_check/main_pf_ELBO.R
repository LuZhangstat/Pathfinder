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
library(dbscan)    #algorithm to cluster log-prob mass
source("../utils/sim_LBFGS_diag_up_ELBO.R")
source("../utils/lp_utils.R")


pd <- pdb_local() # Posterior database connection
pn <- posterior_names(pd)
L_pn = length(pn)

# parameters settings #
alpha = 0.01
mc.cores = parallel::detectCores() - 2
MC = mc.cores    # 10 iterations
MC = 20
init_bound = 2
width = 800; height = 700 # the size of the plot
seed_list = 1:MC 
factr_tol = 1e9

# parameters settings 20 #
# all models
N1 = 200    # Maximum iters in optimization
N_mode_max = 3
N_sam = 1
N_mass = 5
factr_tol = 1e2
MC = 20

load(file = "../results/lp_posteriordb_LBFGS.RData")
set.seed(123)
#model_record <- c(16, 17, 18, 19)

# preallocate results #
# lp_INV <- array(data = NA, dim = c(2, length(model_record)))
# lp_mean <- c()
lp_opath <- c()

which(model_record == 27)
t_0 <- proc.time()
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
  
  ###  run Bob's Phase I  ###
  data <- get_data(po)
  
  # set up lmm #
  posterior <- to_posterior(model, data)
  D <- get_num_upars(posterior)
  lmm = min(max(D, 5), N1)
  cat("No. pars:", D," lmm in L-BFGS: ", lmm, "\n")
  
  # opt_path_stan(model, data, N1 = 30, N_rep = 6, init_bound = 2)
  # t <- proc.time()
  # opath <- opt_path_stan_parallel(seed_list, mc.cores, model, data,
  #                                 N1, N_mode_max, N_sam, N_mass,
  #                                 init_bound = 2.0, factr_tol, lmm)
  # proc.time() - t
  
  t <- proc.time()
  set.seed(123)
  opath <- opt_path_stan_init_parallel(
    initial_ls[[i]], mc.cores, model, data, N1, N_mode_max, N_sam, N_mass,
    init_bound, factr_tol, lmm)
  print(proc.time() - t)
  
  # Check the scaled log prob mass and keep the index with relatively high prob mass
  pick_mode <- filter_mode(opath)
  #pick_mode <- 1:MC
  
  # get an estimate of E(lp)
  pick_samples <- matrix(unlist(lapply(opath[pick_mode], extract_samples)), 
                         nrow = D)
  lps_samples <- unlist(lapply(opath[pick_mode], extract_lps))
  
  # plot the lp of selected plots vs lp_INV 
  lp_opath[[i]] <- list(opath = opath, pick_mode = pick_mode, 
                        pick_samples = pick_samples, lps_samples = lps_samples)
  plot(1:length(lps_samples), lps_samples, 
       ylim = c(min(lp_INV[1, i] - 0.5 * (lp_INV[2, i] - lp_INV[1, i]), 
                    min(lps_samples)), 
                max(lp_INV[2, i] + 0.5 * (lp_INV[2, i] - lp_INV[1, i]), 
                    max(lps_samples))),
       main = paste0(i, "th model"))
  abline(h = lp_INV[, i])
  abline(h = lp_mean[i], col = "red")
  #readline()
}
proc.time() - t_0
# record problematic models 41 94

# save(file = "../results/lp_posteriordb_phI_adapt_set20.RData",
#      list = c("lp_opath", "lp_INV", "model_record"))


# load("../results/lp_posteriordb_phI_adapt_set20.RData")

## check the plots ##

lps_samples_total <- 
  unlist(lapply(1:length(model_record), #1:16, #34:49, #c(17:31,33), #length(model_record), c(1:2, 4:7, 9:14, 16) c(17:20, 22: 31, 33) mode: c(3, 8, 15, 21, 32)
                f <- function(ind){
    lp_opath[[ind]]$lps_samples - lp_mean[ind]}))
index_total <- unlist(lapply(1:length(model_record), #1:16, #34:49, #c(17:31,33), #length(model_record), 
                             f <- function(ind){
                               rep(ind, length(lp_opath[[ind]]$lps_samples))}))
INV_L <-  unlist(lapply(1:length(model_record), #1:16, #34:49, #c(17:31,33), #length(model_record), 
                        f <- function(ind){
                          rep(lp_INV[1, ind] - lp_mean[ind], 
                              length(lp_opath[[ind]]$lps_samples))}))
INV_U <-  unlist(lapply(1:length(model_record), #1:16, #34:49, #c(17:31,33), #length(model_record), 
                        f <- function(ind){
                          rep(lp_INV[2, ind]- lp_mean[ind], 
                              length(lp_opath[[ind]]$lps_samples))}))

pic_array <- list()
pic_array[[1]] <- 1:16;
pic_array[[2]] <- c(17:31,33);
pic_array[[3]] <- 34:49;
for(pic_ind in 1:3){
  cat(pic_ind, "\t")
  illustrate_models = 
    which(sapply(index_total, f <- function(x){any(x == pic_array[[pic_ind]])}))
  df_lps <- data.frame(lps = lps_samples_total[illustrate_models],
                       model = pn[model_record[index_total[illustrate_models]]],
                       ind = index_total[illustrate_models],
                       INV_L = INV_L[illustrate_models],
                       INV_U = INV_U[illustrate_models])
  df_lps$model <- as.factor(df_lps$model)
  
  jpeg(filename = paste0("../pics/phI_adapt_setting20/lps_vs_INV_", 
                         pic_ind, ".jpeg"),
       width = 860, height = 740, units = "px", pointsize = 12)
  p_box_iter <- df_lps %>%
    ggplot(aes(y = model, x = lps)) + 
    geom_point(size = 0.05) + 
    geom_jitter(position=position_jitter(height = 0.1)) +
    geom_errorbar(aes(xmin=INV_L, xmax=INV_U), width = 0.5 ,
                  position = position_dodge(width=0)) +
    ylab("") + 
    xlab("log-density vs 99% CI") + 
    theme_bw(base_size = 26) 
  print(p_box_iter)
  dev.off()
  
}

                          
# 1:16, #34:49, #c(17:31,33), #length(model_record), c(1:2, 4:7, 9:14, 16) 
# c(17:20, 22: 31, 33) mode: c(3, 8, 15, 21, 32)

df_41 <- data.frame(lps = lp_opath[[32]]$lps_samples - lp_mean[32],
                    model = pn[model_record[
                      rep(32, length(lp_opath[[32]]$lps_samples))]],
                    ind = rep(32,length(lp_opath[[32]]$lps_samples)),
                    INV_L = rep(lp_INV[1, 32] - lp_mean[32], 
                                length(lp_opath[[32]]$lps_samples)),
                    INV_U = rep(lp_INV[2, 32] - lp_mean[32], 
                                length(lp_opath[[32]]$lps_samples)))

# plot for model mcycle_gp-accel_gp 
jpeg(filename = paste0("../pics/phI_adapt_setting20/lps_vs_INV_41.jpeg"),
     width = 860, height = 150, #740, 
     units = "px", pointsize = 12)
p_box_iter <- df_41 %>%
  ggplot(aes(y = model, x = lps)) + 
  geom_point(size = 0.05) + 
  geom_jitter(position=position_jitter(height = 0.1)) +
  geom_errorbar(aes(xmin=INV_L, xmax=INV_U), width = 0.5 ,
                position = position_dodge(width=0)) +
  ylab("") + 
  xlab("log-density vs 99% CI") + 
  theme_bw(base_size = 26) 
print(p_box_iter)
dev.off()


## compare the inits, MAPs and picked draws ##
get_init_MAP_lp <- function(ind){
  lp_ind = ncol(lp_opath[[ind]]$opath[[1]]$y)
  inits_lps <- c()
  MAP_lps <- c()
  for(l in 1:length(lp_opath[[ind]]$opath)){
    inits_lps = c(inits_lps, lp_opath[[ind]]$opath[[l]]$y[1, lp_ind] - 
                    lp_mean[ind])
    last_ind <- nrow(lp_opath[[ind]]$opath[[l]]$y)
    MAP_lps = c(MAP_lps, lp_opath[[ind]]$opath[[l]]$y[last_ind, lp_ind] - 
                  lp_mean[ind])
  }
  return(cbind(inits_lps, MAP_lps))
  #return(MAP_lps)
}

#c(3, 8, 15, 21, 32)  multi-model index

lps_inits_MAP_total <- 
  unlist(lapply(1:length(model_record),  #1:16, #34:49, #c(17:31,33), #length(model_record), #c(1:2, 4:7, 9:14, 16) 
                get_init_MAP_lp))

index_inits_MAP <- rep(1:length(model_record), #34:49,  #1:16, c(17:31,33) c(1:2, 4:7, 9:14, 16)
                       each = 20*2)
type_index <- c(rep("picked", length(lps_samples_total)), 
                rep(rep(c("inits", "MAP"), each = 20), length(model_record))) 

INV_L2 <-  unlist(lapply(1:length(model_record), #34:49, #1:16, #34:49, #c(17:31,33), #length(model_record), 
                        f <- function(ind){
                          rep(lp_INV[1, ind] - lp_mean[ind], 20*2)}))
INV_U2 <-  unlist(lapply(1:length(model_record), #34:49, #34:49, #c(17:31,33), #length(model_record), 
                        f <- function(ind){
                          rep(lp_INV[2, ind]- lp_mean[ind], 20*2)}))

### compare MAPs, picked draws and INV ###
pic_comp_array <- list()
pic_comp_array[[1]] <- c(1:2, 4:7, 9:14, 16:17);
pic_comp_array[[2]] <- c(18:20, 22:31, 33:34);
pic_comp_array[[3]] <- 35:49;
pic_comp_array[[4]] <- c(3, 8, 15, 21, 32);  # multi-model problems
illustrate_type = 
  sapply(type_index, f <- function(x){any(x == c("picked", "MAP"))}
)
for(pic_ind in 1:4){
  cat(pic_ind, "\t")
  illustrate_models = 
    sapply(c(index_total, index_inits_MAP), 
           f <- function(x){any(x == pic_comp_array[[pic_ind]])})
  illustrate_ind = which(illustrate_models & illustrate_type)
  df_lps2 <- data.frame(lps = c(lps_samples_total, lps_inits_MAP_total)[illustrate_ind],
                        model = pn[model_record[
                          c(index_total, index_inits_MAP)[illustrate_ind]]],
                        ind = c(index_total, index_inits_MAP)[illustrate_ind],
                        type_index = type_index[illustrate_ind],
                        INV_L = c(INV_L, INV_L2)[illustrate_ind],
                        INV_U = c(INV_U, INV_U2)[illustrate_ind])
  df_lps$model <- as.factor(df_lps$model)
  
  jpeg(filename = paste0("../pics/phI_adapt_setting20/MAP_pick_INV_",
                         pic_ind, ".jpeg"),
       width = 860, height = 740, units = "px", pointsize = 12)
  p_box_iter <- df_lps2 %>%
    ggplot(aes(y = model, x = lps, shape = type_index, colour = type_index)) + 
    geom_point(size = 0.05) + 
    geom_jitter(position=position_jitter(height = 0.1)) +
    geom_errorbar(aes(xmin=INV_L, xmax=INV_U), width = 0.5 ,
                  position = position_dodge(width=0), colour = "black") +
    ylab("") + 
    xlab("log-density vs 99% CI")  + #xlim(-5, 10) +
    theme_bw(base_size = 26) + 
    theme(legend.position = "bottom") +
    theme(legend.title = element_blank()) 
  print(p_box_iter)
  dev.off()
}

### compare inits, picked draws and INV ##
lps = c(lps_samples_total, lps_inits_MAP_total)
models = model_record[c(index_total, index_inits_MAP)]
INV_L3 = c(INV_L, INV_L2)
INV_U3 = c(INV_U, INV_U2)
for(model_ind in model_record){
  data_pick = which(models == model_ind)
  min_lps = min(lps[data_pick], INV_L3[data_pick])
  lps[data_pick] = lps[data_pick] - min_lps + 1.0
  INV_L3[data_pick] = (INV_L3[data_pick] - min_lps + 1.0)
  INV_U3[data_pick] = (INV_U3[data_pick] - min_lps + 1.0)
}

pic_comp_array <- list()
pic_comp_array[[1]] <- c(15:16, 47); # ~50
pic_comp_array[[2]] <- c(1, 6:7, 17:21); #~7500
pic_comp_array[[3]] <- c(9:14);  #2e+07
pic_comp_array[[4]] <- c(22:29); #??e+07
pic_comp_array[[5]] <- c(4, 5, 32); #6e+06
pic_comp_array[[6]] <- c(8); #3e+12
pic_comp_array[[7]] <- c(2, 3, 31, 33:37); #2e+05
pic_comp_array[[8]] <- c(39:46); # 3e+05
pic_comp_array[[9]] <- c(30, 38, 48, 49); #8.2e+09
length(unique(unlist(pic_comp_array)))
#pic_comp_array[[4]] <- c(3, 8, 15, 21, 32);  # multi-model problems
illustrate_type = 
  sapply(type_index, f <- function(x){any(x == c("picked", "inits"))}
  )
for(pic_ind in 1:9){
  cat(pic_ind, "\t")
  illustrate_models = 
    sapply(c(index_total, index_inits_MAP), 
           f <- function(x){any(x == pic_comp_array[[pic_ind]])})
  illustrate_ind = which(illustrate_models & illustrate_type)
  df_lps2 <- data.frame(lps = lps[illustrate_ind],
                        model = pn[model_record[
                          c(index_total, index_inits_MAP)[illustrate_ind]]],
                        ind = c(index_total, index_inits_MAP)[illustrate_ind],
                        type_index = type_index[illustrate_ind],
                        INV_L = INV_L3[illustrate_ind],
                        INV_U = INV_U3[illustrate_ind])
  df_lps$model <- as.factor(df_lps$model)
  
  jpeg(filename = paste0("../pics/phI_adapt_setting20/inits_pick_INV_",
                         pic_ind, ".jpeg"),
       width = 860, height = 740, units = "px", pointsize = 12)
  p_box_iter <- df_lps2 %>%
    ggplot(aes(y = model, x = lps, shape = type_index, colour = type_index)) + 
    geom_point(size = 0.05) + 
    geom_jitter(position=position_jitter(height = 0.1)) +
    geom_errorbar(aes(xmin=INV_L, xmax=INV_U), width = 0.5 ,
                  position = position_dodge(width=0), colour = "black") +
    ylab("") + 
    xlab("log-density vs 99% CI")  + #xlim(-5, 10) +
    scale_x_continuous(trans = "sqrt") +
    theme_bw(base_size = 26) + 
    theme(legend.position = "bottom") +
    theme(legend.title = element_blank()) 
  print(p_box_iter)
  dev.off()
}


## check the efficiency of pathfinder ##
# parameters settings #
alpha = 0.01
L = 1000
M = 20
width = 860; height = 740 # the size of the plot

#takeoff <- c(24)   # model 24 has transformed pars
load("../results/lp_posteriordb_explore4.RData")

pathfinder_fn_call <- 
  sapply(lp_opath, f <- function(x){ 
    sapply(x$opath,  g <- function(z){sum(z$fn_call)}) } )
pathfinder_gr_call <- 
  sapply(lp_opath, f <- function(x){
    sapply(x$opath, g <- function(z){sum(z$gr_call)}) } )

# boxplot of counts
df <- data.frame(fn_call = c(pathfinder_fn_call),
                 gr_call = c(pathfinder_gr_call),
                 n_leapfrogs = c(lp_explore_n_leapfrog),
                 model = rep(pn[model_record], each = M))

jpeg(filename = paste0("../pics/phI_adapt_setting20/box_pathfinder_fn_counts_log.jpeg"),
     width = 860*1.3, height = 740*2, units = "px", pointsize = 12)
p_box_iter <- df %>%
  ggplot(aes(y = reorder(model, n_leapfrogs, FUN = median), 
             x = fn_call)) + 
  geom_boxplot() + #scale_fill_manual(values=c("red", "white")) +
  scale_x_log10() + ylab("") + xlab("No. calls to lp__") +
  scale_alpha_manual(values=c(1,0.1)) +
  theme_bw(base_size = 26 )  +  
  theme(legend.position = "none") 
print(p_box_iter)
dev.off()

jpeg(filename = paste0("../pics/phI_adapt_setting20/box_pathfinder_gr_counts_log.jpeg"),
     width = 860*1.3, height = 740*2, units = "px", pointsize = 12)
p_box_iter <- df %>%
  ggplot(aes(y = reorder(model, n_leapfrogs, FUN = median), 
             x = gr_call)) + 
  geom_boxplot() + #scale_fill_manual(values=c("red", "white")) +
  scale_x_log10() + ylab("") + xlab("No. calls to gradients") +
  scale_alpha_manual(values=c(1,0.1)) + 
  theme_bw(base_size = 26)  +  
  theme(legend.position = "none") 
print(p_box_iter)
dev.off()


## Stan Phase I vs pathfinder ##
df <- data.frame(n_counts = c(c(abs(pathfinder_fn_call)), c(lp_explore_n_leapfrog)),
                 model = rep(rep(pn[model_record], each = M), 2),
                 n_leapfrogs = rep(c(lp_explore_n_leapfrog), 2),
                 multimodel = c(rep(1, M*length(model_record)), 
                                rep(2, M*length(model_record))))

jpeg(filename = paste0("../pics/phI_adapt_setting20/box_compar_fn_log.jpeg"),
     width = 860*2, height = 740*2, units = "px", pointsize = 12)
p_box_compar <- df %>% 
  mutate(type= c("pathfinder ELBO", 
                 "Stan Phase I        ")[multimodel]) %>%
  ggplot(aes(y = reorder(model, n_leapfrogs, FUN = median), 
             x = n_counts, fill = type, color = type)) + 
  geom_boxplot() + 
  scale_fill_manual(values=c("lightseagreen", "lightskyblue")) + 
  scale_color_manual(values=c("darkgreen", "blue")) + 
  scale_x_log10(breaks=c(10, 1e2, 1e3, 1e4, 1e5), 
                labels = c("10", "100", "1000", "10,000", "100,000")) + 
  ylab("") + xlab("calls to lp__") + 
  theme_bw(base_size = 26)  #+
#theme(legend.position = "none") 
print(p_box_compar)
dev.off()


df <- data.frame(n_counts = c(c(abs(pathfinder_gr_call)), c(lp_explore_n_leapfrog)),
                 model = rep(rep(pn[model_record], each = M), 2),
                 n_leapfrogs = rep(c(lp_explore_n_leapfrog), 2),
                 multimodel = c(rep(1, M*length(model_record)), 
                                rep(2, M*length(model_record))))

jpeg(filename = paste0("../pics/phI_adapt_setting20/box_compar_gr_log.jpeg"),
     width = 860*2, height = 740*2, units = "px", pointsize = 12)
p_box_compar <- df %>% 
  mutate(type= c("pathfinder ELBO", 
                 "Stan Phase I        ")[multimodel]) %>%
  ggplot(aes(y = reorder(model, n_leapfrogs, FUN = median), 
             x = n_counts, fill = type, color = type)) + 
  geom_boxplot() + 
  scale_fill_manual(values=c("lightseagreen", "lightskyblue")) + 
  scale_color_manual(values=c("darkgreen", "blue")) + 
  scale_x_log10(breaks=c(10, 1e2, 1e3, 1e4, 1e5), 
                labels = c("10", "100", "1000", "10,000", "100,000")) + 
  ylab("") + xlab("calls to gradient") + 
  theme_bw(base_size = 26)  #+
#theme(legend.position = "none") 
print(p_box_compar)
dev.off()

## ELBO pathfinder vs pathfinder ##
load("../results/lp_posteriordb_phI_adapt_set17_2.RData")
pathfinder_fn_call_pf <- 
  sapply(lp_opath, f <- function(x){ 
    sapply(x$opath,  g <- function(z){sum(z$fn_call)}) } )
pathfinder_gr_call_pf <- 
  sapply(lp_opath, f <- function(x){ 
    sapply(x$opath,  g <- function(z){sum(z$gr_call)}) } )

df <- data.frame(n_counts = c(abs(pathfinder_fn_call), 
                              abs(pathfinder_fn_call_pf)),
                 model = rep(rep(pn[model_record], each = M), 2),
                 n_leapfrogs = rep(abs(pathfinder_fn_call_pf), 2),
                 multimodel = c(rep(1, M*length(model_record)), 
                                rep(2, M*length(model_record))))

jpeg(filename = paste0("../pics/phI_adapt_setting20/box_compar_fn_log_pf.jpeg"),
     width = 860*2, height = 740*2, units = "px", pointsize = 12)
p_box_compar <- df %>% 
  mutate(type= c("pathfinder ELBO", 
                 "pathfinder        ")[multimodel]) %>%
  ggplot(aes(y = reorder(model, n_leapfrogs, FUN = median), 
             x = n_counts, fill = type, color = type)) + 
  geom_boxplot() + 
  scale_fill_manual(values=c("lightseagreen", "lightskyblue")) + 
  scale_color_manual(values=c("darkgreen", "blue")) + 
  scale_x_log10(breaks=c(10, 1e2, 1e3, 1e4), 
                labels = c("10", "100", "1000", "10,000")) + 
  ylab("") + xlab("calls to lp__") + 
  theme_bw(base_size = 26)  #+
#theme(legend.position = "none") 
print(p_box_compar)
dev.off()


df <- data.frame(n_counts = c(c(abs(pathfinder_gr_call)), 
                              c(abs(pathfinder_gr_call_pf))),
                 model = rep(rep(pn[model_record], each = M), 2),
                 n_leapfrogs = rep(c(lp_explore_n_leapfrog), 2),
                 multimodel = c(rep(1, M*length(model_record)), 
                                rep(2, M*length(model_record))))

jpeg(filename = paste0("../pics/phI_adapt_setting20/box_compar_gr_log_pf.jpeg"),
     width = 860*2, height = 740*2, units = "px", pointsize = 12)
p_box_compar <- df %>% 
  mutate(type= c("pathfinder ELBO", 
                 "Stan Phase I        ")[multimodel]) %>%
  ggplot(aes(y = reorder(model, n_leapfrogs, FUN = median), 
             x = n_counts, fill = type, color = type)) + 
  geom_boxplot() + 
  scale_fill_manual(values=c("lightseagreen", "lightskyblue")) + 
  scale_color_manual(values=c("darkgreen", "blue")) + 
  scale_x_log10(breaks=c(10, 1e2, 1e3, 1e4, 1e5), 
                labels = c("10", "100", "1000", "10,000", "100,000")) + 
  ylab("") + xlab("calls to gradient") + 
  theme_bw(base_size = 26)  #+
#theme(legend.position = "none") 
print(p_box_compar)
dev.off()



# check covid model  No. 17
# check ovarian-logistic_regression_rhs No. 62
which(pn == "ovarian-logistic_regression_rhs")
  modelname <- pn[62]
  
  df_ovarian <- data.frame(lps = lps_samples,
                      model = rep(modelname, length(lps_samples)),
                      INV_L = rep(-3100, length(lps_samples)),
                      INV_U = rep(-2850, length(lps_samples)))
  # plot for model mcycle_gp-accel_gp 
  jpeg(filename = paste0("../pics/phI_adapt_setting20/lps_vs_INV_ovarian.jpeg"),
       width = 860, height = 150, #740, 
       units = "px", pointsize = 12)
  p_box_iter <- df_ovarian %>%
    ggplot(aes(y = model, x = lps)) + 
    geom_point(size = 0.05) + 
    geom_jitter(position=position_jitter(height = 0.1)) +
    geom_errorbar(aes(xmin=INV_L, xmax=INV_U), width = 0.5 ,
                  position = position_dodge(width=0)) +
    ylab("") + 
    xlab("log-density vs 99% CI") + 
    theme_bw(base_size = 26) 
  print(p_box_iter)
  dev.off()
  printf("model %s", modelname)
 
  plot(1:length(lps_samples), lps_samples, 
       ylim = c(min(-3200, 
                    min(lps_samples)), 
                max(-2700, 
                    max(lps_samples))),
       main = modelname)
  abline(h = c(-2850, -3100))
  
  # MAP, inits, picked INV
  inits_lps_ovarian <- c()
  MAP_lps_ovarian <- c()
  lp_ind_ovarian = ncol(opath[[1]]$y)
  for(l in 1:length(opath)){
    inits_lps_ovarian = c(inits_lps_ovarian, opath[[l]]$y[1, lp_ind_ovarian])
    last_ind <- nrow(opath[[l]]$y)
    MAP_lps_ovarian = c(MAP_lps_ovarian, opath[[l]]$y[last_ind, lp_ind_ovarian])
  }
  
  df_ovarian <- data.frame(lps = c(lps_samples, inits_lps_ovarian, 
                                   MAP_lps_ovarian),
                           model = rep(modelname, length(lps_samples) + 2 * MC),
                           INV_L = rep(-3100, length(lps_samples) + 2 * MC),
                           INV_U = rep(-2850, length(lps_samples) + 2 * MC),
                           type_index = c(rep("picked", length(lps_samples)),
                                          rep("inits", 20), rep("MAP", 20)))
  df_ovarian$model <- as.factor(df_ovarian$model)
  
  # plot for model mcycle_gp-accel_gp 
  jpeg(filename = paste0("../pics/phI_adapt_setting20/inits_MAP_pick_INV_ovarian.jpeg"),
       width = 860, height = 220, #740, 
       units = "px", pointsize = 12)
  p_box_iter <- df_ovarian %>%
    ggplot(aes(y = model, x = lps, shape = type_index, colour = type_index)) + 
    geom_point(size = 0.05) + 
    geom_jitter(position=position_jitter(height = 0.1)) +
    geom_errorbar(aes(xmin=INV_L, xmax=INV_U), width = 0.5 ,
                  position = position_dodge(width=0), colour = "black") +
    ylab("") + 
    xlab("log-density vs 99% CI") + 
    theme_bw(base_size = 26) +
    theme(legend.position = "bottom") +
    theme(legend.title = element_blank()) 
  print(p_box_iter)
  dev.off()
 
  
  