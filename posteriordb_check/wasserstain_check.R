setwd("./posteriordb") # set working dir to cloned package
library(transport)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(posteriordb)
library(posterior)
source("../utils/lp_utils.R")
source("../utils/sim_pf.R")


load("../results/lp_posteriordb_phI_adapt_default.RData") # Pathfinder #_resam_all #_resam_each
load("../results/lp_posteriordb_LBFGS_h6.RData")
load("../results/PhI_100_h10.RData")
load("../results/ADVI_100.RData")
load("../results/multi_pf_samples_I20.RData")
#load("../results/multi_pf_samples_default.RData")
#load("../results/pf_SIR_WOR_4_draws_default.RData")


pd <- pdb_local() # Posterior database connection
pn <- posterior_names(pd)

## calculate wasserstain distance one from each##
w_d_matrix = matrix(NA, nrow = 49, ncol = 8)
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
  
  ### get reference samples ###
  posterior <- to_posterior(model, data)
  D <- get_num_upars(posterior)
  if(modelname == "eight_schools-eight_schools_noncentered"){
    constrained_draws <- modify_8school_noncen(po)
    unconstrained_draws <-  lapply(constrained_draws, unconstrain_draws, posterior)
  } else if (modelname == "gp_pois_regr-gp_pois_regr") {
    constrained_draws <- modify_draws_gp_pois_regr(po)
    unconstrained_draws <-  lapply(constrained_draws, unconstrain_draws, posterior)
  } else {
    unconstrained_draws <-  lapply(gsd, unconstrain_draws, posterior)
  }
  
  ref_samples = rbind(unconstrained_draws[[1]], unconstrained_draws[[2]], 
                      unconstrained_draws[[3]], unconstrained_draws[[4]],
                      unconstrained_draws[[5]], unconstrained_draws[[6]],
                      unconstrained_draws[[7]], unconstrained_draws[[8]],
                      unconstrained_draws[[9]], unconstrained_draws[[10]])
  
  ### samples from pathfinder ###
  pick_samples <- lp_opath[[i]]$pick_samples
  
  ### inits and optims ###
  inits_optims <- get_init_optim(i)
  
  ### calculate wasserstein distance ###
  # pathfinder #
  if(ncol(pick_samples) == 1){
    a = wpp(rbind(t(pick_samples), t(pick_samples)),
            mass = rep(1 / 2, 2))
  }else{
    a = wpp(t(pick_samples),
            mass = rep(1 / ncol(pick_samples), ncol(pick_samples)))
  }
  b = wpp(ref_samples, mass = rep(1 / nrow(ref_samples), nrow(ref_samples)))
  w_d_pf <- wasserstein(a, b, p = 2); w_d_pf
  
  # multi-path pathfinder #
  set.seed(123)
  pick_samples_IR <- sapply(lp_multi_opath[[i]], f <- function(x){
    x$pick_samples[, 1]}) #sample.int(ncol(x), 1)
  a_IR = wpp(t(pick_samples_IR),
             mass = rep(1 / ncol(pick_samples_IR), ncol(pick_samples_IR)))
  w_d_pf_IR <- wasserstein(a_IR, b, p = 2); w_d_pf_IR

  # optims #
  inits_optims <- get_init_optim(i)
  a_opt = wpp(inits_optims$optims,
              mass = rep(1 / nrow(inits_optims$optims), nrow(inits_optims$optims)))
  w_d_opt <- wasserstein(a_opt, b, p = 2); w_d_opt
  
  # inits #
  a_init = wpp(inits_optims$inits,
               mass = rep(1 / nrow(inits_optims$inits), nrow(inits_optims$inits)))
  w_d_init <- wasserstein(a_init, b, p = 2); w_d_init
  
  # last samples of phase I #
  a_phI = wpp(PhaseI_last_draw[[i]],
               mass = rep(1 / nrow(PhaseI_last_draw[[i]]),
                          nrow(PhaseI_last_draw[[i]])))
  w_d_phI <- wasserstein(a_phI, b, p = 2); w_d_phI

  # ADVI: meanfield #
  a_ADVI_mf = wpp(t(ADVI_meanfield_draw[[i]]),
                  mass = rep(1 / ncol(ADVI_meanfield_draw[[i]]),
                             ncol(ADVI_meanfield_draw[[i]])))
  w_d_ADVI_mf <- wasserstein(a_ADVI_mf, b, p = 2); w_d_ADVI_mf

  # ADVI: meanfield center#
  a_ADVI_mfc = wpp(t(ADVI_meanfield_center[[i]]),
                  mass = rep(1 / ncol(ADVI_meanfield_center[[i]]),
                             ncol(ADVI_meanfield_center[[i]])))
  w_d_ADVI_mfc <- wasserstein(a_ADVI_mfc, b, p = 2); w_d_ADVI_mfc
  
  # ADVI: fullrank #
  a_ADVI_fr = wpp(t(ADVI_fullrank_draw[[i]]),
                  mass = rep(1 / ncol(ADVI_fullrank_draw[[i]]),
                             ncol(ADVI_fullrank_draw[[i]])))
  w_d_ADVI_fr <- wasserstein(a_ADVI_fr, b, p = 2); w_d_ADVI_fr
  
  
  w_d_matrix[i, ] = c(w_d_pf, w_d_pf_IR, w_d_opt, w_d_init, w_d_phI, w_d_ADVI_mf, 
                      w_d_ADVI_mfc, w_d_ADVI_fr)
  
  cat("pf:", w_d_pf, "\t", "pf IR:", w_d_pf_IR, "\t",
      "optims: ", w_d_opt, "\t", "inits: ", w_d_init, "phI:", w_d_phI, "\n", 
      "ADVI meanfield:", w_d_ADVI_mf, "\t",
      "ADVI meanfield center:", w_d_ADVI_mfc, "\t",
      "ADVI fullrank:", w_d_ADVI_fr, "\n")
  
}
proc.time() - t_0

colnames(w_d_matrix) <- c("pf", "pf_20_IR", "max", "random init", "PhI", 
                          "meanfield", "meanfield center", 
                          "fullrank")
rownames(w_d_matrix) <- pn[model_record]

save(file = "../results/wasserstein_phI_adapt_default.RData",
     list = c("w_d_matrix"))

## check the output ##
load("../results/wasserstein_phI_adapt_default.RData") #default; long_hist, short_L, large_K
colMeans(w_d_matrix)
# check pathfinder vs phase I warmup
pf_vs_phI <- (w_d_matrix[, "pf"] / w_d_matrix[, "PhI"])
summary(pf_vs_phI)
quantile(pf_vs_phI, c(0.05, 0.5, 0.95))
pf_IR_vs_phI <- (w_d_matrix[, "pf_20_IR"] / w_d_matrix[, "PhI"])
summary(pf_IR_vs_phI)
quantile(pf_IR_vs_phI, c(0.05, 0.5, 0.95))

## compare pathfinder with ADVI ##
# ADVI meanfield
summary(w_d_matrix[, "pf"] / w_d_matrix[, "meanfield"])
summary(w_d_matrix[, "pf"] / w_d_matrix[, "fullrank"])
summary(w_d_matrix[, "pf"] / w_d_matrix[, "meanfield center"])

## one plot ##
ratios <- c(w_d_matrix[, c("pf", "pf_20_IR", "PhI", "meanfield",
                           "meanfield center", "fullrank")] / 
              w_d_matrix[, "pf_20_IR"])
ratios[ratios > 2^10] <- 2^10
w_d_dat = data.frame(ratios = ratios,
                     label = c(
                       rep("pathfinder", length(model_record)),
                       rep("multi-path pathfinder", length(model_record)),
                       rep("Stan Phase I", length(model_record)),
                       rep("meanfield ADVI", length(model_record)),
                       rep("meanfield center", length(model_record)),
                       rep("dense ADVI", length(model_record))),
                     model = rep(pn[model_record], 6))#6))

p_w_d_compar <- w_d_dat %>%
  ggplot(aes(y = model, #reorder(model, ratios, FUN = mean),
             x = ratios, color = label, shape = label)) +
  geom_point(size = 2) +
  scale_x_continuous(trans = 'log2',
                     limits = c(1/16, 1024),
                     breaks = c(1/16, 1/8, 1/4, 1/2, 1, 2, 4,
                                8, 16, 32, 64, 128, 256, 512, 1024),
                     labels = c("1/16", "1/8", "1/4",
                                "1/2", "1", "2", "4", "8",
                                "16", "32", "64", "128",
                                "256", "512", ">=1024")) +
  ylab("") + xlab("scaled Wasserstein distance") +
  theme_bw(base_size = 12)  +
  theme(legend.position="top", legend.title = element_blank()) #+
#theme(legend.position = "none")

width <- 12.0
height <- 8.0
pointsize <- 16

setEPS()
postscript("../pics/phI_adapt_default/W_d_compar_100_each.eps",  #_resam_each
           width = width, height = height)
print(p_w_d_compar)
dev.off()

## compute the 100 wasserstein distance for each ADVI and pf ##
M = 100
W_d_100_pf <- array(data = NA, dim = c(M, length(model_record)))
W_d_100_pf_IR <- array(data = NA, dim = c(M, length(model_record)))
W_d_100_ADVI_mf <- array(data = NA, dim = c(M, length(model_record)))
W_d_100_ADVI_fr <- array(data = NA, dim = c(M, length(model_record)))

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
  
  ### get reference samples ###
  posterior <- to_posterior(model, data)
  D <- get_num_upars(posterior)
  if(modelname == "eight_schools-eight_schools_noncentered"){
    constrained_draws <- modify_8school_noncen(po)
    unconstrained_draws <-  lapply(constrained_draws, unconstrain_draws, posterior)
  } else if (modelname == "gp_pois_regr-gp_pois_regr") {
    constrained_draws <- modify_draws_gp_pois_regr(po)
    unconstrained_draws <-  lapply(constrained_draws, unconstrain_draws, posterior)
  } else {
    unconstrained_draws <-  lapply(gsd, unconstrain_draws, posterior)
  }
  
  ref_samples = rbind(unconstrained_draws[[1]], unconstrained_draws[[2]], 
                      unconstrained_draws[[3]], unconstrained_draws[[4]],
                      unconstrained_draws[[5]], unconstrained_draws[[6]],
                      unconstrained_draws[[7]], unconstrained_draws[[8]],
                      unconstrained_draws[[9]], unconstrained_draws[[10]])
  
  
  for(j in 1:M){
    cat(j, "\t")
    ### samples from pathfinder ###
    if(lp_opath[[i]]$opath[[j]]$status == "mode"){ # no approximating Normal
      pick_samples <- matrix(lp_opath[[i]]$opath[[j]]$y[
        nrow(lp_opath[[i]]$opath[[j]]$y), -ncol(lp_opath[[i]]$opath[[j]]$y)], ncol = 1)
    } else {
      pick_samples <- extract_samples(lp_opath[[i]]$opath[[j]])
    }
    if(ncol(pick_samples) == 1){
      a = wpp(rbind(t(pick_samples), t(pick_samples)), 
              mass = rep(1 / 2, 2))
    }else{
      a = wpp(t(pick_samples), 
              mass = rep(1 / ncol(pick_samples), ncol(pick_samples)))
    }
    b = wpp(ref_samples, mass = rep(1 / nrow(ref_samples), nrow(ref_samples)))
    w_d_pf <- wasserstein(a, b, p = 2); w_d_pf
    
    W_d_100_pf[j, i] = w_d_pf
    
    ### samples from multi_path pathfinder ###
    pick_samples_M <- lp_multi_opath[[i]][[j]]$pick_samples
    a_M = wpp(t(pick_samples_M), 
            mass = rep(1 / ncol(pick_samples_M), ncol(pick_samples_M)))
    
    w_d_pf_M <- wasserstein(a_M, b, p = 2); w_d_pf_M
    W_d_100_pf_IR[j, i] = w_d_pf_M
    
    # ADVI: meanfield # (for some cases meanfield ADVI returns extreme outcomes 
    # that makes the computation of Wasserstein-2 distance to be extremely slow. 
    # The resulting Wasserstein-2 distances are all very large. We skip those cases
    # and report Inf. We truncate those large W-2 distances in the summary report)
    if( (i == 8 && j == 55) |
        (i == 8 && j == 64) |
        (i == 8 && j == 66) |
        (i == 8 && j == 94) |
        (i == 10 && j == 16) |
        (i == 13 && j == 12) |
        (i == 38 && j == 15) |
        (i == 38 && j == 68) |
        (i == 39 && j == 10) |
        (i == 40 && j == 62)){
      w_d_ADVI_mf = Inf
      W_d_100_ADVI_mf[j, i] = Inf
    }else{
      a_ADVI_mf = wpp(ADVI_meanfield_draw_100[[i]][[j]],
                      mass = rep(1 / nrow(ADVI_meanfield_draw_100[[i]][[j]]),
                                 nrow(ADVI_meanfield_draw_100[[i]][[j]])))
      w_d_ADVI_mf <- wasserstein(a_ADVI_mf, b, p = 2); w_d_ADVI_mf

      W_d_100_ADVI_mf[j, i] <- w_d_ADVI_mf
    }


    # ADVI: fullrank #
    a_ADVI_fr = wpp(ADVI_fullrank_draw_100[[i]][[j]],
                    mass = rep(1 / nrow(ADVI_fullrank_draw_100[[i]][[j]]),
                               nrow(ADVI_fullrank_draw_100[[i]][[j]])))
    w_d_ADVI_fr <- wasserstein(a_ADVI_fr, b, p = 2); w_d_ADVI_fr

    W_d_100_ADVI_fr[j, i] <- w_d_ADVI_fr
    
    cat("pf:", w_d_pf, "\t", "multi pf:", w_d_pf_M, "\t",
        "ADVI meanfield:", W_d_100_ADVI_mf[j, i], "\t",
        "ADVI fullrank:", W_d_100_ADVI_fr[j, i], "\n")
  }
}
proc.time() - t_0

save(file = "../results/wasserstein_100_default.RData",
     list = c("W_d_100_pf", "W_d_100_pf_IR", "W_d_100_ADVI_mf", 
              "W_d_100_ADVI_fr"))



#large_K
load("../results/wasserstein_100_default.RData")
W_d_100_pf[100, ]

## check the output ##
pf_w_d_qtils <- apply(W_d_100_pf, 2, 
                      f <- function(x){quantile(x, c(0.25, 0.5, 0.75))})
pf_M_w_d_qtils <- apply(W_d_100_pf_IR, 2, 
                      f <- function(x){quantile(x, c(0.25, 0.5, 0.75))})
ADVI_mf_w_d_qtils <- apply(W_d_100_ADVI_mf, 2, 
                           f <- function(x){quantile(x, c(0.25, 0.5, 0.75))})
ADVI_fr_w_d_qtils <- apply(W_d_100_ADVI_fr, 2, 
                           f <- function(x){quantile(x, c(0.25, 0.5, 0.75))})

for(i in 1:49){
  modelname <- pn[model_record[i]]
  printf("model %d: %s", model_record[i], modelname)
  print(cbind(pf_w_d_qtils[, i], pf_M_w_d_qtils[, i], ADVI_mf_w_d_qtils[, i], 
              ADVI_fr_w_d_qtils[, i]))
}

## compute wasserstein distance for sensitivity test ##
load("../results/lp_posteriordb_phI_adapt_short_L.RData")
lp_opath_short_L <- lp_opath
load("../results/lp_posteriordb_phI_adapt_large_K.RData")
lp_opath_large_K <- lp_opath
load("../results/lp_posteriordb_phI_adapt_long_hist.RData")
lp_opath_long_hist <- lp_opath
M = 100
W_d_100_pf_short_L <- array(data = NA, dim = c(M, length(model_record)))
W_d_100_pf_large_K <- array(data = NA, dim = c(M, length(model_record)))
W_d_100_pf_long_hist <- array(data = NA, dim = c(M, length(model_record)))

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
  
  ### get reference samples ###
  posterior <- to_posterior(model, data)
  D <- get_num_upars(posterior)
  if(modelname == "eight_schools-eight_schools_noncentered"){
    constrained_draws <- modify_8school_noncen(po)
    unconstrained_draws <-  lapply(constrained_draws, unconstrain_draws, posterior)
  } else if (modelname == "gp_pois_regr-gp_pois_regr") {
    constrained_draws <- modify_draws_gp_pois_regr(po)
    unconstrained_draws <-  lapply(constrained_draws, unconstrain_draws, posterior)
  } else {
    unconstrained_draws <-  lapply(gsd, unconstrain_draws, posterior)
  }
  
  ref_samples = rbind(unconstrained_draws[[1]], unconstrained_draws[[2]], 
                      unconstrained_draws[[3]], unconstrained_draws[[4]],
                      unconstrained_draws[[5]], unconstrained_draws[[6]],
                      unconstrained_draws[[7]], unconstrained_draws[[8]],
                      unconstrained_draws[[9]], unconstrained_draws[[10]])
  
  b = wpp(ref_samples, mass = rep(1 / nrow(ref_samples), nrow(ref_samples)))
  
  for(j in 1:M){
    cat(j, "\t")
    ### samples from pathfinder ###
    # _short_L
    if(lp_opath_short_L[[i]]$opath[[j]]$status == "mode"){ # no approximating Normal
      pick_samples <- matrix(lp_opath_short_L[[i]]$opath[[j]]$y[
        nrow(lp_opath_short_L[[i]]$opath[[j]]$y), 
        -ncol(lp_opath_short_L[[i]]$opath[[j]]$y)], ncol = 1)
    } else {
      pick_samples <- extract_samples(lp_opath_short_L[[i]]$opath[[j]])
    }
    if(ncol(pick_samples) == 1){
      a = wpp(rbind(t(pick_samples), t(pick_samples)), 
              mass = rep(1 / 2, 2))
    }else{
      a = wpp(t(pick_samples), 
              mass = rep(1 / ncol(pick_samples), ncol(pick_samples)))
    }
    w_d_pf1 <- wasserstein(a, b, p = 2); w_d_pf1
    
    W_d_100_pf_short_L[j, i] = w_d_pf1
    
    # _large_K
    if(lp_opath_large_K[[i]]$opath[[j]]$status == "mode"){ # no approximating Normal
      pick_samples <- matrix(lp_opath_large_K[[i]]$opath[[j]]$y[
        nrow(lp_opath_large_K[[i]]$opath[[j]]$y), 
        -ncol(lp_opath_large_K[[i]]$opath[[j]]$y)], ncol = 1)
    } else {
      pick_samples <- extract_samples(lp_opath_large_K[[i]]$opath[[j]])
    }
    if(ncol(pick_samples) == 1){
      a = wpp(rbind(t(pick_samples), t(pick_samples)), 
              mass = rep(1 / 2, 2))
    }else{
      a = wpp(t(pick_samples), 
              mass = rep(1 / ncol(pick_samples), ncol(pick_samples)))
    }
    w_d_pf2 <- wasserstein(a, b, p = 2); w_d_pf2
    
    W_d_100_pf_large_K[j, i] = w_d_pf2
    
    #_long_hist
    if(lp_opath_long_hist[[i]]$opath[[j]]$status == "mode"){ # no approximating Normal
      pick_samples <- matrix(lp_opath_long_hist[[i]]$opath[[j]]$y[
        nrow(lp_opath_long_hist[[i]]$opath[[j]]$y), 
        -ncol(lp_opath_long_hist[[i]]$opath[[j]]$y)], ncol = 1)
    } else {
      pick_samples <- extract_samples(lp_opath_long_hist[[i]]$opath[[j]])
    }
    if(ncol(pick_samples) == 1){
      a = wpp(rbind(t(pick_samples), t(pick_samples)), 
              mass = rep(1 / 2, 2))
    }else{
      a = wpp(t(pick_samples), 
              mass = rep(1 / ncol(pick_samples), ncol(pick_samples)))
    }
    w_d_pf3 <- wasserstein(a, b, p = 2); w_d_pf3
    
    W_d_100_pf_long_hist[j, i] = w_d_pf3
    
    cat("pf:", W_d_100_pf[j, i], "\t",
        "pf_short_L:", w_d_pf1, "\t",
        "pf_large_K:", w_d_pf2, "\t",
        "pf_long_hist:", w_d_pf3, "\n")
  }
}
proc.time() - t_0

save(file = "../results/wasserstein_100_sen.RData",
     list = c("W_d_100_pf_short_L", "W_d_100_pf_large_K", 
              "W_d_100_pf_long_hist"))


## sensitivity test for I ##
load("../results/multi_pf_samples_I5.RData")
lp_multi_opath_I5 <- lp_multi_opath
load("../results/multi_pf_samples_I40.RData")
lp_multi_opath_I40 <- lp_multi_opath
load("../results/wasserstein_100_default.RData")
M = 100
W_d_100_pf_I5 <- array(data = NA, dim = c(M, length(model_record)))
W_d_100_pf_I40 <- array(data = NA, dim = c(M, length(model_record)))

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
  
  ### get reference samples ###
  posterior <- to_posterior(model, data)
  D <- get_num_upars(posterior)
  if(modelname == "eight_schools-eight_schools_noncentered"){
    constrained_draws <- modify_8school_noncen(po)
    unconstrained_draws <-  lapply(constrained_draws, unconstrain_draws, posterior)
  } else if (modelname == "gp_pois_regr-gp_pois_regr") {
    constrained_draws <- modify_draws_gp_pois_regr(po)
    unconstrained_draws <-  lapply(constrained_draws, unconstrain_draws, posterior)
  } else {
    unconstrained_draws <-  lapply(gsd, unconstrain_draws, posterior)
  }
  
  ref_samples = rbind(unconstrained_draws[[1]], unconstrained_draws[[2]], 
                      unconstrained_draws[[3]], unconstrained_draws[[4]],
                      unconstrained_draws[[5]], unconstrained_draws[[6]],
                      unconstrained_draws[[7]], unconstrained_draws[[8]],
                      unconstrained_draws[[9]], unconstrained_draws[[10]])
  
  b = wpp(ref_samples, mass = rep(1 / nrow(ref_samples), nrow(ref_samples)))
  
  for(j in 1:M){
    cat(j, "\t")
    ### samples from pathfinder ###
    # I = 5
    pick_samples_M <- lp_multi_opath_I5[[i]][[j]]$pick_samples
    a_M = wpp(t(pick_samples_M), 
              mass = rep(1 / ncol(pick_samples_M), ncol(pick_samples_M)))
    
    w_d_pf_I5 <- wasserstein(a_M, b, p = 2); w_d_pf_I5
    W_d_100_pf_I5[j, i] = w_d_pf_I5
    
    # I = 40
    pick_samples_M <- lp_multi_opath_I40[[i]][[j]]$pick_samples
    a_M = wpp(t(pick_samples_M), 
              mass = rep(1 / ncol(pick_samples_M), ncol(pick_samples_M)))
    
    w_d_pf_I40 <- wasserstein(a_M, b, p = 2); w_d_pf_I40
    W_d_100_pf_I40[j, i] = w_d_pf_I40
    
    cat("pf I = 5:", w_d_pf_I5, "\t",
        "pf I = 20:", W_d_100_pf_IR[j, i], "\t",
        "pf I = 40:", w_d_pf_I40, "\n")
  }
}
proc.time() - t_0

save(file = "../results/wasserstein_100_sen_I.RData",
     list = c("W_d_100_pf_I5", "W_d_100_pf_I40"))


## check the output ##
pf_w_d_I5_qtils <- apply(W_d_100_pf_I5, 2, 
                          f <- function(x){quantile(x, c(0.25, 0.5, 0.75))})
pf_w_d_I20_qtils <- apply(W_d_100_pf_IR, 2, 
                      f <- function(x){quantile(x, c(0.25, 0.5, 0.75))})
pf_w_d_I40_qtils <- apply(W_d_100_pf_I40, 2, 
                         f <- function(x){quantile(x, c(0.25, 0.5, 0.75))})

for(i in 1:49){
  modelname <- pn[model_record[i]]
  printf("model %d: %s", model_record[i], modelname)
  print(cbind(pf_w_d_I5_qtils[, i], pf_w_d_I20_qtils[, i], 
              pf_w_d_I40_qtils[, i]))
}

