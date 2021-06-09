rm(list = ls())
setwd("./posteriordb") # set working dir to cloned package
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
## install the beta release version of R package posterior
# install.packages("posterior", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# check the dataset in posteriordb #
library(posteriordb)
library(posterior)
library(ggplot2)
library(cmdstanr)
library(bayesplot)
source("../utils/sim_pf.R")
source("../utils/lp_utils.R")

set.seed(123)
pd <- pdb_local() # Posterior database connection
pn <- posterior_names(pd)
L_pn = length(pn)

# parameters settings #
alpha = 0.01
L = 1000
M = 20
width = 600; height = 500 # the size of the plot
mc.cores = parallel::detectCores() - 2
sample_seed = 1234

load(file = "../results/lp_posteriordb_LBFGS_h6.RData")

# preallocate results #
lp_explore_n_iters <- array(data = NA, dim = c(M, length(model_record)))
lp_explore_n_leapfrog <- array(data = NA, dim = c(M, length(model_record)))
lp_data <- c()
PhaseI_last_draw <- list() # Get the last samples of Phase I
PhI_leapfrog_counts <- array(data = NA, dim = c(M, length(model_record)))

i = which(model_record == 27)
for(i in 41:49){ #20 length(model_record)
  modelname <- pn[model_record[i]]
  printf("model %d: %s", model_record[i], modelname)
  # pick model
  po <- posterior(modelname, pdb = pd)
  # get reference posterior samples
  gsd <- reference_posterior_draws(po)
  
  # compile the model
  sc <- stan_code(po)
  model <- stan_model(model_code = sc)
  
  ### setup initials ###
  data = get_data(po)
  posterior <- to_posterior(model, data)
  init = lapply(initial_ls[[i]], f <- function(x){constrain_pars(posterior, x)})
  
  ###  run Stan with a long Phase I warmup time  ###
  # need to run the following line to generate the stan file. no need to run if there is one. 
  # write_stan_file(sc, dir = paste0(getwd(), "/modelcode"),
  #                 basename = paste0(modelname, ".stan")) 
  file <- file.path(getwd(), "modelcode", paste0(modelname, ".stan"))
  mod <- cmdstan_model(file)
  
  ###  random inits  ###
  suppressWarnings(
    fit <- mod$sample(
      data = data,
      seed = 123,
      init = init,
      chains = M,
      parallel_chains = 5,
      refresh = 0,
      save_warmup = TRUE,
      iter_warmup = L,
      iter_sampling = 0,
      init_buffer = L,
      term_buffer = 0,
      show_messages = FALSE,
      sig_figs = 16
    ))
  p1 <- mcmc_trace(fit$draws("lp__", inc_warmup = TRUE)[60:80, ,], iter1 = 60) #c("lp__", "phi[1]", "lambda[1]", "theta1[1]")
  print(p1)
  
  ## record the initial and the last sample of phase I ##
  fit_draws <- fit$draws(inc_warmup = TRUE)
  last_pI_draws <- matrix(fit_draws[75, , ], 
                          nrow = dim(fit_draws)[2])
  colnames(last_pI_draws) <- dimnames(fit_draws)$variable
  unconstrained_draws <- unconstrain_cmd_draws(last_pI_draws, posterior)
  PhaseI_last_draw[[i]] <- unconstrained_draws # record the last sample of phase I
  
  ###  record the number of iterations required to reach INV ###
  # Get the number of iterations and  leapfrogs #
  lp_explore_sum <- lp_explore(fit, lp_INV[, i], L, M, model, data)
  lp_explore_n_iters[, i] = lp_explore_sum$n_iters
  lp_explore_n_leapfrog[, i] = lp_explore_sum$n_sum_leapfrog
  printf("the maximum iter to reach %.1f %% posterior interval of lp__ is %d",
         (1.0 - alpha) * 100, max(lp_explore_sum$n_iters))
  printf("the average leapfrogs is %.2f, sd is %.2f", 
         mean(lp_explore_sum$n_sum_leapfrog), sd(lp_explore_sum$n_sum_leapfrog))
  
  # Get the cost of Phase I warmup
  PhI_leapfrog_counts[, i] <- colSums(
    fit$sampler_diagnostics(inc_warmup = TRUE)[1:75, , "n_leapfrog__"])
  
  # check the trace plot of lp__ #
  pos_d <- list()
  for(ll in 1:dim(fit_draws)[2]){
    pos_d[[ll]] <- as.data.frame(fit_draws[, ll, -1])
  }
  lp_phI <- lp_recover(model, data, pos_d)
  
  L_p = L
  p_lp_trace = data.frame(iter = rep(1:L_p, M), 
                          chain = rep(paste(1:M), each = L_p),
                          lp__ = c(lp_phI[1:L_p, ]))
  
  lp_data[[i]] <- p_lp_trace 
  
  p_lp_s <- ggplot(data = p_lp_trace, 
                   aes(x=iter, y=lp__, group=chain, color=chain)) + geom_line() +
    geom_hline(yintercept = lp_INV[, i]) + 
    geom_hline(yintercept = lp_mean[i], linetype = 2, colour = "blue") + 
    ylim(lp_INV[1, i] - 1.5*(lp_INV[2, i] - lp_INV[1, i]), 
         lp_INV[2, i] + 1*(lp_INV[2, i] - lp_INV[1, i])) + 
    ggtitle(paste("model:", modelname)) + theme_bw() 
  jpeg(filename = paste0("../pics/phI_stan/No", model_record[i], "-", 
                         modelname, "_s.jpeg"),
       width = width, height = height, units = "px", pointsize = 12)
  print(p_lp_s)
  dev.off()
  
  p_lp_L <- ggplot(data = p_lp_trace, 
                   aes(x=iter, y=lp__, group=chain, color=chain)) + 
    geom_line() +
    geom_hline(yintercept = lp_INV[, i]) + 
    geom_hline(yintercept = lp_mean[i], linetype = 2, colour = "blue") + 
    ggtitle(paste("model:", modelname)) + theme_bw() 
  jpeg(filename = paste0("../pics/phI_stan/No", model_record[i], "-", 
                         modelname, "_L.jpeg"),
       width = width, height = height, units = "px", pointsize = 12)
  print(p_lp_L)
  dev.off()
  
  p_lp <- ggplot(data = p_lp_trace, 
                 aes(x=iter, y=lp__, group=chain, color=chain)) + 
    geom_line() +
    geom_hline(yintercept = lp_INV[, i]) + 
    geom_hline(yintercept = lp_mean[i], linetype = 2, colour = "blue") +
    ylim(min(lp_INV[1, i] - 2*(lp_INV[2, i] - lp_INV[1, i]), 
             quantile(p_lp_trace$lp__, 0.05)),
         max(lp_INV[2, i] + 1*(lp_INV[2, i] - lp_INV[1, i]), 
             max(p_lp_trace$lp__))) + 
    ggtitle(paste("model:", modelname)) + theme_bw() 
  jpeg(filename = paste0("../pics/phI_stan/No", model_record[i], "-", 
                         modelname, ".jpeg"),
       width = width, height = height, units = "px", pointsize = 12)
  print(p_lp)
  dev.off()
  
  L_p = ifelse((max(lp_explore_sum$n_iters) <= 40), 50, # pick the x-axis range of the plot
               min(as.integer(1.3*max(lp_explore_sum$n_iters)), L))
  p_lp_trace = data.frame(iter = rep(1:L_p, M), 
                          chain = rep(paste(1:M), each = L_p),
                          lp__ = c(lp_phI[1:L_p, ]))
  
  p_lp <- ggplot(data = p_lp_trace, 
                 aes(x=iter, y=lp__, group=chain, color=chain)) + 
    geom_line() +
    geom_hline(yintercept = lp_INV[, i]) + 
    geom_hline(yintercept = lp_mean[i], linetype = 2, colour = "blue") +
    ylim(min(lp_INV[1, i] - 2*(lp_INV[2, i] - lp_INV[1, i]), 
             quantile(p_lp_trace$lp__, 0.15)),
         max(lp_INV[2, i] + 1*(lp_INV[2, i] - lp_INV[1, i]), 
             max(p_lp_trace$lp__))) + 
    ggtitle(paste("model:", modelname)) + theme_bw() 
  jpeg(filename = paste0("../pics/phI_stan/No", model_record[i], "-", 
                         modelname, "_trun.jpeg"),
       width = width, height = height, units = "px", pointsize = 12)
  print(p_lp)
  dev.off()
}

save(file = "../results/lp_posteriordb_explore_h6.RData",
     list = c("lp_explore_n_iters", "lp_explore_n_leapfrog", "lp_data", 
              "PhaseI_last_draw", "PhI_leapfrog_counts"))

# check the output #
load("../results/lp_posteriordb_explore_h6.RData")

## check the distribution of number of iterations ##
n_iters_mean <- colMeans(lp_explore_n_iters)
mean(n_iters_mean, na.rm = TRUE) # 64.78
median(lp_explore_n_iters, na.rm = TRUE) # 26
sd(n_iters_mean, na.rm = TRUE)   # 141.2
sd(lp_explore_n_iters, na.rm = TRUE) # 142.45
jpeg(filename = paste0("../pics/hist_iters.jpeg"),
     width = width, height = height, units = "px", pointsize = 12)
hist(lp_explore_n_iters, breaks = 100, 
     main = "", ylab = "", axes = TRUE,
     xlab = "iterations")
dev.off()


#' Around 96.5% of phase I MCMC chains reach the target interval within 200 '
#' iterations. 
sum((lp_explore_n_iters <= 200), na.rm = TRUE) / 
  sum(!is.na(lp_explore_n_iters))

#' The 3th, 9th, and 41th model have phase I MCMC chains
#' fail to reach the target interval within 1000 iters
table(model_record[as.integer(which(lp_explore_n_iters == L) / M - 0.5 / M) 
                   + 1])
#9 41 
#20  1 


## check the distribution of sum of leapfrogs ##
n_leapfrog_mean <- colMeans(lp_explore_n_leapfrog)
mean(n_leapfrog_mean, na.rm = TRUE) # 16451.39
median(lp_explore_n_leapfrog, na.rm = TRUE) # 612
sd(n_leapfrog_mean, na.rm = TRUE)   # 91797.28
sd(lp_explore_n_leapfrog, na.rm = TRUE) # 91020.25
jpeg(filename = paste0("../pics/hist_leapfrogs.jpeg"),
     width = width, height = height, units = "px", pointsize = 12)
hist(lp_explore_n_leapfrog, breaks = 200, 
     main = "No. leapfrogs to reach target interval",
     xlab = "leapfrogs")
dev.off()

jpeg(filename = paste0("../pics/hist_leapfrogs_log.jpeg"),
     width = width/2, height = height/2, units = "px", pointsize = 12)
df <- data.frame(sum_leapfrog = c(lp_explore_n_leapfrog[
  !is.na(lp_explore_n_leapfrog)]))
p_leapfrog <- ggplot(data =df , aes(x = sum_leapfrog)) +
  geom_histogram(color="black", fill="white", bins = 60) + scale_x_log10() +
  xlab("No. of leapfrogs") 
#+ labs(title = "No. of leapfrogs to reach target interval")
print(p_leapfrog)
dev.off()


#' Around 79.7% of phase I MCMC chains spend less than 4000
#' leapfrogs for lp__ to reach the 99% posterior interval.
sum((lp_explore_n_leapfrog <= 4000), na.rm = TRUE) / 
  sum(!is.na(lp_explore_n_leapfrog))
mean(lp_explore_n_leapfrog[which(lp_explore_n_leapfrog < 4e3)])
#[1] 789.7
median(lp_explore_n_leapfrog[which(lp_explore_n_leapfrog < 4e3)])
# [1] 336
sd(lp_explore_n_leapfrog[which(lp_explore_n_leapfrog < 4e3)])
#[1] 928.7


#' Around 96.9% of phase I MCMC chains spend less than 30,000
#' leapfrogs for lp__ to reach the 99% posterior interval
sum((lp_explore_n_leapfrog <= 3e4), na.rm = TRUE) / 
  sum(!is.na(lp_explore_n_leapfrog))

mean(lp_explore_n_leapfrog[which(lp_explore_n_iters < 200)])
# [1] 3186.357
median(lp_explore_n_leapfrog[which(lp_explore_n_iters < 200)])
# [1] 475.5
sd(lp_explore_n_leapfrog[which(lp_explore_n_iters < 200)])
# [1] 6669.635

#' The 3th, 9th, 14th and 37th model have phase I MCMC chains
#' fail to reach the target interval within 30,000 leapfrogs
table(model_record[as.integer(which(lp_explore_n_leapfrog > 3e4) / M - 0.5 / M) 
                   + 1])
#9 11 14 37 48 
#20  3  1  5  1


#' check the histogram of number of leapfrogs, distribution of the sum of 
#' leapfrogs before reaching target interval is highly right skewed
lp_explore_n_leapfrog[which(lp_explore_n_leapfrog >= 1e4 & 
                              lp_explore_n_iters < 300)]
table(model_record[as.integer(which(lp_explore_n_leapfrog >= 1e4 & 
                                      lp_explore_n_iters < 300)/M - 0.5/M) + 1])

# 4 10 11 12 14 37 48 
# 2 10 18  9 19 17  7 


