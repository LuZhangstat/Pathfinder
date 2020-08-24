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
source("../utils/sim.R")
source("../utils/lp_utils.R")

set.seed(123)
pd <- pdb_local() # Posterior database connection
pn <- posterior_names(pd)
L_pn = length(pn)

# parameters settings #
alpha = 0.01
L = 1000
M = 20
width = 860; height = 740 # the size of the plot
mc.cores = parallel::detectCores() - 2
sample_seed = 1234

# preallocate results #
lp_explore_n_iters <- array(data = NA, dim = c(M, L_pn))
lp_explore_n_leapfrog <- array(data = NA, dim = c(M, L_pn))
lp_INV <- array(data = NA, dim = c(2, L_pn))

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
  ###  run Stan with a long Phase I warmup time  ###
  suppressWarnings(
    phiI_sample <- sampling(model, data = get_data(po), 
                          seed = sample_seed,
                          iter = L + 1, 
                          warmup = L,
                          chains = M, 
                          cores = mc.cores,
                          algorithm ="NUTS",
                          control = list(adapt_init_buffer = L,
                                         adapt_term_buffer = 0,
                                         adapt_window = 0),
                          save_warmup = TRUE, 
                          refresh = 0))
  
  ###  record the number of iterations required to reach INV ###
  # Get the number of iterations and  leapfrogs #
  lp_explore_sum <- lp_explore(phiI_sample, INV, L, M)
  lp_explore_n_iters[, l] = lp_explore_sum$n_iters
  lp_explore_n_leapfrog[, l] = lp_explore_sum$n_sum_leapfrog
  printf("the maximum iter to reach %.1f %% posterior interval of lp__ is %d",
        (1.0 - alpha) * 100, max(lp_explore_sum$n_iters))
  printf("the average leapfrogs is %.2f, sd is %.2f", 
         mean(lp_explore_sum$n_sum_leapfrog), sd(lp_explore_sum$n_sum_leapfrog))
  
  # check the trace plot of lp__ #
  lp_phI <- ls_lp_phI(phiI_sample, L)
  L_p = ifelse((max(lp_explore_sum$n_iters) <= 40), 50, # pick the x-axis range of the plot
               min(as.integer(1.3*max(lp_explore_sum$n_iters)), L))
  p_lp_trace = data.frame(iter = rep(1:L_p, M), 
                          chain = rep(paste(1:M), each = L_p),
                          lp__ = c(lp_phI[1:L_p, ]))
  p_lp <- ggplot(data = p_lp_trace, 
                 aes(x=iter, y=lp__, group=chain, color=chain)) + geom_line() +
    geom_hline(yintercept = INV)
  jpeg(filename = paste0("../pics/No",l,"-", modelname, ".jpeg"),
       width = width, height = height, units = "px", pointsize = 12)
  print(p_lp)
  dev.off()
}

save(file = "../results/lp_posteriordb_explore.RData", 
     list = c("lp_explore_n_iters", "lp_explore_n_leapfrog",
              "lp_INV"))

## check reference posterior ##
N_models = 0
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
}
N_models
# only 49 out of 97 models have reference posterior samples

## check the distribution of number of iterations ##
n_iters_mean <- colMeans(lp_explore_n_iters)
mean(n_iters_mean, na.rm = TRUE) #86.74592
sd(n_iters_mean, na.rm = TRUE)   #194.8229
hist(lp_explore_n_iters, breaks = 100, 
     main = "No. iters to reach target interval",
     xlab = "No. iters")

#' Around 94% of phase I MCMC chains reach the target interval within 200 '
#' iterations. 
sum((lp_explore_n_iters <= 200), na.rm = TRUE) / 
  sum(!is.na(lp_explore_n_iters))

#' The 9th, 24th, 27th, 40th and 41th model have phase I MCMC chains
#' fail to reach the target interval within 1000 iters
table(as.integer(which(lp_explore_n_iters == L) / M - 0.5 / M) + 1)
# 9 24 27 40 41 
# 20 20  1  1  1 

## check the distribution of sum of leapfrogs ##
n_leapfrog_mean <- colMeans(lp_explore_n_leapfrog)
mean(n_leapfrog_mean, na.rm = TRUE) #21602.9
sd(n_leapfrog_mean, na.rm = TRUE)   #96441.13
hist(lp_explore_n_leapfrog, breaks = 200, 
     main = "No. iters to reach target interval",
     xlab = "No. iters")

#' Around 66% of phase I MCMC chains spend less than 2000
#' leapfrogs for lp__ to reach the 99% posterior interval
sum((lp_explore_n_leapfrog <= 2000), na.rm = TRUE) / 
  sum(!is.na(lp_explore_n_leapfrog))

#' Around 94.8% of phase I MCMC chains spend less than 30,000
#' leapfrogs for lp__ to reach the 99% posterior interval
sum((lp_explore_n_leapfrog <= 3e4), na.rm = TRUE) / 
  sum(!is.na(lp_explore_n_leapfrog))

#' The 3th, 9th, 10th, 14th, 24th, 27th and 37th model have phase I MCMC chains
#' fail to reach the target interval within 30,000 leapfrogs
table(as.integer(which(lp_explore_n_leapfrog > 3e4) / M - 0.5 / M) + 1)

#' In summary, the current Stan algorithm for phase I is slow for model 9 and 24
