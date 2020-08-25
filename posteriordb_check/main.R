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

# check the transform parameter block
# for(id in model_record){
id = 95
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

## check the distribution of number of iterations ##
n_iters_mean <- colMeans(lp_explore_n_iters[, -takeoff])
mean(n_iters_mean, na.rm = TRUE) # 69.11596
sd(n_iters_mean, na.rm = TRUE)   # 144.9198
sd(lp_explore_n_iters[, -takeoff], na.rm = TRUE) # 159.7554
jpeg(filename = paste0("../pics/hist_iters.jpeg"),
     width = width, height = height, units = "px", pointsize = 12)
hist(lp_explore_n_iters[, -takeoff], breaks = 100, 
     main = "No. iters to reach target interval",
     xlab = "No. iters")
dev.off()


#' Around 95.9% of phase I MCMC chains reach the target interval within 200 '
#' iterations. 
sum((lp_explore_n_iters[, -takeoff] <= 200), na.rm = TRUE) / 
  sum(!is.na(lp_explore_n_iters[, -takeoff]))

#' The 9th, 24th, 27th, 40th and 41th model have phase I MCMC chains
#' fail to reach the target interval within 1000 iters
table(as.integer(which(lp_explore_n_iters == L) / M - 0.5 / M) + 1)
# 9 24 27 40 41 
# 20 20  1  1  1 

## check the distribution of sum of leapfrogs ##
n_leapfrog_mean <- colMeans(lp_explore_n_leapfrog[, -takeoff])
mean(n_leapfrog_mean, na.rm = TRUE) # 18013.83
sd(n_leapfrog_mean, na.rm = TRUE)   # 94313.19
sd(lp_explore_n_leapfrog[, -takeoff], na.rm = TRUE) # 98874.89
jpeg(filename = paste0("../pics/hist_leapfrogs.jpeg"),
     width = width, height = height, units = "px", pointsize = 12)
hist(lp_explore_n_leapfrog[, -takeoff], breaks = 200, 
     main = "No. leapfrogs to reach target interval",
     xlab = "No. leapfrogs")
dev.off()

#' Around 80.1% of phase I MCMC chains spend less than 4000
#' leapfrogs for lp__ to reach the 99% posterior interval.
sum((lp_explore_n_leapfrog[, -takeoff] <= 4000), na.rm = TRUE) / 
  sum(!is.na(lp_explore_n_leapfrog[, -takeoff]))
mean(lp_explore_n_leapfrog[, -takeoff][which(lp_explore_n_leapfrog[, -takeoff] < 4e3)])
#[1] 852.5551
sd(lp_explore_n_leapfrog[, -takeoff][which(lp_explore_n_leapfrog[, -takeoff] < 4e3)])
#[1] 958.3538


#' Around 96.7% of phase I MCMC chains spend less than 30,000
#' leapfrogs for lp__ to reach the 99% posterior interval
sum((lp_explore_n_leapfrog[, -takeoff] <= 3e4), na.rm = TRUE) / 
  sum(!is.na(lp_explore_n_leapfrog[, -takeoff]))

mean(lp_explore_n_leapfrog[, -takeoff][which(lp_explore_n_iters[, -takeoff] < 200)])
# [1] 2964.708
sd(lp_explore_n_leapfrog[, -takeoff][which(lp_explore_n_iters[, -takeoff] < 200)])
# [1] 6041.258

#' The 3th, 9th, 10th, 14th, 24th, 27th and 37th model have phase I MCMC chains
#' fail to reach the target interval within 30,000 leapfrogs
table(as.integer(which(lp_explore_n_leapfrog > 3e4) / M - 0.5 / M) + 1)
# 3  9 10 14 24 27 37 
# 1 20  1  2 20  1  6

#' In summary, the current Stan algorithm for phase I is slow for model 9 and 24

#' After checking model 9 and 24, I found that the lp__ posterior interval based
#' on the reference posterior samples seems to be incorrect. See, e.g. 
#' No9-earnings-earn_height_2.jpeg. One chain for model 27 got stuck, 
#' (see No27-hudson_lynx_hare-lotka_volterra.jpeg), but the remaining chains 
#' are fine,(see No27-hudson_lynx_hare-lotka_volterra_no6.jpeg). 
#' Similar for model 40 and 41. After taking off model 9 and 24 and chain 
#' with problem, the hist of sum_leapfrog is:
hist(lp_explore_n_iters[-which(lp_explore_n_iters == 1000)], 
     breaks = 100, 
     main = "No. iters to reach target interval, truncate at 1000",
     xlab = "No. iters")

#' Over 99% of phase I MCMC chains of lp__ reach the target interval 
sum((lp_explore_n_iters[-which(lp_explore_n_iters == 1000)] <= 300), 
    na.rm = TRUE) / 
  sum(!is.na(lp_explore_n_iters[-which(lp_explore_n_iters == 1000)]))

#' The current Stan algorithm for phase I is slow for model 2 and 3 
lp_explore_n_iters[which(lp_explore_n_iters < 1000 & lp_explore_n_iters >= 300)] 
table(as.integer(which(lp_explore_n_iters < 1000 & 
                         lp_explore_n_iters >= 300) / M - 0.5 / M) + 1)
# 2  3 27 
# 3  4  1

#' check the histogram of number of leapfrogs, distribution of the sum of 
#' leapfrogs before reaching target interval is highly right skewed
lp_explore_n_leapfrog[which(lp_explore_n_leapfrog >= 1e4 & 
                              lp_explore_n_iters < 300)]
# [1] 13832 18261 20524 17962 17129 16521 18351 17628 31223 18166 14796 14512 13169 16890
# [15] 27862 14286 21736 13788 23826 20992 15102 10537 17369 20708 22096 19912 14466 13913
# [29] 11544 15548 14767 11165 17868 12032 13170 13517 14020 23162 13882 19798 24102 25818
# [43] 29549 38597 19325 17157 21132 22461 26490 14364 19821 13795 13267 17174 35202 19989
# [57] 12834 18480 34095 26612 60107 21014 26124 12589 10878 16461 12443 24980 38539 33636
# [71] 25617 37588 28496 15163 40488 17515 19067 10190 13982 11008 14914 13890
table(as.integer(which(lp_explore_n_leapfrog >= 1e4 & 
                         lp_explore_n_iters < 300)/M - 0.5/M) + 1)

# 4 10 11 12 14 37 48 
# 2 10 18  9 19 17  7 

