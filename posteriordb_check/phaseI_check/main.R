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
width = 600; height = 500 # the size of the plot
mc.cores = parallel::detectCores() - 2
sample_seed = 1234

# preallocate results #
lp_explore_n_iters <- array(data = NA, dim = c(M, L_pn))
lp_explore_n_leapfrog <- array(data = NA, dim = c(M, L_pn))
lp_INV <- array(data = NA, dim = c(2, L_pn))
lp_mean <- c()
lp_data <- c()

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
  lp_INV[, l] = INV[c(1, 2)]
  lp_mean[l] = INV[3] 
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
  L_p = L
  p_lp_trace = data.frame(iter = rep(1:L_p, M), 
                          chain = rep(paste(1:M), each = L_p),
                          lp__ = c(lp_phI[1:L_p, ]))
  
  lp_data[[l]] <- p_lp_trace 
  
  p_lp_s <- ggplot(data = p_lp_trace, 
                 aes(x=iter, y=lp__, group=chain, color=chain)) + geom_line() +
    geom_hline(yintercept = INV[c(1, 2)]) + 
    geom_hline(yintercept = INV[3], linetype = 2, colour = "blue") + 
    ylim(INV[1] - 1.5*(INV[2] - INV[1]), INV[2] + 1*(INV[2] - INV[1])) + 
    ggtitle(paste("model:", modelname)) + theme_bw() 
  jpeg(filename = paste0("../pics/phI_stan/No",l,"-", modelname, "_s.jpeg"),
       width = width, height = height, units = "px", pointsize = 12)
  print(p_lp_s)
  dev.off()
  
  p_lp_L <- ggplot(data = p_lp_trace, 
                   aes(x=iter, y=lp__, group=chain, color=chain)) + 
    geom_line() +
    geom_hline(yintercept = INV[c(1, 2)]) + 
    geom_hline(yintercept = INV[3], linetype = 2, colour = "blue") + 
    ggtitle(paste("model:", modelname)) + theme_bw() 
  jpeg(filename = paste0("../pics/phI_stan/No",l,"-", modelname, "_L.jpeg"),
       width = width, height = height, units = "px", pointsize = 12)
  print(p_lp_L)
  dev.off()
  
  p_lp <- ggplot(data = p_lp_trace, 
                   aes(x=iter, y=lp__, group=chain, color=chain)) + 
    geom_line() +
    geom_hline(yintercept = INV[c(1, 2)]) + 
    geom_hline(yintercept = INV[3], linetype = 2, colour = "blue") +
    ylim(min(INV[1] - 2*(INV[2] - INV[1]), quantile(p_lp_trace$lp__, 0.05)),
         max(INV[2] + 1*(INV[2] - INV[1]), max(p_lp_trace$lp__))) + 
    ggtitle(paste("model:", modelname)) + theme_bw() 
  jpeg(filename = paste0("../pics/phI_stan/No",l,"-", modelname, ".jpeg"),
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
    geom_hline(yintercept = INV[c(1, 2)]) + 
    geom_hline(yintercept = INV[3], linetype = 2, colour = "blue") +
    ylim(min(INV[1] - 2*(INV[2] - INV[1]), quantile(p_lp_trace$lp__, 0.15)),
         max(INV[2] + 1*(INV[2] - INV[1]), max(p_lp_trace$lp__))) + 
    ggtitle(paste("model:", modelname)) + theme_bw() 
  jpeg(filename = paste0("../pics/phI_stan/No",l,"-", modelname, "_trun.jpeg"),
       width = width, height = height, units = "px", pointsize = 12)
  print(p_lp)
  dev.off()
  
}

# save(file = "../results/lp_posteriordb_explore3.RData",
#      list = c("lp_explore_n_iters", "lp_explore_n_leapfrog",
#               "lp_INV", "lp_mean", "lp_data"))
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

## check the distribution of number of iterations ##
n_iters_mean <- colMeans(lp_explore_n_iters[, -takeoff])
mean(n_iters_mean, na.rm = TRUE) # 69.11596
median(lp_explore_n_iters[, -takeoff], na.rm = TRUE) # 26
sd(n_iters_mean, na.rm = TRUE)   # 144.9198
sd(lp_explore_n_iters[, -takeoff], na.rm = TRUE) # 159.7554
jpeg(filename = paste0("../pics/hist_iters.jpeg"),
     width = width, height = height, units = "px", pointsize = 12)
hist(lp_explore_n_iters[, -takeoff], breaks = 100, 
     main = "", ylab = "", axes = TRUE,
     xlab = "iterations")
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
median(lp_explore_n_leapfrog[, -takeoff], na.rm = TRUE) # 645.5
sd(n_leapfrog_mean, na.rm = TRUE)   # 94313.19
sd(lp_explore_n_leapfrog[, -takeoff], na.rm = TRUE) # 98874.89
jpeg(filename = paste0("../pics/hist_leapfrogs.jpeg"),
     width = width, height = height, units = "px", pointsize = 12)
hist(lp_explore_n_leapfrog[, -takeoff], breaks = 200, 
     main = "No. leapfrogs to reach target interval",
     xlab = "leapfrogs")
dev.off()

jpeg(filename = paste0("../pics/hist_leapfrogs_log.jpeg"),
     width = width/2, height = height/2, units = "px", pointsize = 12)
df <- data.frame(sum_leapfrog = c(lp_explore_n_leapfrog[, -takeoff][
  !is.na(lp_explore_n_leapfrog[, -takeoff])]))
p_leapfrog <- ggplot(data =df , aes(x = sum_leapfrog)) +
  geom_histogram(color="black", fill="white", bins = 60) + scale_x_log10() +
  xlab("No. of leapfrogs") 
#+ labs(title = "No. of leapfrogs to reach target interval")
print(p_leapfrog)
dev.off()


#' Around 80.1% of phase I MCMC chains spend less than 4000
#' leapfrogs for lp__ to reach the 99% posterior interval.
sum((lp_explore_n_leapfrog[, -takeoff] <= 4000), na.rm = TRUE) / 
  sum(!is.na(lp_explore_n_leapfrog[, -takeoff]))
mean(lp_explore_n_leapfrog[, -takeoff][which(lp_explore_n_leapfrog[, -takeoff] < 4e3)])
#[1] 852.5551
median(lp_explore_n_leapfrog[, -takeoff][which(lp_explore_n_leapfrog[, -takeoff] < 4e3)])
# [1] 385
sd(lp_explore_n_leapfrog[, -takeoff][which(lp_explore_n_leapfrog[, -takeoff] < 4e3)])
#[1] 958.3538


#' Around 96.7% of phase I MCMC chains spend less than 30,000
#' leapfrogs for lp__ to reach the 99% posterior interval
sum((lp_explore_n_leapfrog[, -takeoff] <= 3e4), na.rm = TRUE) / 
  sum(!is.na(lp_explore_n_leapfrog[, -takeoff]))

mean(lp_explore_n_leapfrog[, -takeoff][which(lp_explore_n_iters[, -takeoff] < 200)])
# [1] 2964.708
median(lp_explore_n_leapfrog[, -takeoff][which(lp_explore_n_iters[, -takeoff] < 200)])
# [1] 541
sd(lp_explore_n_leapfrog[, -takeoff][which(lp_explore_n_iters[, -takeoff] < 200)])
# [1] 6041.258

#' The 3th, 9th, 10th, 14th, 24th, 27th and 37th model have phase I MCMC chains
#' fail to reach the target interval within 30,000 leapfrogs
table(as.integer(which(lp_explore_n_leapfrog > 3e4) / M - 0.5 / M) + 1)
# 3  9 10 14 24 27 37 
# 1 20  1  2 20  1  6


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

