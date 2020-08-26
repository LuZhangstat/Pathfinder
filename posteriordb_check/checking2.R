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

takeoff <- c(21, 24)   # model 21, 24 has transformed pars
load("../results/lp_posteriordb_explore.RData")

## find the model with fitting results
N_models = 0
model_record = c()
for(l in 1:L_pn){
  modelname <- pn[l]
  # pick model
  po <- posterior(modelname, pdb = pd)
  skip_to_next <- FALSE
  tryCatch(gsd <- reference_posterior_draws(po),
           error = function(e) { skip_to_next <<- TRUE})
  if(any(takeoff == l)){skip_to_next <<- TRUE}
  if(skip_to_next) {
    # print("Error in obtaining reference posterior for this posterior.")
    next }
  N_models = N_models + 1
  model_record = c(model_record, l)
}
N_models

df <- data.frame(n_iters = c(lp_explore_n_iters[, model_record]), 
                 n_leapfrogs = 
                   c(lp_explore_n_leapfrog[, model_record]),
                 model = rep(pn[model_record], each = M))

jpeg(filename = paste0("../pics/box_iters_log.jpeg"),
     width = width*1.3, height = height*2, units = "px", pointsize = 12)
p_box_iter <- ggplot(df, aes(y = reorder(model, n_iters, FUN = median), 
                             x = n_iters)) + geom_boxplot() + 
  scale_x_log10() + ylab("") + xlab("No. of iterations") + 
  theme_grey(base_size = 26)
print(p_box_iter)
dev.off()

jpeg(filename = paste0("../pics/box_leapfrogs_log.jpeg"),
     width = width*1.3, height = height*2, units = "px", pointsize = 12)
p_box_leapfrog <- ggplot(df, aes(y = reorder(model, n_leapfrogs, FUN = median), 
                                 x = n_leapfrogs)) + 
  geom_boxplot() + scale_x_log10() + ylab("") + xlab("No. of leapfrogs") + 
  theme_grey(base_size = 26)
p_box_leapfrog
dev.off()



