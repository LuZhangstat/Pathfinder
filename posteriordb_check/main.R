rm(list = ls())
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
## install the beta release version of R package posterior
# install.packages("posterior", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# check the dataset in posteriordb #
library(posteriordb)
library(posterior)
library(ggplot2)
source("./utils/sim.R")
source("./utils/lp_utils.R")

set.seed(123)
pd <- pdb_default() # Posterior database connection
pn <- posterior_names(pd)
L_pn = length(pn)

# parameters settings #
alpha = 0.01
L = 1000
M = 20
L_p = 50
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
  gsd <- reference_posterior_draws(po)
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
  jpeg(filename = paste0("./pics/No",l,"-", modelname, ".jpeg"),
       width = width, height = height, units = "px", pointsize = 12)
  print(p_lp)
  dev.off()
}

save(file = "./results/lp_posteriordb_explore.RData", 
     list = c("lp_explore_n_iters", "lp_explore_n_leapfrog",
              "lp_INV"))

# check reference posterior
for(l in 1:L_pn){
  modelname <- pn[l]
  printf("model %d: %s", l, modelname)
  
  skip_to_next <- FALSE
  # pick model
  tryCatch(po <- posterior(modelname, pdb = pd), 
           error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next) { 
    print("Error in obtaining model info.")
    next }  
  # get reference posterior samples
  tryCatch(gsd <- reference_posterior_draws(po), 
           error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next) { 
    print("Error in obtaining reference posterior for this posterior.")
    next }  
}
