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
library(dplyr)
source("../utils/sim_pf.R")
source("../utils/lp_utils.R")

set.seed(123)
pd <- pdb_local() # Posterior database connection
pn <- posterior_names(pd)
L_pn = length(pn)

# parameters settings #
alpha = 0.01
L = 1000                # iteration for Stan phase I sampling
N1 = 1000               # maximum iteration for L-BFGS
M = 20                  # No. of chains
width = 860; height = 740 # the size of the plot
mc.cores = parallel::detectCores() - 2
factr_tol = 1e2
init_bound = 2.0

# get model names with reference posterior samples#
N_models = 0
model_record = c()
for(l in 1:L_pn){
  modelname <- pn[l]
  
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
  printf("model %d: %s", l, modelname)
}
N_models


# preallocate results #
lp_LBFGS_n_fn <- array(data = NA, dim = c(M, length(model_record)))
lp_LBFGS_n_gr <- array(data = NA, dim = c(M, length(model_record)))
lp_INV <- array(data = NA, dim = c(2, length(model_record)))
lp_mean <- c()
initial_ls <- list()

#i = which(model_record == 24);i
for(i in 1:length(model_record)){
  modelname <- pn[model_record[i]]
  printf("\n model %d: %s", model_record[i], modelname)
  # pick model
  po <- posterior(modelname, pdb = pd)
  # get reference posterior samples
  gsd <- reference_posterior_draws(po)
  
  # compile the model
  sc <- stan_code(po)
  model <- stan_model(model_code = sc)
  data <- get_data(po)
  
  # obtain posterior interval of lp__
  if(modelname == "eight_schools-eight_schools_noncentered"){
    INV <- lp_Int_q_posteriordb_8school_noncen(po, alpha)
  } else if (modelname == "gp_pois_regr-gp_pois_regr") {
    INV <- lp_Int_q_posteriordb_gp_pois_regr(po, alpha)
  } else {INV <- lp_Int_q_posteriordb(po, alpha)}
  lp_INV[, i] = INV[c(1, 2)]
  lp_mean[i] = INV[3] 
  
  posterior <- to_posterior(model, data)
  D <- get_num_upars(posterior)
  fn <- function(theta) -log_prob(posterior, theta, adjust_transform = TRUE, 
                                  gradient = TRUE)[1] 
  gr <- function(theta) -grad_log_prob(posterior, theta, adjust_transform = TRUE)
  
  # choose the size of history in L-BFGS
  lmm = 6;
  cat("No. pars:", D," lmm in L-BFGS: ", lmm, "\n")
  
  set.seed(1234)
  m = 1
  initial_ls[[i]] = list()
  while (m <= M){
    cat(m, "\t")
    # run M initials and check the optim algorithm
    init <- runif(D, -init_bound, init_bound)
    initial_ls[[i]][[m]] <- init
    restart <- FALSE
    tryCatch(lp_old <- -fn(init), 
             error = function(e) { restart <<- TRUE})
    tryCatch(test_uncon <- constrain_pars(posterior, init), 
             error = function(e) { restart <<- TRUE})
    if(restart){
      print("Error in initialization.")
      next
    }
    if((lp_old >= lp_INV[1, i] && lp_old <= lp_INV[2, i])){
      break_opt <- FALSE
      tryCatch(z <- optim(par = init,
                          fn = fn,  # negate for maximization
                          gr = gr,
                          method = "L-BFGS-B",
                          control = list(maxit = N1, factr = factr_tol,
                                         lmm = lmm #, ndeps = 1e-8 #, 
                                         #trace = 6, REPORT = 1 
                          )), 
               error = function(e) { break_opt <<- TRUE})
      if(break_opt) { 
        print("Error in obtaining optimization path.")
        next
      }
      print("initialized in target region")
      lp_LBFGS_n_fn[m, i] = 1
      lp_LBFGS_n_gr[m, i] = 1
      m = m + 1
      next
    } else if(lp_old > lp_INV[2, i] && is.finite(lp_old)){
      break_opt <- FALSE
      tryCatch(z <- optim(par = init,
                          fn = fn,  # negate for maximization
                          gr = gr,
                          method = "L-BFGS-B",
                          control = list(maxit = N1, factr = factr_tol, 
                                         lmm = lmm #, ndeps = 1e-8 #, 
                                         #trace = 6, REPORT = 1 
                          )), 
               error = function(e) { break_opt <<- TRUE})
      if(break_opt) { 
        print("Error in obtaining optimization path.")
        next
      }
      print("initial lp__ pass the target region")
      lp_LBFGS_n_fn[m, i] = 1
      lp_LBFGS_n_gr[m, i] = 1
      m = m + 1
      next
    } else if(is.infinite(lp_old)){
      print("initial lp__ is Inf")
      next
    }
    
    # run optimization
    for (n in 1:N1) {
      break_opt <- FALSE
      tryCatch(z <- optim(par = init,
                          fn = fn,  # negate for maximization
                          gr = gr,
                          method = "L-BFGS-B",
                          control = list(maxit = n, factr = factr_tol, 
                                         lmm = lmm #, ndeps = 1e-8 #, 
                                         #trace = 6, REPORT = 1 
                          )), 
               error = function(e) { break_opt <<- TRUE})
      if(break_opt) { 
        print("Error in obtaining optimization path.")
        break
      }
      lp_up <- -fn(z$par)
      
      if( (lp_up >= lp_INV[1, i] && lp_up <= lp_INV[2, i]) ){
        # once reach the target region, record the number of calls to fn and gr
        lp_LBFGS_n_fn[m, i] = z$counts[1]
        lp_LBFGS_n_gr[m, i] = z$counts[2]
        m = m + 1
        break
      } else if(lp_up > lp_INV[2, i]){
        print("optimization pass the target region")
        lp_LBFGS_n_fn[m, i] = z$counts[1]
        lp_LBFGS_n_gr[m, i] = z$counts[2]
        m = m + 1
        break
      } else if(lp_old == lp_up){
        lp_LBFGS_n_fn[m, i] = -z$counts[1]
        lp_LBFGS_n_gr[m, i] = -z$counts[1]
        print("converge to a local minimum.")
        m = m + 1
        break
      } else {
        lp_old = lp_up
      }
    }
  }
}


save(file = "../results/lp_posteriordb_LBFGS_h6.RData",
     list = c("lp_LBFGS_n_fn", "lp_LBFGS_n_gr", "initial_ls", "lp_INV",
              "lp_mean", "model_record"))

# check the output #
load(file = "../results/lp_posteriordb_LBFGS_h6.RData")

## check the distribution of number of iterations ##
n_grfn <- colMeans(abs(lp_LBFGS_n_gr))
mean(n_grfn) # 41.67  
median(n_grfn) # 31.05
sd(n_grfn)   # 53.80
sd(abs(lp_LBFGS_n_gr)) # 56.02
jpeg(filename = paste0("../pics/hist_LBFGS_grfn_counts.jpeg"),
     width = width, height = height, units = "px", pointsize = 12)
hist(abs(lp_LBFGS_n_gr), breaks = 100, 
     main = "", ylab = "", axes = TRUE,
     xlab = "iterations")
dev.off()

table(model_record[as.integer(which(abs(lp_LBFGS_n_gr)>200) / M - 0.5 / M) + 1])
#6 27 41 
#20  1  3

# boxplot of counts 
df <- data.frame(n_counts = c(abs(lp_LBFGS_n_gr)),
                 model = rep(pn[model_record], each = M),
                 not_reach_target = 
                   rep(apply(lp_LBFGS_n_gr, 2, 
                             f <- function(x){as.numeric(any(x < 0))}), 
                       each = M))

jpeg(filename = paste0("../pics/box_LBFGS_gr_counts_log.jpeg"),
     width = width*1.3, height = height*2, units = "px", pointsize = 12)
p_box_iter <- df %>% mutate( type=ifelse(not_reach_target == 1, "Highlighted","Normal")) %>%
  ggplot( aes(y = model, 
              x = n_counts, fill=type, alpha=type)) + 
  geom_boxplot() + 
  scale_x_log10() + ylab("") + xlab("No. calls to fn and gr") + 
  theme_grey(base_size = 26)  +
  scale_fill_manual(values=c("red", "white")) +
  theme(legend.position = "none") +
  scale_alpha_manual(values=c(1,0.1))
print(p_box_iter)
dev.off()

jpeg(filename = paste0("../pics/hist_LBFGS_grfn_counts_log.jpeg"),
     width = width, height = height, units = "px", pointsize = 12)
p_hist <- ggplot(df, aes(x = n_counts)) + 
  geom_histogram(bins = 100, color = "black", fill = "gray") +
  theme_bw(base_size = 26) +
  scale_x_log10(breaks=c(1, 10, 100, 300), labels = c("1", "10", "100", "300"))+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), axis.title.y = element_blank())+
  xlab("No. calls to fn and gr") 
print(p_hist)
dev.off()
