rm(list = ls())
library(rstan)
# run the following line to install package posteriordb
# remotes::install_github("MansMeg/posteriordb", subdir = "rpackage")
library(posteriordb)
pd <- pdb_default()  # Posterior database connection
po <- posterior("gp_pois_regr-gp_regr", pdb = pd) # pick model: gp_pois_regr-gp_regr
sc <- stan_code(po) # access data and model
dat <- get_data(po) # Check data 
model <- stan_model(model_code = sc) # compile the model
posterior <- sampling(model, data = dat, chains = 1, iter = 1, refresh = 0,
                      algorithm = "Fixed_param")

# check option adjust_transform on two points: theta1 and theta2 #
theta1 <- c(6.0, 1.4, 1.3); names(theta1) <- c("rho", "alpha", "sigma")
theta2 <- c(5.0, 2.0, 1.0); names(theta2) <- c("rho", "alpha", "sigma")

fn <- function(theta){log_prob(posterior, unconstrain_pars(posterior, theta), 
                               adjust_transform = TRUE)}
fn2 <- function(theta){log_prob(posterior, unconstrain_pars(posterior, theta), 
                                adjust_transform = FALSE)}

fn(theta1); fn(theta2) 
#[1] -29.90545
#[1] -31.21401
fn2(theta1); fn2(theta2)
#[1] 0
#[1] 0

# check the corresponding results with option gradient = TRUE
gn <- function(theta){ log_prob(posterior, unconstrain_pars(posterior, theta), 
                                adjust_transform = TRUE, gradient = TRUE)}
gn2 <- function(theta){log_prob(posterior, unconstrain_pars(posterior, theta), 
                                adjust_transform = FALSE, gradient = TRUE)}

gn(theta1); gn(theta2) 
# [1] 2.861269
# attr(,"gradient")
# [1] 4.593643 6.009790 4.235825
# [1] 1.552708
# attr(,"gradient")
# [1] 8.182420 1.811042 7.207582
gn2(theta1); gn2(theta2)
# [1] 0.4706731
# attr(,"gradient")
# [1] 3.593643 5.009790 3.235825
# [1] -0.7498768
# attr(,"gradient")
# [1] 7.1824202 0.8110424 6.2075820

