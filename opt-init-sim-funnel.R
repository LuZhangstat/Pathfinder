library(rstan)

D <- 8

program_funnel <-
"
data {
  int D;
} parameters {
  vector[D] theta;
  real<lower = 0> sigma;
} model {
  theta ~ normal(0, sigma);
  sigma ~ normal(0, 1);
}
"

model_funnel <- stan_model(model_code = program_funnel)

for (log10_sigma in seq(-2, 2, 0.5)) {
  sigma <- 10^log10_sigma
  print(" ")
  print(c("log10 sigma = ", log10_sigma), quote = FALSE)
  for (k in 1:5) {
    init <- function(n) list(sigma = sigma,
                             theta = rnorm(D, 0, sigma / 4))
    fit <- sampling(model_funnel, data = list(D = D),
                    chains = 1, init = init, refresh = 0)
    buf <- extract(fit, pars = c("lp__"),
                   permute = FALSE, inc_warmup = TRUE)[1:20]
    print(buf, quote = FALSE, digits = 1)
  }
}
