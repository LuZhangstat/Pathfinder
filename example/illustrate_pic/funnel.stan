parameters {
  real log_sigma;
  real alpha;
}
model {
  log_sigma ~ normal(0, 3);
  alpha ~ normal(0, exp(log_sigma / 2));
}
