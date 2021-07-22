parameters {
  real log_sigma;
  vector[99] alpha;
}
model {
  log_sigma ~ normal(0, 3);
  alpha ~ normal(0, exp(log_sigma / 2));
}
