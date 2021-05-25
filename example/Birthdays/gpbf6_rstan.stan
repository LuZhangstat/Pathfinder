functions {
vector diagSPD_EQ(real alpha, real rho, real L, int M) {
  vector[M] one_to_M;
  for (m in 1:M) one_to_M[m] = m^2;
  return sqrt((alpha^2) * sqrt(2*pi()) * rho * exp(-0.5*(rho*pi()/2/L)^2 * one_to_M));
}
/* real spd_Matt(real alpha, real rho, real w) { */
/*   real S = 4*alpha^2 * (sqrt(3)/rho)^3 * 1/((sqrt(3)/rho)^2 + w^2)^2; */
/*   return sqrt(S); */
/* } */
vector diagSPD_periodic(real alpha, real rho, int M) {
  real a = 1/rho^2;
  vector[M] q;
  for (m in 1:M) q[m] = sqrt(alpha^2 * 2 / exp(a) * modified_bessel_first_kind(m, a));
  return append_row(q,q);
}
matrix PHI_EQ(int N, int M, real L, vector x) {
  vector[M] one_to_M;
  for (m in 1:M) one_to_M[m] = m;
  return sin(diag_post_multiply(rep_matrix(pi()/(2*L) * (x+L), M), one_to_M))/sqrt(L);
}
matrix PHI_periodic(int N, int M, real w0, vector x) {
  vector[M] one_to_M;
  matrix[N,M] mw0x;
  for (m in 1:M) one_to_M[m] = m;
  mw0x = diag_post_multiply(rep_matrix(w0*x, M), one_to_M);
  return append_col(cos(mw0x), sin(mw0x));
}
}
data {
  int<lower=1> N;      // number of observations
  vector[N] x;         // univariate covariate
  vector[N] y;         // target variable
  int day_of_week[N];  // 
  int day_of_year[N];  // 
        
  real<lower=0> c_f1;  // factor c to determine the boundary value L
  int<lower=1> M_f1;   // number of basis functions for smooth function
  int<lower=1> J_f2;   // number of cos and sin functions for periodic
}
transformed data {
  // Normalize data
  real xmean = mean(x);
  real ymean = mean(y);
  real xsd = sd(x);
  real ysd = sd(y);
  vector[N] xn = (x - xmean)/xsd;
  vector[N] yn = (y - ymean)/ysd;
  // Basis functions for f1
  real L_f1 = c_f1*max(xn);
  matrix[N,M_f1] PHI_f1 = PHI_EQ(N, M_f1, L_f1, xn);
  // Basis functions for f2
  real period_year = 365.25/xsd;
  matrix[N,2*J_f2] PHI_f2 = PHI_periodic(N, J_f2, 2*pi()/period_year, xn);
  // Concatenated basis functions for f1 and f2
  matrix[N,M_f1+2*J_f2] PHI_f = append_col(PHI_f1, PHI_f2);
}
parameters {
  real intercept0;
  vector[M_f1] beta_f1;         // the basis functions coefficients for f1
  vector[2*J_f2] beta_f2;       // the basis functions coefficients for f2
  vector[6] beta_f3;            // day of week effect
  vector[366] beta_f4;          // day of year effect
  real<lower=0> lengthscale_f1; //
  real<lower=0> lengthscale_f2; //
  real<lower=0> sigma_f1;       // scale of f1
  real<lower=0> sigma_f2;       // scale of f2
  real<lower=0> sigma_f4;       // scale of day of year effect
  real<lower=0> sigma;          // residual scale
}
model {
  // spectral densities for f1 and f2
  vector[M_f1] diagSPD_f1 = diagSPD_EQ(sigma_f1, lengthscale_f1, L_f1, M_f1);
  vector[2*J_f2] diagSPD_f2 = diagSPD_periodic(sigma_f2, lengthscale_f2, J_f2);
  // day of week and day of year effects
  vector[7] f_day_of_week = append_row(0, beta_f3);
  vector[N] intercept = intercept0 + f_day_of_week[day_of_week] + beta_f4[day_of_year];
  // priors
  intercept0 ~ normal(0, 1);
  beta_f1 ~ normal(0, 1);
  beta_f2 ~ normal(0, 1);
  beta_f3 ~ normal(0, 1);
  beta_f4 ~ normal(0, sigma_f4);
  lengthscale_f1 ~ lognormal(log(700/xsd), 1);
  lengthscale_f2 ~ normal(0, .1);
  sigma_f1 ~ normal(0, 1);
  sigma_f2 ~ normal(0, 1);
  sigma_f4 ~ normal(0, 0.1);
  sigma ~ normal(0, 0.5);
  // model
  yn ~ normal_id_glm(PHI_f,
		     intercept,
		     append_row(diagSPD_f1 .* beta_f1, diagSPD_f2 .* beta_f2),
		     sigma);
}

