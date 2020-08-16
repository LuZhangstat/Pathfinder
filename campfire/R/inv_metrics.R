getHessian = function(fit, q) {
  Aqx = function(fit, q, r) {
    dx = 1e-5
    dr = dx * r
    (grad_log_prob(fit, q + dr / 2, adjust_transform = FALSE) -
        grad_log_prob(fit, q - dr / 2, adjust_transform = FALSE)) / dx
  }

  N = length(q)
  A = matrix(0, nrow = N, ncol = N)
  for(i in 1:N) {
    r = rep(0, N)
    r[i] = 1
    A[, i] = Aqx(fit, q, r)
  }
  0.5 * (A + t(A))
}

diag_inv_metric = function(samples) {
  diag(diag(cov(samples)))
}

dense_inv_metric = function(samples, rank_check = TRUE) {
  c = cov(samples)

  e = eigen(c, T)
  nkeep = tail(which(e$values > 1e-10), 1)

  if(nkeep < ncol(samples)) {
    mine = e$values[nkeep]
    c = e$vectors[, 1:nkeep] %*% diag(e$values[1:nkeep] - mine) %*% t(e$vectors[, 1:nkeep])
    c = c + mine * diag(ncol(samples))
  }

  return(c)
}

lw_linear_corr_inv_metric = function(samples) {
  sqrt_D = diag(sqrt(diag(cov(samples))))
  sqrt_Dinv = diag(1 / diag(sqrt_D))
  sqrt_D %*% linshrink_cov(samples %*% sqrt_Dinv) %*% sqrt_D
}

hess_inv_metric = function(stan_fit, Nev, samples) {
  Dsqrt = diag(sqrt(diag(cov(samples))))
  H = Dsqrt %*% getHessian(stan_fit, tail(samples, 1)) %*% Dsqrt
  eh = eigen(H, T)
  sorted = order(eh$values)
  evals = 1.0 / abs(eh$values[sorted])
  evecs = eh$vectors[, sorted]

  etail = evals[Nev + 1]

  Happrox = evecs[, 1:Nev] %*%
    (diag(evals[1:Nev], nrow = Nev) - etail * diag(1.0, nrow = Nev)) %*%
    t(evecs[, 1:Nev]) +
    etail * diag(ncol(samples))

  return(Dsqrt %*% Happrox %*% Dsqrt)
}

hess_wishart_inv_metric = function(stan_fit, Nev, samples) {
  h = hess_inv_metric(stan_fit, Nev, samples)
  c = cov(samples)

  P = ncol(samples)
  N = nrow(samples)

  return((P * h + N * c) / (P + N - 1))
}

compute_inv_metric = function(stan_fit, usamples) {
  Ntest = max(nrow(usamples) / 2, 50)
  Ntrain = nrow(usamples) - Ntest

  Ytrain = head(usamples, Ntrain)
  top_evec = eigen(cov(Ytrain), T)$vectors[, 1]
  Ytest = tail(usamples, Ntest)

  inv_metric_options = list(diag = diag_inv_metric,
                            dense = dense_inv_metric,
                            lw2004 = lw_linear_corr_inv_metric)

  if(ncol(usamples) > 1) {
    inv_metric_options[["hess1"]] = function(samples) hess_inv_metric(stan_fit, 1, samples)
    inv_metric_options[["hess1_wish"]] = function(samples) hess_wishart_inv_metric(stan_fit, 1, samples)
  }

  if(ncol(usamples) > 2) {
    inv_metric_options[["hess2"]] = function(samples) hess_inv_metric(stan_fit, 2, samples)
    inv_metric_options[["hess2_wish"]] = function(samples) hess_wishart_inv_metric(stan_fit, 2, samples)
  }

  perfdf = lapply(names(inv_metric_options), function(inv_metric_name) {
    inv_metric = inv_metric_options[[inv_metric_name]](Ytrain)

    cov_test = cov(Ytest)
    L = t(chol(inv_metric))
    el = eigen(solve(L, t(solve(L, t(cov_test)))), T)
    H = t(L) %*% getHessian(stan_fit, tail(Ytest, 1)) %*% L
    eh = eigen(H, T)

    tibble(name = inv_metric_name,
           c_hybrid = sqrt(max(abs(eh$values)) * max(abs(el$values))))
  }) %>%
    bind_rows %>%
    arrange(-c_hybrid)

  print("Metric calculation info (lower c better, last is always picked, minimum is 1.0):")
  print(perfdf)

  name = perfdf %>% tail(1) %>% pull(name)

  inv_metric = inv_metric_options[[name]](usamples)

  return(inv_metric)
}
