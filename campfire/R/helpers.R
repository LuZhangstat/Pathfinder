getUnconstrainedSamples = function(fit) {
  usamples_list = lapply(fit$diagnostic_files(), function(file) {
    read_csv(file, comment = "#", col_types = cols(.default = col_double())) %>%
      select(-lp__, -accept_stat__, -stepsize__, -treedepth__, -n_leapfrog__, -divergent__, -energy__,
             -starts_with("p_"), -starts_with("g_")) %>%
      as.matrix()
  })

  usamples = array(0, dim = c(nrow(usamples_list[[1]]),
                              length(usamples_list),
                              ncol(usamples_list[[1]])))

  for(i in 1:length(usamples_list)) {
    usamples[, i,] = usamples_list[[i]]
  }

  return(usamples)
}

getExtras = function(fit) {
  lapply(fit$diagnostic_files(), function(file) {
    read_csv(file, comment = "#", col_types = cols(.default = col_double())) %>%
      select(lp__, accept_stat__, stepsize__, treedepth__, n_leapfrog__, divergent__, energy__)
  })
}

getInitFile = function(stan_fit, ldraw) {
  init = constrain_pars(stan_fit, ldraw %>% as.matrix)
  init_file = tempfile("init", fileext = ".dat")
  stan_rdump(names(init), init_file, env = list2env(init))

  return(init_file)
}

getStepsizes = function(fit) {
  sapply(fit$diagnostic_files(), function(file) {
    read_csv(file, comment = "#", col_types = cols(.default = col_double())) %>%
      tail(1) %>%
      pull(stepsize__)
  })
}
