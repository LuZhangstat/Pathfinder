library(transport)

load("../results/lp_posteriordb_phI_adapt_set22.RData") # Pathfinder
load("../results/lp_posteriordb_explore.RData") # phase I output
#load("../results/lp_posteriordb_LBFGS.RData")
load("../results/ADVI_results.RData") # ADVI results

## calculate wasserstain distance ##
w_d_matrix = matrix(NA, nrow = 49, ncol = 6)
t_0 <- proc.time()
for(i in 1:49){ #length(model_record)
  modelname <- pn[model_record[i]]
  printf("model %d: %s", model_record[i], modelname)
  
  # pick model
  po <- posterior(modelname, pdb = pd)
  # get reference posterior samples
  gsd <- reference_posterior_draws(po)
  
  # compile the model
  sc <- stan_code(po)
  model <- stan_model(model_code = sc)
  
  ###  get the data  ###
  data <- get_data(po)
  
  ### get reference samples ###
  posterior <- to_posterior(model, data)
  D <- get_num_upars(posterior)
  if(modelname == "eight_schools-eight_schools_noncentered"){
    constrained_draws <- modify_8school_noncen(po)
    unconstrained_draws <-  lapply(constrained_draws, unconstrain_draws, posterior)
  } else if (modelname == "gp_pois_regr-gp_pois_regr") {
    constrained_draws <- modify_draws_gp_pois_regr(po)
    unconstrained_draws <-  lapply(constrained_draws, unconstrain_draws, posterior)
  } else {
    unconstrained_draws <-  lapply(gsd, unconstrain_draws, posterior)
  }
  ref_samples = rbind(unconstrained_draws[[1]], unconstrained_draws[[2]], 
                      unconstrained_draws[[3]], unconstrained_draws[[4]],
                      unconstrained_draws[[5]], unconstrained_draws[[6]],
                      unconstrained_draws[[7]], unconstrained_draws[[8]],
                      unconstrained_draws[[9]], unconstrained_draws[[10]])
  
  ### samples from pathfinder ###
  pick_samples <- lp_opath[[i]]$pick_samples
  
  ### inits and optims ###
  inits_optims <- get_init_optim(i)
  
  ### calculate wasserstein distance ###
  # pathfinder #
  if(ncol(pick_samples) == 1){
    a = wpp(rbind(t(pick_samples), t(pick_samples)), 
            mass = rep(1 / 2, 2))
  }else{
    a = wpp(t(pick_samples), 
            mass = rep(1 / ncol(pick_samples), ncol(pick_samples)))
  }
  b = wpp(ref_samples, mass = rep(1 / nrow(ref_samples), nrow(ref_samples)))
  w_d_pf <- wasserstein(a, b, p = 2); w_d_pf
  
  
  # optims #
  inits_optims <- get_init_optim(i)
  a_opt = wpp(inits_optims$optims, 
              mass = rep(1 / nrow(inits_optims$optims), nrow(inits_optims$optims)))
  w_d_opt <- wasserstein(a_opt, b, p = 2); w_d_opt
  
  # inits #
  a_init = wpp(inits_optims$inits, 
               mass = rep(1 / nrow(inits_optims$inits), nrow(inits_optims$inits)))
  w_d_init <- wasserstein(a_init, b, p = 2); w_d_init
  
  # last samples of phase I #
  a_phI = wpp(PhaseI_last_draw[[i]],
               mass = rep(1 / nrow(PhaseI_last_draw[[i]]),
                          nrow(PhaseI_last_draw[[i]])))
  w_d_phI <- wasserstein(a_phI, b, p = 2); w_d_phI
  
  # ADVI: meanfield #
  a_ADVI_mf = wpp(ADVI_meanfield_draw[[i]][1:20, ],
              mass = rep(1 / nrow(ADVI_meanfield_draw[[i]][1:20, ]),
                         nrow(ADVI_meanfield_draw[[i]][1:20, ])))
  w_d_ADVI_mf <- wasserstein(a_ADVI_mf, b, p = 2); w_d_ADVI_mf
  
  # ADVI: fullrank #
  a_ADVI_fr = wpp(ADVI_fullrank_draw[[i]][1:20, ],
                  mass = rep(1 / nrow(ADVI_fullrank_draw[[i]][1:20, ]),
                             nrow(ADVI_fullrank_draw[[i]][1:20, ])))
  w_d_ADVI_fr <- wasserstein(a_ADVI_fr, b, p = 2); w_d_ADVI_fr
  
  
  w_d_matrix[i, ] = c(w_d_pf, w_d_opt, w_d_init, w_d_phI, w_d_ADVI_mf, 
                      w_d_ADVI_fr)
  
  cat("pf:", w_d_pf, "\t", "optims: ", w_d_opt, "\t", 
      "inits: ", w_d_init, "phI:", w_d_phI, "\n", 
      "ADVI meanfield:", w_d_ADVI_mf, "\t", "ADVI fullrank:", 
      w_d_ADVI_fr, "\n")
  
}
proc.time() - t_0


which((w_d_matrix[, 2] - w_d_matrix[, 1]) < 0) # 32

which((w_d_matrix[, 3] - w_d_matrix[, 1]) < 0) # 15

# save(file = "../results/wasserstein_phI_adapt_set22.RData",
#      list = c("w_d_matrix"))

load("../results/wasserstein_phI_adapt_set22.RData")
#load("../results/wasserstein_phI_adapt_set22_center.RData")
# wasserstein_phI_adapt_set22.RData
colMeans(w_d_matrix)


w_d_summary <- w_d_matrix[, c(3, 1, 2, 4, 5, 6)] 

colnames(w_d_summary) <- c("random init", "pf", "max", "PhI", "meanfield", "fullrank")
rownames(w_d_summary) <- pn[model_record]

# check pathfinder vs phase I warmup
pf_vs_phI <- (w_d_summary[, 2] / w_d_summary[, 4])
summary(pf_vs_phI)
quantile(pf_vs_phI, c(0.05, 0.5, 0.95))
#write.csv(w_d_summary, file = "../results/wasserstein22.csv")
mean_ratio_orders = order(rowMeans(w_d_summary[, c(2, 3, 4)] / w_d_summary[, 3]))
group = as.integer(order(mean_ratio_orders)/17)
w_d_dat = data.frame(ratios = c(w_d_summary[, c(2, 3, 4)] / w_d_summary[, 3]),
                       label = c(#rep("random", length(model_record)),
                         rep("pathfinder", length(model_record)), 
                         rep("maxima", length(model_record)),
                         rep("warmup", length(model_record))), 
                       model = rep(pn[model_record], 3),
                     group = rep(group, 3))

w_d_dat$label <- as.factor(w_d_dat$label)
w_d_dat$group <- as.factor(w_d_dat$group)

p_w_d_0 <- w_d_dat %>% filter(group == 0) %>%
  ggplot(aes(y = reorder(model, ratios, FUN = mean), 
             x = ratios, color = label, shape = label)) + 
  geom_point(size = 3) + 
  ylab("") + 
  xlab("") + 
  xlim(c(0, 1.5)) +
  theme_bw(base_size = 26) +
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")

p_w_d_1 <- w_d_dat %>% filter(group == 1) %>%
  ggplot(aes(y = reorder(model, ratios, FUN = mean), 
             x = ratios, color = label, shape = label)) + 
  geom_point(size = 3) + 
  ylab("") + 
  xlab("") + #"ratio of Wasserstain distance") + 
  xlim(c(0, 1.5)) +
  theme_bw(base_size = 26) +
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")

p_w_d_2 <- w_d_dat %>% filter(group == 2) %>%
  ggplot(aes(y = reorder(model, ratios, FUN = mean), 
             x = ratios, color = label, shape = label)) + 
  geom_point(size = 3) + 
  ylab("") + 
  xlab("") + 
  #xlim(c(0, 2)) +
  theme_bw(base_size = 26) +
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")

library(gridExtra)
library(grid)

grid_arrange_shared_legend <-
  function(...,
           ncol = length(list(...)),
           nrow = 1,
           position = c("bottom", "right")) {
    
    plots <- list(...)
    position <- match.arg(position)
    g <-
      ggplotGrob(plots[[1]] + theme(legend.position = position,
                                    legend.title = element_blank()))$grobs
    legend <- g[[which(sapply(g, function(x)
      x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    lwidth <- sum(legend$width)
    gl <- lapply(plots, function(x)
      x + theme(legend.position = "none"))
    gl <- c(gl, ncol = ncol, nrow = nrow)
    
    combined <- switch(
      position,
      "bottom" = arrangeGrob(
        do.call(arrangeGrob, gl),
        legend,
        ncol = 1,
        heights = unit.c(unit(1, "npc") - lheight, lheight)
      ),
      "right" = arrangeGrob(
        do.call(arrangeGrob, gl),
        legend,
        ncol = 2,
        widths = unit.c(unit(1, "npc") - lwidth, lwidth)
      )
    )
    
    grid.newpage()
    grid.draw(combined)
    invisible(combined)
    
  }

jpeg(filename = paste0("../pics/phI_adapt_setting22/W_d_compar.jpeg"),
     width = 900, height = 360, #740, 
     units = "px", pointsize = 12)
print(grid_arrange_shared_legend(p_w_d_0, p_w_d_1, p_w_d_2, 
                                 position = "right"))
dev.off()


## compare pathfinder with ADVI ##
# ADVI meanfield 
mean_ratio_orders2 = order(apply(w_d_summary[, c(2, 5, 6)] / w_d_summary[, 5], 
                                 1, max))
mean_ratio_orders2 = order(rowMeans(w_d_summary[, c(2, 5, 6)] / w_d_summary[, 5]))
group2 = as.integer(order(mean_ratio_orders2) / 13)
w_d_dat2 = data.frame(ratios = c(w_d_summary[, c(2, 5, 6)] / w_d_summary[, 5]),
                     label = c(rep("pathfinder", length(model_record)), 
                               rep("meanfield ADVI", length(model_record)),
                               rep("fullrank ADVI", length(model_record))), 
                     model = rep(pn[model_record], 3),
                     group = rep(group2, 3))

w_d_dat2$label <- as.factor(w_d_dat2$label)
w_d_dat2$group <- as.factor(w_d_dat2$group)

p_w_d_0_2 <- w_d_dat2 %>% filter(group == 0) %>%
  ggplot(aes(y = reorder(model, ratios, FUN = max), 
             x = ratios, color = label, shape = label)) + 
  geom_point(size = 3) + 
  ylab("") + 
  xlab("") + 
  xlim(c(0, 2.0)) +
  theme_bw(base_size = 26) +
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")

p_w_d_1_2 <- w_d_dat2 %>% filter(group == 1) %>%
  ggplot(aes(y = reorder(model, ratios, FUN = max), 
             x = ratios, color = label, shape = label)) + 
  geom_point(size = 3) + 
  ylab("") + 
  xlab("") + #"ratio of Wasserstain distance") + 
  xlim(c(0, 2.0)) +
  theme_bw(base_size = 26) +
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")

p_w_d_2_2 <- w_d_dat2 %>% filter(group == 2) %>%
  ggplot(aes(y = reorder(model, ratios, FUN = max), 
             x = ratios, color = label, shape = label)) + 
  geom_point(size = 3) + 
  ylab("") + 
  xlab("") + 
  xlim(c(0, 2.0)) +
  theme_bw(base_size = 26) +
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")

p_w_d_2_3 <- w_d_dat2 %>% filter(group == 3) %>%
  ggplot(aes(y = reorder(model, ratios, FUN = max), 
             x = ratios, color = label, shape = label)) + 
  geom_point(size = 3) + 
  ylab("") + 
  xlab("") + 
  theme_bw(base_size = 26) +
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")


jpeg(filename = paste0("../pics/phI_adapt_setting22/W_d_pf_ADVI.jpeg"),
     width = 900, height = 300, #740, 
     units = "px", pointsize = 12)
print(grid_arrange_shared_legend(p_w_d_0_2, p_w_d_1_2, p_w_d_2_2, p_w_d_2_3,
                                 position = "right"))
dev.off()

