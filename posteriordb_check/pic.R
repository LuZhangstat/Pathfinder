## plots: Stan Phase I vs L-BFGS ##
load("../results/lp_posteriordb_LBFGS.RData")
load("../results/lp_posteriordb_explore.RData")

df <- data.frame(n_counts = c(c(abs(lp_LBFGS_n_gr)), c(lp_explore_n_leapfrog)),
                 model = rep(rep(pn[model_record], each = M), 2),
                 n_leapfrogs = rep(c(lp_explore_n_leapfrog), 2),
                 not_reach_target = 
                   c(rep(apply(lp_LBFGS_n_gr, 2, 
                               f <- function(x){as.numeric(any(x < 0))}), 
                         each = M), rep(2, M*length(model_record))))

jpeg(filename = paste0("../pics/box_compar_LBFGS_log.jpeg"),
     width = 860*2, height = 740*2, units = "px", pointsize = 12)
p_box_compar <- df %>% mutate(type= c("L-BFGS", "L-BFGS(multimodal)", "Stan Phase I ")
                              [not_reach_target + 1]) %>%
  ggplot(aes(y =  reorder(model, n_leapfrogs, FUN = median), 
             x = n_counts, fill = type, color = type)) + 
  geom_boxplot() + 
  scale_fill_manual(values=c("lightseagreen", "pink", "lightskyblue")) + 
  scale_color_manual(values=c("darkgreen", "red", "blue")) + 
  scale_x_log10(breaks=c(10, 1e2, 1e3, 1e4, 1e5), 
                labels = c("10", "100", "1000", "10,000", "100,000")) + 
  ylab("") + xlab("calls to lp__ + gradient") + 
  theme_bw(base_size = 26)  #+
#theme(legend.position = "none") 
print(p_box_compar)
dev.off()


## function for extract optims and inits ##
get_init_optim <- function(ind){
  lp_ind = ncol(lp_opath[[ind]]$opath[[1]]$y)
  inits <- c()
  optims <- c()
  for(l in 1:length(lp_opath[[ind]]$opath)){
    inits = rbind(inits, lp_opath[[ind]]$opath[[l]]$y[1, 1:(lp_ind - 1)])
    last_ind <- nrow(lp_opath[[ind]]$opath[[l]]$y)
    optims = rbind(optims, 
                   lp_opath[[ind]]$opath[[l]]$y[last_ind, 1:(lp_ind - 1)])
  }
  return(list(inits = inits, optims = optims))
}


## plots for section 3.2 ##
# 8 school centered #
get_opt_tr <- function(opath){
  
  ###
  #' function for retreveing optimization trajectories
  #' 

  lp_ind = ncol(opath[[1]]$y)
  opt_tr <- c()
  ind_tr <- c()
  tr_id <- c()
  for(l in 1:length(opath)){
    opt_tr = rbind(opt_tr, opath[[l]]$y[, 1:(lp_ind - 1)])
    ind_tr = c(ind_tr, 
               1:nrow(opath[[l]]$y[, 1:(lp_ind - 1)]))
    tr_id = c(tr_id, 
               rep(l, nrow(opath[[l]]$y[, 1:(lp_ind - 1)])))
    
  }
  return(list(opt_tr = opt_tr, ind_tr = ind_tr, tr_id = tr_id))
}

opt_tr_res <- get_opt_tr(opath)#opath #lp_opath[[15]]$opath
check_dim <- c(8, 10)

## run wasserstain_check with i = 15 and generate the plots
# then
## run main_pf_variational.R with init_bound = 15.0 and generate the plots
dta_sam <- data.frame(sam_x = pick_samples[check_dim[1], ],
                      sam_y = pick_samples[check_dim[2], ])

dta_check <- data.frame(x = ref_samples[, check_dim[1]],
                        y = ref_samples[, check_dim[2]])

dta_opt <- data.frame(
  optim_x = opt_tr_res$opt_tr[, check_dim[1]],
  optim_y = opt_tr_res$opt_tr[, check_dim[2]],
  optim_ind = opt_tr_res$ind_tr,
  tr_id = opt_tr_res$tr_id
)
dta_opt$tr_id = factor(dta_opt$tr_id)

p_check <- ggplot(dta_check, aes(x=x, y=y) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=-1) +
  scale_x_continuous(expand = c(0, 0), limits = c(-18, 22)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-20, 4)) +
  geom_point(data = dta_sam, aes(x=sam_x, y=sam_y), colour="red", size = 3) +
  xlab("") + ylab("") +
  theme(
    legend.position='none'
  )

p_check
jpeg(filename = paste0("../pics/phI_adapt_setting22/8-school_points.jpeg"),
     width = 500, height = 400, #740, 
     units = "px", pointsize = 12)
print(p_check)
dev.off()

p_check1 <- ggplot(dta_check, aes(x=x, y=y) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=-1) +
  scale_x_continuous(expand = c(0, 0), limits = c(-18, 22)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-20, 4)) +
  geom_point(data = dta_opt, 
            aes(x = optim_x, y = optim_y, group = tr_id, 
                color = optim_ind, alpha = 0.4), size = 1) +
  geom_line(data = dta_opt, 
            aes(x = optim_x, y = optim_y, group = tr_id, 
                color = optim_ind, alpha = 0.4)) +
  scale_color_gradient(low="white", high="orange") +
  xlab("") + ylab("") +
  theme(
    legend.position='none'
  ) 
p_check1

jpeg(filename = paste0("../pics/phI_adapt_setting22/8-school_opt_tr.jpeg"),
     width = 500, height = 400, #740, 
     units = "px", pointsize = 12)
print(p_check1)
dev.off()

# multimodality #
## run wasserstain_check with i = 3 and generate the plots
opt_tr_res <- get_opt_tr(lp_opath[[3]]$opath)
check_dim <- c(4, 6)  #24
dta_sam <- data.frame(sam_x = pick_samples[check_dim[1], ],
                      sam_y = pick_samples[check_dim[2], ])

dta_check <- data.frame(x = ref_samples[, check_dim[1]],
                        y = ref_samples[, check_dim[2]])



dta_phI <- data.frame(x = PhaseI_last_draw[[3]][, check_dim[1]],
                      y = PhaseI_last_draw[[3]][, check_dim[2]])

dta_opt <- data.frame(
  optim_x = opt_tr_res$opt_tr[, check_dim[1]],
  optim_y = opt_tr_res$opt_tr[, check_dim[2]],
  optim_ind = opt_tr_res$ind_tr,
  tr_id = opt_tr_res$tr_id
)
dta_opt$tr_id = factor(dta_opt$tr_id)

p_check <- ggplot(dta_check, aes(x=x, y=y) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=-1) +
  scale_x_continuous(expand = c(0, 0), limits = c(-1.5, 2.5)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-5, 3)) +
  geom_point(data = dta_sam, aes(x=sam_x, y=sam_y), colour="red",
             size = 3) +
  xlab("") + ylab("") +
  theme(
    legend.position='none'
  )

p_check
jpeg(filename = paste0("../pics/phI_adapt_setting22/3_points.jpeg"),
     width = 400, height = 400, #740, 
     units = "px", pointsize = 12)
print(p_check)
dev.off()

p_check1 <- ggplot(dta_check, aes(x=x, y=y) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=-1) +
  scale_x_continuous(expand = c(0, 0), limits = c(-1.5, 2.5)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-5, 3)) +
  geom_point(data = dta_opt, 
             aes(x = optim_x, y = optim_y, group = tr_id, 
                 color = optim_ind, alpha = 0.4), size = 1) +
  geom_line(data = dta_opt, 
            aes(x = optim_x, y = optim_y, group = tr_id, 
                color = optim_ind, alpha = 0.4)) +
  scale_color_gradient(low="white", high="orange") +
  xlab("") + ylab("") +
  theme(
    legend.position='none'
  ) 
p_check1

jpeg(filename = paste0("../pics/phI_adapt_setting22/3_opt_tr.jpeg"),
     width = 400, height = 400, #740, 
     units = "px", pointsize = 12)
print(p_check1)
dev.off()

p_phI <- ggplot(dta_check, aes(x=x, y=y) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=-1) +
  scale_x_continuous(expand = c(0, 0), limits = c(-1.5, 2.5)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-5, 3)) +
  geom_point(data = dta_phI, aes(x=x, y=y), colour="red", size = 3) +
  xlab("") + ylab("") +
  theme(
    legend.position='none'
  )

p_phI
jpeg(filename = paste0("../pics/phI_adapt_setting22/3_phI.jpeg"),
     width = 400, height = 400, #740, 
     units = "px", pointsize = 12)
print(p_phI)
dev.off()

# varying curvature #
## run wasserstain_check with i = 32 and generate the plots
opt_tr_res <- get_opt_tr(lp_opath[[32]]$opath)
check_dim <- c(45, 46)  #24
dta_sam <- data.frame(sam_x = pick_samples[check_dim[1], ],
                      sam_y = pick_samples[check_dim[2], ])

dta_check <- data.frame(x = ref_samples[, check_dim[1]],
                        y = ref_samples[, check_dim[2]])



dta_phI <- data.frame(x = PhaseI_last_draw[[32]][, check_dim[1]],
                      y = PhaseI_last_draw[[32]][, check_dim[2]])

dta_opt <- data.frame(
  optim_x = opt_tr_res$opt_tr[, check_dim[1]],
  optim_y = opt_tr_res$opt_tr[, check_dim[2]],
  optim_ind = opt_tr_res$ind_tr,
  tr_id = opt_tr_res$tr_id
)
dta_opt$tr_id = factor(dta_opt$tr_id)

p_check <- ggplot(dta_check, aes(x=x, y=y) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=-1) +
  scale_x_continuous(expand = c(0, 0), limits = c(-2.2, 4.6)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-7, 1)) +
  geom_point(data = dta_sam, aes(x=sam_x, y=sam_y), colour="red",
             size = 3) +
  xlab("") + ylab("") +
  theme(
    legend.position='none'
  )

p_check
jpeg(filename = paste0("../pics/phI_adapt_setting22/32_points.jpeg"),
     width = 400, height = 400, #740, 
     units = "px", pointsize = 12)
print(p_check)
dev.off()

p_check1 <- ggplot(dta_check, aes(x=x, y=y) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=-1) +
  scale_x_continuous(expand = c(0, 0), limits = c(-2.2, 4.6)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-7, 1)) +
  geom_point(data = dta_opt, 
             aes(x = optim_x, y = optim_y, group = tr_id, 
                 color = optim_ind, alpha = 0.4), size = 1) +
  geom_line(data = dta_opt, 
            aes(x = optim_x, y = optim_y, group = tr_id, 
                color = optim_ind, alpha = 0.4)) +
  scale_color_gradient(low="white", high="orange") +
  xlab("") + ylab("") +
  theme(
    legend.position='none'
  ) 
p_check1

jpeg(filename = paste0("../pics/phI_adapt_setting22/32_opt_tr.jpeg"),
     width = 400, height = 400, #740, 
     units = "px", pointsize = 12)
print(p_check1)
dev.off()

p_phI <- ggplot(dta_check, aes(x=x, y=y) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=-1) +
  scale_x_continuous(expand = c(0, 0), limits = c(-2.2, 4.6)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-7, 1)) +
  geom_point(data = dta_phI, aes(x=x, y=y), colour="red", size = 3) +
  xlab("") + ylab("") +
  theme(
    legend.position='none'
  )

p_phI
jpeg(filename = paste0("../pics/phI_adapt_setting22/32_phI.jpeg"),
     width = 400, height = 400, #740, 
     units = "px", pointsize = 12)
print(p_phI)
dev.off()

## plots for section 2.4 ##
## run wasserstain_check with i = 1 and generate the plots
# 0.1806857 vs  0.1276343
pick_samples_center <- 
  sapply(lp_opath[[i]]$opath[lp_opath[[i]]$pick_mode], extract_samples_center)
a_center = wpp(t(pick_samples_center), 
               mass = rep(1 / ncol(pick_samples_center), ncol(pick_samples_center)))
b = wpp(ref_samples, mass = rep(1 / nrow(ref_samples), nrow(ref_samples)))
w_d_c <- wasserstein(a_center, b, p = 2); w_d_c #0.191

pick_samples <- lp_opath[[i]]$pick_samples
a = wpp(t(pick_samples), 
        mass = rep(1 / ncol(pick_samples), ncol(pick_samples)))
b = wpp(ref_samples, mass = rep(1 / nrow(ref_samples), nrow(ref_samples)))
w_d_pf <- wasserstein(a, b, p = 2); w_d_pf # 0.122

check_dim <- c(1, 2)
dta_sam <- data.frame(sam_x = pick_samples[check_dim[1], ],
                      sam_y = pick_samples[check_dim[2], ])
dta_sam_center <- data.frame(sam_x = pick_samples_center[check_dim[1], ],
                             sam_y = pick_samples_center[check_dim[2], ])

dta_check <- data.frame(x = ref_samples[, check_dim[1]],
                        y = ref_samples[, check_dim[2]])


p_check <- ggplot(dta_check, aes(x=x, y=y) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=-1) +
  scale_x_continuous(expand = c(0, 0)) + #, limits = c(-18, 22)) +
  scale_y_continuous(expand = c(0, 0)) + #, limits = c(-20, 4)) +
  geom_point(data = dta_sam, aes(x=sam_x, y=sam_y), colour="red", size = 3) +
  xlab("") + ylab("") +
  theme(
    legend.position='none'
  )

p_check
jpeg(filename = paste0("../pics/phI_adapt_setting22/1-ark_pf.jpeg"),
     width = 400, height = 400, #740, 
     units = "px", pointsize = 12)
print(p_check)
dev.off()


p_check_c <- ggplot(dta_check, aes(x=x, y=y) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=-1) +
  scale_x_continuous(expand = c(0, 0)) + #, limits = c(-18, 22)) +
  scale_y_continuous(expand = c(0, 0)) + #, limits = c(-20, 4)) +
  geom_point(data = dta_sam_center, aes(x=sam_x, y=sam_y),
             colour="red", size = 3) +
  xlab("") + ylab("") +
  theme(
    legend.position='none'
  )

p_check_c

jpeg(filename = paste0("../pics/phI_adapt_setting22/1-ark_c.jpeg"),
     width = 400, height = 400, #740, 
     units = "px", pointsize = 12)
print(p_check_c)
dev.off()





