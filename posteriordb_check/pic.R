setwd("./posteriordb")
library(ggplot2)
library(posteriordb)
library(RColorBrewer)
source("../utils/lp_utils.R")
pd <- pdb_local() # Posterior database connection
pn <- posterior_names(pd)


# colorblind-friendly palettes
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


load("../results/lp_posteriordb_LBFGS_h10.RData")
#load("../results/lp_posteriordb_explore_h10.RData")
load("../results/lp_posteriordb_phI_adapt_default.RData")

### Wasserstein distance check for 100 repeats ###
## PSIS one sample from each ##
M = 100
load("../results/wasserstein_phI_adapt_default.RData")
# load("../results/W_d_IR_4_30.RData")
# w_d_summary <- cbind(w_d_matrix[, c(3, 1, 2, 4, 5, 6, 7)], w_d_IR)
# colnames(w_d_summary) <- c("random init", "pf", "max", "PhI", "meanfield",
#                            "meanfield center", "fullrank", "pf_4_IR")
# rownames(w_d_summary) <- pn[model_record]

# check pathfinder vs phase I warmup
pf_vs_phI <- (w_d_matrix[, "pf"] / w_d_matrix[, "PhI"])
summary(pf_vs_phI)
quantile(pf_vs_phI, c(0.05, 0.5, 0.95))

pf_vs_phI <- (w_d_matrix[, "pf_4_IR"] / w_d_matrix[, "PhI"])
summary(pf_vs_phI)
quantile(pf_vs_phI, c(0.05, 0.5, 0.95))
table(pf_vs_phI<1.2)
max(pf_vs_phI)

## rank methods ##
rank_score <- apply(w_d_matrix[, c("pf", "pf_4_IR", "PhI", "meanfield", 
                                   "meanfield center", "fullrank")], 1, 
                    f <- function(x){order(order(x))})
table(rank_score[2, ]==1)
table(rank_score[2, ]==2)
table(rank_score[2, ]<4)
table(rank_score[1, ]==1)
table(rank_score[3, ]<3)
rowSums(rank_score)

# compare pathfinder with ADVI #
summary(w_d_matrix[, "pf"] / w_d_matrix[, "meanfield"])
summary(w_d_matrix[, "pf"] / w_d_matrix[, "fullrank"])
summary(w_d_matrix[, "pf"] / w_d_matrix[, "meanfield center"])
summary(w_d_matrix[, "pf_4_IR"] / w_d_matrix[, "meanfield"])
summary(w_d_matrix[, "pf_4_IR"] / w_d_matrix[, "fullrank"])
summary(w_d_matrix[, "pf_4_IR"] / w_d_matrix[, "meanfield center"])

## one plot ##
summary(w_d_matrix[, "pf_4_IR"] / w_d_matrix[, "PhI"])
ratios <- c(w_d_matrix[, c("pf", "pf_4_IR", "PhI", "meanfield", 
                           "meanfield center", "fullrank")] / 
              w_d_matrix[, "pf_4_IR"])
ratios[ratios > 2^10] <- 2^10
w_d_dat = data.frame(ratios = ratios,
                     label = c(
                       rep("pathfinder", length(model_record)), 
                       rep("pathfinder IR with 4 runs", length(model_record)),
                       rep("Stan Phase I", length(model_record)),
                       rep("mean-field ADVI", length(model_record)),
                       rep("mean-field center", length(model_record)),
                       rep("dense ADVI", length(model_record))),
                     model = rep(pn[model_record], 6))#6))

p_w_d_compar <- w_d_dat %>% 
  ggplot(aes(y = model, #reorder(model, ratios, FUN = mean), 
             x = ratios, color = label, shape = label)) + 
  geom_point(size = 2) +
  scale_x_continuous(trans = 'log2',
                     limits = c(1/4, 1024),
                     breaks = c(1/4, 1/2, 1, 2, 4, 
                                8, 16, 32, 64, 128, 256, 512, 1024),
                     labels = c("1/4", 
                                "1/2", "1", "2", "4", "8", 
                                "16", "32", "64", "128", 
                                "256", "512", ">=1024")) + 
  ylab("") + xlab("scaled Wasserstein distance") + 
  theme_bw(base_size = 12)  +
  theme(legend.position="top", legend.title = element_blank()) #+
#theme(legend.position = "none") 

width <- 12.0
height <- 8.0
pointsize <- 16

setEPS()
postscript("../pics/phI_adapt_default/W_d_compar_100_each_2.eps",  #_resam_each
           width = width, height = height)
print(p_w_d_compar)
dev.off()

## 100 for each ##
load("../results/wasserstein_100_default.RData")

ratio_M_wd_mf <- (apply(W_d_100_ADVI_mf, 2, f <- function(x){quantile(x, 0.5)})/
                    apply(W_d_100_pf, 2, f <- function(x){quantile(x, 0.5)}))

mean(ratio_M_wd_mf)
table(ratio_M_wd_mf > 2)
table(ratio_M_wd_mf < 0.5)

mean_M_wd_mf <- colMeans(W_d_100_ADVI_mf) / colMeans(W_d_100_pf)
table(mean_M_wd_mf > 2)
table(mean_M_wd_mf < 0.5)

ratio_M_wd_fr <- (apply(W_d_100_ADVI_fr, 2, f <- function(x){quantile(x, 0.5)})/
                    apply(W_d_100_pf, 2, f <- function(x){quantile(x, 0.5)}))

table(ratio_M_wd_fr > 2)
table(ratio_M_wd_fr < 0.5)

mean_M_wd_fr <- colMeans(W_d_100_ADVI_fr) / colMeans(W_d_100_pf)
table(mean_M_wd_fr > 2)
table(mean_M_wd_fr < 0.5)

table((ratio_M_wd_fr > 2) & (ratio_M_wd_mf > 2))
table((ratio_M_wd_fr < 0.5) & (ratio_M_wd_mf < 0.5))

table((mean_M_wd_fr > 2) & (mean_M_wd_mf > 2))
table((mean_M_wd_fr < 0.5) & (mean_M_wd_mf < 0.5))


w_d_median_pf <- rep(apply(W_d_100_pf, 2, median), each = M)
apply(W_d_100_ADVI_mf / w_d_median_pf, 2, median)
apply(W_d_100_ADVI_fr / w_d_median_pf, 2, median)
w_d_scaled = c(c(W_d_100_pf / w_d_median_pf), 
               c(W_d_100_ADVI_mf / w_d_median_pf), 
               c(W_d_100_ADVI_fr / w_d_median_pf))
w_d_scaled[w_d_scaled >= 2^12] <- 2^12
df <- data.frame(w_d = w_d_scaled,
                 model = rep(rep(pn[model_record], each = M), 3),
                 type = rep(c("pathfinder", "mean-field ADVI", "dense ADVI"), 
                            each = length(model_record)*M))


# w_d_PhI_scaled <- c(w_d_summary[, "PhI"] / apply(W_d_100_pf, 2, median))
# w_d_point = data.frame(w_d_PhI_scaled = w_d_PhI_scaled,
#                        model = pn[model_record], 
#                        point = rep("Stan Phase I", length(model_record)))

width <- 12.0
height <- 8.0
setEPS()
postscript("../pics/phI_adapt_default/W_d_box_compar_100_each.eps",  #_resam_each
           width = width, height = height)
p_box_compar <- df %>%
  ggplot(aes(y = model, #reorder(model, n_leapfrogs, FUN = median), 
             x = w_d, color = type)) + 
  geom_boxplot(outlier.size = 0.1) + 
  # geom_point(aes(y = model, x = w_d_PhI_scaled, shape = point), 
  #            data = w_d_point, size = 2, inherit.aes=FALSE) +
  # scale_shape_manual(values=c(23)) + 
  #scale_fill_manual(values = cbbPalette) +
  scale_colour_manual(values=cbbPalette) + 
  scale_x_continuous(trans = 'log2',
                     limits = c(1/32, 2^12),
                     breaks = c(1/32, 1/16, 1/8, 1/4, 1/2, 1, 2, 4, 
                                8, 16, 32, 64, 128, 256, 512, 1024,
                                2^11, 2^12),
                     labels = c("1/32", "1/16", "1/8", "1/4", 
                                "1/2", "1", "2", "4", "8", 
                                "16", "32", "64", "128", 
                                "256", "512", "1024", "2048", ">4096")) + 
  ylab("") + xlab("") + #xlab("calls to log density and gradient") + 
  theme_bw(base_size = 12 )+
  theme(legend.position="top", legend.title = element_blank()) 
print(p_box_compar)
dev.off()


## computational cost comparision ##
#load("../results/lp_posteriordb_phI_adapt_set30.RData") # Pathfinder #_resam_all #_resam_each
load("../results/lp_posteriordb_phI_adapt_default.RData") # Pathfinder #_resam_all #_resam_each
load("../results/PhI_100_h10.RData")
load("../results/ADVI_100.RData")


pathfinder_fn_call <- 
  sapply(lp_opath, f <- function(x){
    sapply(x$opath, g <- function(z){sum(z$fn_call)}) } )
pathfinder_gr_call <- 
  sapply(lp_opath, f <- function(x){
    sapply(x$opath, g <- function(z){sum(z$gr_call)}) } )

summary(colSums(PhI_leapfrog_counts) / colSums(pathfinder_fn_call))
summary(colSums(PhI_leapfrog_counts) / colSums(pathfinder_gr_call))
table((colSums(PhI_leapfrog_counts) / colSums(pathfinder_gr_call)) > 100)
summary(colSums(calls_lp_mean) / colSums(pathfinder_fn_call))
summary(colSums(calls_gr_mean) / colSums(pathfinder_gr_call))
summary(colSums(calls_lp_full) / colSums(pathfinder_fn_call))
summary(colSums(calls_gr_full) / colSums(pathfinder_gr_call))

summary(calls_lp_mean / colMeans(pathfinder_fn_call))
summary(calls_gr_mean / colMeans(pathfinder_gr_call))
summary(calls_lp_full / colMeans(pathfinder_fn_call))
summary(calls_gr_full / colMeans(pathfinder_gr_call))


df <- data.frame(n_counts = c(c(pathfinder_fn_call), 
                              c(PhI_leapfrog_counts),
                              c(calls_lp_mean), c(calls_lp_full)),
                 model = rep(rep(pn[model_record], each = M), 4),
                 type = rep(c("Pathfinder", "Stan Phase I", 
                              "mean-field ADVI", "dense ADVI"), 
                            each = M * length(model_record)))

p_lp_compar <- df %>%
  ggplot(aes(y = model, #reorder(model, n_leapfrogs, FUN = median), 
             x = n_counts, color = type)) + 
  geom_boxplot(outlier.size = 0.1) +
  scale_colour_manual(values = cbbPalette) +
  scale_x_log10(breaks=c(10, 1e2, 1e3, 1e4, 5e4), 
                labels = c("10", "100", "1000", "10,000", "50,000")) + 
  ylab("") + xlab("") + #xlab("calls to lp__") + 
  theme_bw(base_size = 12)+
  theme(legend.position="top", legend.title = element_blank())   

width <- 12.0
height <- 8.0

setEPS()
postscript("../pics/phI_adapt_default/cost_lp_100_each.eps",  #_resam_each
           width = width, height = height)
print(p_lp_compar)
dev.off()

df <- data.frame(n_counts = c(c(pathfinder_gr_call), 
                              c(PhI_leapfrog_counts),
                              c(calls_gr_mean), c(calls_gr_full)),
                 model = rep(rep(pn[model_record], each = M), 4),
                 type = rep(c("Pathfinder", "Stan Phase I", 
                              "mean-field ADVI", "dense ADVI"), 
                            each = M * length(model_record)))

p_gr_compar <- df %>%
  ggplot(aes(y = model, #reorder(model, n_leapfrogs, FUN = median), 
             x = n_counts, color = type)) + 
  geom_boxplot(outlier.size = 0.1) +
  scale_colour_manual(values = cbbPalette) +
  scale_x_log10(breaks=c(10, 1e2, 1e3, 1e4, 5e4), 
                labels = c("10", "100", "1000", "10,000", "50,000")) + 
  ylab("") + xlab("") + #xlab("calls to lp__") + 
  theme_bw(base_size = 12)+
  theme(legend.position="top", legend.title = element_blank())   

width <- 12.0
height <- 8.0

setEPS()
postscript("../pics/phI_adapt_default/cost_gr_100_each.eps",  #_resam_each
           width = width, height = height)
print(p_gr_compar)
dev.off()

## sensitivity test ##
## 100 for each ##
load("../results/wasserstein_100_default.RData")
W_d_100_pf_default <- W_d_100_pf
load("../results/wasserstein_100_short_L.RData")
W_d_100_pf_short_L <- W_d_100_pf
load("../results/wasserstein_100_large_K.RData")
W_d_100_pf_large_K <- W_d_100_pf
load("../results/wasserstein_100_long_hist.RData")
W_d_100_pf_long_hist <- W_d_100_pf


## Number of monte carlo samples in ELBO estimation ##
mean(apply(W_d_100_pf_large_K, 2, f <- function(x){quantile(x, 0.5)}) / 
       apply(W_d_100_pf_default, 2, f <- function(x){quantile(x, 0.5)}))

min(apply(W_d_100_pf_large_K, 2, f <- function(x){quantile(x, 0.5)}) / 
      apply(W_d_100_pf_default, 2, f <- function(x){quantile(x, 0.5)}))

load("../results/lp_posteriordb_phI_adapt_large_K.RData") #
pathfinder_large_K_fn_call <- 
  sapply(lp_opath, f <- function(x){
    sapply(x$opath, g <- function(z){sum(z$fn_call)}) } )
pathfinder_large_K_gr_call <- 
  sapply(lp_opath, f <- function(x){
    sapply(x$opath, g <- function(z){sum(z$gr_call)}) } )
sum(pathfinder_large_K_fn_call)/sum(pathfinder_fn_call)
load("../results/lp_posteriordb_phI_adapt_default.RData") 


w_d_median_pf <- rep(apply(W_d_100_pf_default, 2, median), each = M)
apply(W_d_100_pf_large_K / w_d_median_pf, 2, median)
apply(W_d_100_pf_short_L / w_d_median_pf, 2, median)
w_d_scaled = c(#c(W_d_100_pf_28 / w_d_median_pf), 
  c(W_d_100_pf_large_K / w_d_median_pf), 
  c(W_d_100_pf_default / w_d_median_pf))
range(w_d_scaled)
w_d_scaled[w_d_scaled >= 2^12] <- 2^12
df <- data.frame(w_d = w_d_scaled,
                 model = rep(rep(pn[model_record], each = M), 2),
                 type = rep(c("Lmax = 1000, tol = 1e-13, K = 30, J = 6", 
                              "Lmax = 1000, tol = 1e-13, K = 5, J = 6"), 
                   each = length(model_record)*M))


width <- 12.0
height <- 8.0
setEPS()
postscript("../pics/phI_adapt_default/W_d_box_sens_K.eps",  #_resam_each
           width = width, height = height)
p_box_compar <- df %>%
  ggplot(aes(y = model, #reorder(model, n_leapfrogs, FUN = median), 
             x = w_d, color = type)) + 
  geom_boxplot(outlier.size = 0.1) +
  scale_colour_manual(values=cbbPalette) + 
  scale_x_continuous(trans = 'log2',
                     limits = c(1/16, 128),
                     breaks = c(1/16, 1/8, 1/4, 1/2, 1, 2, 4, 
                                8, 16, 32, 64, 128),
                     labels = c("1/16", "1/8", "1/4", 
                                "1/2", "1", "2", "4", "8", 
                                "16", "32", "64", "128")) + 
  ylab("") + xlab("") + #xlab("calls to log density and gradient") + 
  theme_bw(base_size = 12 )+
  theme(legend.position="top", legend.title = element_blank()) 
print(p_box_compar)
dev.off()

## length and convergence tolerance ##
w_d_scaled = c(c(W_d_100_pf_short_L / w_d_median_pf), 
               c(W_d_100_pf_default / w_d_median_pf))
range(w_d_scaled)
df <- data.frame(w_d = w_d_scaled,
                 model = rep(rep(pn[model_record], each = M), 2),
                 type = rep(c("Lmax = 200, tol = 1e-8, K = 5, J = 6",  
                              "Lmax = 1000, tol = 1e-13, K = 5, J = 6"), 
                            each = length(model_record)*M))


width <- 12.0
height <- 8.0
setEPS()
postscript("../pics/phI_adapt_default/W_d_box_sens_L.eps",  #_resam_each
           width = width, height = height)
p_box_compar <- df %>%
  ggplot(aes(y = model, #reorder(model, n_leapfrogs, FUN = median), 
             x = w_d, color = type)) + 
  geom_boxplot(outlier.size = 0.1) +
  scale_colour_manual(values=cbbPalette) + 
  scale_x_continuous(trans = 'log2',
                     limits = c(1/16, 128),
                     breaks = c(1/16, 1/8, 1/4, 1/2, 1, 2, 4, 
                                8, 16, 32, 64, 128),
                     labels = c("1/16", "1/8", "1/4", 
                                "1/2", "1", "2", "4", "8", 
                                "16", "32", "64", "128")) + 
  ylab("") + xlab("") + #xlab("calls to log density and gradient") + 
  theme_bw(base_size = 12 )+
  theme(legend.position="top", legend.title = element_blank()) 
print(p_box_compar)
dev.off()

## history size ##
mean((apply(W_d_100_pf_long_hist, 2, f <- function(x){quantile(x, 0.5)}) / 
        apply(W_d_100_pf_default, 2, f <- function(x){quantile(x, 0.5)}))[c(-5, -8)])
# 0.9763456
median((apply(W_d_100_pf_long_hist, 2, f <- function(x){quantile(x, 0.5)}) / 
          apply(W_d_100_pf_default, 2, f <- function(x){quantile(x, 0.5)}))[c(-5, -8)])
# 0.9560745
range((apply(W_d_100_pf_long_hist, 2, f <- function(x){quantile(x, 0.5)}) / 
         apply(W_d_100_pf_default, 2, f <- function(x){quantile(x, 0.5)}))[c(-5, -8)])
# 0.8065852 1.4026032

w_d_scaled = c(c(W_d_100_pf_long_hist / w_d_median_pf), 
               c(W_d_100_pf_default / w_d_median_pf))
range(w_d_scaled)
df <- data.frame(w_d = w_d_scaled,
                 model = rep(rep(pn[model_record], each = M), 2),
                 type = rep(c("Lmax = 1000, tol = 1e-13, K = 5, J = 60",  
                              " Lmax = 1000, tol = 1e-13, K = 5, J = 6"), 
                            each = length(model_record)*M))

width <- 12.0
height <- 8.0
setEPS()
postscript("../pics/phI_adapt_default/W_d_box_sens_J.eps",  #_resam_each
           width = width, height = height)
p_box_compar <- df %>%
  ggplot(aes(y = model, #reorder(model, n_leapfrogs, FUN = median), 
             x = w_d, color = type)) + 
  geom_boxplot(outlier.size = 0.1) +
  scale_colour_manual(values=cbbPalette) + 
  scale_x_continuous(trans = 'log2',
                     limits = c(1/16, 128),
                     breaks = c(1/16, 1/8, 1/4, 1/2, 1, 2, 4, 
                                8, 16, 32, 64, 128),
                     labels = c("1/16", "1/8", "1/4", 
                                "1/2", "1", "2", "4", "8", 
                                "16", "32", "64", "128")) + 
  ylab("") + xlab("") + #xlab("calls to log density and gradient") + 
  theme_bw(base_size = 12 )+
  theme(legend.position="top", legend.title = element_blank()) 
print(p_box_compar)
dev.off()


## computational cost comparision ##
load("../results/lp_posteriordb_phI_adapt_default.RData") # Pathfinder #_resam_all #_resam_each
pathfinder_fn_call_default <- 
  sapply(lp_opath, f <- function(x){
    sapply(x$opath, g <- function(z){sum(z$fn_call)}) } )
pathfinder_gr_call_default <- 
  sapply(lp_opath, f <- function(x){
    sapply(x$opath, g <- function(z){sum(z$gr_call)}) } )

load("../results/lp_posteriordb_phI_adapt_short_L.RData") # Pathfinder #_resam_all #_resam_each
pathfinder_fn_call_short_L <- 
  sapply(lp_opath, f <- function(x){
    sapply(x$opath, g <- function(z){sum(z$fn_call)}) } )
pathfinder_gr_call_short_L <- 
  sapply(lp_opath, f <- function(x){
    sapply(x$opath, g <- function(z){sum(z$gr_call)}) } )
load("../results/lp_posteriordb_phI_adapt_large_K.RData") # Pathfinder #_resam_all #_resam_each
pathfinder_fn_call_large_K <- 
  sapply(lp_opath, f <- function(x){
    sapply(x$opath, g <- function(z){sum(z$fn_call)}) } )
pathfinder_gr_call_large_K <- 
  sapply(lp_opath, f <- function(x){
    sapply(x$opath, g <- function(z){sum(z$gr_call)}) } )

df <- data.frame(n_counts = c(c(pathfinder_fn_call_short_L), 
                              c(pathfinder_fn_call_large_K),
                              c(pathfinder_fn_call_default)),
                 model = rep(rep(pn[model_record], each = M), 3),
                 type = rep(c("K = 5, Lmax = 200, tol = 1e-8", 
                              "K = 30, Lmax = 1000, tol = 1e-13", 
                              "K = 5, Lmax = 1000, tol = 1e-13"), 
                            each = M * length(model_record)))

p_lp_compar <- df %>%
  ggplot(aes(y = model, #reorder(model, n_leapfrogs, FUN = median), 
             x = n_counts, color = type)) + 
  geom_boxplot(outlier.size = 0.1) +
  scale_colour_manual(values = cbbPalette) +
  scale_x_log10(breaks=c(10, 1e2, 1e3, 1e4, 5e4), 
                labels = c("10", "100", "1000", "10,000", "50,000")) + 
  ylab("") + xlab("") + #xlab("calls to lp__") + 
  theme_bw(base_size = 12)+
  theme(legend.position="top", legend.title = element_blank())   

width <- 12.0
height <- 8.0

setEPS()
postscript("../pics/phI_adapt_default/cost_lp_sens.eps",  #_resam_each
           width = width, height = height)
print(p_lp_compar)
dev.off()

df <- data.frame(n_counts = c(c(pathfinder_gr_call_short_L), 
                              c(pathfinder_gr_call_default)),
                 model = rep(rep(pn[model_record], each = M), 2),
                 type = rep(c("K = 5, Lmax = 200, tol = 1e-8", 
                              "K = 5, Lmax = 1000, tol = 1e-13"), 
                            each = M * length(model_record)))

p_gr_compar <- df %>%
  ggplot(aes(y = model, #reorder(model, n_leapfrogs, FUN = median), 
             x = n_counts, color = type)) + 
  geom_boxplot(outlier.size = 0.1) +
  scale_colour_manual(values = cbbPalette) +
  scale_x_log10(breaks=c(10, 1e2, 1e3, 1e4, 5e4), 
                labels = c("10", "100", "1000", "10,000", "50,000")) + 
  ylab("") + xlab("") + #xlab("calls to lp__") + 
  theme_bw(base_size = 12)+
  theme(legend.position="top", legend.title = element_blank())   

width <- 12.0
height <- 8.0

setEPS()
postscript("../pics/phI_adapt_default/cost_gr_sens.eps",  #_resam_each
           width = width, height = height)
print(p_gr_compar)
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

## run main_pf with i = 15, seed_list = 1:20, sample with PSIS WOR for 20 inits
# run wasserstain_check with i = 15 and generate the plots
# then
## run main_pf.R with init_bound = 15.0 and generate the plots
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
setEPS()
postscript("../pics/phI_adapt_default/8-school_points22.eps",  #8-school_points22.eps
           width = 5.0, height = 5.0)
print(p_check)
dev.off()

p_check1 <- ggplot(dta_check, aes(x=x, y=y) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=-1) +
  scale_x_continuous(expand = c(0, 0), limits = c(-18, 22)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-20, 4)) +
  geom_point(data = dta_opt, 
             aes(x = optim_x, y = optim_y, group = tr_id, 
                 color = optim_ind, alpha = 0.5), size = 1) +
  geom_line(data = dta_opt, 
            aes(x = optim_x, y = optim_y, group = tr_id, 
                color = optim_ind, alpha = 0.5)) +
  scale_color_gradient(low="white", high="orange") +
  xlab("") + ylab("") +
  theme(
    legend.position='none'
  ) 
p_check1

ggsave("8-school_opt_tr22.eps", #"8-school_opt_tr22.eps"
       plot = p_check1,
       device = cairo_ps,
       path = "../pics/phI_adapt_default/",
       width = 5.0, height = 5.0, units = "in")


# multimodality #
## run main_pf with i = 3, seed_list = 1:20, sample with PSIS WOR for 20 inits
# run wasserstain_check with i = 3 and generate the plots
opt_tr_res <- get_opt_tr(opath)
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
setEPS()
postscript("../pics/phI_adapt_default/3_points.eps",  #_resam_each
           width = 5.0, height = 5.0)
print(p_check)
dev.off()

p_check1 <- ggplot(dta_check, aes(x=x, y=y) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=-1) +
  scale_x_continuous(expand = c(0, 0), limits = c(-1.5, 2.5)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-5, 3)) +
  geom_point(data = dta_opt, 
             aes(x = optim_x, y = optim_y, group = tr_id, 
                 color = optim_ind), size = 1) +
  geom_line(data = dta_opt, 
            aes(x = optim_x, y = optim_y, group = tr_id, 
                color = optim_ind)) +
  scale_color_gradient(low="white", high="orange") +
  xlab("") + ylab("") +
  theme(
    legend.position='none'
  ) 
p_check1

ggsave("3_opt_tr.eps",
       plot = p_check1,
       device = cairo_ps,
       path = "../pics/phI_adapt_default/",
       width = 5.0, height = 5.0, units = "in")


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
setEPS()
postscript("../pics/phI_adapt_default/3_phI.eps",  #_resam_each
           width = 5.0, height = 5.0)
print(p_phI)
dev.off()

# non-identifiable problem #
## run main_pf with i = 32, seed_list = 1:20, sample with PSIS WOR for 20 inits
## run wasserstain_check with i = 32 
opt_tr_res <- get_opt_tr(opath)
check_dim <- c(45, 46)# c(1, 2) #c(43, 44) #c(45, 46) #  c(3, 4)
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

print(p_check)
setEPS()
postscript("../pics/phI_adapt_default/32_points.eps",  #_resam_each
           width = 5.0, height = 5.0)
print(p_check)
dev.off()

p_check1 <- ggplot(dta_check, aes(x=x, y=y) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=-1) +
  scale_x_continuous(expand = c(0, 0), limits = c(-2.2, 4.6)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-7, 1)) +
  geom_point(data = dta_opt,
             aes(x = optim_x, y = optim_y, group = tr_id,
                 color = optim_ind), size = 1) +
  geom_line(data = dta_opt,
            aes(x = optim_x, y = optim_y, group = tr_id,
                color = optim_ind)) +
  scale_color_gradient(low="white", high="orange") +
  xlab("") + ylab("") +
  theme(
    legend.position='none'
  )
print(p_check1)


ggsave("32_opt_tr.eps",
       plot = p_check1,
       device = cairo_ps,
       path = "../pics/phI_adapt_default/",
       width = 5.0, height = 5.0, units = "in")

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
setEPS()
postscript("../pics/phI_adapt_default/32_phI.eps",  #_resam_each
           width = 5.0, height = 5.0)
print(p_phI)
dev.off()

opt_tr_res <- get_opt_tr(opath)
check_dim <- c(1, 2)# c(1, 2) #c(43, 44) #c(45, 46) #  c(3, 4)
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

p_check1 <- ggplot(dta_check, aes(x=x, y=y) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=-1) +
  scale_x_continuous(expand = c(0, 0), limits = c(-50, 35))+#, limits = c(-4, 1.5)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-4.5, 6))+#, limits = c(-2.5, 6)) +
  geom_point(data = dta_opt, 
             aes(x = optim_x, y = optim_y, group = tr_id, 
                 color = optim_ind), size = 1) +
  geom_line(data = dta_opt, 
            aes(x = optim_x, y = optim_y, group = tr_id, 
                color = optim_ind)) +
  scale_color_gradient(low="white", high="orange") +
  geom_point(data = dta_sam, aes(x=sam_x, y=sam_y), colour="red", size = 3) +
  xlab("") + ylab("") +
  theme(
    legend.position='none'
  ) 
print(p_check1)
ggsave("32_opt_tr_12.eps",
       plot = p_check1,
       device = cairo_ps,
       path = "../pics/phI_adapt_default/",
       width = 5.0, height = 5.0, units = "in")

p_phI <- ggplot(dta_check, aes(x=x, y=y) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=-1) +
  scale_x_continuous(expand = c(0, 0), limits = c(-50, 35)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-4.5, 6)) +
  geom_point(data = dta_phI, aes(x=x, y=y), colour="red", size = 3) +
  xlab("") + ylab("") +
  theme(
    legend.position='none'
  )

p_phI
ggsave("32_phI_12.eps",
       plot = p_phI,
       device = cairo_ps,
       path = "../pics/phI_adapt_default/",
       width = 5.0, height = 5.0, units = "in")

# amplified plots
p_check1 <- ggplot(dta_check, aes(x=x, y=y) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=-1) +
  scale_x_continuous(expand = c(0, 0), limits = c(-4, 2))+#, limits = c(-4, 2)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-4.5, 6))+#, limits = c(-2.5, 6)) +
  geom_point(data = dta_opt, 
             aes(x = optim_x, y = optim_y, group = tr_id, 
                 color = optim_ind), size = 1) +
  geom_line(data = dta_opt, 
            aes(x = optim_x, y = optim_y, group = tr_id, 
                color = optim_ind)) +
  scale_color_gradient(low="white", high="orange") +
  geom_point(data = dta_sam, aes(x=sam_x, y=sam_y), colour="red",
             size = 3)+
  xlab("") + ylab("") +
  theme(
    legend.position='none'
  ) 
print(p_check1)
ggsave("32_opt_tr_12_s.eps",
       plot = p_check1,
       device = cairo_ps,
       path = "../pics/phI_adapt_default/",
       width = 5.0, height = 5.0, units = "in")

# special initials #
init_ls <- list()
for (s in 1:length(seed_list)){
  init_ls[[s]] <- ref_samples[s, ]
  init_ls[[s]][2] <- init_ls[[s]][2] - 2
}
t <- proc.time()
opath <- opt_path_stan_init_parallel(
  init_ls, mc.cores, model, data, init_bound = 2.0, 
  N1, N_sam_DIV, N_sam, factr_tol, lmm, seed_list)
print(proc.time() - t)
pick_samples <- Imp_Resam_WOR(opath, n_inits = 20, seed = 1)

opt_tr_res <- get_opt_tr(opath)
check_dim <- c(1, 2)# c(1, 2) #c(43, 44) #c(45, 46) #  c(3, 4)
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

p_check1 <- ggplot(dta_check, aes(x=x, y=y) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=-1) +
  scale_x_continuous(expand = c(0, 0), limits = c(-50, 35))+#, limits = c(-4, 1.5)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-4.5, 6))+#, limits = c(-2.5, 6)) +
  geom_point(data = dta_opt, 
             aes(x = optim_x, y = optim_y, group = tr_id, 
                 color = optim_ind), size = 1) +
  geom_line(data = dta_opt, 
            aes(x = optim_x, y = optim_y, group = tr_id, 
                color = optim_ind)) +
  scale_color_gradient(low="white", high="orange") +
  geom_point(data = dta_sam, aes(x=sam_x, y=sam_y), colour="red", size = 3) +
  xlab("") + ylab("") +
  theme(
    legend.position='none'
  ) 
print(p_check1)
ggsave("32_opt_tr_12_2.eps",
       plot = p_check1,
       device = cairo_ps,
       path = "../pics/phI_adapt_default/",
       width = 5.0, height = 5.0, units = "in")

# non-Gaussian problem #
# i = 18
opt_tr_res <- get_opt_tr(opath)
check_dim <- c(2, 3)#c(2, 3)
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

p_check1 <- ggplot(dta_check, aes(x=x, y=y) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=-1) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 2))+#, limits = c(-4, 1.5)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 3.5))+#, limits = c(-2.5, 6)) +
  geom_point(data = dta_sam, aes(x=sam_x, y=sam_y), colour="red",
             size = 3) +
  xlab("") + ylab("") +
  theme(
    legend.position='none'
  )
print(p_check1)
ggsave("18_opt_points.eps",
       plot = p_check1,
       device = cairo_ps,
       path = "../pics/phI_adapt_default/",
       width = 5.0, height = 5.0, units = "in")


p_check2 <- ggplot(dta_check, aes(x=x, y=y) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=-1) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 2))+#, limits = c(-4, 1.5)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 3.5))+#, limits = c(-2.5, 6)) +
  geom_point(data = dta_opt,
             aes(x = optim_x, y = optim_y, group = tr_id,
                 color = optim_ind), size = 1) +
  geom_line(data = dta_opt,
            aes(x = optim_x, y = optim_y, group = tr_id,
                color = optim_ind)) +
  scale_color_gradient(low="white", high="orange") +
  xlab("") + ylab("") +
  theme(
    legend.position='none'
  )
print(p_check2)
ggsave("18_opt_tr.eps",
       plot = p_check2,
       device = cairo_ps,
       path = "../pics/phI_adapt_default/",
       width = 5.0, height = 5.0, units = "in")


## check plots ##
opt_tr_res <- get_opt_tr(opath)
for(first_check in 1:6){
  check_dim <- c(2 * first_check - 1, 2 * first_check)
  check_dim <- c(2, 3)#c(2, 3)
  cat(check_dim, "\n")
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
  
  p_check1 <- ggplot(dta_check, aes(x=x, y=y) ) +
    stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
    scale_fill_distiller(palette=4, direction=-1) +
    scale_x_continuous(expand = c(0, 0))+#, limits = c(-2.5, 6))+#, limits = c(-4, 1.5)) +
    scale_y_continuous(expand = c(0, 0))+#, limits = c(-5, 1))+#, limits = c(-2.5, 6)) +
    geom_point(data = dta_opt,
               aes(x = optim_x, y = optim_y, group = tr_id,
                   color = optim_ind), size = 1) +
    geom_line(data = dta_opt,
              aes(x = optim_x, y = optim_y, group = tr_id,
                  color = optim_ind)) +
    scale_color_gradient(low="white", high="orange") +
    geom_point(data = dta_sam, aes(x=sam_x, y=sam_y), colour="red",
               size = 3) +
    xlab("") + ylab("") +
    theme(
      legend.position='none'
    )
  print(p_check1)
  readline(prompt="Press [enter] to continue:")
  
  p_check2 <- ggplot(dta_check, aes(x=x, y=y) ) +
    geom_point(size = 3) + 
    geom_point(data = dta_opt,
               aes(x = optim_x, y = optim_y, group = tr_id,
                   color = optim_ind), size = 1) +
    geom_line(data = dta_opt,
              aes(x = optim_x, y = optim_y, group = tr_id,
                  color = optim_ind)) +
    scale_color_gradient(low="white", high="orange") +
    geom_point(data = dta_sam, aes(x=sam_x, y=sam_y), colour="red",
               size = 3) +
    xlab("") + ylab("") +
    theme(
      legend.position='none'
    ) 
  print(p_check2)
  readline(prompt="Press [enter] to continue:")
  
  p_phI <- ggplot(dta_check, aes(x=x, y=y) ) +
    stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
    scale_fill_distiller(palette=4, direction=-1) +
    scale_x_continuous(expand = c(0, 0))+ #, limits = c(-2.5, 6)) +
    scale_y_continuous(expand = c(0, 0))+ #, limits = c(-5, 1)) +
    geom_point(data = dta_phI, aes(x=x, y=y), colour="red", size = 3) +
    xlab("") + ylab("") +
    theme(
      legend.position='none'
    )
  
  p_phI
}
#1,2, 3,4, 6
## plots for section 2.4 ##
## run wasserstain_check with i = 6 
# 0.03386512 vs  # 0.01100075
pick_samples_center <- ADVI_meanfield_center[[i]]
a_center = wpp(t(pick_samples_center), 
               mass = rep(1 / ncol(pick_samples_center), ncol(pick_samples_center)))
b = wpp(ref_samples, mass = rep(1 / nrow(ref_samples), nrow(ref_samples)))
w_d_c <- wasserstein(a_center, b, p = 2); w_d_c #0.03386512

pick_samples <- lp_opath[[i]]$pick_samples
a = wpp(t(pick_samples), 
        mass = rep(1 / ncol(pick_samples), ncol(pick_samples)))
b = wpp(ref_samples, mass = rep(1 / nrow(ref_samples), nrow(ref_samples)))
w_d_pf <- wasserstein(a, b, p = 2); w_d_pf #0.01100075

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
  scale_x_continuous(expand = c(0, 0), limits = c(-0.45, -0.26)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.03, 0.16)) +
  geom_point(data = dta_sam, aes(x=sam_x, y=sam_y), colour="red", size = 3) +
  xlab("") + ylab("") +
  theme(
    legend.position='none'
  )

p_check
ggsave("6-dog_pf.eps",
       plot = p_check,
       device = cairo_ps,
       path = "../pics/phI_adapt_default/",
       width = 5.0, height = 5.0, units = "in")


p_check_c <- ggplot(dta_check, aes(x=x, y=y) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=-1) +
  scale_x_continuous(expand = c(0, 0), limits = c(-0.45, -0.26)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-0.03, 0.16)) + 
  geom_point(data = dta_sam_center, aes(x=sam_x, y=sam_y),
             colour="red", size = 3) +
  xlab("") + ylab("") +
  theme(
    legend.position='none'
  )

p_check_c
ggsave("6-dog_c.eps",
       plot = p_check_c,
       device = cairo_ps,
       path = "../pics/phI_adapt_default/",
       width = 5.0, height = 5.0, units = "in")


