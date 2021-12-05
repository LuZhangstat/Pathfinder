setwd("./posteriordb")
library(ggplot2)
library(posteriordb)
library(RColorBrewer)
source("../utils/lp_utils.R")
source("../utils/sim_pf.R")
pd <- pdb_local("./") # Posterior database connection
pn <- posterior_names(pd)


# colorblind-friendly palettes
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


## HMC vs L-BFGS plot in Introduction (old test)##
# load("../results/lp_posteriordb_LBFGS_h6.RData")
# load("../results/lp_posteriordb_explore_h6.RData")
# M = 20
# df <- data.frame(n_counts = c(c(abs(lp_LBFGS_n_gr)), c(lp_explore_n_leapfrog)),
#                  model = rep(rep(pn[model_record], each = M), 2),
#                  n_leapfrogs = rep(c(lp_explore_n_leapfrog), 2),
#                  not_reach_target = 
#                    c(rep(apply(lp_LBFGS_n_gr, 2, 
#                                f <- function(x){as.numeric(any(x < 0))}), 
#                          each = M), rep(2, M*length(model_record))))
# 
# width <- 12.0
# height <- 8.0
# setEPS()
# postscript("../pics/box_compar_LBFGS_log.eps",  #_resam_each
#            width = width, height = height)
# p_box_compar <- df %>% mutate(type= c("L-BFGS", "L-BFGS(have minor modes)", 
#                                       "Stan Phase I")
#                               [not_reach_target + 1]) %>%
#   ggplot(aes(y = model, #reorder(model, n_leapfrogs, FUN = median), 
#              x = n_counts, color = type)) + 
#   geom_boxplot() + 
#   scale_colour_manual(values=cbbPalette) + 
#   scale_x_log10(breaks=c(10, 1e2, 1e3, 1e4, 1e5), 
#                 labels = c("10", "100", "1000", "10,000", "100,000")) + 
#   ylab("") + xlab("count of log density and its gradient evaluation ") + #xlab("calls to log density and gradient") + 
#   theme_bw(base_size = 12 )+
#   theme(legend.position="top", legend.title = element_blank()) 
# print(p_box_compar)
# dev.off()

## Plots in Section 3 ##
load("../results/lp_posteriordb_LBFGS_h6.RData")
load("../results/lp_posteriordb_phI_adapt_default.RData")

# ### Wasserstein distance check for 100 repeats ###
# M = 100
# load("../results/wasserstein_phI_adapt_default_RVI.RData")
# 
# # check pathfinder vs phase I warmup
# pf_vs_phI <- (w_d_matrix[, "pf"] / w_d_matrix[, "PhI"])
# summary(pf_vs_phI)
# quantile(pf_vs_phI, c(0.05, 0.5, 0.95))
# 
# pf_vs_phI <- (w_d_matrix[, "pf_20_IR"] / w_d_matrix[, "PhI"])
# summary(pf_vs_phI)
# quantile(pf_vs_phI, c(0.05, 0.5, 0.95))
# table(pf_vs_phI<1.2)
# max(pf_vs_phI)
# 
# ## rank methods ##
# rank_score <- apply(w_d_matrix[, c("pf", "pf_20_IR", "PhI", "meanfield", 
#                                    "meanfield center", "fullrank")], 1, 
#                     f <- function(x){order(order(x))})
# table(rank_score[2, ]==1)
# table(rank_score[2, ]==2)
# table(rank_score[2, ]<4)
# table(rank_score[1, ]==1)
# table(rank_score[3, ]<3)
# rowSums(rank_score)
# round(rowMeans(rank_score), 1)
# 
# # faceted rank histograms
# rank_dat = data.frame(rank = c(rank_score),
#                       label = rep(c("pathfinder", 
#                                     "multi-path pathfinder", 
#                                     "Stan Phase I",
#                                     "mean-field ADVI",
#                                     "mean-field center",
#                                     "dense ADVI"), length(model_record)))
# 
# p_rank <- ggplot(rank_dat, aes(x = rank)) + geom_histogram(bins = 6) + 
#   facet_wrap(~label)
# p_rank
# width <- 8.0
# height <- 4.0
# setEPS()
# postscript("../pics/phI_adapt_default/rank_test.eps",  #_resam_each
#            width = width, height = height)
# print(p_rank)
# dev.off()
# 
# ## rank comparision ##
# rank_M <- matrix(NA, 6, 6)
# method_names <- c("pf", "pf_20_IR", "PhI", "meanfield", 
#                   "meanfield center", "fullrank")
# for(A in 1:6){
#   for(B in 1:6){
#     rank_M[A, B] <- sum(w_d_matrix[, method_names[A]] < 
#                           w_d_matrix[, method_names[B]])
#   }
# }
# rownames(rank_M) <- c("pf <", "pf_20_IR <", "PhI <", "meanfield <", 
#                       "meanfield center <", "fullrank <")
# colnames(rank_M) <- c("pf", "pf_20_IR", "PhI", "meanfield", 
#                       "meanfield center", "fullrank")
# rank_M <- round(rank_M/49, 2)
# rank_M
# 
# # compare pathfinder with ADVI #
# summary(w_d_matrix[, "pf"] / w_d_matrix[, "meanfield"])
# summary(w_d_matrix[, "pf"] / w_d_matrix[, "fullrank"])
# summary(w_d_matrix[, "pf"] / w_d_matrix[, "meanfield center"])
# summary(w_d_matrix[, "pf_20_IR"] / w_d_matrix[, "meanfield"])
# summary(w_d_matrix[, "pf_20_IR"] / w_d_matrix[, "fullrank"])
# summary(w_d_matrix[, "pf_20_IR"] / w_d_matrix[, "meanfield center"])
# 
# ## one plot ##
# summary(w_d_matrix[, "pf_20_IR"] / w_d_matrix[, "PhI"])
# ratios <- c(w_d_matrix[, c("pf", "pf_20_IR", "PhI", "meanfield", 
#                            "meanfield center", "fullrank")] / 
#               w_d_matrix[, "pf_20_IR"])
# ratios[ratios > 2^10] <- 2^10
# w_d_dat = data.frame(ratios = ratios,
#                      label = c(
#                        rep("pathfinder", length(model_record)), 
#                        rep("multi-path pathfinder", length(model_record)),
#                        rep("Stan Phase I", length(model_record)),
#                        rep("mean-field ADVI", length(model_record)),
#                        rep("mean-field center", length(model_record)),
#                        rep("dense ADVI", length(model_record))),
#                      model = rep(pn[model_record], 6))#6))
# 
# p_w_d_compar <- w_d_dat %>% 
#   ggplot(aes(y = model, #reorder(model, ratios, FUN = mean), 
#              x = ratios, color = label, shape = label)) + 
#   geom_point(size = 2) +
#   scale_x_continuous(trans = 'log2',
#                      limits = c(1/4, 1024),
#                      breaks = c(1/4, 1/2, 1, 2, 4, 
#                                 8, 16, 32, 64, 128, 256, 512, 1024),
#                      labels = c("1/4", 
#                                 "1/2", "1", "2", "4", "8", 
#                                 "16", "32", "64", "128", 
#                                 "256", "512", ">=1024")) + 
#   ylab("") + xlab("scaled Wasserstein distance") + 
#   theme_bw(base_size = 12)  +
#   theme(legend.position="top", legend.title = element_blank()) #+
# #theme(legend.position = "none") 
# 
# width <- 12.0
# height <- 8.0
# pointsize <- 16
# 
# setEPS()
# postscript("../pics/phI_adapt_default/W_d_compar_100_each_2.eps",  #_resam_each
#            width = width, height = height)
# print(p_w_d_compar)
# dev.off()


## 100 for each ##
load("../results/wasserstein_phI.RData")
load("../results/wasserstein_100_default_W1_WR_updat.RData")
load("../results/Laplace_results.RData")
M = 100

#pick_ind <- c(1:7, 13, 15:21, 26, 31:32, 37, 46, 47:49) # index to select the models in the summary
pick_ind <- c(1:7, 13, 15:21, 31:32, 46, 47:48) # index to select the models in the summary


# meanfield ADVI #
ratio_M_wd_mf <- (apply(W_d_100_ADVI_mf[, pick_ind], 2, f <- function(x){quantile(x, 0.5)})/
                    apply(W_d_100_pf[, pick_ind], 2, f <- function(x){quantile(x, 0.5)}))

mean(ratio_M_wd_mf)
table(ratio_M_wd_mf > 2)
# FALSE  TRUE 
# 12     8
table(ratio_M_wd_mf < 0.5)

mean_M_wd_mf <- colMeans(W_d_100_ADVI_mf[, pick_ind]) / colMeans(W_d_100_pf[, pick_ind])
table(mean_M_wd_mf > 2)
table(mean_M_wd_mf < 0.5)

ratio_M_wd_mf_IR <- (apply(W_d_100_ADVI_mf[, pick_ind], 2, f <- function(x){quantile(x, 0.5)})/
                       apply(W_d_100_pf_IR[, pick_ind], 2, f <- function(x){quantile(x, 0.5)}))
mean(ratio_M_wd_mf_IR)
mean_M_wd_mf_IR <- colMeans(W_d_100_ADVI_mf[, pick_ind]) / colMeans(W_d_100_pf_IR[, pick_ind])

# dense ADVI #
ratio_M_wd_fr <- (apply(W_d_100_ADVI_fr[, pick_ind], 2, f <- function(x){quantile(x, 0.5)})/
                    apply(W_d_100_pf[, pick_ind], 2, f <- function(x){quantile(x, 0.5)}))
mean(ratio_M_wd_fr)
table(ratio_M_wd_fr > 2)
# FALSE  TRUE 
# 11     9
table(ratio_M_wd_fr < 0.5)

mean_M_wd_fr <- colMeans(W_d_100_ADVI_fr[, pick_ind]) / colMeans(W_d_100_pf[, pick_ind])
table(mean_M_wd_fr > 2)
table(mean_M_wd_fr < 0.5)

ratio_M_wd_fr_IR <- (apply(W_d_100_ADVI_fr[, pick_ind], 2, f <- function(x){quantile(x, 0.5)})/
                       apply(W_d_100_pf_IR[, pick_ind], 2, f <- function(x){quantile(x, 0.5)}))
mean(ratio_M_wd_fr_IR)
mean_M_wd_fr_IR <- colMeans(W_d_100_ADVI_fr[, pick_ind]) / colMeans(W_d_100_pf_IR[, pick_ind])

# check both #
table((ratio_M_wd_fr > 2) & (ratio_M_wd_mf > 2))
table((ratio_M_wd_fr < 0.5) & (ratio_M_wd_mf < 0.5))

table((mean_M_wd_fr > 2) & (mean_M_wd_mf > 2) & (mean_M_wd_fr_IR > 2) & 
        (mean_M_wd_mf_IR > 2)) # 34
table((mean_M_wd_fr < 0.5) & (mean_M_wd_mf < 0.5))

# check Stan Phase I #
table((w_d_phI[pick_ind] / apply(W_d_100_pf[, pick_ind], 2, f <- function(x){quantile(x, 0.5)})) > 2)
# FALSE  TRUE 
# 13     7 
table((w_d_phI[pick_ind] / apply(W_d_100_pf[, pick_ind], 2, f <- function(x){quantile(x, 0.5)})) < 0.5)
# FALSE  TRUE 
# 19     1 

# check model bball_drive_event_0-hmm_drive_0 #
median(W_d_100_ADVI_mf[, 3]) / median(W_d_100_pf[, 3])
mean(W_d_100_pf_IR[, 3]) / mean(W_d_100_pf[, 3])

w_d_median_pf_IR <- rep(apply(W_d_100_pf_IR[, pick_ind], 2, median), each = M)
apply(W_d_100_pf[, pick_ind] / w_d_median_pf_IR, 2, median)
apply(W_d_100_ADVI_mf[, pick_ind] / w_d_median_pf_IR, 2, median)
apply(W_d_100_ADVI_fr[, pick_ind] / w_d_median_pf_IR, 2, median)

w_d_median_pf <- rep(apply(W_d_100_pf[, pick_ind], 2, median), each = M)
w_d_median_phI <- rep(w_d_phI[pick_ind], each = M)
w_d_scaled = c(#c(W_d_100_pf_IR / w_d_median_pf),
               c(W_d_100_pf[, pick_ind] / w_d_median_pf), 
               c(W_d_100_ADVI_mf[, pick_ind] / w_d_median_pf), 
               c(W_d_100_ADVI_fr[, pick_ind] / w_d_median_pf))
w_d_scaled[w_d_scaled >= 2^9] <- 2^9


w_d_PhI_scaled <- c(w_d_phI[pick_ind] / apply(W_d_100_pf[, pick_ind], 2, median))

modelnamels <- pn[model_record[pick_ind]]
modelnamels[9] <- "eight_schools_centered" 
modelnamels[10] <- "eight_schools_noncentered" 
modelnamels[14] <- "hmm_example"
modelnamels[16] <- "low_dim_gauss_mix"
modelnamels[19] <- "one_comp_mm_elim_abs"
w_d_point = data.frame(w_d_PhI_scaled = w_d_PhI_scaled,
                       model = modelnamels,
                       point = rep("Stan Phase I", length(model_record[pick_ind])))
df <- data.frame(w_d = w_d_scaled,
                 model = rep(rep(modelnamels, each = M), 3),
                 type = rep(c(#"multi-path pathfinder", 
                              1, #"pathfinder", 
                              2, #"mean-field ADVI", 
                              3 #"dense ADVI"
                              ), 
                            each = length(modelnamels)*M),
                 npar = rep(rep(dimension_record, each = M), 3), 
                 w_d_phI_scaled = rep(rep(w_d_PhI_scaled, each = M), 3))
df$type <- factor(df$type, levels = 1:3,
                     labels = c("pathfinder", "mean-field ADVI",
                                         "dense ADVI"))


width <- 12.0
height <- 10.0
setEPS()
# postscript("../pics/phI_adapt_default/W_d_box_compar_100_each.eps",  #_resam_each
#            width = width, height = height)
postscript("../pics/phI_adapt_default/W_d_box_compar_100_VI.eps",  #_resam_each
           width = width, height = height)
p_box_compar <- df %>%
  ggplot(aes(y = reorder(model, w_d_phI_scaled, FUN = median), #model, #reorder(model, n_leapfrogs, FUN = median), 
             x = w_d, color = type)) + 
  geom_boxplot(outlier.size = 0.1) + 
  geom_point(aes(y = model, x = w_d_PhI_scaled, shape = point), 
             color = "#009E73",
             data = w_d_point, size = 3, inherit.aes=FALSE) +
  scale_shape_manual(values=c(17)) +
  scale_x_continuous(trans = 'log2',
                     limits = c(1/32, 2^9),
                     breaks = c(1/32, 1/16, 1/8, 1/4,
                                1/2, 1, 2, 4, 8, 16, 32, 64, 128, 256,
                                512),
                     labels = c("1/32", "1/16",
                                "1/8", "1/4", "1/2", "1", "2", "4", "8",
                                "16", "32", "64", "128",
                                "256", ">512")) +
  scale_colour_manual(values=c("#56B4E9", "#E69F00", "#000000")) +  #cbbPalette
  ylab("") + xlab("scaled 1-Wasserstein distance") + #xlab("calls to log density and gradient") + 
  theme_bw(base_size = 20)+
  theme(legend.position="top", legend.title = element_blank()) 
print(p_box_compar)
dev.off()

## compare single- & multi-path Pathfinder and Stan phase I ##
w_d_scaled = c(c(W_d_100_pf_IR[, pick_ind] / w_d_median_pf),
  c(W_d_100_pf[, pick_ind] / w_d_median_pf))

w_d_scaled[w_d_scaled >= 2^9] <- 2^9

w_d_PhI_scaled <- c(w_d_phI[pick_ind] / apply(W_d_100_pf[, pick_ind], 2, median))
w_d_point = data.frame(w_d_PhI_scaled = w_d_PhI_scaled,
                       model = modelnamels,
                       point = rep("Stan Phase I", 
                                   length(modelnamels)))

df <- data.frame(w_d = w_d_scaled,
                 model = rep(rep(modelnamels, each = M), 2),
                 type = rep(c("multi-path pathfinder", "pathfinder"), 
                   each = length(modelnamels)*M), 
                 npar = rep(rep(dimension_record, each = M), 2),
                 w_d_phI_scaled = rep(rep(w_d_PhI_scaled, each = M), 2))



w_d_PhI_scaled_2 <- c(w_d_phI[pick_ind] / apply(W_d_100_pf_IR[, pick_ind], 2, median))
table(w_d_PhI_scaled > 2)
table(w_d_PhI_scaled < 0.5)

width <- 12.0
height <- 10.0
setEPS()
postscript("../pics/phI_adapt_default/W_d_box_compar_100_pf.eps",  #_resam_each
           width = width, height = height)
p_box_compar <- df %>%
  ggplot(aes(y = reorder(model, w_d_phI_scaled, FUN = median), #model, #reorder(model, n_leapfrogs, FUN = median), 
             x = w_d, color = type)) + 
  geom_boxplot(outlier.size = 0.1) + 
  geom_point(aes(y = model, x = w_d_PhI_scaled, shape = point), 
             color = "#009E73",
             data = w_d_point, size = 3, inherit.aes=FALSE) +
  scale_shape_manual(values=c(17)) +
  #geom_vline(xintercept=1, color = "red", size=0.5)+
  scale_colour_manual(values=c("#CC79A7", "#56B4E9")) + 
  scale_x_continuous(trans = 'log2',
                     limits = c(1/32, 128),
                     breaks = c(1/32, 1/16, 1/8, 1/4,
                                1/2, 1, 2, 4, 8, 16, 32, 64, 128),
                     labels = c("1/32", "1/16",
                                "1/8", "1/4", "1/2", "1", "2", "4", "8",
                                "16", "32", "64", "128")) + 
  ylab("") + xlab("scaled 1-Wasserstein distance") + #xlab("calls to log density and gradient") + 
  theme_bw(base_size = 20)+
  theme(legend.position="top", legend.title = element_blank()) 
print(p_box_compar)
dev.off()


## computational cost comparision ##
load("../results/lp_posteriordb_phI_adapt_default.RData") # Pathfinder #_resam_all #_resam_each
load("../results/PhI_100_h10.RData")
load("../results/ADVI_100_updat.RData")
load("../results/multi_pf_samples_default.RData")

pathfinder_fn_call <- 
  sapply(lp_opath[pick_ind], f <- function(x){
    sapply(x$opath, g <- function(z){z$fn_call})})
pathfinder_gr_call <- 
  sapply(lp_opath[pick_ind], f <- function(x){
    sapply(x$opath, g <- function(z){sum(z$gr_call)}) } )

summary(colSums(PhI_leapfrog_counts[, pick_ind]) / colSums(pathfinder_fn_call))
summary(colSums(PhI_leapfrog_counts[, pick_ind]) / colSums(pathfinder_gr_call))
table((colSums(PhI_leapfrog_counts[, pick_ind]) / colSums(pathfinder_gr_call)) > 100)
summary(colSums(calls_lp_mean[, pick_ind] + 
                  calls_gr_mean[, pick_ind]) / colSums(pathfinder_fn_call))
summary(colSums(calls_gr_mean[, pick_ind]) / colSums(pathfinder_gr_call))
summary(colSums(calls_lp_full[, pick_ind] + calls_gr_full[, pick_ind]) / 
          colSums(pathfinder_fn_call))
summary(colSums(calls_gr_full[, pick_ind]) / colSums(pathfinder_gr_call))

# summary(calls_lp_mean[, pick_ind] / colMeans(pathfinder_fn_call))
# summary(calls_gr_mean[, pick_ind] / colMeans(pathfinder_gr_call))
# summary(calls_lp_full[, pick_ind] / colMeans(pathfinder_fn_call))
# summary(calls_gr_full[, pick_ind] / colMeans(pathfinder_gr_call))

df <- data.frame(n_counts = c(c(pathfinder_fn_call), 
                              c(calls_lp_mean[, pick_ind] + calls_gr_mean[, pick_ind]), 
                              c(calls_lp_full[, pick_ind] + calls_gr_full[, pick_ind]),
                              c(PhI_leapfrog_counts[, pick_ind])),
                 model = rep(rep(modelnamels, each = M), 4),
                 type = rep(1:4, 
                            each = M * length(model_record[pick_ind])),
                 w_d_phI_scaled = rep(rep(w_d_PhI_scaled, each = M), 4))

df$type <- factor(df$type, levels = 1:4, 
                  labels = c("Pathfinder", "mean-field ADVI", 
                             "dense ADVI", "Stan Phase I"))
p_lp_compar <- df %>%
  ggplot(aes(y = reorder(model, w_d_phI_scaled, FUN = median), #model, #reorder(model, n_leapfrogs, FUN = median), 
             x = n_counts, color = type)) + 
  geom_boxplot(outlier.size = 0.1) +
  scale_colour_manual(values = c("#56B4E9", "#E69F00", "#000000", "#009E73")) + #cbbPalette
  scale_x_log10(breaks=c(10, 1e2, 1e3, 1e4, 5e4), 
                labels = c("10", "100", "1000", "10,000", "50,000")) + 
  ylab("") + xlab("count of log density evaluation ") + #xlab("calls to lp__") + 
  theme_bw(base_size = 20)+
  theme(legend.position="top", legend.title = element_blank())   

width <- 12.0
height <- 8.0

setEPS()
postscript("../pics/phI_adapt_default/cost_lp_100_each.eps",  #_resam_each
           width = width, height = height)
print(p_lp_compar)
dev.off()

df <- data.frame(n_counts = c(c(pathfinder_gr_call), 
                              c(calls_gr_mean[, pick_ind]), 
                              c(calls_gr_full[, pick_ind]),
                              c(PhI_leapfrog_counts[, pick_ind])),
                 model = rep(rep(modelnamels, each = M), 4),
                 type = rep(1:4, 
                            each = M * length(model_record[pick_ind])),
                 w_d_phI_scaled = rep(rep(w_d_PhI_scaled, each = M), 4))

df$type <- factor(df$type, levels = 1:4, 
                  labels = c("Pathfinder", "mean-field ADVI", 
                             "dense ADVI", "Stan Phase I"))

p_gr_compar <- df %>%
  ggplot(aes(y = reorder(model, w_d_phI_scaled, FUN = median), #model, 
             x = n_counts, color = type)) + 
  geom_boxplot(outlier.size = 0.1) +
  scale_colour_manual(values = c("#56B4E9", "#E69F00", "#000000", "#009E73")) +
  scale_x_log10(breaks=c(10, 1e2, 1e3, 1e4, 5e4), 
                labels = c("10", "100", "1000", "10,000", "50,000")) + 
  ylab("") + xlab("count of log density gradient evaluation") + #xlab("calls to lp__") + 
  theme_bw(base_size = 20)+
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
load("../results/wasserstein_100_default_W1_WR.RData")
W_d_100_pf_default <- W_d_100_pf
load("../results/wasserstein_100_sen_W1.RData")
load("../results/wasserstein_100_sen_W1_L_revise.RData") # add the sensitive test for revision


## Number of monte carlo samples in ELBO estimation ##
mean(apply(W_d_100_pf_large_K[, pick_ind], 2, f <- function(x){quantile(x, 0.5)}) / 
       apply(W_d_100_pf_default[, pick_ind], 2, f <- function(x){quantile(x, 0.5)}))
# 1 - 0.9677836
# [1] 0.0322164
min(apply(W_d_100_pf_large_K[, pick_ind], 2, f <- function(x){quantile(x, 0.5)}) / 
      apply(W_d_100_pf_default[, pick_ind], 2, f <- function(x){quantile(x, 0.5)}))
# 1 - 0.8489219
# [1] 0.1510781

load("../results/lp_posteriordb_phI_adapt_large_K.RData") #
pathfinder_large_K_fn_call <- 
  sapply(lp_opath[pick_ind], f <- function(x){
    sapply(x$opath, g <- function(z){sum(z$fn_call)}) } )
pathfinder_large_K_gr_call <- 
  sapply(lp_opath[pick_ind], f <- function(x){
    sapply(x$opath, g <- function(z){sum(z$gr_call)}) } )
sum(pathfinder_large_K_fn_call)/sum(pathfinder_fn_call)
load("../results/lp_posteriordb_phI_adapt_default.RData") 


w_d_median_pf <- rep(apply(W_d_100_pf_default[, pick_ind], 2, median), each = M)
apply(W_d_100_pf_large_K[, pick_ind] / w_d_median_pf, 2, median)
apply(W_d_100_pf_short_L[, pick_ind] / w_d_median_pf, 2, median)
w_d_scaled = c(
  c(W_d_100_pf_large_K[, pick_ind] / w_d_median_pf), 
  c(W_d_100_pf_default[, pick_ind] / w_d_median_pf))
range(w_d_scaled)
w_d_scaled[w_d_scaled >= 2^7] <- 2^7
df <- data.frame(w_d = w_d_scaled,
                 model = rep(rep(modelnamels, each = M), 2),
                 type = rep(1:2, 
                            each = length(model_record[pick_ind])*M),
                 w_d_phI_scaled = rep(rep(w_d_PhI_scaled, each = M), 2))
df$type <- factor(df$type, levels = 1:2, 
                  labels = c("Lmax = 1000, tol = 1e-13, K = 30, J = 6", 
                             "Lmax = 1000, tol = 1e-13, K = 5, J = 6"))

width <- 12.0
height <- 10.0
setEPS()
postscript("../pics/phI_adapt_default/W_d_box_sens_K.eps",  #_resam_each
           width = width, height = height)
p_box_compar <- df %>%
  ggplot(aes(y = reorder(model, w_d_phI_scaled, FUN = median), #model, # 
             x = w_d, color = type)) + 
  geom_boxplot(outlier.size = 0.1) +
  scale_colour_manual(values = c("#000000", "#E69F00")) +
  scale_x_continuous(trans = 'log2',
                     limits = c(1/16, 2^7),
                     breaks = c(1/16, 1/8, 1/4, 1/2, 1, 2, 4, 
                                8, 16, 32, 64, 128),
                     labels = c("1/16", "1/8", "1/4", 
                                "1/2", "1", "2", "4", "8", 
                                "16", "32", "64", "128")) + 
  ylab("") + xlab("scaled 1-Wasserstein distance") + #xlab("calls to log density and gradient") + 
  theme_bw(base_size = 19)+
  theme(legend.position="top", legend.title = element_blank()) 
print(p_box_compar)
dev.off()

## length and convergence tolerance ##
w_d_scaled = c(c(W_d_100_pf_short_L[, pick_ind] / w_d_median_pf), 
               c(W_d_100_pf_default[, pick_ind] / w_d_median_pf))
range(w_d_scaled)
df <- data.frame(w_d = w_d_scaled,
                 model = rep(rep(modelnamels, each = M), 2),
                 type = rep(1:2, 
                            each = length(model_record[pick_ind])*M),
                 w_d_phI_scaled = rep(rep(w_d_PhI_scaled, each = M), 2))
df$type <- factor(df$type, levels = 1:2, 
                  labels = c("Lmax = 200, tol = 1e-8, K = 5, J = 6",  
                             "Lmax = 1000, tol = 1e-13, K = 5, J = 6"))

width <- 12.0
height <- 10.0
setEPS()
postscript("../pics/phI_adapt_default/W_d_box_sens_L.eps",  #_resam_each
           width = width, height = height)
p_box_compar <- df %>%
  ggplot(aes(y = reorder(model, w_d_phI_scaled, FUN = median), #model, # 
             x = w_d, color = type)) + 
  geom_boxplot(outlier.size = 0.1) +
  scale_colour_manual(values = c("#000000", "#E69F00")) +
  scale_x_continuous(trans = 'log2',
                     limits = c(1/16, 128),
                     breaks = c(1/16, 1/8, 1/4, 1/2, 1, 2, 4, 
                                8, 16, 32, 64, 128),
                     labels = c("1/16", "1/8", "1/4", 
                                "1/2", "1", "2", "4", "8", 
                                "16", "32", "64", "128")) + 
  ylab("") + xlab("scaled 1-Wasserstein distance") + #xlab("calls to log density and gradient") + 
  theme_bw(base_size = 19)+
  theme(legend.position="top", legend.title = element_blank()) 
print(p_box_compar)
dev.off()

## length and convergence tolerance revised ##
w_d_scaled = c(c(W_d_100_pf_short_L_revise[, pick_ind] / w_d_median_pf), 
               c(W_d_100_pf_default[, pick_ind] / w_d_median_pf))
range(w_d_scaled)
df <- data.frame(w_d = w_d_scaled,
                 model = rep(rep(modelnamels, each = M), 2),
                 type = rep(1:2, 
                            each = length(model_record[pick_ind])*M),
                 w_d_phI_scaled = rep(rep(w_d_PhI_scaled, each = M), 2))
df$type <- factor(df$type, levels = 1:2, 
                  labels = c("Lmax = 100, tol = 1e-3, K = 5, J = 6",  
                             "Lmax = 1000, tol = 1e-13, K = 5, J = 6"))
width <- 12.0
height <- 10.0
setEPS()
postscript("../pics/phI_adapt_default/W_d_box_sens_L_revise.eps",  #_resam_each
           width = width, height = height)
p_box_compar <- df %>%
  ggplot(aes(y = reorder(model, w_d_phI_scaled, FUN = median), #model, # 
             x = w_d, color = type)) + 
  geom_boxplot(outlier.size = 0.1) +
  scale_colour_manual(values = c("#000000", "#E69F00")) +
  scale_x_continuous(trans = 'log2',
                     limits = c(1/16, 128),
                     breaks = c(1/16, 1/8, 1/4, 1/2, 1, 2, 4, 
                                8, 16, 32, 64, 128),
                     labels = c("1/16", "1/8", "1/4", 
                                "1/2", "1", "2", "4", "8", 
                                "16", "32", "64", "128")) + 
  ylab("") + xlab("scaled 1-Wasserstein distance") + #xlab("calls to log density and gradient") + 
  theme_bw(base_size = 19)+
  theme(legend.position="top", legend.title = element_blank()) 
print(p_box_compar)
dev.off()

## history size ##
mean((apply(W_d_100_pf_long_hist[, pick_ind], 2, f <- function(x){quantile(x, 0.5)}) / 
        apply(W_d_100_pf_default[, pick_ind], 2, f <- function(x){quantile(x, 0.5)}))[-5])
# 0.9735699
median((apply(W_d_100_pf_long_hist[, pick_ind], 2, f <- function(x){quantile(x, 0.5)}) / 
          apply(W_d_100_pf_default[, pick_ind], 2, f <- function(x){quantile(x, 0.5)}))[c(-5)])
#1- 0.9926707
#[1] 0.0073293
range((apply(W_d_100_pf_long_hist[, pick_ind], 2, f <- function(x){quantile(x, 0.5)}) / 
         apply(W_d_100_pf_default[, pick_ind], 2, f <- function(x){quantile(x, 0.5)}))[c(-5)])
# 0.8367152 1.0984029

w_d_scaled = c(c(W_d_100_pf_long_hist[, pick_ind] / w_d_median_pf), 
               c(W_d_100_pf_default[, pick_ind] / w_d_median_pf))
range(w_d_scaled)
df <- data.frame(w_d = w_d_scaled,
                 model = rep(rep(modelnamels, each = M), 2),
                 type = rep(1:2, each = length(model_record[pick_ind])*M),
                 w_d_phI_scaled = rep(rep(w_d_PhI_scaled, each = M), 2))

df$type <- factor(df$type, levels = 1:2, 
                  labels = c("Lmax = 1000, tol = 1e-13, K = 5, J = 60",  
                             "Lmax = 1000, tol = 1e-13, K = 5, J = 6"))



width <- 12.0
height <- 10.0
setEPS()
postscript("../pics/phI_adapt_default/W_d_box_sens_J.eps",  #_resam_each
           width = width, height = height)
p_box_compar <- df %>%
  ggplot(aes(y =reorder(model, w_d_phI_scaled, FUN = median), # model, #
             x = w_d, color = type)) + 
  geom_boxplot(outlier.size = 0.1) +
  scale_colour_manual(values = c("#000000", "#E69F00")) +
  scale_x_continuous(trans = 'log2',
                     limits = c(1/16, 128),
                     breaks = c(1/16, 1/8, 1/4, 1/2, 1, 2, 4, 
                                8, 16, 32, 64, 128),
                     labels = c("1/16", "1/8", "1/4", 
                                "1/2", "1", "2", "4", "8", 
                                "16", "32", "64", "128")) + 
  ylab("") + xlab("scaled 1-Wasserstein distance") + #xlab("calls to log density and gradient") + 
  theme_bw(base_size = 19)+
  theme(legend.position="top", legend.title = element_blank()) 
print(p_box_compar)
dev.off()

## sensitivity test for I ##
load("../results/wasserstein_100_default_W1_WR.RData")
load("../results/wasserstein_100_sen_I_W1_WR.RData")

ratio_M_wd_I5 <- (apply(W_d_100_pf_I5[, pick_ind], 2, f <- function(x){quantile(x, 0.5)})/
                    apply(W_d_100_pf_IR[, pick_ind], 2, f <- function(x){quantile(x, 0.5)}))

mean(ratio_M_wd_I5)
range(ratio_M_wd_I5)
table(ratio_M_wd_I5 > 2)
table(ratio_M_wd_I5 < 0.5)

mean_M_wd_I5 <- colMeans(W_d_100_pf_I5[, pick_ind]) / 
  colMeans(W_d_100_pf_IR[, pick_ind])
table(mean_M_wd_I5 > 2)
table(mean_M_wd_I5 < 0.5)

ratio_M_wd_I40 <- (apply(W_d_100_pf_I40[, pick_ind], 2, f <- function(x){quantile(x, 0.5)})/
                     apply(W_d_100_pf_IR[, pick_ind], 2, f <- function(x){quantile(x, 0.5)}))

mean(ratio_M_wd_I40)
range(ratio_M_wd_I40)
table(ratio_M_wd_I40 > 2)
table(ratio_M_wd_I40 < 0.5)

mean_M_wd_I40 <- colMeans(W_d_100_pf_I40[, pick_ind]) / 
  colMeans(W_d_100_pf_IR[, pick_ind])
table(mean_M_wd_I40 > 2)
table(mean_M_wd_I40 < 0.5)

w_d_median_pf_IR <- rep(apply(W_d_100_pf_IR[, pick_ind], 2, median), each = M)
apply(W_d_100_pf_I5[, pick_ind] / w_d_median_pf_IR, 2, median)
apply(W_d_100_pf_I40[, pick_ind] / w_d_median_pf_IR, 2, median)

# check variance #
apply(W_d_100_pf_I5[, pick_ind], 2, var) / apply(W_d_100_pf_IR[, pick_ind], 2, var)
mean((apply(W_d_100_pf_I5[, pick_ind], 2, var) / apply(W_d_100_pf_IR[, pick_ind], 2, var))[-3])
mean((apply(W_d_100_pf_I40[, pick_ind], 2, var) / apply(W_d_100_pf_IR[, pick_ind], 2, var))[-3])
which((apply(W_d_100_pf_I5[, pick_ind], 2, var) / apply(W_d_100_pf_I40[, pick_ind], 2, var)) > 4)
summary(apply(W_d_100_pf_I40[, pick_ind], 2, var) / apply(W_d_100_pf_IR[, pick_ind], 2, var))
range(apply(W_d_100_pf_IR[, pick_ind], 2, var)/apply(W_d_100_pf_I40[, pick_ind], 2, var))


w_d_scaled = c(c(W_d_100_pf_I40[, pick_ind] / w_d_median_pf_IR), 
               c(W_d_100_pf_IR[, pick_ind] / w_d_median_pf_IR),
               c(W_d_100_pf_I5[, pick_ind] / w_d_median_pf_IR))
w_d_scaled[w_d_scaled >= 2^11] <- 2^11
df <- data.frame(w_d = w_d_scaled,
                 model = rep(rep(modelnamels, each = M), 3),
                 type = rep(1:3, each = length(model_record[pick_ind])*M),
                 w_d_phI_scaled = rep(rep(w_d_PhI_scaled, each = M), 3))
df$type <- factor(df$type, levels = 1:3, 
                  labels = c("I = 40", "I = 20", "I = 5"))



width <- 12.0
height <- 10.0
setEPS()
postscript("../pics/phI_adapt_default/W_d_box_sens_I.eps",  #_resam_each
           width = width, height = height)
p_box_compar <- df %>%
  ggplot(aes(y = reorder(model, w_d_phI_scaled, FUN = median), #model, # 
             x = w_d, color = type)) + 
  geom_boxplot(outlier.size = 0.1) + 
  scale_colour_manual(values = c("#000000", "#E69F00", "#009E73")) +
  #scale_colour_manual(values=cbbPalette) + 
  scale_x_continuous(trans = 'log2',
                     limits = c(1/4, 16),
                     breaks = c(1/4, 1/2, 1, 2, 4, 
                                8, 16),
                     labels = c("1/4", 
                                "1/2", "1", "2", "4", "8", 
                                "16")) + 
  ylab("") + xlab("scaled 1-Wasserstein distance") + #xlab("calls to log density and gradient") + 
  theme_bw(base_size = 19)+
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

## case studies from posteriordb in Section 3.4 ##
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


## plots for case studies ##
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

# 8 school centered #
## run main_pf with i = 15, seed_list = 1:20, sample with PSIS WR for 100 samples
# run wasserstain_check with i = 15 and generate the plots
# then
## run main_pf.R with init_bound = 15.0 and generate the plots
opt_tr_res <- get_opt_tr(opath)#opath #lp_opath[[15]]$opath
check_dim <- c(8, 10)

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
ggsave("8-school_points.eps", #"8-school_opt_tr22.eps"
       plot = p_check,
       device = cairo_ps,
       path = "../pics/phI_adapt_default/",
       width = 5.0, height = 5.0, units = "in")

p_check1 <- ggplot(dta_check, aes(x=x, y=y) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=-1) +
  scale_x_continuous(expand = c(0, 0), limits = c(-18, 22)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-20, 4)) +
  geom_point(data = dta_opt, 
             aes(x = optim_x, y = optim_y, group = tr_id, 
                 color = optim_ind, alpha = 0.5), size = 1) +
  geom_path(data = dta_opt, 
            aes(x = optim_x, y = optim_y, group = tr_id, 
                color = optim_ind, alpha = 0.5)) +
  scale_color_gradient(low="white", high="orange") +
  xlab("") + ylab("") +
  theme(
    legend.position='none'
  ) 
p_check1

ggsave("8-school_opt_tr.eps", #"8-school_opt_tr22.eps"
       plot = p_check1,
       device = cairo_ps,
       path = "../pics/phI_adapt_default/",
       width = 5.0, height = 5.0, units = "in")


# multimodality #
## run main_pf with i = 3, seed_list = 1:20, sample with PSIS WR for 100 approximate draws
# run wasserstain_check with i = 3 and generate the plots
opt_tr_res <- get_opt_tr(opath)
check_dim <- c(4, 6)  
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
ggsave("3_points.eps",
       plot = p_check,
       device = cairo_ps,
       path = "../pics/phI_adapt_default/",
       width = 5.0, height = 5.0, units = "in")


p_check1 <- ggplot(dta_check, aes(x=x, y=y) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=-1) +
  scale_x_continuous(expand = c(0, 0), limits = c(-1.5, 2.5)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-5, 3)) +
  geom_point(data = dta_opt, 
             aes(x = optim_x, y = optim_y, group = tr_id, 
                 color = optim_ind), size = 1) +
  geom_path(data = dta_opt, 
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
ggsave("3_phI.eps",
       plot = p_phI,
       device = cairo_ps,
       path = "../pics/phI_adapt_default/",
       width = 5.0, height = 5.0, units = "in")


# non-identifiable problem #
## run main_pf with i = 32, seed_list = 1:20, sample with PSIS WOR for 20 inits
## run wasserstain_check with i = 32 
opt_tr_res <- get_opt_tr(opath)
check_dim <- c(45, 46) # c(1, 2) #c(43, 44) #c(45, 46) #  c(3, 4)
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
ggsave("32_points.eps",
       plot = p_check,
       device = cairo_ps,
       path = "../pics/phI_adapt_default/",
       width = 5.0, height = 5.0, units = "in")

p_check1 <- ggplot(dta_check, aes(x=x, y=y) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=-1) +
  scale_x_continuous(expand = c(0, 0), limits = c(-2.2, 4.6)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-7, 1)) +
  geom_point(data = dta_opt,
             aes(x = optim_x, y = optim_y, group = tr_id,
                 color = optim_ind), size = 1) +
  geom_path(data = dta_opt,
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
ggsave("32_phI.eps",
       plot = p_phI,
       device = cairo_ps,
       path = "../pics/phI_adapt_default/",
       width = 5.0, height = 5.0, units = "in")


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
  geom_path(data = dta_opt, 
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
  geom_path(data = dta_opt, 
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
pick_samples <- Imp_Resam_WR(opath, n_sam = 100, seed = 1)

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
  geom_path(data = dta_opt, 
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

# check local hessian of points in high probability mass region
library(numDeriv)
posterior <- to_posterior(model, data)
D <- get_num_upars(posterior)
lp <- function(theta) -log_prob(posterior, theta, adjust_transform = TRUE, 
                                gradient = TRUE)[1]

H_l <- hessian(lp, ref_samples[10, ], method="Richardson")
inv_H_l <- solve(H_l)
diag(inv_H_l)[1:10]
apply(ref_samples[, 1:10], 2, var)
apply(pick_samples[1:10, ], 1, var)
apply((ADVI_meanfield_draw_100[[32]][[10]][, 1:20]), 2, var)

# non-Gaussian problem #
## run main_pf with i = 18, seed_list = 1:20, sample with PSIS WR for 100 samples
## run wasserstain_check with i = 18 

opt_tr_res <- get_opt_tr(opath)
check_dim <- c(2, 3)
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
  geom_path(data = dta_opt,
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


## codes for checking ##
opt_tr_res <- get_opt_tr(opath)
for(first_check in 1:33){
  check_dim <- c(2 * first_check - 1, 2 * first_check)
  cat(check_dim, "\n")
  dta_sam <- data.frame(sam_x = pick_samples[check_dim[1], ],
                        sam_y = pick_samples[check_dim[2], ])
  
  dta_check <- data.frame(x = ref_samples[, check_dim[1]],
                          y = ref_samples[, check_dim[2]])
  
  # dta_phI <- data.frame(x = PhaseI_last_draw[[32]][, check_dim[1]],
  #                       y = PhaseI_last_draw[[32]][, check_dim[2]])
  
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
    geom_path(data = dta_opt,
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
    geom_path(data = dta_opt,
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
  
  # p_phI <- ggplot(dta_check, aes(x=x, y=y) ) +
  #   stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  #   scale_fill_distiller(palette=4, direction=-1) +
  #   scale_x_continuous(expand = c(0, 0))+ #, limits = c(-2.5, 6)) +
  #   scale_y_continuous(expand = c(0, 0))+ #, limits = c(-5, 1)) +
  #   geom_point(data = dta_phI, aes(x=x, y=y), colour="red", size = 3) +
  #   xlab("") + ylab("") +
  #   theme(
  #     legend.position='none'
  #   )
  # p_phI
}


## plots for wasserstein distance illustration ##
## run wasserstain_check with i = 6 
# 0.03430658 vs  vs 0.00878406
pick_samples_center <- ADVI_meanfield_center[[i]]
#pick_samples_center <- t(ADVI_meanfield_draw_100[[i]][[1]])
a_center = wpp(t(pick_samples_center), 
               mass = rep(1 / ncol(pick_samples_center), ncol(pick_samples_center)))
b = wpp(ref_samples, mass = rep(1 / nrow(ref_samples), nrow(ref_samples)))
w_d_c <- wasserstein(a_center, b, p = 1); w_d_c #0.03430658


pick_samples_mf <- t(ADVI_meanfield_draw_100[[i]][[1]])
a_mf = wpp(t(pick_samples_mf), 
               mass = rep(1 / ncol(pick_samples_mf), ncol(pick_samples_mf)))
w_d_mf <- wasserstein(a_mf, b, p = 1); w_d_mf #0.01676967


pick_samples <- lp_opath[[i]]$opath[[1]]$DIV_save$repeat_draws   # pick_samples
a = wpp(t(pick_samples), 
        mass = rep(1 / ncol(pick_samples), ncol(pick_samples)))
b = wpp(ref_samples, mass = rep(1 / nrow(ref_samples), nrow(ref_samples)))
w_d_pf <- wasserstein(a, b, p = 1); w_d_pf #0.00878406

check_dim <- c(1, 2)
dta_sam <- data.frame(sam_x = pick_samples[check_dim[1], ],
                      sam_y = pick_samples[check_dim[2], ])
dta_sam_center <- data.frame(sam_x = pick_samples_center[check_dim[1], ],
                             sam_y = pick_samples_center[check_dim[2], ])

dta_sam_mf <- data.frame(sam_x = pick_samples_mf[check_dim[1], ],
                             sam_y = pick_samples_mf[check_dim[2], ])

dta_check <- data.frame(x = ref_samples[, check_dim[1]],
                        y = ref_samples[, check_dim[2]])


p_check <- ggplot(dta_check, aes(x=x, y=y) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=-1) +
  scale_x_continuous(expand = c(0, 0), limits = c(-0.45, -0.25)) +
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
  scale_x_continuous(expand = c(0, 0), limits = c(-0.45, -0.25)) + 
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

p_check_mf <- ggplot(dta_check, aes(x=x, y=y) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=-1) +
  scale_x_continuous(expand = c(0, 0), limits = c(-0.45, -0.25)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-0.03, 0.16)) + 
  geom_point(data = dta_sam_mf, aes(x=sam_x, y=sam_y),
             colour="red", size = 3) +
  xlab("") + ylab("") +
  theme(
    legend.position='none'
  )

p_check_mf
ggsave("6-dog_mf.eps",
       plot = p_check_mf,
       device = cairo_ps,
       path = "../pics/phI_adapt_default/",
       width = 5.0, height = 5.0, units = "in")


### plot of ELBO ###
library(ggallin)
load("../results/ADVI_100_updat.RData")
load("../results/ELBO_pf.RData")
M = 100

# ELBO #

ELBO_pf_median = apply(ELBO_pf[, pick_ind], 2, f <- function(x)(median(x, na.rm = TRUE)))
ELBO_pf_median[6] = median(ELBO_pf[, 6][is.finite(ELBO_pf[, 6])]) # correct the -Inf bug in Stan code for model dogs-dogs_log
ELBO_ADVI_mf_median = apply(ELBO_ADVI_mf[, pick_ind], 2, 
                           f <- function(x)(median(x, na.rm = TRUE)))
ELBO_ADVI_fr_median = apply(ELBO_ADVI_fr[, pick_ind], 2, 
                           f <- function(x)(median(x, na.rm = TRUE)))

plot(ELBO_pf_median[-2], ELBO_ADVI_mf_median[-2])
abline(a = 0, b = 1)

ELBO_pf_shift = ELBO_pf[, pick_ind] - rep(ELBO_pf_median, each = M)
ELBO_ADVI_mf_shift = ELBO_ADVI_mf[, pick_ind] - rep(ELBO_pf_median, each = M)
ELBO_ADVI_fr_shift = ELBO_ADVI_fr[, pick_ind] - rep(ELBO_pf_median, each = M)
ELBO = c(c(ELBO_ADVI_fr_shift), c(ELBO_ADVI_mf_shift), c(ELBO_pf_shift))
ELBO[ELBO >= 100] <- 100
ELBO[ELBO <= -10000] <- -10000
na_id <- which(is.na(c(ELBO)))

df <- data.frame(ELBO = ELBO[-na_id],
                 model = rep(rep(modelnamels, each = M), 3)[-na_id],
                 type = rep(c(1, 2, 3), each = length(model_record[pick_ind])*M)[-na_id],
                 w_d_phI_scaled = rep(rep(w_d_PhI_scaled, each = M), 3)[-na_id])
df$type <- factor(df$type, levels = 1:3,
                  labels = c("dense ADVI", "mean-field ADVI", "pathfinder"))

width <- 12.0
height <- 10.0
setEPS()
postscript("../pics/phI_adapt_default/ELBO_compar_100_VI.eps",  #_resam_each
           width = width, height = height)
p_box_compar <- df %>%
  ggplot(aes(y = reorder(model, w_d_phI_scaled, FUN = median), #reorder(model, n_leapfrogs, FUN = median), 
             x = ELBO, color = type)) + 
  geom_boxplot(outlier.size = 0.1) +
  scale_x_continuous(trans = pseudolog10_trans,
    limits = c(-10000, 100),
                     breaks = c(-10000, -1000, -100, -10,  0, 10, 100),
                     labels = c("<-10000", "-1000", "-100", "-10", "0", 
                                "10", ">100")) +
  scale_colour_manual(values=c("#000000", "#E69F00", "#56B4E9")) +  #cbbPalette
  ylab("") + xlab("shifted ELBO") + #xlab("calls to log density and gradient") + 
  theme_bw(base_size = 20)+
  theme(legend.position="top", legend.title = element_blank()) 
print(p_box_compar)
dev.off()


#### Laplace approximation ####
w_d_median_pf <- rep(apply(W_d_100_pf[, pick_ind], 2, median), each = 100)

w_d_scaled = c(c(W_d_100_pf[, pick_ind] / w_d_median_pf))

range(w_d_scaled)
df <- data.frame(w_d = w_d_scaled,
                 model = rep(rep(modelnamels, each = 100), 1),
                 type = rep(c("pathfinder"),
                            each = length(model_record[pick_ind])*100),
                 npar = rep(dimension_record, each = 100),
                 w_d_phI_scaled = rep(rep(w_d_PhI_scaled, each = M), 1))


w_d_lp_scaled <- c(w_d_laplace / apply(W_d_100_pf[, pick_ind], 2, median))
range(w_d_lp_scaled, na.rm = T)
a <- ifelse(is.na(w_d_lp_scaled), "red", "black")[order(w_d_PhI_scaled)]


w_d_point = data.frame(w_d_lp_scaled = c(w_d_lp_scaled),
                       model = modelnamels,
                       point = c(#rep("Stan Phase I", length(model_record)),
                         rep("Laplace", length(model_record[pick_ind]))),
                       npar = rep(dimension_record, each = 100))

width <- 12.0
height <- 10.0
setEPS()
postscript("../pics/phI_adapt_default/W_d_box_laplace_pf.eps",  #_resam_each
           width = width, height = height)
p_box_compar <- df %>%
  ggplot(aes(y = reorder(model, w_d_phI_scaled, FUN = median),
             x = w_d, color = type)) +
  geom_boxplot(outlier.size = 0.1) +
  geom_point(aes(y = model, x = w_d_lp_scaled, shape = point),
             color = c("#0072B2"),
             data = w_d_point, size = 3, inherit.aes=FALSE) +
  scale_shape_manual(values=c(18)) +
  #geom_vline(xintercept=1, color = "red", size=0.5)+
  scale_colour_manual(values=c("#56B4E9")) +
  scale_x_continuous(trans = 'log2',
                     limits = c(1/4, 16),
                     breaks = c(1/4, 1/2, 1, 2, 4, 8, 16),
                     labels = c("1/4", "1/2", "1", "2", "4", "8", "16")) +
  ylab("") + xlab("scaled 1-Wasserstein distance") + #xlab("calls to log density and gradient") +
  theme_bw(base_size = 20)+
  theme(axis.text.y = element_text(colour = a)) + 
  theme(legend.position="top", legend.title = element_blank())
print(p_box_compar)
dev.off()

w_d_median_pf_IR <- rep(apply(W_d_100_pf_IR[, pick_ind], 2, median), each = 100)

w_d_scaled = c(c(W_d_100_pf_IR[, pick_ind] / w_d_median_pf_IR))

df <- data.frame(w_d = w_d_scaled,
                 model = rep(rep(modelnamels, each = 100), 1),
                 type = rep(c("multi-path pathfinder"), 
                            each = length(model_record[pick_ind])*100),
                 npar = rep(dimension_record, each = 100),
                 w_d_phI_scaled = rep(rep(w_d_PhI_scaled, each = M), 1))

w_d_lp_scaled <- c(w_d_laplace / apply(W_d_100_pf_IR[, pick_ind], 2, median))
range(w_d_lp_scaled, na.rm = TRUE)
w_d_point = data.frame(w_d_lp_scaled = c(w_d_lp_scaled),
                       model = modelnamels,
                       point = c(rep("Laplace", length(model_record[pick_ind]))),
                       npar = rep(dimension_record, each = 100))

a <- ifelse(is.na(w_d_lp_scaled), "red", "black")[order(w_d_PhI_scaled)]

width <-12.0
height <- 10.0
setEPS()
postscript("../pics/phI_adapt_default/W_d_box_compar_laplace_mpf.eps",  #_resam_each
           width = width, height = height)
p_box_compar <- df %>%
  ggplot(aes(y = reorder(model, w_d_phI_scaled, FUN = median), # model 
             x = w_d, color = type)) + 
  geom_boxplot(outlier.size = 0.1) + 
  geom_point(aes(y = model, #reorder(model, npar, FUN = median), 
                 x = w_d_lp_scaled, shape = point), 
             color = c("#0072B2"),
             data = w_d_point, size = 3, inherit.aes=FALSE) +
  scale_shape_manual(values=c(18)) +
  #geom_vline(xintercept=1, color = "red", size=0.5)+
  scale_colour_manual(values=c("#CC79A7")) + 
  scale_x_continuous(trans = 'log2',
                     limits = c(1/4, 16),
                     breaks = c(1/4, 1/2, 1, 2, 4, 8, 16),
                     labels = c("1/4", "1/2", "1", "2", "4", "8", "16")) + 
  ylab("") + xlab("scaled 1-Wasserstein distance") + #xlab("calls to log density and gradient") + 
  theme_bw(base_size = 20)+
  theme(axis.text.y = element_text(colour = a)) + 
  theme(legend.position="top", legend.title = element_blank()) 
print(p_box_compar)
dev.off()



