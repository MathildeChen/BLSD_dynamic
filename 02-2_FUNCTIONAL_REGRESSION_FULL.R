
# Black Leaf Streak Disease (BLSD) dynamic on banana leaves

# Authors: M. Seidel, M. Chen, J. Avelino, C. Ababdie CIRAD, 2025

# ----------------------------------
# Packages

# > General data management 
library(tidyverse) 

# > Functional analyses
library(refund) 

# -----------------------------------------
# STEP 1
# Data for functional regression

# > Functional sparse response in long format
ydata2 <- ydata %>% rename(".obs"=".id")
unique(ydata2$.obs) # should be 192

# > Temporal index on which disease is measured
yindex2 <- sort(unique(ydata2$.index)) ; yindex2
# -4 -3  0  3  4  7 10 11 14 17 18 21 24 25 28 31 32 35 38 39 42 45 46 49 52 56 59 60 63

# > Functional sparse response in wide format
Y_t_pivot <- ydata %>% 
  arrange(.index) %>% 
  pivot_wider(names_from = ".index", values_from = ".value") %>% 
  arrange(.id)

# > Merge with X to have the same order between X and Y_t
X <- X_init %>% 
  left_join(., Y_t_pivot, by=c("ID"=".id")) %>%
  # ZONE (ref = ZONE 6, ie the zone with the lowest damaged surface on average)
  mutate(COMP_1 = if_else(compartiment == "1", 1, 0),
         COMP_2 = if_else(compartiment == "2", 1, 0),
         COMP_3 = if_else(compartiment == "3", 1, 0),
         COMP_4 = if_else(compartiment == "4", 1, 0),
         COMP_5 = if_else(compartiment == "5", 1, 0),
         COMP_7 = if_else(compartiment == "7", 1, 0),
         COMP_8 = if_else(compartiment == "8", 1, 0)) %>% 
  mutate(compartment = as.numeric(as.character(compartiment)) - 1) %>% 
  # TREE*LEAF
  mutate(T1L1 = if_else(feuille_unique == "B1F1", 1, 0),
         T1L2 = if_else(feuille_unique == "B1F2", 1, 0),
         T1L3 = if_else(feuille_unique == "B1F3", 1, 0),
         T2L1 = if_else(feuille_unique == "B2F1", 1, 0),
         T2L2 = if_else(feuille_unique == "B2F2", 1, 0),
         T2L3 = if_else(feuille_unique == "B2F3", 1, 0),
         T3L1 = if_else(feuille_unique == "B3F1", 1, 0),
         T3L2 = if_else(feuille_unique == "B3F2", 1, 0),
         T3L3 = if_else(feuille_unique == "B3F3", 1, 0)) %>%
  # CLUSTERS FROM KMEANS
  left_join(., X_scores1 %>% 
              dplyr::select(fragment_unique, starts_with("km")), 
            by="fragment_unique") %>%
  mutate(km_4 = factor(km_4, levels = c("2", "1", "3", "4"))) %>% 
  mutate(KM4_1 = if_else(km_4==1, 1, 0),
         KM4_2 = if_else(km_4==2, 1, 0),
         KM4_3 = if_else(km_4==3, 1, 0),
         KM4_4 = if_else(km_4==4, 1, 0)) 

# > Merge both in a list 
data  <- NULL

# functional observation, but with missing data 
data$X <- X %>% dplyr::select(-bananier, -fragment_unique, -fragment, -feuille, -compartiment, -feuille_unique, -starts_with("T"), -ID, -nblesions1, -compartment, -starts_with("COMP_"), -starts_with("km_"), -starts_with("KM"))

# other covariates
# compartment
data$COMP_1 <- X$COMP_1
data$COMP_2 <- X$COMP_2
data$COMP_3 <- X$COMP_3
data$COMP_4 <- X$COMP_4
data$COMP_5 <- X$COMP_5
data$COMP_7 <- X$COMP_7
data$COMP_8 <- X$COMP_8
data$compartment <- X$compartment
data$compartment <- factor(data$compartment)

# tree leaf 
data$T1L1 <- X$T1L1
data$T1L2 <- X$T1L2
data$T1L3 <- X$T1L3
data$T2L1 <- X$T2L1
data$T2L2 <- X$T2L2
data$T2L3 <- X$T2L3
data$T3L1 <- X$T3L1
data$T3L2 <- X$T3L2
data$T3L3 <- X$T3L3
data$tree_leaf <- X$feuille_unique
data$tree_leaf <- factor(data$tree_leaf)

data$KM4_1 <- X$KM4_1
data$KM4_2 <- X$KM4_2
data$KM4_3 <- X$KM4_3
data$KM4_4 <- X$KM4_4
data$km_4 <- factor(X$km_4)

# nb of lesions at the first date
data$nblesions1 <- as.numeric(as.character(X$nblesions1))

# nb de lesions
#data$nblesions <- as.numeric(as.character(X$nblesions))

data <- as.data.frame(data)

str(data)

# -----------------------------------------
# STEP 2
# Function-on-scalar models fitting 

# Model with no adjustment 
mod_0 <- pffr(Y_t ~ 1, 
              data = data, 
              ydata = ydata2,
              bs.yindex = list(bs = "ps", k=27), # k=27 -> length(yindex2)-2
              bs.int = list(bs = "ps", k=27))    # k=27 -> length(yindex2)-2

summary(mod_0) ; AIC(mod_0)
# R-sq.(adj) = 0.729  ; REML score: 21735 ; AIC = 43483.38  

# adjusted for tree * leaf (as fixed effects)
mod_0_bl <- pffr(Y_t ~ T1L2 + T1L3 + T2L1 + T2L2 + T2L3 + T3L1 + T3L2 + T3L3, 
                 data = data, 
                 ydata = ydata2,
                 bs.yindex = list(bs = "ps", k=27),
                 bs.int = list(bs = "ps", k=27))

summary(mod_0_bl) ; AIC(mod_0_bl)
# R-sq.(adj) =  0.825 ; REML score =  21059 ; AIC = 42349.05

# adjusted for cluster and tree * leaf (as fixed effects) 
mod_0_blc4 <- pffr(Y_t ~ KM4_4 + KM4_2 + KM4_3 + 
                     T1L2 + T1L3 + T2L1 + T2L2 + T2L3 + T3L1 + T3L2 + T3L3, 
                   data = data, 
                   ydata = ydata2,
                   bs.yindex = list(bs = "ps", k=27),
                   bs.int = list(bs = "ps", k=27))

summary(mod_0_blc4) ; AIC(mod_0_blc4)
# R-sq.(adj) =  0.952 ; REML score =  19351 ; AIC = 38898.56

# adjusted for cluster, nb of lesions at the first date, and tree * leaf (as fixed effects) 
mod_0_blc4l <- pffr(Y_t ~ KM4_4 + KM4_2 + KM4_3 + 
                      nblesions1 + 
                      T1L2 + T1L3 + T2L1 + T2L2 + T2L3 + T3L1 + T3L2 + T3L3, 
                    data = data, 
                    ydata = ydata2,
                    bs.yindex = list(bs = "ps", k=27),
                    bs.int = list(bs = "ps", k=27))

summary(mod_0_blc4l) ; AIC(mod_0_blc4l)
# R-sq.(adj) =  0.955  ; REML score =  19257 ; AIC = 38709.55

# Similar models but with tree*leaf as random effects
mod_1_bl <- pffr(Y_t ~ 0 + c(0) +      # dropping the constant and time-varying intercepts via 0 + c(0) to estimate zone + c(zone) 
                   s(tree_leaf, bs="re"),   # leaf*tree-specific smooth residuals
                 data = data, 
                 ydata = ydata2,
                 bs.yindex = list(bs = "ps", k=5),
                 bs.int = list(bs = "ps", k=20)) 

summary(mod_1_bl) ; AIC(mod_1_bl)
# R-sq. (adj) = 0.824 ; REML score =  21257 ; AIC = 42340.86

# including clusters 
mod_1_blc4 <- pffr(Y_t ~ 0 + c(0) +             # dropping the constant and time-varying intercepts via 0 + c(0) to estimate zone + c(zone) 
                     km_4 + c(km_4) +           # time-varying cluster effects not centered at zero 
                     s(tree_leaf, bs="re"),     # leaf*tree-specific smooth residuals
                   data = data, 
                   ydata = ydata2,
                   bs.yindex = list(bs = "ps", k=5),
                   bs.int = list(bs = "ps", k=20)) 

summary(mod_1_blc4) ; AIC(mod_1_blc4)
# R-sq.(adj) =  0.946 ; REML score =  19582 ; AIC = 39140.37

# additionally adjusted for the initial number of lesions
mod_1_blc4l <- pffr(Y_t ~ 0 + c(0) +             # dropping the constant and time-varying intercepts via 0 + c(0) to estimate zone + c(zone) 
                      km_4 + c(km_4) +           # time-varying cluster effects not centered at zero 
                      nblesions1 +               # number of initial lesions as potential confounding factor
                      s(tree_leaf, bs="re"),     # leaf*tree-specific smooth residuals
                    data = data, 
                    ydata = ydata2,
                    bs.yindex = list(bs = "ps", k=5),
                    bs.int = list(bs = "ps", k=20)) 

summary(mod_1_blc4l) ; AIC(mod_1_blc4l)
# R-sq.(adj) =  0.949 ; REML score =  19501 ; AIC = 38976.56 

# Extraction of the coefficients for the clusters
pa <- list("mod_1_blc4" = mod_1_blc4,
           "mod_1_blc4l"= mod_1_blc4l) %>% 
  map_dfr(., ~{
    
    coef(.x, seWithMean = FALSE, useVc = FALSE)$smterms %>%
      map_dfr(., ~{
        
        .x$coef %>% 
          dplyr::select("xi" = "yindex.vec", "mean"="value", se) %>%
          # > Compute 95% IC (under normality hypothese)
          mutate(ic_up   = mean + 1.96*se,
                 ic_down = mean - 1.96*se) %>% 
          # > Compute significativity range of the coefficient
          mutate(sign = if_else(ic_up*ic_down > 0, "Significativ", "No significativ"))
        
        
      }, .id = "variable") %>% 
      mutate(AIC = as.numeric(as.character(.x$aic)),
             REML = as.numeric(as.character(.x$gcv.ubre[[1]])),
             R2 = as.numeric(as.character(summary(.x)[10]))) %>% 
      mutate(lab = paste0("AIC: ", round(AIC,digits = 1), ";\nREML: ", round(REML, digits = 1), ";\nR²: ", round(R2, digits = 3)))
    
  }, .id = "model") %>% 
  filter(variable %in% c("km_31(yindex)","km_32(yindex)","km_33(yindex)",
                         "km_41(yindex)","km_42(yindex)","km_43(yindex)","km_44(yindex)",
                         "km_51(yindex)","km_52(yindex)","km_53(yindex)","km_54(yindex)","km_55(yindex)")) %>%
  mutate(adjustment_lesions = if_else(model %in% c("mod_1_blc3", "mod_1_blc4", "mod_1_blc5"), 
                                      "Model accounting for\nthe initial number of lesions", 
                                      "Model not accounting for\nthe initial number of lesions")) %>%
  mutate(nb_cluster = "Fragments partitioned into\n4 clusters") %>% 
  mutate(cluster = case_when(
    variable %in% c("km_41(yindex)") ~ "Cluster 1",
    variable %in% c("km_42(yindex)") ~ "Cluster 2",
    variable %in% c("km_43(yindex)") ~ "Cluster 3",
    variable %in% c("km_44(yindex)") ~ "Cluster 4"
  )) %>%
  #mutate(variable = substr(variable, 1, 5)) %>% 
  ggplot(., aes(x = xi, y = mean, color=cluster)) +
  geom_hline(yintercept = 0, color="black", lty=2) +
  geom_text(aes(x = 0, y = 2000, label = lab), check_overlap = T, color = "black", nudge_x = -4, hjust = 0) +
  geom_path(linewidth = 1) +
  geom_ribbon(aes(ymin=ic_down, ymax=ic_up, fill=cluster), color="transparent", alpha=0.3) +
  facet_wrap(adjustment_lesions~.) +
  scale_color_manual(values=pal_cluster, name="") +
  scale_fill_manual(values=pal_cluster, name="") +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Days after rank 4", y = "Mean functional coefficient\n(95% confidence interval)")



pred_mod_1_blc4 <- predict(mod_1_blc4, se.fit = TRUE, type = "response",
                           newdata = data.frame(km_4 = c(1,2,3, 4), 
                                                tree_leaf=c("B1F1","B1F1","B1F1","B1F1")))

pred_mod_1_blc4l <- predict(mod_1_blc4l, se.fit = TRUE, type = "response",
                            newdata = data.frame(km_4 = c(1,2,3, 4), 
                                                 nblesions1 = rep(mean(data$nblesions1), 4),
                                                 tree_leaf=c("B1F1","B1F1","B1F1","B1F1")))


estimated_mean <- data.frame(
  model = c(rep("Fragments partitioned into\n4 clusters", 29*4), 
            rep("Fragments partitioned into\n4 clusters", 29*4)),
  adjustment = c(rep("Model not accounting for\nthe initial number of  lesions", 29*4),
                 rep("Model accounting for\nthe initial number of  lesions", 29*4)),
  cluster = c(rep("1", 29),  rep("2", 29),  rep("3", 29), rep("4", 29),
              rep("1", 29),  rep("2", 29),  rep("3", 29), rep("4", 29)), 
  x_i     = rep(yindex2, 8),
  f_i     = c(pred_mod_1_blc4$fit[1,],  pred_mod_1_blc4$fit[2,],  pred_mod_1_blc4$fit[3,],  pred_mod_1_blc4$fit[4,],
              pred_mod_1_blc4l$fit[1,],  pred_mod_1_blc4l$fit[2,],  pred_mod_1_blc4l$fit[3,],  pred_mod_1_blc4l$fit[4,]),
  se_i    = c(pred_mod_1_blc4$se.fit[1,],  pred_mod_1_blc4$se.fit[2,],  pred_mod_1_blc4$se.fit[3,],  pred_mod_1_blc4$se.fit[4,],
              pred_mod_1_blc4l$se.fit[1,],  pred_mod_1_blc4l$se.fit[2,],  pred_mod_1_blc4l$se.fit[3,],  pred_mod_1_blc4l$se.fit[4,])) %>% 
  # Compute 95% CI
  mutate(
    ic_up   = f_i + 1.96*se_i,
    ic_down = f_i - 1.96*se_i) %>% 
  mutate(ic_up = if_else(ic_up>3600, 3600, ic_up),
         ic_up = if_else(ic_up<0, 0, ic_up),
         ic_down = if_else(ic_down>3600, 3600, ic_down),
         ic_down = if_else(ic_down<0, 0, ic_down)) 



pb <- estimated_mean %>%
  ggplot(., aes(x = x_i, color=cluster)) +
  geom_hline(yintercept = c(0, 3512), color="black", lty=2) +
  # > mean 
  geom_path(data = Fit.mu1, 
            aes(x=d, y=mu), color="black", lwd=1) +
  geom_line(aes(y=f_i), linewidth=1.2) +
  geom_ribbon(aes(ymin=ic_down, ymax=ic_up, fill=cluster), color="transparent", alpha=0.3) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x="Days after rank 4", y="Mean infected surface in mm²\n95% confidence interval) ") +
  facet_wrap(adjustment~.) +
  guides(color = guide_legend(nrow=1)) +
  lims(y=c(-100, 3600)) +
  scale_color_manual(values=pal_cluster, name="") +
  scale_fill_manual(values=pal_cluster, name="") ; pb


plot_grid(pb+theme(legend.position = "none")+ggtitle("a."), pa+ggtitle("b."), nrow=2, rel_heights = c(0.43, 0.57))






