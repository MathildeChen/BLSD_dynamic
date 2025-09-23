
# Black Leaf Streak Disease (BLSD) dynamic on banana leaves

# Authors: M. Seidel, M. Chen, J. Avelino, C. Ababdie CIRAD, 2025

# ----------------------------------
# Packages

# > General data management 
library(tidyverse) 

# > Functional analyses
library(refund)  # functional PCA and regressions

# > Clustering
library(cluster)
library(FactoMineR)
library(factoextra)
library(corrplot)

# > Mixed models 
library(emmeans)
library(ggpubr)
library(lme4)

# > Graphics
library(cowplot)                   # add several panels in a plot
library(ggrepel) ; library(ggtext) # deal with text 
library(wesanderson)               # generate palettes for clustering 

# > Tables 
library(arsenal)    # contingency tables 
library(multcomp)   # pairwise comparison, tukey tests
library(kableExtra) # print tables in htlm

# ----------------------------------
# Path to project
path_project <- "D:/Mes Donnees/CIRAD/COLLABORATIONS/2024_MARINE/"


# ----------------------------------
# Load home-made functions
source(paste0(path_project, "01_SCRIPT/00_FUNCTIONS.R"))

# ----------------------------------
# Data loading 

# to produce the data, run the following script:
#source(paste0(path_project, "01_SCRIPT/00_DATA_IMPUTATION.R"))

# or directly load: 
load(paste0(path_project, "00_DONNEES/data_fragment_for_analyses.rda"))

# It contains a dataset ready for functional PCA:
ydata <- data_fragment_for_analyses$ydata ; str(ydata)

#tibble [2,730 × 3] (S3: tbl_df/tbl/data.frame)
# $ .id   : int [1:2730] 1 1 1 1 1 1 1 1 1 1 ...
# $ .index: num [1:2730] 3 7 14 17 21 28 31 35 38 42 ...
# $ .value: num [1:2730] 0 42 603 622 2533 ...

yindex <- data_fragment_for_analyses$yindex ; yindex
# [1] -4 -3  0  3  4  7 10 11 14 17 18 21 24 25 28 31 32 35 38 39 42 45 46 49 52 56 59 60 63

# It also contains a dataset with the different cofactors: 
X_init <- data_fragment_for_analyses$X_init ; str(X_init)
#tibble [192 × 9] (S3: tbl_df/tbl/data.frame)
# $ bananier       : Factor w/ 3 levels "Tree 1","Tree 2",..: 1 1 1 1 1 1 1 1 1 1 ...
# $ feuille_unique : Factor w/ 9 levels "B1F1","B1F2",..: 1 1 1 1 1 1 1 1 1 1 ...
# $ compartiment   : Factor w/ 8 levels "1","2","3","4",..: 1 1 2 2 2 2 3 3 4 4 ...
# $ fragment       : num [1:192] 1 2 3 4 5 6 7 8 10 11 ...
# $ fragment_unique: Factor w/ 192 levels "P1H1z1","P1H1z10",..: 1 12 16 17 18 19 20 21 2 3 ...
# $ ID             : int [1:192] 1 2 3 4 5 6 7 8 9 10 ...
# $ feuille        : Factor w/ 3 levels "Leaf 1","Leaf 2",..: 1 1 1 1 1 1 1 1 1 1 ...
# $ date1          : Date[1:192], format: "2012-10-08" "2012-10-08" "2012-10-08" ...
# $ nblesions1     : int [1:192] 21 60 20 9 46 36 42 79 20 12 ...

# ----------------------------------
# a. Functional principal component analysis (FPCA)

# > prior the analysis, it is required to optimally set the FPCA
#source(paste0(path_project, "01_SCRIPT/01_OPTIM_FUNCTIONAL_PCA.R"))

# > perform the functional PCA
Fit.MM1 = fpca.sc(ydata  = ydata,  # data
                  nbasis = 9,
                  npc    = 3, 
                  simul  = T,
                  var    = TRUE,   # model-based estimates for the variance should be computed
                  center = TRUE)   # centering the data by retrieving an estimated mean function

# > estimated mean function (overal dynamic over the fragments)
Fit.mu1 = data.frame(mu = Fit.MM1$mu,
                     d  = sort(unique(ydata$.index)))

# > corresponding eigenfunctions 
Fit.basis1 = data.frame(phi = Fit.MM1$efunctions,
                        d   = sort(unique(ydata$.index)))


# > scores derivated from the functional PCs
X_scores1 <- Fit.MM1$scores %>% 
  as.data.frame(.) %>%
  rename("pca_score1"=1, "pca_score2"=2, "pca_score3"=3) %>%
  cbind(X_init, .) 


# ------------------------------
# b. Clustering on the scores derived from FPCA

# > Scores from FPCA
fpca_scores <- X_scores1 %>% 
  dplyr::select(starts_with("pca_score"))

summary(fpca_scores) # not centered and standardized

z_fpca_scores <- fpca_scores %>% 
  scale(center=T) %>% 
  as.data.frame(.)

summary(z_fpca_scores) # centered (mean = 0) and standardized (sd = 1)
sd(z_fpca_scores$pca_score1)
sd(z_fpca_scores$pca_score2)
sd(z_fpca_scores$pca_score3)

# > Perform clustering for K=4
set.seed(123)
KM4 <- kmeans(z_fpca_scores,  centers = 4, nstart = 25) 

# > Add in initial data
X_scores1$km_4 <- KM4$cluster
X_scores1$km_4 <- as.factor(X_scores1$km_4)

# Change labels 
X_scores1 <- X_scores1 %>%
  mutate(km_4 = case_when(
    km_4 == "1" ~ "2",
    km_4 == "2" ~ "1",
    TRUE ~ km_4
  ))

# Chracteristics of the clusters
extract_km_results(KM4)

table(X_scores1$km_4)
#  1  2  3  4 
# 30 62 72 28 

table(X_scores1$bananier, X_scores1$km_4)
#           1  2  3  4
#   Tree 1 26 12 23  5
#   Tree 2 30 13 19  0
#   Tree 3  6  5 30 23

table(X_scores1$feuille_unique, X_scores1$km_4)
#           1  2  3  4
#   Tree 1 12 26 23  5
#   Tree 2 13 30 19  0
#   Tree 3  5  6 30 23

tableby(km_4 ~ pca_score1 + pca_score2 +pca_score3 + bananier + feuille + compartiment + nblesions1, data = X_scores1 %>% mutate(compartiment=if_else(compartiment=="8", "5", compartiment))) %>% 
  summary(total=F, text=T,  pfootnote = T, digits = 2)


# -----------------------------------------
# c1. Description in terms of FPCA scores and number of inital lesiosn

# by FPCA scores
# Tukey tests
tab_tukey <- X_scores1 %>% 
  gather(key=FPC, value=score, starts_with("pca_score")) %>% 
  group_by(km_4, FPC) %>% 
  mutate(mean_score = round(mean(score),1),
         sd_score = round(sd(score), 1)) %>% 
  left_join(., 
            rbind(
              lazy_tukey(pca_score = "pca_score1") %>% mutate(FPC="pca_score1"),
              lazy_tukey(pca_score = "pca_score2") %>% mutate(FPC="pca_score2"),
              lazy_tukey(pca_score = "pca_score3") %>% mutate(FPC="pca_score3")), 
            by=c("km_4", "FPC")) %>% 
  mutate(lab = paste0(mean_score, "\n(", sd_score, ")")) %>%
  mutate(FPC=recode(FPC, "pca_score1"="FPC score 1", "pca_score2"="FPC score 2", "pca_score3"="FPC score 3")) 

# by the number of first lesions
# Table format
data_fragment_first_lesions <- X_scores1
data_fragment_first_lesions <- data_fragment_first_lesions %>% relocate("feuille", .after = "bananier")

# Examine the number of first lesions VS data-driven clusters using a negative binomial glmm

# > more details on the choice of the model are provided in the script: 02-1_GLMM_FULL.R

mod_clust6 <- glmer.nb(nblesions1 ~ km_4 + (1|bananier/feuille), 
                       data = data_fragment_first_lesions)

# AIC 
AIC(mod_clust6)

# Check qqplot
ggqqplot(residuals(mod_clust6))

# Check over dispersion
chi2_6 <- sum(residuals(mod_clust6, type = "pearson")^2)
pvalue_6 <- 1 - pchisq(chi2_6, df = df.residual(mod_clust6)) 
coeff_6 <- chi2_6 / df.residual(mod_clust6)

chi2_6
pvalue_6 # p >0.05 = no over dispersion of the residuals of the model
coeff_6 # coefficient of dispersion ~ 1

# Estimated marginal means on the model 6
# check the difference between the clusters 
EMMmod_clust6 <- emmeans(mod_clust6, ~ km_4)
pairs(EMMmod_clust6)

tuk.cld <- cld(EMMmod_clust6, Letters = letters, reversed = TRUE)
tuk.cld
# Results : 
# km_4 emmean    SE  df asymp.LCL asymp.UCL .group
# 2      4.14 0.209 Inf      3.73      4.55  a    
# 1      3.49 0.187 Inf      3.12      3.86   b   
# 4      2.88 0.239 Inf      2.41      3.35    c  
# 3      2.83 0.183 Inf      2.47      3.19    c 

# Attribute to each fragment / cluster the corresponding letter 
data_fragment_first_lesions <- fonction_plot_ameans()

# -----------------------------------------
# c2. Functional regression 

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

# STEP 2
# Function-on-scalar models fitting 

# Fit a function on scalar regression to examine the differences in terms of dynamics between 
# the four clusters of fragment
# the model accounts for leaf (as random effect) 

# > Check the 02-2_FUNCTIONAL_REGRESSION_FULL.R script to see details on other functional models fit 

mod_1_blc4 <- pffr(Y_t ~ 0 + c(0) +             # dropping the constant and time-varying intercepts via 0 + c(0) to estimate zone + c(zone) 
                     km_4 + c(km_4) +           # time-varying cluster effects not centered at zero 
                     s(tree_leaf, bs="re"),     # leaf*tree-specific smooth residuals
                   data = data, 
                   ydata = ydata2,
                   bs.yindex = list(bs = "ps", k=5),
                   bs.int = list(bs = "ps", k=20)) 

summary(mod_1_blc4) ; AIC(mod_1_blc4)
# R-sq.(adj) =  0.946 ; REML score =  19582 ; AIC = 39140.37

# The estimated mean dynamic per cluster, according to the model 
pred_mod_1_blc4 <- predict(mod_1_blc4, se.fit = TRUE, type = "response",
                           newdata = data.frame(km_4 = c(1,2,3, 4), 
                                                tree_leaf=c("B1F1","B1F1","B1F1","B1F1")))

estimated_mean <- data.frame(
  model = rep("Fragments partitioned into\n4 clusters", 29*4),
  adjustment = rep("Model not accounting for\nthe initial number of  lesions", 29*4),
  cluster = c(rep("1", 29),  rep("2", 29),  rep("3", 29), rep("4", 29)), 
  x_i     = rep(yindex2, 4),
  f_i     = c(pred_mod_1_blc4$fit[1,],  pred_mod_1_blc4$fit[2,],  pred_mod_1_blc4$fit[3,],  pred_mod_1_blc4$fit[4,]),
  se_i    = c(pred_mod_1_blc4$se.fit[1,],  pred_mod_1_blc4$se.fit[2,],  pred_mod_1_blc4$se.fit[3,],  pred_mod_1_blc4$se.fit[4,])) %>% 
  # Compute 95% CI
  mutate(
    ic_up   = f_i + 1.96*se_i,
    ic_down = f_i - 1.96*se_i) %>% 
  mutate(ic_up = if_else(ic_up>3600, 3600, ic_up),
         ic_up = if_else(ic_up<0, 0, ic_up),
         ic_down = if_else(ic_down>3600, 3600, ic_down),
         ic_down = if_else(ic_down<0, 0, ic_down)) 

# Extraction of the functional coefficients for the clusters
tab_coef <- list("mod_1_blc4"=  mod_1_blc4) %>% 
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
  filter(variable %in% c("km_41(yindex)","km_42(yindex)","km_43(yindex)","km_44(yindex)")) %>%
  mutate(adjustment_lesions = "Model accounting for\nthe initial number of lesions") %>%
  mutate(nb_cluster = "Fragments partitioned into\n4 clusters") %>% 
  mutate(cluster = case_when(
    variable %in% c("km_41(yindex)") ~ "Cluster 1",
    variable %in% c("km_42(yindex)") ~ "Cluster 2",
    variable %in% c("km_43(yindex)") ~ "Cluster 3",
    variable %in% c("km_44(yindex)") ~ "Cluster 4"
  )) 



# -----------------------------------------
# c3. Correspondance analyses 

# Contingency tables using Fisher exact test 
tableby(km_4 ~ bananier + feuille + compartiment, 
        data=X_scores1, 
        numeric.stats=c("meansd")) %>% 
  summary(., text=TRUE, pfootnote=TRUE, digits = 1, total=F) %>% 
  kbl(.) %>%
  kable_styling(bootstrap_options = "striped", full_width = F)


X_scores1$compartiment_old <- if_else(X_scores1$compartiment=="8", "5", X_scores1$compartiment)
tableby(compartiment_old ~ km_4, 
        data=X_scores1, 
        numeric.stats=c("meansd")) %>% 
  summary(., text=TRUE, pfootnote=TRUE, digits = 1, total=F) %>% 
  kbl(.) %>%
  kable_styling(bootstrap_options = "striped", full_width = F)


# > Fisher exact tests
t <- X_scores1 %>% 
  dplyr::select(bananier, feuille, compartiment_old, fragment, km_4) %>% 
  mutate(compartiment_old = paste0("Comp. ", compartiment_old)) %>% 
  mutate(km_4 = paste0("Cluster ", km_4)) %>% 
  mutate(fragment = paste0("F.", fragment))

t$compartiment_old <- as.factor(t$compartiment_old)
t$km_4 <- as.factor(t$km_4)
t$fragment <- as.factor(t$fragment)
summary(t)

# Test de Fisher
fisher.test(x = X_scores1$km_4, y=X_scores1$compartiment_old, simulate.p.value = T) 
# p-value = 0.0004998

# Correspondance analyses on the predefined compartiments vs clusters

chi2 <- chisq.test(table(t$km_4, t$compartiment_old))

chi2$observed
#          Comp. 1 Comp. 2 Comp. 3 Comp. 4 Comp. 5 Comp. 6 Comp. 7
#Cluster 1      11       2       9       4       1       1       2
#Cluster 2       6      15       2      10      18       5       6
#Cluster 3       0      11       3      19      26       6       7
#Cluster 4       1       8       2       3       7       4       3

chi2$expected
#           Comp. 1 Comp. 2  Comp. 3 Comp. 4   Comp. 5  Comp. 6 Comp. 7
#Cluster 1  2.8125   5.625 2.500000   5.625  8.125000 2.500000  2.8125
#Cluster 2  5.8125  11.625 5.166667  11.625 16.791667 5.166667  5.8125
#Cluster 3  6.7500  13.500 6.000000  13.500 19.500000 6.000000  6.7500
#Cluster 4  2.6250   5.250 2.333333   5.250  7.583333 2.333333  2.6250

#AFC sur compartiments
afc_compartiments <- CA(table(t$km_4, t$compartiment_old)) ; summary(afc_compartiments)
#The chi square of independence between the two variables is equal to 72.05223 (p-value =  2.025249e-08 ).

#Eigenvalues
#                       Dim.1   Dim.2   Dim.3
#Variance               0.324   0.039   0.011
#% of var.             86.470  10.496   3.035
#Cumulative % of var.  86.470  96.965 100.000
#
#Rows
#    Iner*1000     Dim.1     ctr    cos2     Dim.2     ctr    cos2     Dim.3     ctr    cos2  
#1 |    16.939 |  -0.118   1.394   0.267 |   0.146  17.371   0.404 |  -0.131  48.944   0.329 |
#2 |   265.224 |   1.300  81.428   0.996 |  -0.077   2.345   0.003 |   0.021   0.601   0.000 |
#3 |    68.383 |  -0.368  15.631   0.742 |  -0.215  44.017   0.254 |   0.029   2.852   0.005 |
#4 |    24.725 |  -0.186   1.547   0.203 |   0.313  36.267   0.578 |   0.193  47.603   0.219 |
#
#Columns
#    Iner*1000     Dim.1     ctr    cos2     Dim.2     ctr    cos2     Dim.3     ctr    cos2  
#1 |   164.566 |   1.308  49.408   0.974 |   0.095   2.163   0.005 |  -0.190  29.751   0.021 |
#2 |    27.184 |  -0.229   3.041   0.363 |   0.303  43.842   0.635 |  -0.016   0.438   0.002 |
#3 |   106.190 |   1.096  30.869   0.943 |  -0.132   3.701   0.014 |   0.234  40.074   0.043 |
#4 |    20.321 |  -0.172   1.708   0.273 |  -0.280  37.242   0.722 |  -0.024   0.954   0.005 |
#5 |    44.513 |  -0.395  13.001   0.948 |  -0.083   4.736   0.042 |  -0.041   4.054   0.010 |
#6 |    10.916 |  -0.246   1.551   0.461 |   0.193   7.877   0.284 |   0.183  24.407   0.255 |
#7 |     1.581 |  -0.121   0.423   0.867 |   0.043   0.439   0.109 |   0.020   0.322   0.023 |

afc_compartiments$eig
#       eigenvalue percentage of variance cumulative percentage of variance
# dim 1 0.32449688              86.469770                          86.46977
# dim 2 0.03938704              10.495596                          96.96537
# dim 3 0.01138813               3.034634                         100.00000

# Check the other axes (2, 3) instead of (1,2)
plot(afc_compartiments, c(2,3))

# Coordinates
rbind(data.frame(afc_compartiments$col$cos2),
      data.frame(afc_compartiments$row$cos2)) %>% 
  View(.)

# clusters on col (fragments) coordinates 
coords_cols <- rbind(data.frame(afc_compartiments$col$coord),
                     data.frame(afc_compartiments$row$coord))
summary(coords_cols)

# Scale prior clustering
z_coords_cols <- as.data.frame(scale(coords_cols[,1:2]))
summary(z_coords_cols)

# Choice of the number of clusters for the clustering on CA coodinates
plot_grid(
  fviz_nbclust(z_coords_cols, kmeans, method = "wss", k.max = 4)        + ggtitle(label = paste0("kmeans"), subtitle = "Elbow method"),
  fviz_nbclust(z_coords_cols, kmeans, method = "silhouette", k.max = 4) + ggtitle(label = "",               subtitle = "Silhouette method"),
  fviz_nbclust(z_coords_cols, kmeans, method = "gap_stat", k.max = 4, 
               nstart = 5, nboot = 50, verbose = FALSE)     + ggtitle(label = "",               subtitle = "Gap statistic method"),
  ncol = 3)

set.seed(123)
clusters_afc_comprtiments <- kmeans(z_coords_cols, centers = 3, iter.max = 25)

coords_cols$km <- as.vector(clusters_afc_comprtiments$cluster)
coords_cols$lab <- rownames(coords_cols)

coords_cols %>% 
  arrange(km)

coords_row <- data.frame(afc_compartiments$row$coord)
coords_row$lab <- rownames(afc_compartiments$row$coord)

