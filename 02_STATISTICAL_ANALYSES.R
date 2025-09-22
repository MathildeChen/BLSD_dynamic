
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

# > Graphics
library(cowplot)                   # add several panels in a plot
library(ggrepel) ; library(ggtext) # deal with text 
library(wesanderson)               # generate palettes for clustering 

# > Tables 
library(arsenal)    # contingency tables 
library(multcomp)   # pairwise comparison, tukey tests
library(kableExtra) # print tables in htlm

# ----------------------------------
path_project <- "D:/Mes Donnees/CIRAD/COLLABORATIONS/2024_MARINE/"

# Data loading 
load(paste0(path_project, "00_DONNEES/data_fragment_for_analyses.rda"))
ydata <- data_fragment_for_analyses$ydata ; str(ydata)

#tibble [2,730 × 3] (S3: tbl_df/tbl/data.frame)
# $ .id   : int [1:2730] 1 1 1 1 1 1 1 1 1 1 ...
# $ .index: num [1:2730] 3 7 14 17 21 28 31 35 38 42 ...
# $ .value: num [1:2730] 0 42 603 622 2533 ...

yindex <- data_fragment_for_analyses$yindex ; yindex
# [1] -4 -3  0  3  4  7 10 11 14 17 18 21 24 25 28 31 32 35 38 39 42 45 46 49 52 56 59 60 63

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


# Variation in the disease progressions

ydata %>%
  left_join(., X_init, by=c(".id"="ID")) %>% 
  group_by(feuille_unique, .index) %>% 
  mutate(mean = mean(.value, na.rm=T),
         sd = sd(.value, na.rm=T)) %>% 
  mutate(bananier = recode(bananier, "Tree 1"="Plant 1", "Tree 2"="Plant 2", "Tree 3"="Plant 3")) %>%
  mutate(lab = case_when(
    bananier == "Plant 1" & feuille == "Leaf 1" ~ "a. Leaf 1 - Plant 1",
    bananier == "Plant 1" & feuille == "Leaf 2" ~ "b. Leaf 2 - Plant 1",
    bananier == "Plant 1" & feuille == "Leaf 3" ~ "c. Leaf 3 - Plant 1",
    bananier == "Plant 2" & feuille == "Leaf 1" ~ "d. Leaf 1 - Plant 2",
    bananier == "Plant 2" & feuille == "Leaf 2" ~ "e. Leaf 2 - Plant 2",
    bananier == "Plant 2" & feuille == "Leaf 3" ~ "f. Leaf 3 - Plant 2",
    bananier == "Plant 3" & feuille == "Leaf 1" ~ "g. Leaf 1 - Plant 3",
    bananier == "Plant 3" & feuille == "Leaf 2" ~ "h. Leaf 2 - Plant 3",
    bananier == "Plant 3" & feuille == "Leaf 3" ~ "i. Leaf 3 - Plant 3"
  )) %>% 
  ggplot() +
  geom_line(aes(x=.index, y=.value, group=.id), linetype=2) +
  geom_line(aes(x=.index, y=mean, group=feuille_unique), linewidth=1.25, color="blue")+
  geom_ribbon(aes(x=.index, ymin=mean-sd, ymax=mean+sd, group=feuille_unique), fill="blue", alpha=0.2) +
  facet_wrap(. ~ lab, ncol=3) +
  #geom_text(aes(x = -2, y = 3500, label = lab), check_overlap = T, size=5) + 
  labs(x = "Days after leaf reached rank 4", y = "Diseased surface (mm²)") +
  theme_cowplot() +
  theme(panel.border = element_rect(color="black"),
        strip.background = element_blank(),
        panel.grid.major = element_line(color="grey87"))

ggsave(filename = paste0(path_project, "/04_ARTICLE/FIGURES_R1/FIGURE_2.tiff"),
       bg="white", dpi = 600, width=10, height=10)

ggsave(filename = paste0(path_project, "/04_ARTICLE/FIGURES_R1/FIGURE_2.png"),
       bg="white", dpi = 300, width=10, height=10)

# ----------------------------------
# a. Functional principal component analysis (FPCA)
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


# > Plot mean and basis functions
# Which source of variability are represented by each function
Fit.basis1 %>% 
  left_join(Fit.mu1, by="d") %>% 
  gather(key="phi", value=value, -d, -mu) %>% 
  mutate(mu_plus_phi  = mu+5000*value,
         mu_minus_phi = mu-5000*value) %>% 
  gather(key=value_type, value=value, mu_plus_phi, mu_minus_phi) %>% 
  mutate(value_type = dplyr::recode(value_type, 
                             'mu_plus_phi'="Mean + eigenfunction", 
                             'mu_minus_phi' = "Mean - eigenfunction"), 
         value_type = factor(value_type, levels=c("Mean + eigenfunction", "Mean - eigenfunction"))) %>%
  mutate(phi=dplyr::recode(phi, 
                    "phi.1"=paste0("1st functional component (",  round((Fit.MM1$evalues[1]/sum(Fit.MM1$evalues))*100, 1), "%)"),
                    "phi.2"=paste0("2nd functional component (", round((Fit.MM1$evalues[2]/sum(Fit.MM1$evalues))*100, 1), "%)"),
                    "phi.3"=paste0("3rd functional component (",  round((Fit.MM1$evalues[3]/sum(Fit.MM1$evalues))*100, 1), "%)"))) %>%
  ggplot(.) + 
  geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
  geom_hline(yintercept = 0, color="darkgrey", lty=2) +
  geom_path(aes(x = d, y = mu, group=interaction(phi, value_type)), color="black") +
  geom_path(aes(x = d, y = value, color = value_type, group=interaction(phi, value_type))) +
  ggtitle("Which pattern of variability is explained by each principal component", 
          "(and how much does it explain the variability in the full dataset)") +
  facet_wrap(.~phi, nrow=1) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Days after rank 4", 
       y = "Infected area",
       color="") +
  scale_color_manual(values=c("red", "blue"))

ggsave(filename = paste0(path_project, "/04_ARTICLE/FIGURES/FPCA_interpretation_supp_fig_1.png"),
       bg="white", dpi = 300, width=12, height=5)


# Function for plotting for a set of individual fragments
plot_dyn_id <- function(frags=NULL, lab=NULL, 
                        score="pca_score1", score_lab = "1",
                        xlab='Day after rank 4', color_curves="blue")
{
  
  # > if you want all plots
  if(is.null(frags))
  {
    
    frags <- unique(X_init$ID)
    
  }
  
  # > IDs associated 
  frag_fitted <- X_init %>% 
    filter(ID %in% frags) %>% 
    split(.$ID) %>% 
    map_dfr(., ~{
      
      # > Retrieve data for this fragment in each tree and leaf
      # Unique fragment
      EX <- unique(.x$ID)
      
      # Fitted data
      #Yhat_ex     <- retransf(Fit.MM1$Yhat[EX,])
      #diag.var_EX <- retransf(Fit.MM1$diag.var[EX,])
      
      Yhat_ex     <- Fit.MM1$Yhat[EX,]
      diag.var_EX <- Fit.MM1$diag.var[EX,]
      
      data.frame(fitted    = Yhat_ex,
                 ptwise.UB = Yhat_ex + 1.96 * sqrt(diag.var_EX),
                 ptwise.LB = Yhat_ex - 1.96 * sqrt(diag.var_EX),
                 d = sort(unique(ydata$.index)))
      
      
    }, .id=".id")
  
  # Plot
  if(is.null(score))
  {
    
    
    # Obs
    ydata_frag <- ydata %>% 
      filter(.id %in% frags) %>%
      mutate(.id=factor(.id))
    
    # Plot
    plot_frags <- frag_fitted %>% 
      mutate(.id=factor(.id)) %>%
      ggplot(.) +
      # > fitted
      geom_path(aes(x = d, y = fitted, 
                    group=.id), color=color_curves) +
      # > observations
      geom_point(data = ydata_frag, 
                 aes(x = .index, y =.value), 
                 color = color_curves,
                 size=0.75) +
      # > mean function
      geom_path(data = Fit.mu1, 
                aes(x=d, y=mu), color="black", lwd=1.2) +
      labs(x = paste0(xlab), y = 'Diseased surface (mm²)') +
      theme_bw() +
      theme(legend.position="bottom", 
            legend.title.position = "top") +
      ggtitle(paste0(lab))
    
    
  }
  
  if(is.null(score)==F)
  {
    
    # PCA scores
    X_scores_frags <- X_scores1 %>% 
      filter(ID %in% frags) %>% 
      dplyr::select(ID, "pca_score"=starts_with(paste0(score))) %>%
      mutate(ID=as.character(ID)) %>% 
      arrange(pca_score) %>% 
      mutate(order_score = 1:n())
    
    # Obs
    ydata_frag <- ydata %>% 
      filter(.id %in% frags) %>%
      mutate(.id=factor(.id)) %>% 
      #mutate(.value = (3512*.value)/(1-.value)) %>% 
      left_join(X_scores_frags, by=c(".id"="ID"))
    
    # Plot
    plot_frags <- frag_fitted %>% 
      left_join(X_scores_frags, by=c(".id"="ID")) %>% 
      mutate(.id=factor(.id)) %>%
      ggplot(.) +
      # > fitted
      geom_path(aes(x = d, y = fitted, 
                    color = order_score, group=.id)) +
      # > observations
      geom_point(data = ydata_frag, 
                 aes(x = .index, y =.value, 
                     color = order_score),
                 size=0.75) +
      # > mean 
      geom_path(data = Fit.mu1, 
                aes(x=d, y=mu), color="black", lwd=1.2) +
      scale_color_gradientn(#colors=c(viridis::rocket(15)[3:12], 
                            #         viridis::mako(15, direction=-1)[3:12]), 
                            breaks = c(1, 40),
                            colors = viridis::inferno(100, direction = -1),
                            labels = c("Lower score", "Higher score"),
                            guide = guide_colorbar(title = paste0("Order in FPC score ", score_lab),
                                                   barwidth = 7, barheight = 0.5)) +
      labs(x = paste0(xlab), y = 'Diseased surface (mm²)') +
      theme_bw() +
      theme(legend.position="bottom", 
            legend.title.position = "top") +
      ggtitle(paste0(lab))
    
  }
  
  return(plot_frags)
  
}


# Curves for all fragments
#p1 <- plot_dyn_id(frags = unique(X_scores1$ID), score="pca_score1") +
#  geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
#  geom_hline(yintercept = 0, color="darkgrey", lty=2) ; p1
#
#p2 <- plot_dyn_id(frags = unique(X_scores1$ID), score="pca_score2", score_lab = "2") +
#  geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
#  geom_hline(yintercept = 0, color="darkgrey", lty=2) ; p2
#
#p3 <- plot_dyn_id(frags = unique(X_scores1$ID), score="pca_score3", score_lab = "3") +
#  geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
#  geom_hline(yintercept = 0, color="darkgrey", lty=2) ; p3
#
#plot_grid(p1,p2,p3,nrow=1)


# Fragments with the highest and lowest values 
# for each score 

n_frag <- 20 

# > FPC 1
high_pca_score1 <- X_scores1 %>% ungroup() %>% 
  arrange(desc(pca_score1)) %>% 
  head(., n_frag) %>% 
  pull(ID) 

low_pca_score1 <- X_scores1 %>% 
  arrange(pca_score1) %>% 
  head(., n_frag) %>% 
  pull(ID) 

frags_pca_score1 <- c(high_pca_score1, low_pca_score1)

# > FPC 2
high_pca_score2 <- X_scores1 %>% 
  arrange(desc(pca_score2)) %>% 
  head(., n_frag) %>% 
  pull(ID) 

low_pca_score2 <- X_scores1 %>% 
  arrange(pca_score2) %>% 
  head(., n_frag) %>% 
  pull(ID) 

frags_pca_score2 <- c(high_pca_score2, low_pca_score2)

# > FPC 3
high_pca_score3 <- X_scores1 %>% 
  arrange(desc(pca_score3)) %>% 
  head(., n_frag) %>% 
  pull(ID) 

# Unique fragments with lowest score
low_pca_score3 <- X_scores1 %>% 
  arrange(pca_score3) %>% 
  head(., n_frag) %>% 
  pull(ID)

frags_pca_score3 <- c(high_pca_score3, low_pca_score3)

# Plots 
p1 <- plot_dyn_id(frags = frags_pca_score1, score="pca_score1") +
  geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
  geom_hline(yintercept = 0, color="darkgrey", lty=2)

p2 <- plot_dyn_id(frags = frags_pca_score2, 
            score="pca_score2", score_lab = "2") +
  geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
  geom_hline(yintercept = 0, color="darkgrey", lty=2)

p3 <- plot_dyn_id(frags = frags_pca_score3, 
            score="pca_score3", score_lab = "3") +
  geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
  geom_hline(yintercept = 0, color="darkgrey", lty=2)

plot_grid(p1 + ggtitle("a."), 
          p2 + ggtitle("b."),
          p3 + ggtitle("c."),
          nrow=1)

ggsave(filename = paste0(path_project, "/04_ARTICLE/FIGURES/FPCA_interpretation.png"),
       bg="white", dpi = 300, width=12, height=5)

ggsave(filename = paste0(path_project, "/04_ARTICLE/FIGURES_R1/FIGURE_3.tiff"),
       bg="white", dpi = 600, width=12, height=5)

ggsave(filename = paste0(path_project, "/04_ARTICLE/FIGURES_R1/FIGURE_3.png"),
       bg="white", dpi = 300, width=12, height=5)

# plot other scores according to the 1rst one
p2_1 <- plot_dyn_id(frags = frags_pca_score2, 
            score="pca_score1", score_lab = "1") +
  geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
  geom_hline(yintercept = 0, color="darkgrey", lty=2)

p3_1 <- plot_dyn_id(frags = frags_pca_score3, 
            score="pca_score1", score_lab = "1") +
  geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
  geom_hline(yintercept = 0, color="darkgrey", lty=2)

# plot 3rd scores according to the 2nd one
p3_2 <- plot_dyn_id(frags = frags_pca_score3, 
            score="pca_score2", score_lab = "2") +
  geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
  geom_hline(yintercept = 0, color="darkgrey", lty=2)

plot_grid(
  plot_grid(p1 + ggtitle("a."), ggplot()+theme_void(), ggplot()+theme_void(), nrow=1, align="hv"),
  plot_grid(p2, p2_1 + ggtitle("b."), ggplot()+theme_void(), nrow=1, align="hv"),
  plot_grid(p3, p3_1, p3_2 + ggtitle("c."), nrow=1, align="hv"),
  ncol=1)

ggsave(filename = paste0(path_project, "/04_ARTICLE/FIGURES/FPCA_interpretation_supp_fig_2.png"),
       bg="white", dpi = 300, width=12, height=12)



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

# > Correlations between the standardized scores
fpca_scores %>% 
  cor(.) %>% 
  corrplot(., type = "lower", diag = F, method = 'color', addCoef.col = 'black')

# > Choose the optimal number of clusters 
plot_grid(
  fviz_nbclust(z_fpca_scores, kmeans, method = "wss")        + ggtitle(label = paste0("kmeans"), subtitle = "Elbow method"),
  fviz_nbclust(z_fpca_scores, kmeans, method = "silhouette") + ggtitle(label = "",               subtitle = "Silhouette method"),
  fviz_nbclust(z_fpca_scores, kmeans, method = "gap_stat", 
               nstart = 25, nboot = 50, verbose = FALSE)     + ggtitle(label = "",               subtitle = "Gap statistic method"),
  ncol = 3)
# > K=4

# > Perform clustering for K=4, or 5
# 4 functional components
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

extract_km_results <- function(KM)
{
  
  results <- data.frame(
    tot_sum_squares = KM$totss, 
    tot_within_sum_squares = KM$tot.withinss,
    between_sum_squares = KM$betweenss
  )
  
  return(results)
  
}

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

# > Cluster visualization
# > per cluster
pal_cluster <- rev(wesanderson::wes_palette("Zissou1", 4, "continuous"))

# Examine the scores
plot_grid(
  
  # normal
  ggplot(X_scores1, aes(x=pca_score1, y = pca_score2)) + geom_point() + theme_bw() +
    labs(x="FPC score 1", y = "FPC score 2") + ggtitle("a. "),
  
  ggplot(X_scores1, aes(x=pca_score1, y = pca_score3)) + geom_point() + theme_bw() +
    labs(x="FPC score 1", y = "FPC score 3"),
  
  ggplot(X_scores1, aes(x=pca_score2, y = pca_score3)) + geom_point() + theme_bw() +
    labs(x="FPC score 2", y = "FPC score 3"),
  
  # > per tree
  ggplot(X_scores1) +
    geom_point(aes(x=pca_score1, y=pca_score2, color=bananier)) + 
    labs(x="FPC score 1", y="FPC score 2") + ggtitle("b. ") +
    scale_color_manual(values = pals::parula(3)) +
    theme_bw() + theme(legend.position = "none"),
  
  ggplot(X_scores1) +
    geom_point(aes(x=pca_score1, y=pca_score3, color=bananier )) + 
    scale_color_manual(values = pals::parula(3), name="") +
    labs(x="FPC score 1", y="FPC score 3") +
    theme_bw() + theme(legend.position = "none"),
  
  ggplot(X_scores1 %>% 
           mutate(bananier = recode(bananier, "Tree 1"="Plant 1", "Tree 2"="Plant 2", "Tree 3"="Plant 3"))) +
    geom_point(aes(x=pca_score2, y=pca_score3, color=bananier )) + 
    scale_color_manual(values = pals::parula(3), name="") +
    labs(x="FPC score 2", y="FPC score 3") +
    theme_bw(), 
  
  # > per leaf
  ggplot(X_scores1) +
    geom_point(aes(x=pca_score1, y=pca_score2, color=feuille)) + 
    labs(x="FPC score 1", y="FPC score 2") + ggtitle("c. ") +
    scale_color_manual(values = pals::inferno(5)[2:4]) +
    theme_bw() + theme(legend.position = "none"),
  
  ggplot(X_scores1) +
    geom_point(aes(x=pca_score1, y=pca_score3, color=feuille )) + 
    scale_color_manual(values = pals::inferno(5)[2:4]) +
    labs(x="FPC score 1", y="FPC score 3") +
    theme_bw() + theme(legend.position = "none"),
  
  ggplot(X_scores1) +
    geom_point(aes(x=pca_score2, y=pca_score3, color=feuille )) + 
    scale_color_manual(values = pals::inferno(5)[2:4], name="") +
    labs(x="FPC score 2", y="FPC score 3") +
    theme_bw(), 
  
  # per cluster
  ggplot(X_scores1) +
    geom_point(aes(x=pca_score1, y=pca_score2, color=km_4)) + 
    scale_color_manual(values = pal_cluster) +
    labs(x="FPC score 1", y="FPC score 2") + ggtitle("d.") +
    theme_bw() + theme(legend.position = "none"),
  
  ggplot(X_scores1) +
    geom_point(aes(x=pca_score1, y=pca_score3, color=km_4)) + 
    scale_color_manual(values = pal_cluster) +
    labs(x="FPC score 1", y="FPC score 3") + 
    theme_bw() + theme(legend.position = "none"),
  
  ggplot(X_scores1 %>% 
           mutate(km_4 = paste0("Cluster ", km_4))) +
    geom_point(aes(x=pca_score2, y=pca_score3, color=km_4)) + 
    scale_color_manual(values = pal_cluster, name = "") +
    labs(x="FPC score 2", y="FPC score 3") + theme_bw(),
  
  nrow=4, axis = "tbrl", align = "hv"
  
)

ggsave(filename = paste0(path_project, "/04_ARTICLE/FIGURES/Supp_Figures_Scores_tree_leaf.png"),
       bg="white", dpi = 300, width=14, height=10)


# -----------------------------------------
# c1. Description in terms of FPCA scores and number of inital lesiosn

# Home made function to perfom Tuckey tests
lazy_tukey <- function(pca_score=NULL, data=X_scores1)
{
  
  data <- data %>% mutate(km_4 = factor(km_4, levels = c("1", "2", "3", "4")))
  anova_score <- aov(data[,colnames(data)==pca_score] ~ km_4, data=data)
  p.val.anova_score <- summary(anova_score)[[1]][[5]][1]
  aov_sorties <- as.data.frame(predict(anova_score, 
                                       newdata = data.frame(km_4 = levels(data$km_4)), se.fit = T)) %>% 
    dplyr::select(fit, se.fit) %>%
    mutate(km_4 = levels(data$km_4))
  
  require(multcomp)
  tukey_score <- glht(anova_score, linfct = mcp(km_4 = "Tukey"))
  cld.var_score <- cld(tukey_score, decreasing = F)
  tab.mult.comp_score <- data.frame(km_4 = levels(data$km_4),
                                    tukey = cld.var_score$mcletters$Letters)
  
  
  return(left_join(aov_sorties, tab.mult.comp_score, by="km_4"))
  
}


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
  

p1 <- ggplot(tab_tukey) +
  geom_text(aes(x = km_4, y=0, label = lab),size=3, check_overlap = T) +
  facet_grid(.~FPC)+
  theme_cowplot() +
  labs(y="Mean\n(standard\ndeviation)") +
  theme(legend.position = "none",
        panel.border = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(angle=0,size=10, vjust = 0.5),
        axis.text = element_blank(),
        plot.margin = unit(c(0, 1, 0, 0), "cm"),
        axis.ticks = element_blank(),
        axis.line = element_blank()) ; p1
  
p2 <- tab_tukey %>% 
  mutate(lab_y = case_when(
    FPC == "FPC score 1" ~ 11500,
    FPC == "FPC score 2" ~ 4500,
    FPC == "FPC score 3" ~ 2500
  )) %>% 
  ggplot(., aes(x = km_4, y = score)) +
  geom_text(aes(x=km_4, y=lab_y, label = tukey, color=tukey), size=3, check_overlap = T)+
  #geom_boxplot(aes(fill=km_4), width=0.75, outlier.shape = NA) + 
  #scale_fill_manual(values=rev(wesanderson::wes_palette("Zissou1", 4, "continuous"))) +
  geom_boxplot(aes(fill=tukey, color=tukey), width=0.75, outlier.shape = NA, alpha=0.5) + 
  scale_color_manual(values =rev(c("#575C6DFF", "#4576EF", "#730000", "#2D0F38"))) +
  scale_fill_manual(values =rev(c("#575C6DFF", "#4576EF", "#730000", "#2D0F38"))) +
  geom_jitter(aes(color=tukey), width=0.25, size=1) +
  facet_wrap(.~FPC, scales="free", nrow=1) +
  theme_cowplot() +
  theme(legend.position = "none",
        panel.border = element_rect(color="black", linewidth = 1),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title = element_text(size=10),
        axis.text =  element_text(size=10),
        plot.margin = unit(c(0, 1, 0, 0), "cm")) +
  labs(x = "Clusters of fragments", y = "FPC score value") ; p2

plot_grid(p1,p2, rel_heights = c(0.3, 0.7), ncol=1, axis = "rl")

ggsave(filename = paste0(path_project, "/04_ARTICLE/FIGURES/clusters_caracteristics_panel_a.png"),
       bg="white", dpi = 300, width=10, height=5)

ggsave(filename = paste0(path_project, "/04_ARTICLE/FIGURES_R1/FIGURE_4.tiff"),
       bg="white", dpi = 600, width=10, height=5)



# Revised version
tab_tukey %>% 
  mutate(lab_y = case_when(
    FPC == "FPC score 1" ~ 11500,
    FPC == "FPC score 2" ~ 4500,
    FPC == "FPC score 3" ~ 2500
  )) %>% 
  ggplot(., aes(x = km_4, y = score)) +
  geom_text(aes(x=km_4, y=lab_y, label = tukey, color=tukey), size=3, check_overlap = T)+
  #geom_boxplot(aes(fill=km_4), width=0.75, outlier.shape = NA) + 
  #scale_fill_manual(values=rev(wesanderson::wes_palette("Zissou1", 4, "continuous"))) +
  geom_boxplot(aes(fill=tukey, color=tukey), width=0.75, outlier.shape = NA, alpha=0.5) + 
  scale_color_manual(values =rev(c("#575C6DFF", "#4576EF", "#730000", "#2D0F38"))) +
  scale_fill_manual(values =rev(c("#575C6DFF", "#4576EF", "#730000", "#2D0F38"))) +
  geom_jitter(aes(color=tukey), width=0.25, size=1) +
  facet_wrap(.~FPC, scales="free", nrow=1) +
  theme_cowplot() +
  theme(legend.position = "none",
        panel.border = element_rect(color="black", linewidth = 1),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        axis.title = element_text(size=10),
        axis.text =  element_text(size=10)) +
  labs(x = "Clusters of fragments", y = "FPC score value")


ggsave(filename = paste0(path_project, "/04_ARTICLE/FIGURES/clusters_caracteristics_panel_a_R1.png"),
       bg="white", dpi = 300, width=10, height=3.5)

ggsave(filename = paste0(path_project, "/04_ARTICLE/FIGURES_R1/FIGURE_4_R1.tiff"),
       bg="white", dpi = 600, width=10, height=3.5)


# Plot the curves for each clusters
p3 <- plot_grid(
  plot_dyn_id(frags = X_scores1 %>% filter(km_4==1) %>% pull(ID), score = NULL, color = pal_cluster[1]) +
    geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
    geom_hline(yintercept = 0, color="darkgrey", lty=2) +
    ggtitle("a. Cluster 1", subtitle = paste0("(N=", paste0(KM4$size[2]), ")")) +
    theme_cowplot() +
    theme(legend.position = "none",
          panel.border = element_rect(color="black", linewidth = 1),
          strip.background = element_blank(),
          strip.text = element_blank(),
          axis.title = element_text(size=10),
          plot.subtitle = element_text(size=11),
          axis.text =  element_text(size=10),
          plot.margin = unit(c(0, 0, 0, 0), "cm")),
  
  plot_dyn_id(frags = X_scores1 %>% filter(km_4==2) %>% pull(ID), score = NULL, color = pal_cluster[2]) +
    geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
    geom_hline(yintercept = 0, color="darkgrey", lty=2) +
    ggtitle("b. Cluster 2", subtitle = paste0("(N=", paste0(KM4$size[1]), ")")) +
    theme_cowplot() +
    theme(legend.position = "none",
          panel.border = element_rect(color="black", linewidth = 1),
          strip.background = element_blank(),
          strip.text = element_blank(),
          axis.title.x = element_text(size=10),
          axis.text =  element_text(size=10),
          plot.margin = unit(c(0, 0, 0, 0), "cm"),
          plot.subtitle = element_text(size=11),
          axis.title.y = element_blank(),
          axis.text.y = element_blank()),
  
  plot_dyn_id(frags = X_scores1 %>% filter(km_4==3) %>% pull(ID), score = NULL, color = pal_cluster[3]) +
    geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
    geom_hline(yintercept = 0, color="darkgrey", lty=2) +
    ggtitle("c. Cluster 3", subtitle = paste0("(N=", paste0(KM4$size[3]), ")")) +
    theme_cowplot() +
    theme(legend.position = "none",
          panel.border = element_rect(color="black", linewidth = 1),
          strip.background = element_blank(),
          strip.text = element_blank(),
          axis.title.x = element_text(size=10),
          axis.text =  element_text(size=10),
          plot.subtitle = element_text(size=11),
          plot.margin = unit(c(0, 0, 0, 0), "cm"),
          axis.title.y = element_blank(),
          axis.text.y = element_blank()),
  
  plot_dyn_id(frags = X_scores1 %>% filter(km_4==4) %>% pull(ID), score = NULL, color = pal_cluster[4]) +
    geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
    geom_hline(yintercept = 0, color="darkgrey", lty=2) +
    ggtitle("d. Cluster 4", subtitle = paste0("(N=", paste0(KM4$size[4]), ")")) +
    theme_cowplot() +
    theme(legend.position = "none",
          panel.border = element_rect(color="black", linewidth = 1),
          strip.background = element_blank(),
          strip.text = element_blank(),
          plot.subtitle = element_text(size=11),
          axis.title.x = element_text(size=10),
          axis.text.x =  element_text(size=10),
          plot.margin = unit(c(0, 1, 0, 0), "cm"),
          axis.title.y = element_blank(),
          axis.text.y = element_blank()),
  
  nrow=1, align="hv", axis="tbrl", rel_widths = c(0.31, 0.29, 0.29, 0.29)
) ; p3


#plot_grid(plot_grid(p1,p2, rel_heights = c(0.3, 0.7), ncol=1, axis = "btrl", align = "v"), 
#          p3, 
#          ncol=1, axis="tbrl", align = "hv", labels = c("a.", "b."))
#
#
#ggsave(filename = paste0(path_project, "/04_ARTICLE/FIGURES/clusters_caracteristics.png"),
#       bg="white", dpi = 300, width=10, height=6)
#

p3
ggsave(filename = paste0(path_project, "/04_ARTICLE/FIGURES/Dynamic_Clusters_4.png"),
       bg="white", dpi = 300, width=14, height=4)

ggsave(filename = paste0(path_project, "/04_ARTICLE/FIGURES_R1/FIGURE_6.tiff"),
       bg="white", dpi = 600, width=14, height=4)

ggsave(filename = paste0(path_project, "/04_ARTICLE/FIGURES_R1/FIGURE_6.png"),
       bg="white", dpi = 300, width=14, height=4)


# Number of first lesions

library(agricolae)
library(cowplot)
library(DHARMa)
library(dtplyr)
library(emmeans)
library(GGally)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(kableExtra)
library(lme4)
library(MASS)
library(multcomp)
library(multcompView)
library(MuMIn)
library(plotly)
library(pROC)
library(readr)
library(rstatix)
library(scales)
library(table1)
library(tidyverse)


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
mod_1_blc4 <- pffr(Y_t ~ 0 + c(0) +           # dropping the constant and time-varying intercepts via 0 + c(0) to estimate zone + c(zone) 
                     km_4 + c(km_4) +           # time-varying cluster effects not centered at zero 
                     s(tree_leaf, bs="re"),     # leaf*tree-specific smooth residuals
                   data = data, 
                   ydata = ydata2,
                   bs.yindex = list(bs = "ps", k=5),
                   bs.int = list(bs = "ps", k=20)) 

summary(mod_1_blc4) ; AIC(mod_1_blc4)
# R-sq.(adj) =  0.946 ; REML score =  19582 ; AIC = 39140.37

# additionally adjusted for the initial number of lesions
mod_1_blc4l <- pffr(Y_t ~ 0 + c(0) +           # dropping the constant and time-varying intercepts via 0 + c(0) to estimate zone + c(zone) 
                      km_4 + c(km_4) +           # time-varying cluster effects not centered at zero 
                      nblesions1 + 
                      s(tree_leaf, bs="re"),     # leaf*tree-specific smooth residuals
                    data = data, 
                    ydata = ydata2,
                    bs.yindex = list(bs = "ps", k=5),
                    bs.int = list(bs = "ps", k=20)) 

summary(mod_1_blc4l) ; AIC(mod_1_blc4l)
# R-sq.(adj) =  0.949 ; REML score =  19501 ; AIC = 38976.56 

# Extraction of the coefficients for the clusters
pa <- list("mod_1_blc4"=  mod_1_blc4,
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
  mutate(nb_cluster = case_when(
    model %in% c("mod_1_blc3", "mod_1_blc3l") ~ "Fragments partitioned into\n3 clusters",
    model %in% c("mod_1_blc4", "mod_1_blc4l") ~ "Fragments partitioned into\n4 clusters",
    model %in% c("mod_1_blc5", "mod_1_blc5l") ~ "Fragments partitioned into\n5 clusters"
  )) %>% 
  mutate(cluster = case_when(
    variable %in% c("km_31(yindex)", "km_41(yindex)","km_51(yindex)") ~ "Cluster 1",
    variable %in% c("km_32(yindex)", "km_42(yindex)","km_52(yindex)") ~ "Cluster 2",
    variable %in% c("km_33(yindex)", "km_43(yindex)","km_53(yindex)") ~ "Cluster 3",
    variable %in% c("km_44(yindex)","km_54(yindex)") ~ "Cluster 4",
    variable %in% c("km_55(yindex)") ~ "Cluster 5"
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

ggsave(filename = paste0(path_project, "/04_ARTICLE/FIGURES/Dynamic_Clusters_4_regression.png"),
       bg="white", dpi = 300, width=10, height=7)

estimated_mean %>%
  filter(adjustment=="Model not accounting for\nthe initial number of  lesions")%>% 
  ggplot(., aes(x = x_i, color=cluster)) +
  geom_hline(yintercept = c(0, 3512), color="black", lty=2) +
  # > mean 
  geom_path(data = Fit.mu1, 
            aes(x=d, y=mu), color="black", lwd=1) +
  geom_line(aes(y=f_i), linewidth=1.2) +
  geom_ribbon(aes(ymin=ic_down, ymax=ic_up, fill=cluster), color="transparent", alpha=0.3) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x="Days after rank 4", y="Mean diseased surface in mm²\n95% confidence interval) ") +
  guides(color = guide_legend(nrow=1)) +
  lims(y=c(-100, 3600)) +
  scale_color_manual(values=pal_cluster, name="Cluster") +
  scale_fill_manual(values=pal_cluster, name="Cluster")

ggsave(filename = paste0(path_project, "/04_ARTICLE/FIGURES/Dynamic_Clusters_4_regression_final.png"),
       bg="white", dpi = 300, width=5, height=5)

ggsave(filename = paste0(path_project, "/04_ARTICLE/FIGURES_R1/FIGURE_7.tiff"),
       bg="white", dpi = 600, width=5, height=5)

ggsave(filename = paste0(path_project, "/04_ARTICLE/FIGURES_R1/FIGURE_7.png"),
       bg="white", dpi = 300, width=5, height=5)

plot_dyn_id(frags = X_scores1 %>% filter(km_4==2) %>% pull(ID), score = NULL, color = "grey") +
  geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
  geom_hline(yintercept = 0, color="darkgrey", lty=2) +
  ggtitle("a.", subtitle = paste0("Cluster 2 (N=", paste0(KM4$size[2]), ")")) +
  theme(legend.position = "none") +
  geom_line(data = estimated_mean %>% 
              filter(cluster=="2", adjustment=="Model not accounting for\nthe initial number of  lesions"),
            aes(x = x_i, y=f_i, color=pal_cluster[2])) +
  geom_ribbon(data = estimated_mean %>% 
                filter(cluster=="2", adjustment=="Model not accounting for\nthe initial number of  lesions"),
              aes(x = x_i, ymin=ic_down, ymax=ic_up, fill=pal_cluster[2]), color="transparent", alpha=0.3) +
  lims(y=c(-100, 3600))



# Plot the curves for each clusters
plot_grid(
  
  plot_dyn_id(frags = X_scores1 %>% filter(km_4==2) %>% pull(ID), score = NULL, color = "grey") +
    geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
    geom_hline(yintercept = 0, color="darkgrey", lty=2) +
    ggtitle("a.", subtitle = paste0("Cluster 2 (N=", paste0(KM4$size[2]), ")")) +
    theme(legend.position = "none") +
    geom_line(data = estimated_mean %>% 
                filter(cluster=="2", adjustment=="Model not accounting for\nthe initial number of  lesions"),
              aes(x = x_i, y=f_i), color="darkred", size=1.5) +
    geom_ribbon(data = estimated_mean %>% 
                  filter(cluster=="2", adjustment=="Model not accounting for\nthe initial number of  lesions"),
                aes(x = x_i, ymin=ic_down, ymax=ic_up), fill=pal_cluster[2], color="transparent", alpha=0.5) +
    lims(y=c(-100, 3600)),
  
  plot_dyn_id(frags = X_scores1 %>% filter(km_4==1) %>% pull(ID), score = NULL, color = "grey") +
    geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
    geom_hline(yintercept = 0, color="darkgrey", lty=2) +
    ggtitle("", subtitle = paste0("Cluster 1 (N=", paste0(KM4$size[1]), ")")) +
    theme(legend.position = "none") +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank()) +
    geom_line(data = estimated_mean %>% 
                filter(cluster=="1", adjustment=="Model not accounting for\nthe initial number of  lesions"),
              aes(x = x_i, y=f_i), color="darkorange", size=1.5) +
    geom_ribbon(data = estimated_mean %>% 
                  filter(cluster=="1", adjustment=="Model not accounting for\nthe initial number of  lesions"),
                aes(x = x_i, ymin=ic_down, ymax=ic_up), fill=pal_cluster[1], color="transparent", alpha=0.5) +
    lims(y=c(-100, 3600)),
  
  plot_dyn_id(frags = X_scores1 %>% filter(km_4==3) %>% pull(ID), score = NULL, color = "grey") +
    geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
    geom_hline(yintercept = 0, color="darkgrey", lty=2) +
    ggtitle("", subtitle = paste0("Cluster 3 (N=", paste0(KM4$size[3]), ")")) +
    theme(legend.position = "none") +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank()) +
    geom_line(data = estimated_mean %>% 
                filter(cluster=="3", adjustment=="Model not accounting for\nthe initial number of  lesions"),
              aes(x = x_i, y=f_i), color="darkgreen", size=1.5) +
    geom_ribbon(data = estimated_mean %>% 
                  filter(cluster=="3", adjustment=="Model not accounting for\nthe initial number of  lesions"),
                aes(x = x_i, ymin=ic_down, ymax=ic_up), fill=pal_cluster[3], color="transparent", alpha=0.5) +
    lims(y=c(-100, 3600)),
  
  plot_dyn_id(frags = X_scores1 %>% filter(km_4==4) %>% pull(ID), score = NULL, color = "grey") +
    geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
    geom_hline(yintercept = 0, color="darkgrey", lty=2) +
    ggtitle("", subtitle = paste0("Cluster 4 (N=", paste0(KM4$size[4]), ")")) +
    theme(legend.position = "none") +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank()) +
    geom_line(data = estimated_mean %>% 
                filter(cluster=="4", adjustment=="Model not accounting for\nthe initial number of  lesions"),
              aes(x = x_i, y=f_i), color="darkblue", size=1.5) +
    geom_ribbon(data = estimated_mean %>% 
                  filter(cluster=="4", adjustment=="Model not accounting for\nthe initial number of  lesions"),
                aes(x = x_i, ymin=ic_down, ymax=ic_up), fill=pal_cluster[4], color="transparent", alpha=0.5) +
    lims(y=c(-100, 3600)),
  
  plot_dyn_id(frags = X_scores1 %>% filter(km_4==2) %>% pull(ID), score = NULL, color = "grey") +
    geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
    geom_hline(yintercept = 0, color="darkgrey", lty=2) +
    ggtitle("b.", subtitle = paste0("Cluster 2 (N=", paste0(KM4$size[2]), ")")) +
    theme(legend.position = "none") +
    geom_line(data = estimated_mean %>% 
                filter(cluster=="2", adjustment!="Model not accounting for\nthe initial number of  lesions"),
              aes(x = x_i, y=f_i), color="darkred", size=1.5) +
    geom_ribbon(data = estimated_mean %>% 
                  filter(cluster=="2", adjustment=="Model not accounting for\nthe initial number of  lesions"),
                aes(x = x_i, ymin=ic_down, ymax=ic_up), fill=pal_cluster[2], color="transparent", alpha=0.5) +
    lims(y=c(-100, 3600)),
  
  plot_dyn_id(frags = X_scores1 %>% filter(km_4==1) %>% pull(ID), score = NULL, color = "grey") +
    geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
    geom_hline(yintercept = 0, color="darkgrey", lty=2) +
    ggtitle("", subtitle = paste0("Cluster 1 (N=", paste0(KM4$size[1]), ")")) +
    theme(legend.position = "none") +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank()) +
    geom_line(data = estimated_mean %>% 
                filter(cluster=="1", adjustment!="Model not accounting for\nthe initial number of  lesions"),
              aes(x = x_i, y=f_i), color="darkorange", size=1.5) +
    geom_ribbon(data = estimated_mean %>% 
                  filter(cluster=="1", adjustment=="Model not accounting for\nthe initial number of  lesions"),
                aes(x = x_i, ymin=ic_down, ymax=ic_up), fill=pal_cluster[1], color="transparent", alpha=0.5) +
    lims(y=c(-100, 3600)),
  
  plot_dyn_id(frags = X_scores1 %>% filter(km_4==3) %>% pull(ID), score = NULL, color = "grey") +
    geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
    geom_hline(yintercept = 0, color="darkgrey", lty=2) +
    ggtitle("", subtitle = paste0("Cluster 3 (N=", paste0(KM4$size[3]), ")")) +
    theme(legend.position = "none") +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank()) +
    geom_line(data = estimated_mean %>% 
                filter(cluster=="3", adjustment!="Model not accounting for\nthe initial number of  lesions"),
              aes(x = x_i, y=f_i), color="darkgreen", size=1.5) +
    geom_ribbon(data = estimated_mean %>% 
                  filter(cluster=="3", adjustment=="Model not accounting for\nthe initial number of  lesions"),
                aes(x = x_i, ymin=ic_down, ymax=ic_up), fill=pal_cluster[3], color="transparent", alpha=0.5) +
    lims(y=c(-100, 3600)),
  
  
  # b
  plot_dyn_id(frags = X_scores1 %>% filter(km_4==4) %>% pull(ID), score = NULL, color = "grey") +
    geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
    geom_hline(yintercept = 0, color="darkgrey", lty=2) +
    ggtitle("", subtitle = paste0("Cluster 4 (N=", paste0(KM4$size[4]), ")")) +
    theme(legend.position = "none") +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank()) +
    geom_line(data = estimated_mean %>% 
                filter(cluster=="4", adjustment!="Model not accounting for\nthe initial number of  lesions"),
              aes(x = x_i, y=f_i), color="darkblue", size=1.5) +
    geom_ribbon(data = estimated_mean %>% 
                  filter(cluster=="4", adjustment=="Model not accounting for\nthe initial number of  lesions"),
                aes(x = x_i, ymin=ic_down, ymax=ic_up), fill=pal_cluster[4], color="transparent", alpha=0.5) +
    lims(y=c(-100, 3600)),
  
  nrow=2, align="hv", axis="tbrl", rel_widths = c(0.31, 0.29, 0.29, 0.29)
)

# 
ggsave(filename = paste0(path_project, "/04_ARTICLE/FIGURES/Dynamic_Clusters_model.png"),
       bg="white", dpi = 300, width=14, height=8)


p1 <- list("mod_1_blc4l"= mod_1_blc4l) %>% 
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
  mutate(nb_cluster = case_when(
    model %in% c("mod_1_blc3", "mod_1_blc3l") ~ "Fragments partitioned into\n3 clusters",
    model %in% c("mod_1_blc4", "mod_1_blc4l") ~ "Fragments partitioned into\n4 clusters",
    model %in% c("mod_1_blc5", "mod_1_blc5l") ~ "Fragments partitioned into\n5 clusters"
  )) %>% 
  mutate(cluster = case_when(
    variable %in% c("km_31(yindex)", "km_41(yindex)","km_51(yindex)") ~ "Cluster 1",
    variable %in% c("km_32(yindex)", "km_42(yindex)","km_52(yindex)") ~ "Cluster 2",
    variable %in% c("km_33(yindex)", "km_43(yindex)","km_53(yindex)") ~ "Cluster 3",
    variable %in% c("km_44(yindex)","km_54(yindex)") ~ "Cluster 4",
    variable %in% c("km_55(yindex)") ~ "Cluster 5"
  )) %>%
  #mutate(variable = substr(variable, 1, 5)) %>% 
  ggplot(., aes(x = xi, y = mean, color=cluster)) +
  geom_hline(yintercept = 0, color="black", lty=2) +
  geom_text(aes(x = 0, y = 2000, label = lab), check_overlap = T, color = "black", nudge_x = -4, hjust = 0) +
  geom_path(linewidth = 1) +
  geom_ribbon(aes(ymin=ic_down, ymax=ic_up, fill=cluster), color="transparent", alpha=0.3) +
  #facet_grid(adjustment_lesions~nb_cluster) +
  scale_color_manual(values=pal_cluster, name="") +
  scale_fill_manual(values=pal_cluster, name="") +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Days after rank 4", y = "Mean functional coefficient (95% confidence interval)") + 
  ggtitle("a.")

p2 <- data.frame(
  model = c(rep("Fragments partitioned into\n3 clusters", 29*3), rep("Fragments partitioned into\n4 clusters", 29*4), rep("Fragments partitioned into\n5 clusters", 29*5), 
            rep("Fragments partitioned into\n3 clusters", 29*3), rep("Fragments partitioned into\n4 clusters", 29*4), rep("Fragments partitioned into\n5 clusters", 29*5)),
  adjustment = c(rep("Model not accounting for\nthe initial number of  lesions", 29*12),
                 rep("Model accounting for\nthe initial number of  lesions", 29*12)),
  cluster = c(rep("Cluster 1", 29),  rep("Cluster 2", 29),  rep("Cluster 3", 29),
              rep("Cluster 1", 29),  rep("Cluster 2", 29),  rep("Cluster 3", 29), rep("Cluster 4", 29),
              rep("Cluster 1", 29),  rep("Cluster 2", 29),  rep("Cluster 3", 29), rep("Cluster 4", 29), rep("Cluster 5", 29),
              rep("Cluster 1", 29),  rep("Cluster 2", 29),  rep("Cluster 3", 29),
              rep("Cluster 1", 29),  rep("Cluster 2", 29),  rep("Cluster 3", 29), rep("Cluster 4", 29),
              rep("Cluster 1", 29),  rep("Cluster 2", 29),  rep("Cluster 3", 29), rep("Cluster 4", 29), rep("Cluster 5", 29)), 
  x_i     = c(rep(yindex2, 3*24)),
  f_i     = c(pred_mod_1_blc3$fit[1,],  pred_mod_1_blc3$fit[2,],  pred_mod_1_blc3$fit[3,],
              pred_mod_1_blc4$fit[1,],  pred_mod_1_blc4$fit[2,],  pred_mod_1_blc4$fit[3,],  pred_mod_1_blc4$fit[4,],
              pred_mod_1_blc5$fit[1,],  pred_mod_1_blc5$fit[2,],  pred_mod_1_blc5$fit[3,],  pred_mod_1_blc5$fit[4,],  pred_mod_1_blc5$fit[5,],
              pred_mod_1_blc3l$fit[1,],  pred_mod_1_blc3l$fit[2,],  pred_mod_1_blc3l$fit[3,],
              pred_mod_1_blc4l$fit[1,],  pred_mod_1_blc4l$fit[2,],  pred_mod_1_blc4l$fit[3,],  pred_mod_1_blc4l$fit[4,],
              pred_mod_1_blc5l$fit[1,],  pred_mod_1_blc5l$fit[2,],  pred_mod_1_blc5l$fit[3,],  pred_mod_1_blc5l$fit[4,],  pred_mod_1_blc5l$fit[5,]),
  se_i    = c(pred_mod_1_blc3$se.fit[1,],  pred_mod_1_blc3$se.fit[2,],    pred_mod_1_blc3$se.fit[3,],
              pred_mod_1_blc4$se.fit[1,],  pred_mod_1_blc4$se.fit[2,],  pred_mod_1_blc4$se.fit[3,],  pred_mod_1_blc4$se.fit[4,],
              pred_mod_1_blc5$se.fit[1,],  pred_mod_1_blc5$se.fit[2,],  pred_mod_1_blc5$se.fit[3,],  pred_mod_1_blc5$se.fit[4,],  pred_mod_1_blc5$se.fit[5,],
              pred_mod_1_blc3l$se.fit[1,],  pred_mod_1_blc3l$se.fit[2,],  pred_mod_1_blc3l$se.fit[3,],
              pred_mod_1_blc4l$se.fit[1,],  pred_mod_1_blc4l$se.fit[2,],  pred_mod_1_blc4l$se.fit[3,],  pred_mod_1_blc4l$se.fit[4,],
              pred_mod_1_blc5l$se.fit[1,],  pred_mod_1_blc5l$se.fit[2,],  pred_mod_1_blc5l$se.fit[3,],  pred_mod_1_blc5l$se.fit[4,],  pred_mod_1_blc5l$se.fit[5,])) %>% 
  # Compute 95% CI
  mutate(
    ic_up   = f_i + 1.96*se_i,
    ic_down = f_i - 1.96*se_i) %>% 
  mutate(ic_up = if_else(ic_up>3600, 3600, ic_up),
         ic_up = if_else(ic_up<0, 0, ic_up),
         ic_down = if_else(ic_down>3600, 3600, ic_down),
         ic_down = if_else(ic_down<0, 0, ic_down)) %>%
  filter(adjustment == "Model accounting for\nthe initial number of  lesions",
         model == "Fragments partitioned into\n4 clusters") %>% 
  ggplot(., aes(x = x_i, color=cluster)) +
  geom_hline(yintercept = c(0, 3512), color="black", lty=2) +
  # > mean 
  geom_path(data = Fit.mu1, 
            aes(x=d, y=mu), color="black", lwd=1) +
  geom_line(aes(y=f_i), linewidth=1.2) +
  geom_ribbon(aes(ymin=ic_down, ymax=ic_up, fill=cluster), color="transparent", alpha=0.3) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x="Days after rank 4", y="Mean infected surface in mm² (and)95% confidence interval) ") +
  #facet_grid(adjustment~model) +
  guides(color = guide_legend(nrow=1)) +
  lims(y=c(-100, 3600)) +
  scale_color_manual(values=pal_cluster, name="") +
  scale_fill_manual(values=pal_cluster, name="")+ 
  ggtitle("b.") 

plot_grid(p1, plot_grid(p2, ggplot()+theme_void(), ncol=1, rel_heights = c(0.91, 0.09)), nrow=1)

ggsave(filename = paste0(path_project, "/04_ARTICLE/FIGURES/Regression_functional.png"),
       bg="white", dpi = 300, width=12, height=6)


# > AFC
# > Contingency tables using Fisher exact test 

tableby(km_4 ~ bananier + feuille + compartiment, 
        data=X_scores1, 
        numeric.stats=c("meansd")) %>% 
  summary(., text=TRUE, pfootnote=TRUE, digits = 1, total=F) %>% 
  kbl(.) %>%
  kable_styling(bootstrap_options = "striped", full_width = F)


tableby(compartiment_old ~ km_4, 
        data=X_scores1, 
        numeric.stats=c("meansd")) %>% 
  summary(., text=TRUE, pfootnote=TRUE, digits = 1, total=F) %>% 
  kbl(.) %>%
  kable_styling(bootstrap_options = "striped", full_width = F)



# > Fisher exact tests

X_scores1$compartiment_old <- if_else(X_scores1$compartiment=="8", "5", X_scores1$compartiment)

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

# AFC sur fragments
afc_fragments <- CA(table(t$km_4, t$fragment)) ; summary(afc_fragments)

# clusters on col (fragments) coordinates 
coords_cols <- data.frame(afc_fragments$col$coord)
summary(coords_cols)

z_coords_cols <- scale(coords_cols)

plot_grid(
  fviz_nbclust(z_coords_cols, kmeans, method = "wss")        + ggtitle(label = paste0("kmeans"), subtitle = "Elbow method"),
  fviz_nbclust(z_coords_cols, kmeans, method = "silhouette") + ggtitle(label = "",               subtitle = "Silhouette method"),
  fviz_nbclust(z_coords_cols, kmeans, method = "gap_stat", 
               nstart = 25, nboot = 50, verbose = FALSE)     + ggtitle(label = "",               subtitle = "Gap statistic method"),
  ncol = 3)

set.seed(123)
clusters_afc_fragments <- kmeans(z_coords_cols, centers = 3, iter.max = 25)

coords_cols$km <- as.vector(clusters_afc_fragments$cluster)
coords_cols$lab <- rownames(afc_fragments$col$coord)

coords_row <- data.frame(afc_fragments$row$coord)
coords_row$lab <- rownames(afc_fragments$row$coord)

afc_fragments$eig

# Plot

p1 <- ggplot() + 
  geom_hline(yintercept = 0, lty=2, color="darkgrey") +
  geom_vline(xintercept = 0, lty=2, color="darkgrey") +
  # Clusters (rows)
  geom_point(data = coords_row, aes(x = Dim.1, y = Dim.2), pch=8, size=5) +
  geom_text_repel(data = coords_row, 
                  aes(x = Dim.1, y = Dim.2, label = lab), 
                  color="black", box.padding = 0.75, nudge_y = 0.05, size= 4) +
  # Fragments (columns)
  geom_point(data = coords_cols, aes(x = Dim.1, y = Dim.2), color="darkblue", 
             size=2, pch=18) + 
  geom_text_repel(data = coords_cols, 
                  aes(x = Dim.1, y = Dim.2, label = lab), color="darkblue") +
  theme_cowplot() +
  theme(legend.position = "none") +
  ggtitle("a.") +
  labs(x = paste0("Dim. 1 (", round(afc_fragments$eig[1,2], 1), "%)"), 
       y = paste0("Dim. 2 (", round(afc_fragments$eig[2,2], 1), "%)")) +
  scale_x_continuous(breaks = c(-1, -0.5, 0, 0.5, 1, 1.5, 2)) +
  scale_y_continuous(breaks = c(-1, -0.5, 0, 0.5, 1)) ; p1

p3 <- ggplot() + 
  geom_hline(yintercept = 0, lty=2, color="darkgrey") +
  geom_vline(xintercept = 0, lty=2, color="darkgrey") +
  # Clusters (rows)
  geom_point(data = coords_row, aes(x = Dim.1, y = Dim.3), pch=8, size=5) +
  geom_text_repel(data = coords_row, 
                  aes(x = Dim.1, y = Dim.3, label = lab), 
                  color="black", box.padding = 0.75, nudge_y = 0.05, size= 4) +
  # Fragments (columns)
  geom_point(data = coords_cols, aes(x = Dim.3, y = Dim.2), color="darkblue", 
             size=2, pch=18) + 
  geom_text_repel(data = coords_cols, 
                  aes(x = Dim.1, y = Dim.3, label = lab), color="darkblue") +
  theme_cowplot() +
  theme(legend.position = "none") +
  ggtitle(" ") +
  labs(x = paste0("Dim. 1 (", round(afc_fragments$eig[1,2], 1), "%)"), 
       y = paste0("Dim. 3 (", round(afc_fragments$eig[3,2], 1), "%)")) +
  scale_x_continuous(breaks = c(-1, -0.5, 0, 0.5, 1, 1.5, 2)) +
  scale_y_continuous(breaks = c(-1, -0.5, 0, 0.5, 1)) ; p3

p2 <- ggplot() + 
  geom_hline(yintercept = 0, lty=2, color="darkgrey") +
  geom_vline(xintercept = 0, lty=2, color="darkgrey") +
  geom_point(data = coords_row, aes(x = Dim.1, y = Dim.2), pch=8, size=5, color="darkgrey") +
  geom_point(data = coords_cols, aes(x = Dim.1, y = Dim.2, color = as.factor(km)), 
             size=2, pch=18) + 
  geom_text_repel(data = coords_cols, 
            aes(x = Dim.1, y = Dim.2, label = lab, color = as.factor(km))) +
  scale_color_viridis_d(option = "C", direction = -1, end=0.8) +
  #stat_ellipse(data = coords_cols, aes(x = Dim.1, y = Dim.2, color = as.factor(km))) +
  theme_cowplot() +
  theme(legend.position = "none") +
  ggtitle("b.") +
  labs(x = paste0("Dim. 1 (", round(afc_fragments$eig[1,2], 1), "%)"), 
       y = paste0("Dim. 2 (", round(afc_fragments$eig[2,2], 1), "%)")) +
  scale_x_continuous(breaks = c(-1, -0.5, 0, 0.5, 1, 1.5, 2)) +
  scale_y_continuous(breaks = c(-1, -0.5, 0, 0.5, 1)) ; p2

p4 <- ggplot() + 
  geom_hline(yintercept = 0, lty=2, color="darkgrey") +
  geom_vline(xintercept = 0, lty=2, color="darkgrey") +
  geom_point(data = coords_row, aes(x = Dim.1, y = Dim.3), pch=8, size=5, color="darkgrey") +
  geom_point(data = coords_cols, aes(x = Dim.1, y = Dim.3, color = as.factor(km)), 
             size=2, pch=18) + 
  geom_text_repel(data = coords_cols,
            aes(x = Dim.1, y = Dim.3, label = lab, color = as.factor(km))) +
  scale_color_viridis_d(option = "C", direction = -1, end=0.8) +
  #stat_ellipse(data = coords_cols, aes(x = Dim.1, y = Dim.2, color = as.factor(km))) +
  theme_cowplot() +
  theme(legend.position = "none") +
  ggtitle(" ") +
  labs(x = paste0("Dim. 1 (", round(afc_fragments$eig[1,2], 1), "%)"), 
       y = paste0("Dim. 3 (", round(afc_fragments$eig[3,2], 1), "%)")) +
  scale_x_continuous(breaks = c(-1, -0.5, 0, 0.5, 1, 1.5, 2)) +
  scale_y_continuous(breaks = c(-1, -0.5, 0, 0.5, 1)) ; p4

plot_grid(p1, p3, p2, p4, ncol=2, rel_heights = c(0.525, 0.475), rel_widths = c(0.475, 0.525))

ggsave(filename = paste0(path_project, "/04_ARTICLE/FIGURES/CA_analysis_ellipse.png"),
       bg="white", dpi = 300, width=10, height=10)




# AFC sur compartiments

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

plot(afc_compartiments, c(2,3))
  
  
  
View(afc_compartiments$eig)

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

rbind(data.frame(afc_compartiments$col$cos2),
      data.frame(afc_compartiments$row$cos2)) %>% 
  View(.)

# clusters on col (fragments) coordinates 
coords_cols <- rbind(data.frame(afc_compartiments$col$coord),
                     data.frame(afc_compartiments$row$coord))
summary(coords_cols)

z_coords_cols <- as.data.frame(scale(coords_cols[,1:2]))
summary(z_coords_cols)

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

#coords_row <- data.frame(afc_compartiments$row$coord)
#coords_row$lab <- rownames(afc_compartiments$row$coord)

# Plot
f <- function(x) {
  res <- scales::percent(x, accuracy = .01)
  res
}


p0 <- X_scores1 %>% 
  group_by(compartiment_old, km_4) %>% 
  count() %>% 
  group_by(compartiment_old) %>% 
  mutate(freq2=f(n / sum(n)), 
         freq = round((n / sum(n)) * 100, 1),
         cum_freq = cumsum(freq)) %>% 
  #mutate(km_4 =paste0("Cluster ", km_4),
  #       compartiment_old = paste0("Comp. ", compartiment_old)) %>% 
  mutate(km_4 = factor(km_4, levels = c("4", "3", "2", "1"))) %>% 
  ggplot(.) +
  geom_col(aes(x=compartiment_old, y = freq, fill=km_4)) +
  geom_text(aes(x=compartiment_old, y = cum_freq - freq/2, label = freq2), 
            size=5)+
  #coord_flip() + 
  scale_fill_manual(values=rev(pal_cluster), 
                    name = expression(paste(italic("A posteriori"), plain("clusters"))), 
                    guide = guide_legend(reverse = T)) +
  facet_grid(.~compartiment_old, scales="free", space="free") +
  theme_cowplot() +
  theme(strip.background = element_blank(), 
        strip.text = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(y = "Frequence (%)", x=expression(paste(italic("A priori"), plain("compartments")))) +
  ggtitle("a.") ; p0

p1 <- ggplot() + 
  geom_hline(yintercept = 0, lty=2, color="darkgrey") +
  geom_vline(xintercept = 0, lty=2, color="darkgrey") +
  # clusters
  #geom_point(data = coords_row, aes(x = Dim.1, y = Dim.2), pch=8, size=5, color="darkgrey") +
  #geom_text(data = coords_row, 
  #          aes(x = Dim.1, y = Dim.2, label = lab), 
  #          color="darkgrey", nudge_x = 0.125, nudge_y = 0.025, size= 4) +
  # compartiments
  geom_point(data = coords_cols, aes(x = Dim.1, y = Dim.2, color = as.factor(km), pch=as.factor(km)), 
             size=3) + 
  geom_text(data = coords_cols, 
                  aes(x = Dim.1, y = Dim.2, label = lab, color = as.factor(km)),
            nudge_y = -0.025, size=5) +
  scale_color_viridis_d(option = "C", direction = -1, end=0.8) +
  #stat_ellipse(data = coords_cols, aes(x = Dim.1, y = Dim.2, color = as.factor(km))) +
  theme_cowplot() +
  theme(legend.position = "none") +
  labs(x = paste0("Axis 1 (explained variability: ", round(afc_compartiments$eig[1,2], 1), "%)"), 
       y = paste0("Axis 2 (explained variability: ", round(afc_compartiments$eig[2,2], 1), "%)")) +
  lims(x=c(-0.45,1.5)) +
  ggtitle("b.")

p2 <- ggplot() + 
  geom_hline(yintercept = 0, lty=2, color="darkgrey") +
  geom_vline(xintercept = 0, lty=2, color="darkgrey") +
  # clusters
  geom_point(data = coords_row, aes(x = Dim.3, y = Dim.2), pch=8, size=5, color="darkgrey") +
  geom_text(data = coords_row, 
            aes(x = Dim.3, y = Dim.2, label = lab), 
            color="darkgrey", nudge_x = 0.025, nudge_y = 0.025, size= 4) +
  ## compartiments
  geom_point(data = coords_cols, aes(x = Dim.3, y = Dim.2, color = as.factor(km)), 
             size=2, pch=18) + 
  geom_text(data = coords_cols, 
            aes(x = Dim.3, y = Dim.2, label = lab, color = as.factor(km)),
            nudge_y = -0.025) +
  scale_color_viridis_d(option = "C", direction = -1, end=0.8) +
  #stat_ellipse(data = coords_cols, aes(x = Dim.1, y = Dim.2, color = as.factor(km))) +
  theme_cowplot() +
  theme(legend.position = "none") +
  labs(x = paste0("Axis 3 (", round(afc_compartiments$eig[3,2], 1), "%)"), 
       y = paste0("Axis 2 (", round(afc_compartiments$eig[2,2], 1), "%)")) +
  lims(x=c(-0.25,0.25)) +
  ggtitle("c.")



p3 <- ggplot() + 
  geom_hline(yintercept = 0, lty=2, color="darkgrey") +
  geom_vline(xintercept = 0, lty=2, color="darkgrey") +
  # clusters
  geom_point(data = coords_row, aes(x = Dim.3, y = Dim.1), pch=8, size=5, color="darkgrey") +
  geom_text(data = coords_row, 
            aes(x = Dim.3, y = Dim.1, label = lab), 
            color="darkgrey", nudge_x = 0.025, nudge_y = 0.025, size= 4) +
  ## compartiments
  geom_point(data = coords_cols, aes(x = Dim.3, y = Dim.1, color = as.factor(km)), 
             size=2, pch=18) + 
  geom_text(data = coords_cols, 
            aes(x = Dim.3, y = Dim.1, label = lab, color = as.factor(km)),
            nudge_y = -0.025) +
  scale_color_viridis_d(option = "C", direction = -1, end=0.8) +
  #stat_ellipse(data = coords_cols, aes(x = Dim.1, y = Dim.2, color = as.factor(km))) +
  theme_cowplot() +
  theme(legend.position = "none") +
  labs(x = paste0("Axis 3 (", round(afc_compartiments$eig[3,2], 1), "%)"), 
       y = paste0("Axis 1 (", round(afc_compartiments$eig[1,2], 1), "%)")) +
  lims(x=c(-0.25,0.25)) +
  ggtitle("c.")

#plot_grid(plot_grid(p0, ggplot() +theme_void(), rel_widths = c(0.85, 0.15)), 
#          plot_grid(p1,p2, p3), nrow=2)
#
#ggsave(filename = paste0(path_project, "/04_ARTICLE/FIGURES/CA_analysis_compartiments_ellipse.png"),
#       bg="white", dpi = 300, width=10, height=10)

plot_grid(p0, p1, nrow=2)

ggsave(filename = paste0(path_project, "/04_ARTICLE/FIGURES/CA_analysis_compartiments_ellipse.png"),
       bg="white", dpi = 300, width=10, height=12)

ggsave(filename = paste0(path_project, "/04_ARTICLE/FIGURES_R1/FIGURE_8.tiff"),
       bg="white", dpi = 600, width=10, height=12)

ggsave(filename = paste0(path_project, "/04_ARTICLE/FIGURES_R1/FIGURE_8.png"),
       bg="white", dpi = 300, width=10, height=12)

t <- t %>% 
  mutate(compartiment_old_2 = if_else(compartiment_old=="Comp. 7", "Comp. 6", compartiment_old))

chi2_2 <- chisq.test(table(t$km_4, t$compartiment_old_2))
chi2_2$observed

# Save 
data_clean_processed <- X_scores1 %>% 
  dplyr::select(-starts_with("pca_score"), -km_4)

save(data_clean_processed, 
     file = paste0(path_project, "/04_ARTICLE/DATA_R1/00_Data_pre_processed.rda"))

results_fpca_kmeans <- X_scores1 %>% 
  dplyr::select(bananier, feuille, compartiment_old, fragment, starts_with("pca_score"), km_4) %>% 
  distinct()

save(results_fpca_kmeans, 
     file = paste0(path_project, "/04_ARTICLE/DATA_R1/01_Results_FPCA_kmeans.rda"))

results_CA <- coords_cols
save(results_CA, file = paste0(path_project, "/04_ARTICLE/DATA_R1/02_Results_CA.rda"))
