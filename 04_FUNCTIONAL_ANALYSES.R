
# Black Leaf Streak Disease (BLSD) dynamic on banana leaves

# Authors: M. Seidel, M. Chen, J. Avelino, C. Ababdie CIRAD, 2025

# ----------------------------------
# Packages

# > General data management 
library(tidyverse) # 

# > Functional analyses
library(tidyfun) # functional data visualization
library(refund)  # functional PCA and regressions

# > Graphics
library(cowplot) # add several panels in a plot

# ----------------------------------
path_project <- "D:/Mes Donnees/CIRAD/COLLABORATIONS/2024_MARINE/"

# Data loading 
load(paste0(path_project, "00_DONNEES/data_fragment.RData"))

# ----------------------------------
# 0. Preparing data for statistical analyses

# Get the first obs per fragment 
first_obs <- data_fragment %>% 
  group_by(bananier, feuille, fragment) %>% 
  summarise(first_date = first(date))

# Merge with the data to keep only the first obs
data_first_date <- data_fragment %>% 
  left_join(first_obs, by = c("bananier", "feuille", "fragment")) %>% 
  filter(date==first_date)

# Format data for functional analysis using refund

# > Deal with ID obs 
tab_ID <- data.frame(fragment_unique = unique(data_fragment$fragment_unique),
                     ID              = 1:length(unique(data_fragment$fragment_unique)))

# > Scalar variables 
X_init <- data_fragment %>% 
  dplyr::select(bananier, fragment_unique, "feuille_unique"="feuille", fragment, compartiment) %>% 
  # Fragment
  mutate(fragment = substr(fragment, 2, 4)) %>%
  mutate(fragment = as.numeric(as.character(fragment))) %>% 
  # > Leaf 
  mutate(feuille = substr(feuille_unique, 3,4),
         feuille = recode(feuille, "F1"="Leaf 1", "F2" = "Leaf 2", "F3"= "Leaf 3"),
         feuille = factor(feuille, levels = c("Leaf 1", "Leaf 2", "Leaf 3"))) %>% 
  # > Unique leaf
  mutate(feuille_unique = as.factor(feuille_unique)) %>% 
  # > Tree 
  mutate(bananier = recode(bananier, "B1"="Tree 1", "B2" = "Tree 2", "B3"= "Tree 3"),
         bananier = factor(bananier, levels = c("Tree 1", "Tree 2", "Tree 3"))) %>% 
  arrange(fragment_unique) %>% 
  distinct(.) %>% 
  # > ID of fragment
  left_join(., tab_ID) %>%
  # > nb of lesions at the 1rst date
  left_join(., data_first_date %>% 
              dplyr::select(fragment_unique, "nblesions1"="nblesions"))

str(X_init)
#'data.frame':	192 obs. of  7 variables:
#$ bananier       : Factor w/ 3 levels "Tree 1","Tree 2",..: 1 1 1 1 1 1 1 1 1 1 ...
#$ fragment_unique: Factor w/ 192 levels "P1H1z1","P1H1z10",..: 1 2 3 4 5 6 7 8 9 10 ...
#$ feuille_unique : Factor w/ 9 levels "B1F1","B1F2",..: 1 1 1 1 1 1 1 1 1 1 ...
#$ fragment       : Factor w/ 22 levels "z1","z10","z11",..: 1 2 3 4 5 6 7 8 9 10 ...
#$ compartiment   : Factor w/ 8 levels "1","2","3","4",..: 1 4 4 4 5 5 8 8 8 8 ...
#$ feuille        : Factor w/ 3 levels "Leaf 1","Leaf 2",..: 1 1 1 1 1 1 1 1 1 1 ...
#$ ID             : int  1 2 3 4 5 6 7 8 9 10 ...
#$ nblesions1     : int  21 20 12 25 8 8 4 10 18 16 ...

# > Functional sparse response in long format
Y_t <- data_fragment %>% 
  left_join(., tab_ID) %>% 
  dplyr::select(".id"= ID, 
                #".index"=jour_depuis_app_feuille, 
                ".index"=jour_depuis_rang4, 
                ".value"=surface) 

# > Create a tdf (tidy functional object) to ease the visualization using tidyfun package
Y_t_df <- Y_t |> 
  left_join(X_init, by=c(".id"="ID")) |>
  tf_nest(.value, .id = .id, .arg = .index)

class(Y_t_df) # data.frame 
class(Y_t_df$.value) # tdf data contained in the dataframe

# > Transformation to make sure that the data are not <0 or >3512
Y_t_trans <- data_fragment %>% 
  left_join(., tab_ID) %>%
  mutate(surface_trans=surface / (3512 + surface)) %>% 
  dplyr::select(".id"= ID, 
                #".index"=jour_depuis_app_feuille, 
                ".index"=jour_depuis_rang4, 
                ".value"=surface_trans) 

# > Functional sparse response in long format for analyses using the refund package
ydata <- Y_t
unique(ydata$.id) # should be numbers, length = number of rows in X
dim(X_init)        # 192 unique fragments

dim(ydata)         # 2332 observations in total

# > Temporal index on which disease is measured
#   in days after the appearance of the leaf
yindex <- sort(unique(ydata$.index)) ; yindex
#  -4 -3  0  3  4  7 10 11 14 17 18 21 24 25 28 31 32 35 38 39 42 45 46 49 52 56 59 60 63

# ----------------------------------
# 1. Difference in number of lesions at the first observation 
# between the compartments

# Plotting the data 
data_first_date  %>%
  ggplot(., aes(x = reorder(compartiment, -nblesions, median), y = nblesions)) + 
  geom_boxplot(aes(fill = compartiment), outlier.shape = NA) +
  geom_jitter(aes(#color = as.factor(fragment), 
    shape = fragment),
    width = 0.25,
    shape = 19) +
  theme_cowplot() +
  theme(legend.position = "bottom", 
        strip.text = element_text(angle=0),panel.grid = element_blank()) + 
  labs(x = "Compartment", 
       y = "Number of lesions") +
  scale_fill_manual(values = wesanderson::wes_palette("Zissou1", 8, "continuous"), name="Compartment")


# Save the plot
#ggsave(filename = paste0(path_project, "/04_ARTICLE/FIGURES/Nb_lesions_date_1_col.png"), 
#       width = 7, height = 5, bg = "white", dpi = 300)

# Model 



# ----------------------------------
# 2. Temporal patterns of BLSD progression on banana leaves

# 


# Plot of the raw data and the mean by compartment
plot_grid(
  
  # Raw data 
  ggplot() + 
    geom_hline(yintercept = c(0, 3512/2, 3512), color = "lightgrey", linetype = 2) +
    geom_meatballs(data = Y_t_df, aes(y = .value, color=compartiment)) + 
    facet_grid(.~bananier) +
    theme_cowplot() +
    theme(panel.grid.major.x = element_line(color = "lightgrey", linetype = 2),
          legend.position = "none") +
    ggtitle("a.") +
    labs(x = "Days after rank 4", y = "Infected surface (mm²)") +
    #scale_color_manual(values = wesanderson::wes_palette("Zissou1", 8, "continuous"), name="Compartment")
    scale_color_viridis_d(name="Compartment"),
  
  # Mean by leaves
  Y_t_df |>
    group_by(bananier, feuille, compartiment) |> 
    summarise(mean_surface = mean(.value, na.rm=T)) |>
    mutate(smooth_mean = tf_smooth(mean_surface)) |>
    ggplot() + 
    geom_hline(yintercept = c(0, 3512/2, 3512), color = "lightgrey", linetype = 2) +
    geom_spaghetti(aes(y = smooth_mean, color = compartiment), linewidth = 1, alpha=1) +
    facet_grid(.~bananier) +
    theme_cowplot() +
    theme(panel.grid.major.x = element_line(color = "lightgrey", linetype = 2),
          legend.position = "bottom") +
    ggtitle("b.") +
    labs(x = "Days after rank 4", y = "Infected surface (mm²)") +
    #scale_color_manual(values = wesanderson::wes_palette("Zissou1", 8, "continuous"), name="Compartment"),
    scale_color_viridis_d(name="Compartment"),
  
  ncol=1, rel_heights = c(0.45, 0.55)
)

# Save
ggsave(filename = paste0(path_project, "/04_ARTICLE/FIGURES/Dynamic_per_compartment.png"),
       bg="white", dpi = 300, width=9, height=9)


# Functional principal component analysis (FPCA)
# > perform the functional PCA
Fit.MM1 = fpca.sc(ydata  = ydata,  # data
                  nbasis = 9,
                  npc    = 4, 
                  simul  = T,
                  var    = TRUE,   # model-based estimates for the variance should be computed
                  center = TRUE)   # centering the data by retrieving an estimated mean function

# > estimated mean function (overal dynamic over the fragments)
#Fit.mu1 = data.frame(mu = (3512*Fit.MM1$mu)/(1-Fit.MM1$mu),
#                     d  = sort(unique(ydata$.index)))
Fit.mu1 = data.frame(mu = Fit.MM1$mu,
                     d  = sort(unique(ydata$.index)))
# > corresponding eigenfunctions 
#Fit.basis1 = data.frame(phi = (3512*Fit.MM1$efunctions)/(1-Fit.MM1$efunctions),
#                        d   = sort(unique(ydata$.index)))
#
Fit.basis1 = data.frame(phi = Fit.MM1$efunctions,
                        d   = sort(unique(ydata$.index)))

# > scores derivated from the functional PCs
X_scores1 <- Fit.MM1$scores %>% 
  as.data.frame(.) %>%
  rename("pca_score1"=1, "pca_score2"=2, "pca_score3"=3, "pca_score4"=4) %>%
  cbind(X_init, .) 

# > explained variability (in %)
round((Fit.MM1$evalues/sum(Fit.MM1$evalues))*100, 2)
# 82.15 14.81  2.18  0.86
round(cumsum(Fit.MM1$evalues/sum(Fit.MM1$evalues))*100, 2)
# 82.15  96.96  99.14 100.00
# the 2 first PCs explain 97% of the variability in the dataset

# > Plot explained variance (screeplot)
data.frame(phi = paste0("Composante ", 1:Fit.MM1$npc), 
           evalues = Fit.MM1$evalues,
           exp_var = round((Fit.MM1$evalues/sum(Fit.MM1$evalues))*100, 2),
           cum_exp_var = round(cumsum(Fit.MM1$evalues/sum(Fit.MM1$evalues))*100, 2)) %>% 
  ggplot(., aes(x=phi)) +
  geom_col(aes(y=exp_var)) +
  geom_point(aes(y=cum_exp_var), color="red") +
  geom_path(aes(y=cum_exp_var, group=1), color="red") +
  theme_bw() +
  labs(x="\nFunctional principal components", y="Explained variability in the data set (%)")

# > Plot mean and basis functions
# Which source of variability are represented by each function
Fit.basis1 %>% 
  left_join(Fit.mu1, by="d") %>% 
  gather(key="phi", value=value, -d, -mu) %>% 
  mutate(mu_plus_phi  = mu+5000*value,
         mu_minus_phi = mu-5000*value) %>% 
  gather(key=value_type, value=value, mu_plus_phi, mu_minus_phi) %>% 
  mutate(value_type = recode(value_type, 
                             'mu_plus_phi'="Mean + eigenfunction", 
                             'mu_minus_phi' = "Mean - eigenfunction"), 
         value_type = factor(value_type, levels=c("Mean + eigenfunction", "Mean - eigenfunction"))) %>%
  mutate(phi=recode(phi, 
                    "phi.1"=paste0("1st functional component (",  round((Fit.MM1$evalues[1]/sum(Fit.MM1$evalues))*100, 1), "%)"),
                    "phi.2"=paste0("2nd functional component (", round((Fit.MM1$evalues[2]/sum(Fit.MM1$evalues))*100, 1), "%)"),
                    "phi.3"=paste0("3rd functional component (",  round((Fit.MM1$evalues[3]/sum(Fit.MM1$evalues))*100, 1), "%)"),
                    "phi.4"=paste0("4th functional component (", round((Fit.MM1$evalues[4]/sum(Fit.MM1$evalues))*100, 1), "%)"))) %>%
  ggplot(.) + 
  geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
  geom_hline(yintercept = 0, color="darkgrey", lty=2) +
  geom_path(aes(x = d, y = mu, group=interaction(phi, value_type)), color="black") +
  geom_path(aes(x = d, y = value, color = value_type, group=interaction(phi, value_type))) +
  ggtitle("Which pattern of variability is explained by each principal component", 
          "(and how much does it explain the variability in the full dataset)") +
  facet_wrap(.~phi, ncol=2) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Days after rank 4", 
       y = "Infected area",
       color="") +
  scale_color_manual(values=c("red", "blue"))


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
      #Yhat_ex     <- (3512*Fit.MM1$Yhat[EX,])/(1-Fit.MM1$Yhat[EX,])
      #diag.var_EX <- (3512*Fit.MM1$diag.var[EX,])/(1-Fit.MM1$diag.var[EX,])
      
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
      labs(x = paste0(xlab), y = 'Infected surface') +
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
      arrange(desc(pca_score)) %>% 
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
      scale_color_gradientn(colors=c(viridis::rocket(15)[3:12], 
                                     viridis::mako(15, direction=-1)[3:12]), 
                            breaks = c(1, 20),
                            labels = c("Higher score", "Lower score"),
                            guide = guide_colorbar(title = paste0("Order in FPCA score ", score_lab),
                                                   barwidth = 15, barheight = 0.5)) +
      labs(x = paste0(xlab), y = 'Infected surface (mm²)') +
      theme_bw() +
      theme(legend.position="bottom", 
            legend.title.position = "top") +
      ggtitle(paste0(lab))
    
  }
  
  return(plot_frags)
  
}



# > Pca score 1
# Unique fragments with higher score
high_pca_score1 <- X_scores1 %>% 
  arrange(desc(pca_score1)) %>% 
  head(., 10) %>% 
  pull(ID) 

# Unique fragments with lowest score
low_pca_score1 <- X_scores1 %>% 
  arrange(pca_score1) %>% 
  head(., 10) %>% 
  pull(ID) 

# Plot
plot_dyn_id(frags = c(low_pca_score1, high_pca_score1), score="pca_score1") +
  geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
  geom_hline(yintercept = 0, color="darkgrey", lty=2)

# > Pca score 2
# Unique fragments with higher score
high_pca_score2 <- X_scores1 %>% 
  arrange(desc(pca_score2)) %>% 
  head(., 10) %>% 
  pull(ID) 

# Unique fragments with lowest scores
low_pca_score2 <- X_scores1 %>% 
  arrange(pca_score2) %>% 
  head(., 10) %>% 
  pull(ID) 

# Plot
plot_grid(
  plot_dyn_id(frags = c(low_pca_score2, high_pca_score2), 
              score="pca_score2", score_lab = "2") +
    geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
    geom_hline(yintercept = 0, color="darkgrey", lty=2),
  
  plot_dyn_id(frags = c(low_pca_score2, high_pca_score2), 
              score="pca_score1") +
    geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
    geom_hline(yintercept = 0, color="darkgrey", lty=2),
  
  nrow=1
)

# > Pca score 3
# Unique fragments with higher score
high_pca_score3 <- X_scores1 %>% 
  arrange(desc(pca_score3)) %>% 
  head(., 10) %>% 
  pull(ID) 

# Unique fragments with lowest score
low_pca_score3 <- X_scores1 %>% 
  arrange(pca_score3) %>% 
  head(., 10) %>% 
  pull(ID) 

# Plot
plot_dyn_id(frags = c(low_pca_score3, high_pca_score3), 
            score="pca_score3", score_lab = "3") +
  geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
  geom_hline(yintercept = 0, color="darkgrey", lty=2)

plot_grid(
  plot_dyn_id(frags = c(low_pca_score3, high_pca_score3), 
              score="pca_score3", score_lab = "3") +
    geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
    geom_hline(yintercept = 0, color="darkgrey", lty=2),
  
  plot_dyn_id(frags = c(low_pca_score3, high_pca_score3), 
              score="pca_score1") +
    geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
    geom_hline(yintercept = 0, color="darkgrey", lty=2),
  
  plot_dyn_id(frags = c(low_pca_score3, high_pca_score3), 
              score="pca_score2", score_lab = "2") +
    geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
    geom_hline(yintercept = 0, color="darkgrey", lty=2),
  
  nrow=1
)

# > Pca score 4
# Unique fragments with higher score
high_pca_score4 <- X_scores1 %>% 
  arrange(desc(pca_score4)) %>% 
  head(., 10) %>% 
  pull(ID) 

# Unique fragments with lowest score
low_pca_score4 <- X_scores1 %>% 
  arrange(pca_score4) %>% 
  head(., 10) %>% 
  pull(ID) 

# Plot
plot_dyn_id(frags = c(low_pca_score4, high_pca_score4), 
            score="pca_score4", score_lab = "4") +
  geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
  geom_hline(yintercept = 0, color="darkgrey", lty=2)



plot_grid(
  plot_dyn_id(frags = c(low_pca_score4, high_pca_score4), 
              score="pca_score4", score_lab = "4") +
    geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
    geom_hline(yintercept = 0, color="darkgrey", lty=2),
  
  plot_dyn_id(frags = c(low_pca_score4, high_pca_score4), 
              score="pca_score1") +
    geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
    geom_hline(yintercept = 0, color="darkgrey", lty=2),
  
  plot_dyn_id(frags = c(low_pca_score4, high_pca_score4), 
              score="pca_score2", score_lab = "2") +
    geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
    geom_hline(yintercept = 0, color="darkgrey", lty=2),
  
  plot_dyn_id(frags = c(low_pca_score4, high_pca_score4), 
              score="pca_score3", score_lab = "3") +
    geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
    geom_hline(yintercept = 0, color="darkgrey", lty=2),
  
  nrow=1
)



plot_grid(
  plot_dyn_id(frags = c(low_pca_score1, high_pca_score1), score="pca_score1") +
    geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
    geom_hline(yintercept = 0, color="darkgrey", lty=2) +
    ggtitle("a."),
  
  plot_dyn_id(frags = c(low_pca_score2, high_pca_score2), 
              score="pca_score2", score_lab = "2") +
    geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
    geom_hline(yintercept = 0, color="darkgrey", lty=2) +
    ggtitle("b."),
  
  plot_dyn_id(frags = c(low_pca_score3, high_pca_score3), 
              score="pca_score3", score_lab = "3") +
    geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
    geom_hline(yintercept = 0, color="darkgrey", lty=2) +
    ggtitle("c."),
  
  
  plot_dyn_id(frags = c(low_pca_score4, high_pca_score4), 
              score="pca_score4", score_lab = "4") +
    geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
    geom_hline(yintercept = 0, color="darkgrey", lty=2) +
    ggtitle("d."),
  
  nrow=2
  
)

ggsave(filename = paste0(path_project, "/04_ARTICLE/FIGURES/FPCA_interpretation.png"),
       bg="white", dpi = 300, width=9, height=9)

# ------------------------------
# K-means clustering for infected surface from FPCA
# packages for PCA and clustering
library(cluster)
library(FactoMineR)
library(factoextra)
library(corrplot)

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
sd(z_fpca_scores$pca_score4)


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
# > K= 3 or 4 or 5, depending on the methods

# > Perform clustering for K=4, or 5
# 4 functional components
set.seed(123)
KM3 <- kmeans(z_fpca_scores,  centers = 3, nstart = 25) 
set.seed(123)
KM4 <- kmeans(z_fpca_scores,  centers = 4, nstart = 25) 
set.seed(123)
KM5 <- kmeans(z_fpca_scores,  centers = 5, nstart = 25) 
# 2 first functional components

# > Add in initial data
X_scores1$km_3 <- KM3$cluster
X_scores1$km_3 <- as.factor(X_scores1$km_3)

X_scores1$km_4 <- KM4$cluster
X_scores1$km_4 <- as.factor(X_scores1$km_4)

X_scores1$km_5 <- KM5$cluster
X_scores1$km_5 <- as.factor(X_scores1$km_5)


extract_km_results <- function(KM)
{
  
  results <- data.frame(
    tot_sum_squares = KM$totss, 
    tot_within_sum_squares = KM$tot.withinss,
    between_sum_squares = KM$betweenss
  )
  
  return(results)
  
}

list(KM3, KM4, KM5) %>% 
  map_dfr(., extract_km_results, .id = "km")

KM3$size

# > Cluster visualization
# > per cluster
pal_cluster <- wesanderson::wes_palette("Zissou1", 5, "continuous")

plot_grid(
  
  ggplot(X_scores1) +
    geom_point(aes(x=pca_score1, y=pca_score2, color=km_3)) + 
    ggtitle("3 clusters") +
    scale_color_manual(values = pal_cluster, name = "Cluster") +
    labs(x="Score 1", y="Score 2") + theme_bw() + theme(legend.position = "none"),
  
  ggplot(X_scores1) +
    geom_point(aes(x=pca_score1, y=pca_score3, color=km_3)) + 
    scale_color_manual(values = pal_cluster, name = "Cluster") +
    labs(x="Score 1", y="Score 3") + theme_bw() + theme(legend.position = "bottom"),
  
  ggplot(X_scores1) +
    geom_point(aes(x=pca_score1, y=pca_score4, color=km_3)) + 
    scale_color_manual(values = pal_cluster, name = "Cluster") +
    labs(x="Score 1", y="Score 4") + theme_bw() + theme(legend.position = "none"),
  
  
  ggplot(X_scores1) +
    geom_point(aes(x=pca_score1, y=pca_score2, color=km_4)) + 
    scale_color_manual(values = pal_cluster, name = "Cluster") +
    ggtitle("4 clusters") +
    labs(x="Score 1", y="Score 2") + theme_bw() + theme(legend.position = "none"),
  
  ggplot(X_scores1) +
    geom_point(aes(x=pca_score1, y=pca_score3, color=km_4)) + 
    scale_color_manual(values = pal_cluster, name = "Cluster") +
    labs(x="Score 1", y="Score 3") + theme_bw() + theme(legend.position = "bottom"),
  
  ggplot(X_scores1) +
    geom_point(aes(x=pca_score1, y=pca_score4, color=km_4)) + 
    scale_color_manual(values = pal_cluster, name = "Cluster") +
    labs(x="Score 1", y="Score 4") + theme_bw() + theme(legend.position = "none"),
  
  ggplot(X_scores1) +
    geom_point(aes(x=pca_score1, y=pca_score2, color=km_5)) + 
    ggtitle("5 clusters") +
    scale_color_manual(values = pal_cluster, name = "Cluster") +
    labs(x="Score 1", y="Score 2") + theme_bw() +  theme(legend.position = "none"),
  
  ggplot(X_scores1) +
    geom_point(aes(x=pca_score1, y=pca_score3, color=km_5)) + 
    scale_color_manual(values = pal_cluster, name = "Cluster") +
    labs(x="Score 1", y="Score 3") + theme_bw() + theme(legend.position = "bottom"),
  
  ggplot(X_scores1) +
    geom_point(aes(x=pca_score1, y=pca_score4, color=km_5)) + 
    scale_color_manual(values = pal_cluster, name = "Cluster") +
    labs(x="Score 1", y="Score 4") + theme_bw() + theme(legend.position = "none"),
  nrow=3, axis = "tbrl", align = "hv"
  
)



# > per tree

plot_grid(
  
  ggplot(X_scores1) +
    geom_point(aes(x=pca_score1, y=pca_score2, color=bananier)) + 
    labs(x="Score 1", y="Score 2") + ggtitle("b. ") +
    scale_color_manual(values = pals::parula(3)) +
    theme_bw() + theme(legend.position = "none"),
  
  ggplot(X_scores1) +
    geom_point(aes(x=pca_score1, y=pca_score3, color=bananier )) + 
    scale_color_manual(values = pals::parula(3)) +
    labs(x="Score 1", y="Score 3") +
    theme_bw() + theme(legend.position = "none"),
  
  ggplot(X_scores1 %>%
           rename("Tree"="bananier")) +
    geom_point(aes(x=pca_score1, y=pca_score4, color=Tree )) + 
    scale_color_manual(values = pals::parula(3)) +
    labs(x="Score 1", y="Score 4") +
    theme_bw(), 
  
  ggplot(X_scores1) +
    geom_point(aes(x=pca_score2, y=pca_score3, color=bananier )) + 
    scale_color_manual(values = pals::parula(3)) +
    labs(x="Score 2", y="Score 3") +
    theme_bw() + theme(legend.position = "none"), 
  
  ggplot(X_scores1) +
    geom_point(aes(x=pca_score2, y=pca_score4, color=bananier )) + 
    scale_color_manual(values = pals::parula(3)) +
    labs(x="Score 2", y="Score 4") +
    theme_bw() + theme(legend.position = "none"), 
  
  ggplot(X_scores1) +
    geom_point(aes(x=pca_score3, y=pca_score4, color=bananier )) + 
    scale_color_manual(values = pals::parula(3)) +
    labs(x="Score 3", y="Score 4") +
    theme_bw() + theme(legend.position = "none"), 
  
  nrow=2, axis = "tbrl", align = "hv"
  
)



plot_grid(
  
  ggplot(X_scores1) +
    geom_point(aes(x=pca_score1, y=pca_score2, color=km_3)) + 
    ggtitle("3 clusters") +
    scale_color_manual(values = pal_cluster, name = "Cluster") +
    labs(x="Score 1", y="Score 2") + theme_bw() + theme(legend.position = "none"),
  
  ggplot(X_scores1) +
    geom_point(aes(x=pca_score1, y=pca_score3, color=km_3)) + 
    scale_color_manual(values = pal_cluster, name = "Cluster") +
    labs(x="Score 1", y="Score 3") + theme_bw() + theme(legend.position = "bottom"),
  
  ggplot(X_scores1) +
    geom_point(aes(x=pca_score1, y=pca_score4, color=km_3)) + 
    scale_color_manual(values = pal_cluster, name = "Cluster") +
    labs(x="Score 1", y="Score 4") + theme_bw() + theme(legend.position = "none"),
  
  
  ggplot(X_scores1) +
    geom_point(aes(x=pca_score1, y=pca_score2, color=km_4)) + 
    scale_color_manual(values = pal_cluster, name = "Cluster") +
    ggtitle("4 clusters") +
    labs(x="Score 1", y="Score 2") + theme_bw() + theme(legend.position = "none"),
  
  ggplot(X_scores1) +
    geom_point(aes(x=pca_score1, y=pca_score3, color=km_4)) + 
    scale_color_manual(values = pal_cluster, name = "Cluster") +
    labs(x="Score 1", y="Score 3") + theme_bw() + theme(legend.position = "bottom"),
  
  ggplot(X_scores1) +
    geom_point(aes(x=pca_score1, y=pca_score4, color=km_4)) + 
    scale_color_manual(values = pal_cluster, name = "Cluster") +
    labs(x="Score 1", y="Score 4") + theme_bw() + theme(legend.position = "none"),
  
  ggplot(X_scores1) +
    geom_point(aes(x=pca_score1, y=pca_score2, color=km_5)) + 
    ggtitle("5 clusters") +
    scale_color_manual(values = pal_cluster, name = "Cluster") +
    labs(x="Score 1", y="Score 2") + theme_bw() +  theme(legend.position = "none"),
  
  ggplot(X_scores1) +
    geom_point(aes(x=pca_score1, y=pca_score3, color=km_5)) + 
    scale_color_manual(values = pal_cluster, name = "Cluster") +
    labs(x="Score 1", y="Score 3") + theme_bw() + theme(legend.position = "bottom"),
  
  ggplot(X_scores1) +
    geom_point(aes(x=pca_score1, y=pca_score4, color=km_5)) + 
    scale_color_manual(values = pal_cluster, name = "Cluster") +
    labs(x="Score 1", y="Score 4") + theme_bw() + theme(legend.position = "none"),
  nrow=3, axis = "tbrl", align = "hv"
  
)





# Plot the curves for each clusters
# K=3
plot_grid(
  plot_dyn_id(frags = X_scores1 %>% filter(km_3==1) %>% pull(ID), score = NULL, color = pal_cluster[1]) +
    geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
    geom_hline(yintercept = 0, color="darkgrey", lty=2) +
    ggtitle(paste0("Cluster 1 (N=", paste0(KM3$size[1]), ")")) +
    theme(legend.position = "none"),
  
  plot_dyn_id(frags = X_scores1 %>% filter(km_3==3) %>% pull(ID), score = NULL, color = pal_cluster[3]) +
    geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
    geom_hline(yintercept = 0, color="darkgrey", lty=2) +
    ggtitle(paste0("Cluster 3 (N=", paste0(KM3$size[3]), ")")) +
    theme(legend.position = "none"),
  
  plot_dyn_id(frags = X_scores1 %>% filter(km_3==2) %>% pull(ID), score = NULL, color = pal_cluster[2]) +
    geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
    geom_hline(yintercept = 0, color="darkgrey", lty=2) +
    ggtitle(paste0("Cluster 2 (N=", paste0(KM3$size[2]), ")")) +
    theme(legend.position = "none"),
  
  nrow=1
)

# K=4
plot_grid(
  plot_dyn_id(frags = X_scores1 %>% filter(km_4==1) %>% pull(ID), score = NULL, color = pal_cluster[1]) +
    geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
    geom_hline(yintercept = 0, color="darkgrey", lty=2) +
    ggtitle(paste0("Cluster 1 (N=", paste0(KM4$size[1]), ")")) +
    theme(legend.position = "none"),
  
  plot_dyn_id(frags = X_scores1 %>% filter(km_4==4) %>% pull(ID), score = NULL, color = pal_cluster[4]) +
    geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
    geom_hline(yintercept = 0, color="darkgrey", lty=2) +
    ggtitle(paste0("Cluster 4 (N=", paste0(KM4$size[4]), ")")) +
    theme(legend.position = "none"),
  
  plot_dyn_id(frags = X_scores1 %>% filter(km_4==2) %>% pull(ID), score = NULL, color = pal_cluster[2]) +
    geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
    geom_hline(yintercept = 0, color="darkgrey", lty=2) +
    ggtitle(paste0("Cluster 2 (N=", paste0(KM4$size[2]), ")")) +
    theme(legend.position = "none"),
  
  plot_dyn_id(frags = X_scores1 %>% filter(km_4==3) %>% pull(ID), score = NULL, color = pal_cluster[3]) +
    geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
    geom_hline(yintercept = 0, color="darkgrey", lty=2) +
    ggtitle(paste0("Cluster 3 (N", paste0(KM4$size[3]), ")")) +
    theme(legend.position = "none"),
  
  nrow=1
)

ggsave(filename = paste0(path_project, "/04_ARTICLE/FIGURES/Dynamic_Clusters_4.png"),
       bg="white", dpi = 300, width=18, height=5)

# K=5
plot_grid(
  plot_dyn_id(frags = X_scores1 %>% filter(km_5==1) %>% pull(ID), score = NULL, color = pal_cluster[1]) +
    geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
    geom_hline(yintercept = 0, color="darkgrey", lty=2) +
    ggtitle("Cluster 1") +
    theme(legend.position = "none"),
  
  plot_dyn_id(frags = X_scores1 %>% filter(km_5==5) %>% pull(ID), score = NULL, color = pal_cluster[5]) +
    geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
    geom_hline(yintercept = 0, color="darkgrey", lty=2) +
    ggtitle("Cluster 5") +
    theme(legend.position = "none"),
  
  plot_dyn_id(frags = X_scores1 %>% filter(km_5==4) %>% pull(ID), score = NULL, color = pal_cluster[4]) +
    geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
    geom_hline(yintercept = 0, color="darkgrey", lty=2) +
    ggtitle("Cluster 4") +
    theme(legend.position = "none"),
  
  plot_dyn_id(frags = X_scores1 %>% filter(km_5==3) %>% pull(ID), score = NULL, color = pal_cluster[3]) +
    geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
    geom_hline(yintercept = 0, color="darkgrey", lty=2) +
    ggtitle("Cluster 3") +
    theme(legend.position = "none"),
  
  plot_dyn_id(frags = X_scores1 %>% filter(km_5==2) %>% pull(ID), score = NULL, color = pal_cluster[2]) +
    geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
    geom_hline(yintercept = 0, color="darkgrey", lty=2) +
    ggtitle("Cluster 2") +
    theme(legend.position = "none"),
  
 nrow=1
)


# -----------------------------------------
# Functional regression 

# STEP 1
# Data for functional regression

# > Functional sparse response in long format
ydata2 <- Y_t %>% rename(".obs"=".id")
unique(ydata2$.obs) # should be numbers, length = number of rows in X

# > Temporal index on which disease is measured
yindex2 <- sort(unique(ydata2$.index)) ; yindex2
# -4 -3  0  3  4  7 10 11 14 17 18 21 24 25 28 31 32 35 38 39 42 45 46 49 52 56 59 60 63

# > Functional sparse response in wide format
Y_t_pivot <- Y_t %>% 
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
  mutate(km_3 = factor(km_3, levels = c("1", "3", "2"))) %>% 
  mutate(KM3_1 = if_else(km_3==1, 1, 0),
         KM3_2 = if_else(km_3==2, 1, 0),
         KM3_3 = if_else(km_3==3, 1, 0)) %>% 
  mutate(km_4 = factor(km_4, levels = c("1", "4", "2", "3"))) %>% 
  mutate(KM4_1 = if_else(km_4==1, 1, 0),
         KM4_2 = if_else(km_4==2, 1, 0),
         KM4_3 = if_else(km_4==3, 1, 0),
         KM4_4 = if_else(km_4==4, 1, 0)) %>% 
  mutate(km_5 = factor(km_5, levels = c("1", "5", "4", "3", "2"))) %>% 
  mutate(KM5_1 = if_else(km_5==1, 1, 0),
         KM5_2 = if_else(km_5==2, 1, 0),
         KM5_3 = if_else(km_5==3, 1, 0),
         KM5_4 = if_else(km_5==4, 1, 0),
         KM5_5 = if_else(km_5==5, 1, 0))

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

data$KM3_1 <- X$KM3_1
data$KM3_2 <- X$KM3_2
data$KM3_3 <- X$KM3_3
data$km_3  <- factor(X$km_3)

data$KM4_1 <- X$KM4_1
data$KM4_2 <- X$KM4_2
data$KM4_3 <- X$KM4_3
data$KM4_4 <- X$KM4_4
data$km_4 <- factor(X$km_4)

data$KM5_1 <- X$KM5_1
data$KM5_2 <- X$KM5_2
data$KM5_3 <- X$KM5_3
data$KM5_4 <- X$KM5_4
data$KM5_5 <- X$KM5_5
data$km_5 <- X$km_5

# nb of lesions at the first date
data$nblesions1 <- as.numeric(as.character(X$nblesions1))

# nb de lesions
data$nblesions <- as.numeric(as.character(X$nblesions))

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
# R-sq.(adj) = 0.663  ; REML score: 18838 ; AIC = 37693.91  

# adjusted for tree * leaf (as fixed effects)
mod_0_bl <- pffr(Y_t ~ T1L2 + T1L3 + T2L1 + T2L2 + T2L3 + T3L1 + T3L2 + T3L3, 
                 data = data, 
                 ydata = ydata2,
                 bs.yindex = list(bs = "ps", k=27),
                 bs.int = list(bs = "ps", k=27))

summary(mod_0_bl) ; AIC(mod_0_bl)
# R-sq.(adj) =  0.781 ; REML score =  18251 ; AIC = 36748.26

# adjusted for cluster and tree * leaf (as fixed effects) 
mod_0_blc3 <- pffr(Y_t ~ KM3_2 + KM3_3 +
                    T1L2 + T1L3 + T2L1 + T2L2 + T2L3 + T3L1 + T3L2 + T3L3, 
                 data = data, 
                 ydata = ydata2,
                 bs.yindex = list(bs = "ps", k=27),
                 bs.int = list(bs = "ps", k=27))

summary(mod_0_blc3) ; AIC(mod_0_blc3)
# R-sq.(adj) =  0.924 ; REML score =  17039 ; AIC = 34314.83

mod_0_blc4 <- pffr(Y_t ~ KM4_4 + KM4_2 + KM4_3 + 
                     T1L2 + T1L3 + T2L1 + T2L2 + T2L3 + T3L1 + T3L2 + T3L3, 
                   data = data, 
                   ydata = ydata2,
                   bs.yindex = list(bs = "ps", k=27),
                   bs.int = list(bs = "ps", k=27))

summary(mod_0_blc4) ; AIC(mod_0_blc4)
# R-sq.(adj) =  0.93 ; REML score =  16921 ; AIC = 34109.3

mod_0_blc5 <- pffr(Y_t ~ KM5_5 + KM5_4 + KM5_3 + KM5_2 + 
                     T1L2 + T1L3 + T2L1 + T2L2 + T2L3 + T3L1 + T3L2 + T3L3, 
                   data = data, 
                   ydata = ydata2,
                   bs.yindex = list(bs = "ps", k=27),
                   bs.int = list(bs = "ps", k=27))

summary(mod_0_blc5) ; AIC(mod_0_blc5)
# R-sq.(adj) =  0.942 ; REML score =  16709 ; AIC = 33689.85

# adjusted for cluster, nb of lesions at the first date, and tree * leaf (as fixed effects) 
mod_0_blc3l <- pffr(Y_t ~ KM3_2 + KM3_3 +
                    nblesions1 + 
                    T1L2 + T1L3 + T2L1 + T2L2 + T2L3 + T3L1 + T3L2 + T3L3, 
                  data = data, 
                  ydata = ydata2,
                  bs.yindex = list(bs = "ps", k=27),
                  bs.int = list(bs = "ps", k=27))

summary(mod_0_blcl) ; AIC(mod_0_blcl)
# R-sq.(adj) =  0.928 ; REML score =  16968 ; AIC = 34181.8

mod_0_blc4l <- pffr(Y_t ~ KM4_4 + KM4_2 + KM4_3 + 
                      nblesions1 + 
                     T1L2 + T1L3 + T2L1 + T2L2 + T2L3 + T3L1 + T3L2 + T3L3, 
                   data = data, 
                   ydata = ydata2,
                   bs.yindex = list(bs = "ps", k=27),
                   bs.int = list(bs = "ps", k=27))

summary(mod_0_blc4l) ; AIC(mod_0_blc4l)
# R-sq.(adj) =  0.938  ; REML score =  16791   ; AIC = 33848.2

mod_0_blc5l <- pffr(Y_t ~ KM5_5 + KM5_4 + KM5_3 + KM5_2 + 
                      nblesions1 + 
                     T1L2 + T1L3 + T2L1 + T2L2 + T2L3 + T3L1 + T3L2 + T3L3, 
                   data = data, 
                   ydata = ydata2,
                   bs.yindex = list(bs = "ps", k=27),
                   bs.int = list(bs = "ps", k=27))

summary(mod_0_blc5l) ; AIC(mod_0_blc5l)
# R-sq.(adj) =  0.947 ; REML score =  16603 ; AIC = 33478.07

# Similar models but with tree*leaf as random effects
mod_1_bl <- pffr(Y_t ~ 0 + c(0) +      # dropping the constant and time-varying intercepts via 0 + c(0) to estimate zone + c(zone) 
                   s(tree_leaf, bs="re"),   # leaf*tree-specific smooth residuals
                 data = data, 
                 ydata = ydata2,
                 bs.yindex = list(bs = "ps", k=5),
                 bs.int = list(bs = "ps", k=20)) 

summary(mod_1_bl) ; AIC(mod_1_bl)

# including clusters 
mod_1_blc3 <- pffr(Y_t ~ 0 + c(0) +           # dropping the constant and time-varying intercepts via 0 + c(0) to estimate zone + c(zone) 
                   km_3 + c(km_3) +           # time-varying cluster effects not centered at zero 
                   s(tree_leaf, bs="re"),     # leaf*tree-specific smooth residuals
                 data = data, 
                 ydata = ydata2,
                 bs.yindex = list(bs = "ps", k=5),
                 bs.int = list(bs = "ps", k=20)) 

summary(mod_1_blc3) ; AIC(mod_1_blc3)
# R-sq.(adj) =  0.919 ; REML score =  17215 ; AIC = 34417.85 

mod_1_blc4 <- pffr(Y_t ~ 0 + c(0) +           # dropping the constant and time-varying intercepts via 0 + c(0) to estimate zone + c(zone) 
                     km_4 + c(km_4) +           # time-varying cluster effects not centered at zero 
                     s(tree_leaf, bs="re"),     # leaf*tree-specific smooth residuals
                   data = data, 
                   ydata = ydata2,
                   bs.yindex = list(bs = "ps", k=5),
                   bs.int = list(bs = "ps", k=20)) 

summary(mod_1_blc4) ; AIC(mod_1_blc4)
# R-sq.(adj) =  0.925 ; REML score =  17114 ; AIC = 34240.57

mod_1_blc5 <- pffr(Y_t ~ 0 + c(0) +           # dropping the constant and time-varying intercepts via 0 + c(0) to estimate zone + c(zone) 
                     km_5 + c(km_5) +           # time-varying cluster effects not centered at zero 
                     s(tree_leaf, bs="re"),     # leaf*tree-specific smooth residuals
                   data = data, 
                   ydata = ydata2,
                   bs.yindex = list(bs = "ps", k=5),
                   bs.int = list(bs = "ps", k=20)) 

summary(mod_1_blc5) ; AIC(mod_1_blc5)
# R-sq.(adj) =  0.936 ; REML score =  16922 ; AIC = 33868.92

# additionally adjusted for the initial number of lesions
mod_1_blc3l <- pffr(Y_t ~ 0 + c(0) +           # dropping the constant and time-varying intercepts via 0 + c(0) to estimate zone + c(zone) 
                      km_3 + c(km_3) +           # time-varying cluster effects not centered at zero 
                      nblesions1 + 
                      s(tree_leaf, bs="re"),     # leaf*tree-specific smooth residuals
                   data = data, 
                   ydata = ydata2,
                   bs.yindex = list(bs = "ps", k=5),
                   bs.int = list(bs = "ps", k=20)) 

summary(mod_1_blc3l) ; AIC(mod_1_blc3l)
# R-sq.(adj) =  0.923 ; REML score =  17151 ; AIC = 34292.2 

mod_1_blc4l <- pffr(Y_t ~ 0 + c(0) +           # dropping the constant and time-varying intercepts via 0 + c(0) to estimate zone + c(zone) 
                      km_4 + c(km_4) +           # time-varying cluster effects not centered at zero 
                      nblesions1 + 
                      s(tree_leaf, bs="re"),     # leaf*tree-specific smooth residuals
                    data = data, 
                    ydata = ydata2,
                    bs.yindex = list(bs = "ps", k=5),
                    bs.int = list(bs = "ps", k=20)) 

summary(mod_1_blc4l) ; AIC(mod_1_blc4l)
# R-sq.(adj) =  0.932 ; REML score =  16995 ; AIC = 34001.02 

mod_1_blc5l <- pffr(Y_t ~ 0 + c(0) +           # dropping the constant and time-varying intercepts via 0 + c(0) to estimate zone + c(zone) 
                      km_5 + c(km_5) +           # time-varying cluster effects not centered at zero 
                      nblesions1 + 
                      s(tree_leaf, bs="re"),     # leaf*tree-specific smooth residuals
                    data = data, 
                    ydata = ydata2,
                    bs.yindex = list(bs = "ps", k=5),
                    bs.int = list(bs = "ps", k=20)) 

summary(mod_1_blc5l) ; AIC(mod_1_blc5l)
# R-sq.(adj) =  0.941 ; REML score =  16829 ; AIC = 33685.95 

# Extraction of the coefficients for the clusters
list("mod_1_blc3"=  mod_1_blc3,
     "mod_1_blc4"=  mod_1_blc4,
     "mod_1_blc5"=  mod_1_blc5,
     "mod_1_blc3l"= mod_1_blc3l,
     "mod_1_blc4l"= mod_1_blc4l,
     "mod_1_blc5l"= mod_1_blc5l) %>% 
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
    facet_grid(adjustment_lesions~nb_cluster) +
    scale_color_manual(values=pal_cluster, name="") +
    scale_fill_manual(values=pal_cluster, name="") +
    theme_bw() +
    theme(legend.position = "bottom") +
    labs(x = "Days after rank 4", y = "Mean functional coefficient (95% confidence interval)")




list("mod_1_blc3"=  mod_1_blc3,
     "mod_1_blc4"=  mod_1_blc4,
     "mod_1_blc5"=  mod_1_blc5,
     "mod_1_blc3l"= mod_1_blc3l,
     "mod_1_blc4l"= mod_1_blc4l,
     "mod_1_blc5l"= mod_1_blc5l) %>% 
  map(., ~{
    
    pred <- 
    
    
    
  })

pred_mod_1_blc3 <- predict(mod_1_blc3, se.fit = TRUE, type = "response",
        newdata = data.frame(km_3 = c(1,2,3), 
                             tree_leaf=c("B1F1","B1F1","B1F1")))

pred_mod_1_blc3l <- predict(mod_1_blc3l, se.fit = TRUE, type = "response",
                           newdata = data.frame(km_3 = c(1,2,3), 
                                                nblesions1 = rep(mean(data_first_date$nblesions), 3),
                                                tree_leaf=c("B1F1","B1F1","B1F1")))


pred_mod_1_blc4 <- predict(mod_1_blc4, se.fit = TRUE, type = "response",
                           newdata = data.frame(km_4 = c(1,2,3, 4), 
                                                tree_leaf=c("B1F1","B1F1","B1F1","B1F1")))

pred_mod_1_blc4l <- predict(mod_1_blc4l, se.fit = TRUE, type = "response",
                            newdata = data.frame(km_4 = c(1,2,3, 4), 
                                                 nblesions1 = rep(mean(data_first_date$nblesions), 4),
                                                 tree_leaf=c("B1F1","B1F1","B1F1","B1F1")))


pred_mod_1_blc5 <- predict(mod_1_blc5, se.fit = TRUE, type = "response",
                           newdata = data.frame(km_5 = c(1,2,3, 4, 5), 
                                                tree_leaf=c("B1F1","B1F1","B1F1","B1F1","B1F1")))

pred_mod_1_blc5l <- predict(mod_1_blc5l, se.fit = TRUE, type = "response",
                            newdata = data.frame(km_5 = c(1,2,3, 4, 5), 
                                                 nblesions1 = rep(mean(data_first_date$nblesions), 5),
                                                 tree_leaf=c("B1F1","B1F1","B1F1","B1F1","B1F1")))

data.frame(
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
  ggplot(., aes(x = x_i, color=cluster)) +
  geom_hline(yintercept = c(0, 3512), color="black", lty=2) +
  # > mean 
  geom_path(data = Fit.mu1, 
            aes(x=d, y=mu), color="black", lwd=1) +
  geom_line(aes(y=f_i), linewidth=1.2) +
  geom_ribbon(aes(ymin=ic_down, ymax=ic_up, fill=cluster), color="transparent", alpha=0.3) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x="Days after rank 4", y="Mean infected surface in mm² (and)95% confidence interval) ") +
  facet_grid(adjustment~model) +
  guides(color = guide_legend(nrow=1)) +
  lims(y=c(-100, 3600)) +
  scale_color_manual(values=pal_cluster, name="") +
  scale_fill_manual(values=pal_cluster, name="")




# 


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





mod_2_blc4l <- pffr(Y_t ~ 0 + c(0)   +           # dropping the constant and time-varying intercepts via 0 + c(0) to estimate zone + c(zone) 
                      km_4 + c(km_4) +           # time-varying cluster effects not centered at zero 
                      s(nblesions)   + 
                      s(tree_leaf, bs="re"),     # leaf*tree-specific smooth residuals
                    data = data, 
                    ydata = ydata2,
                    bs.yindex = list(bs = "ps", k=5),
                    bs.int = list(bs = "ps", k=20)) 

summary(mod_2_blc4l) ; AIC(mod_2_blc4l)



mod_2_blc4l$coefficients



