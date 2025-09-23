
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
# Function for plotting for a set of individual fragments
# the continuous dynamic reconstructed from observations 
# using the FPCA 

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
      theme_bw(base_size = 15) +
      theme(legend.position="bottom", 
            legend.title.position = "top",
            axis.text = element_text(color="black"),
            title = element_text(size=12)) +
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
      theme_bw(base_size = 15) +
      theme(legend.position="bottom", 
            legend.title.position = "top",
            axis.text = element_text(color="black")) +
      ggtitle(paste0(lab))
    
  }
  
  return(plot_frags)
  
}

# ----------------------------------
# Clusters' characteristics from K-means
# extract metrics about the clustering

extract_km_results <- function(KM)
{
  
  results <- data.frame(
    tot_sum_squares = KM$totss, 
    tot_within_sum_squares = KM$tot.withinss,
    between_sum_squares = KM$betweenss
  )
  
  return(results)
  
}


# ----------------------------------
# Home made function to perfom Tuckey tests
# used to examine the differences between clusters in terms of FPCA scores 
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


# ----------------------------------
# Home made function to attribute correct letter to each fragment based on the GLMM fitting
# to examine the differences between clusters in terms of number of first lesions
fonction_plot_ameans <- function(){
  df.2 <- data_fragment_first_lesions
  ameans <- tuk.cld
  df_tukey <- data.frame(ameans$.group)
  rownames(df_tukey) <- ameans$km_4
  df.2$Groupes=rep(NA)
  I <- length(rownames(df_tukey))
  K <- length(df.2$km_4)
  for ( k in 1:K) {
    for (i in 1 :I) {
      if ( df.2$km_4[k] == rownames(df_tukey)[i]) {
        df.2$Groupes[k]<- df_tukey[i,1] }
    }
  }
  return(df.2)}

