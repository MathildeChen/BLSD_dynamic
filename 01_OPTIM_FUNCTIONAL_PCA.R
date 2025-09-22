
# Black Leaf Streak Disease (BLSD) dynamic on banana leaves

# Authors: M. Seidel, M. Chen, J. Avelino, C. Ababdie CIRAD, 2025

# ----------------------------------
# Packages

# > General data management 
library(tidyverse) # 

# > Functional analyses
library(caret)   # measure of fit quality 
library(refund)  # functional PCA and regressions

# > Graphics
library(cowplot) # add several panels in a plot

# ----------------------------------
path_project <- "D:/Mes Donnees/CIRAD/COLLABORATIONS/2024_MARINE/"

# Data loading 
load(paste0(path_project, "00_DONNEES/data_fragment_for_analyses.rda"))
ydata <- data_fragment_for_analyses$ydata

# ----------------------------------
# > Choose the optimal number of basis and the optimal number of functions
#   i.e. the combination giving the lowest error when re-estimating the initial data set from the compressed one
list_error <- NULL

# Test of the number of basis and the number of functions used
# > select a value of number of functions
for(n_functions in 1:6)
{
  # > select a value of number of basis (B-splines)
  for(n_basis in seq(5, 25, by=1))
  {
    
    # > perform the functional PCA
    Fit.MM_test = fpca.sc(ydata  = ydata,       # data
                          nbasis = n_basis,     # number of basis 
                          npc    = n_functions, # force the number of functions
                          #var    = TRUE,        # model-based estimates for the variance should be computed
                          center = TRUE)        # centering the data by retrieving an estimated mean function
    
    # > reconstruction of the curves for all individuals from the FPCA 
    #   to estimate the error of reconstruction
    error <- NULL
    
    for(EX in unique(ydata$.id))
    {
      
      # > Observed data points
      ydata1EX <- ydata %>% filter(.id==EX)
      
      # > Fitted data from FPCA
      EX.MM = data.frame(fitted    = Fit.MM_test$Yhat[EX,],
                         d         = sort(unique(ydata$.index)))
      
      # > Merge
      errorEX <- ydata1EX %>% 
        left_join(EX.MM, by=c(".index"="d")) %>% 
        mutate(.value_est = fitted)
      
      error <- rbind(error, errorEX)
      
    }
    
    # > Add data to other simulations results
    list_error[[paste0("test_", n_functions, "_", n_basis)]] <- error
    
  }
  
}

# > Compute error of reconstruction for each combination of n_basis and n_functions and 
list_error %>% 
  map_dfr(.,~{ 
    
    r2   <- R2(obs = .x$.value, pred = .x$.value_est, na.rm=T, )
    RMSE <- RMSE(obs = .x$.value, pred = .x$.value_est, na.rm=T)
    
    data.frame(r2, RMSE)
    
    
  }, .id="test") %>% 
  separate(test, c("to_remove", "n_functions", "n_basis"), sep="_") %>%
  dplyr::select(-to_remove) %>% 
  gather(key=perf_metric, value=perf_value,  r2, RMSE) %>% 
  mutate(n_basis = factor(n_basis, levels=1:25)) %>% 
  mutate(line_to_plot = if_else(perf_metric=="r2", 0.98, 10*(3512/100))) %>% 
  ggplot(., aes(x=n_functions, y=perf_value, color=n_basis, group=n_basis)) +
  geom_hline(aes(yintercept = line_to_plot)) +
  geom_point() + geom_path() +
  facet_wrap(.~perf_metric,scales = "free") +
  theme_bw() +
  labs(x="Number of functional principal component", y="Quality of fit", color="Number of B-splines")

# > Compute error of reconstruction for each combination of n_basis and n_functions and 
R2 <- list_error %>% 
  map_dfr(.,~{ 
    
    r2   <- R2(obs = .x$.value, pred = .x$.value_est, na.rm=T, )
    RMSE <- RMSE(obs = .x$.value, pred = .x$.value_est, na.rm=T)
    
    data.frame(r2, RMSE)
    
    
  }, .id="test") %>% 
  separate(test, c("to_remove", "n_functions", "n_basis"), sep="_") %>%
  dplyr::select(-to_remove) %>% 
  mutate(n_basis = factor(n_basis, levels=1:25))

p1 <- R2 %>%
  ggplot(., aes(x=n_functions, y=r2, color=n_basis, group=n_basis)) +
  geom_hline(yintercept = 0.98, color = "darkgrey", linetype=2) +
  geom_point() + geom_path() +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x="Number of functional principal component", y="Quality of curve re-estimation by the FPCA", color="Number of B-splines")

p2 <- R2 %>% 
  filter(n_functions %in% c(2,3),
         r2 > 0.98) %>% 
  ggplot(., aes(x=n_functions, y=r2, color=n_basis, group=n_basis)) +
  geom_point() + geom_path(linetype=2) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x="Number of functional principal component", y="Quality of curve re-estimation\nby the FPCA", color="Number of B-splines")

plot_grid(p1,p2)


ggsave(p1,
       filename = paste0(path_project, "/04_ARTICLE/FIGURES/FPCA_optt_supp_figures.png"),
       bg="white", dpi = 300, width=6, height=6)


ggsave(p2,
       filename = paste0(path_project, "/04_ARTICLE/FIGURES/FPCA_optt_supp_figures_2.png"),
       bg="white", dpi = 300, width=6, height=4)


# ----------------------------------
# Functional principal component analysis (FPCA)
# > perform the functional PCA with optimal setting
Fit.MM1 = fpca.sc(ydata  = ydata,  # data
                  nbasis = 9,
                  npc    = 2, 
                  simul  = T,
                  var    = TRUE,   # model-based estimates for the variance should be computed
                  center = TRUE)   # centering the data by retrieving an estimated mean function

# > explained variability (in %)
round((Fit.MM1$evalues/sum(Fit.MM1$evalues))*100, 2)
# 85.49 12.22  2.29

round(cumsum(Fit.MM1$evalues/sum(Fit.MM1$evalues))*100, 2)
# 85.49  97.71 100.00
# the 2 first PCs explain 97.7% of the variability in the dataset

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

ggsave(filename = paste0(path_project, "/04_ARTICLE/FIGURES/FPCA_screeplot_supp_figures.png"),
       bg="white", dpi = 300, width=5, height=5)









