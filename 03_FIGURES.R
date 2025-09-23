
# Black Leaf Streak Disease (BLSD) dynamic on banana leaves

# 03. Script to produce the figures of the paper

# Authors: M. Seidel, M. Chen, J. Avelino, C. Ababdie CIRAD, 2025

# ----------------------------------
# Path to project
path_project <- "D:/Mes Donnees/CIRAD/COLLABORATIONS/2024_MARINE/"
path_figures_project <- "D:/Mes Donnees/CIRAD/COLLABORATIONS/2024_MARINE/04_ARTICLE/FIGURES_R1/"

# ----------------------------------
# Run the analyses to generate input tables and objects for the figures
source(paste0(path_project, "01_SCRIPT/SCRIPTS_FOR_Github/02_STATISTICAL_ANALYSES.R"))

# Color palette
pal_cluster <- wesanderson::wes_palette("Zissou1", 5, "discrete")[c(1,3,4,5)]

# ----------------------------------
# Figure 2
# Variation in the disease progressions

Figure2 <- ydata %>%
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
  labs(x = "Days after leaf reached rank 4", y = "Diseased surface (mm²)") +
  theme_cowplot() +
  theme(panel.border = element_rect(color="black"),
        strip.background = element_blank(),
        panel.grid.major = element_line(color="grey87")) ; Figure2

ggsave(plot = Figure2,
       filename = paste0(path_figures_project, "FIGURE_2.tiff"),
       bg="white", dpi = 600, width=10, height=10)


# ----------------------------------
# Figure 3: 
# Fragments with the highest and lowest values for each score 

# Choose the number of fragments to plot
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
p3a <- plot_dyn_id(frags = frags_pca_score1, score="pca_score1") +
  geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
  geom_hline(yintercept = 0, color="darkgrey", lty=2)

p3b <- plot_dyn_id(frags = frags_pca_score2, 
                  score="pca_score2", score_lab = "2") +
  geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
  geom_hline(yintercept = 0, color="darkgrey", lty=2) 

p3c <- plot_dyn_id(frags = frags_pca_score3, 
                  score="pca_score3", score_lab = "3") +
  geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
  geom_hline(yintercept = 0, color="darkgrey", lty=2)

Figure3 <- plot_grid(p3a + ggtitle("a."), 
          p3b + ggtitle("b."),
          p3c + ggtitle("c."),
          nrow=1)

ggsave(plot = Figure3,
       filename = paste0(path_figures_project, "FIGURE_3.tiff"),
       bg="white", dpi = 600, width=12, height=5)

#--------------------
# Figure 4
# Description of the clusters in terms of the FPCA scores 
# based on the results of Tukey tests

Figure4 <- tab_tukey %>% 
  mutate(lab_y = case_when(
    FPC == "FPC score 1" ~ 11500,
    FPC == "FPC score 2" ~ 4500,
    FPC == "FPC score 3" ~ 2500
  )) %>% 
  ggplot(., aes(x = km_4, y = score)) +
  geom_text(aes(x=km_4, y=lab_y, label = tukey, color=tukey), size=5, check_overlap = T)+
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
        axis.title = element_text(size=15),
        axis.text =  element_text(size=15)) +
  labs(x = "Clusters of fragments", y = "FPC score value") ; Figure4


ggsave(plot = Figure4, 
       filename = paste0(path_figures_project, "FIGURE_4.tiff"),
       bg="white", dpi = 600, width=10, height=3.5)

#--------------------
# Figure 5
# Description of the clusters in terms of the number of first lesions
# based on the results of GLM models

Figure5 <- ggplot(data_fragment_first_lesions, 
                  aes(x=fct_rev(fct_reorder(km_4,nblesions1,.fun="mean")), 
                      y=nblesions1, colour = Groupes, fill=Groupes))+
  geom_boxplot(outlier.alpha = 0, alpha=0.5)+
  geom_jitter(width=0.25)+  
  stat_summary(fun=mean, colour="black", geom="point", shape=4, size=3)+
  theme_bw(base_size = 15)+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(hjust=1, vjust=1))+
  geom_text(data = data_fragment_first_lesions, aes(label = Groupes),
            y=150, size=5, check_overlap = TRUE) + ylim(c(0,150))+
  theme(plot.caption = element_text(hjust = 0.5), 
        axis.text.x = element_text(hjust = 0.5, color="black")) +
  ylab("Number of initial lesions")+
  xlab("Clusters")+
  scale_fill_manual(values = c("#4777EFFF","#7A0403FF","#30123BFF"))+
  scale_colour_manual(values = c("#4777EFFF","#7A0403FF","#30123BFF"))

Figure5
ggsave(plot = Figure5, 
       paste0(path_figures_project, "FIGURE_5.tiff"), 
       dpi = 600, width=7, height=5)

#--------------------
# Figure 6
# Plot the curves for each clusters

Figure6 <- plot_grid(
  plot_dyn_id(frags = X_scores1 %>% filter(km_4==1) %>% pull(ID), score = NULL, color = pal_cluster[1]) +
    geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
    geom_hline(yintercept = 0, color="darkgrey", lty=2) +
    ggtitle(paste0("a. Cluster 1 (N=", paste0(KM4$size[2]), ")")),
  
  plot_dyn_id(frags = X_scores1 %>% filter(km_4==2) %>% pull(ID), score = NULL, color = pal_cluster[2]) +
    geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
    geom_hline(yintercept = 0, color="darkgrey", lty=2) +
    ggtitle(paste0("b. Cluster 2 (N=", paste0(KM4$size[1]), ")")),
  
  plot_dyn_id(frags = X_scores1 %>% filter(km_4==3) %>% pull(ID), score = NULL, color = pal_cluster[3]) +
    geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
    geom_hline(yintercept = 0, color="darkgrey", lty=2) +
    ggtitle(paste0("c. Cluster 3 (N=", paste0(KM4$size[3]), ")")),
  
  plot_dyn_id(frags = X_scores1 %>% filter(km_4==4) %>% pull(ID), score = NULL, color = pal_cluster[4]) +
    geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
    geom_hline(yintercept = 0, color="darkgrey", lty=2) +
    ggtitle(paste0("d. Cluster 4 (N=", paste0(KM4$size[4]), ")")),
  
  nrow=1, align="hv", axis="tbrl", rel_widths = c(0.31, 0.29, 0.29, 0.29)
) ; Figure6


ggsave(plot = Figure6, 
       filename = paste0(path_project, "/04_ARTICLE/FIGURES_R1/FIGURE_6.tiff"),
       bg="white", dpi = 600, width=14, height=4)

# ----------------------------------
# Figure 7

Figure7 <- estimated_mean %>%
  ggplot(., aes(x = x_i, color=cluster)) +
  geom_hline(yintercept = c(0, 3512), color="black", lty=2) +
  # > mean 
  geom_path(data = Fit.mu1, 
            aes(x=d, y=mu), color="black", lwd=1) +
  geom_line(aes(y=f_i), linewidth=1.2) +
  geom_ribbon(aes(ymin=ic_down, ymax=ic_up, fill=cluster), color="transparent", alpha=0.3) +
  theme_bw(base_size = 15) +
  theme(legend.position = "bottom", 
        axis.text = element_text(color="black")) +
  labs(x="Days after rank 4", y="Mean diseased surface in mm²\n95% confidence interval) ") +
  guides(color = guide_legend(nrow=1)) +
  lims(y=c(-100, 3600)) +
  scale_color_manual(values=pal_cluster, name="Cluster") +
  scale_fill_manual(values=pal_cluster, name="Cluster"); Figure7


ggsave(plot = Figure7, 
       filename = paste0(path_figures_project, "FIGURE_7.tiff"), 
       bg="white", dpi = 600, width=5, height=5)

# ----------------------------------
# Correspondance between a-priori and a-posteriori clusters 
# Results from the correspondance analysis

# Function to rescale the % 
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

Figure8 <- plot_grid(p0, p1, nrow=2)

ggsave(plot = Figure8,
       filename = paste0(path_figures_project, "FIGURE_8.tiff"), 
       bg="white", dpi = 600, width=10, height=12)

# ----------------------------------
# Supplementary Figures X

# > plot other scores according to the 1rst one
p3b_2 <- plot_dyn_id(frags = frags_pca_score2, 
                     score="pca_score1", score_lab = "1") +
  geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
  geom_hline(yintercept = 0, color="darkgrey", lty=2)

p3c_2 <- plot_dyn_id(frags = frags_pca_score3, 
                     score="pca_score1", score_lab = "1") +
  geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
  geom_hline(yintercept = 0, color="darkgrey", lty=2)

# plot 3rd scores according to the 2nd one
p3c_3 <- plot_dyn_id(frags = frags_pca_score3, 
                     score="pca_score2", score_lab = "2") +
  geom_hline(yintercept = 3512, color="darkgrey", lty=2) +
  geom_hline(yintercept = 0, color="darkgrey", lty=2)

plot_grid(
  plot_grid(p3a + ggtitle("a."), ggplot()+theme_void(), ggplot()+theme_void(), nrow=1, align="hv"),
  plot_grid(p3b, p3b_2 + ggtitle("b."), ggplot()+theme_void(), nrow=1, align="hv"),
  plot_grid(p3c, p3c_2, p3c_3 + ggtitle("c."), nrow=1, align="hv"),
  ncol=1)

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


# -----

# Plot the curves for each clusters
plot_grid(
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
) 



# ---------------------------
# Represent the mean dynamic + the functional coefficient for each cluster
# estimated from a function-on-scalar regression with leaf as random effect

pa <- tab_coef %>%
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


plot_grid(pb+theme(legend.position = "none")+ggtitle("a."), 
          pa+ggtitle("b."), nrow=2, rel_heights = c(0.43, 0.57))


