
# Black Leaf Streak Disease (BLSD) dynamic on banana leaves

# Authors: M. Seidel, M. Chen, J. Avelino, C. Ababdie CIRAD, 2025

# ----------------------------------
# Packages

# > General data management 
library(tidyverse) 

# > Mixed models 
library(emmeans)
library(ggpubr)
library(lme4)

# ----------------------------------
# Path to project
path_project <- "D:/Mes Donnees/CIRAD/COLLABORATIONS/2024_MARINE/"

# ----------------------------------
# Run the analyses to generate input tables and objects for the figures
source(paste0(path_project, "01_SCRIPT/02_STATISTICAL_ANALYSES.R"))

# Examine the number of first lesions VS data-driven clusters using GLM and GLMM
# fit several models, check residuals, overdispersion and choose the one with the lowest AIC

# > family = Poisson
mod_clust1 <- glm (nblesions1  ~ km_4, family = "poisson", 
                   data = data_fragment_first_lesions)
mod_clust2 <- glmer(nblesions1 ~ km_4 + (1|bananier) + (1|feuille), family = "poisson", 
                    data = data_fragment_first_lesions)
mod_clust3 <- glmer(nblesions1 ~ km_4 + (1|bananier/feuille), family = "poisson", 
                    data = data_fragment_first_lesions) #isSingular

# > family = negative binomial 
mod_clust4 <- glm.nb(nblesions1   ~ km_4, 
                     data = data_fragment_first_lesions)
mod_clust5 <- glmer.nb(nblesions1 ~ km_4 + (1|bananier) + (1|feuille), 
                       data = data_fragment_first_lesions)
mod_clust6 <- glmer.nb(nblesions1 ~ km_4 + (1|bananier/feuille), 
                       data = data_fragment_first_lesions)

# Check AIC
AIC_total <- AIC(mod_clust1, mod_clust2, mod_clust3, mod_clust4, mod_clust5, mod_clust6)
AIC_total

AIC_1 <- AIC(mod_clust1)
AIC_2 <- AIC(mod_clust2)
AIC_3 <- AIC(mod_clust3)
AIC_4 <- AIC(mod_clust4)
AIC_5 <- AIC(mod_clust5)
AIC_6 <- AIC(mod_clust6)


# Check qqplot
qqplot_mod_clust1 <- ggqqplot(residuals(mod_clust1))
qqplot_mod_clust1

qqplot_mod_clust2 <- ggqqplot(residuals(mod_clust2))
qqplot_mod_clust2

qqplot_mod_clust3 <- ggqqplot(residuals(mod_clust3))
qqplot_mod_clust3

qqplot_mod_clust4 <- ggqqplot(residuals(mod_clust4))
qqplot_mod_clust4

qqplot_mod_clust5 <- ggqqplot(residuals(mod_clust5))
qqplot_mod_clust5

qqplot_mod_clust6 <- ggqqplot(residuals(mod_clust6))
qqplot_mod_clust6


# Over dispersion
# Mod_clust1
chi2_1 <- sum(residuals(mod_clust1, type = "pearson")^2)
pvalue_1 <- 1 - pchisq(chi2_1, df = df.residual(mod_clust1)) 
coeff_1 <- chi2_1 / df.residual(mod_clust1)

chi2_1
pvalue_1 # p < 0.05 = over dispersion of the residuals of the model
coeff_1 # coefficient of dispersion >> 1

# Mod_clust2
chi2_2 <- sum(residuals(mod_clust2, type = "pearson")^2)
pvalue_2 <- 1 - pchisq(chi2_2, df = df.residual(mod_clust2)) 
coeff_2 <- chi2_2 / df.residual(mod_clust2)

chi2_2
pvalue_2 # p < 0.05 = over dispersion of the residuals of the model
coeff_2 # coefficient of dispersion >> 1

# Mod_clust3
chi2_3 <- sum(residuals(mod_clust3, type = "pearson")^2)
pvalue_3 <- 1 - pchisq(chi2_3, df = df.residual(mod_clust3)) 
coeff_3 <- chi2_3 / df.residual(mod_clust3)

chi2_3
pvalue_3 # p < 0.05 = over dispersion of the residuals of the model
coeff_3 # coefficient of dispersion >> 1

# All models are over dispersed

# Mod_clust4
chi2_4 <- sum(residuals(mod_clust4, type = "pearson")^2)
pvalue_4 <- 1 - pchisq(chi2_4, df = df.residual(mod_clust4)) 
coeff_4 <- chi2_4 / df.residual(mod_clust4)

chi2_4
pvalue_4 # p >0.05 = no over dispersion of the residuals of the model
coeff_4 # coefficient of dispersion ~ 1

# Mod_clust5
chi2_5 <- sum(residuals(mod_clust5, type = "pearson")^2)
pvalue_5 <- 1 - pchisq(chi2_5, df = df.residual(mod_clust5)) 
coeff_5 <- chi2_5 / df.residual(mod_clust5)

chi2_5
pvalue_5 # p >0.05 = no over dispersion of the residuals of the model
coeff_5 # coefficient of dispersion ~ 1

# Mod_clust6
chi2_6 <- sum(residuals(mod_clust6, type = "pearson")^2)
pvalue_6 <- 1 - pchisq(chi2_6, df = df.residual(mod_clust6)) 
coeff_6 <- chi2_6 / df.residual(mod_clust6)

chi2_6
pvalue_6 # p >0.05 = no over dispersion of the residuals of the model
coeff_6 # coefficient of dispersion ~ 1

# MOD CLUST 6 : BETTER AIC + BETTER PVALUE AND COEFFICIENT OF DISPERSION
# > Choice of the model 6
