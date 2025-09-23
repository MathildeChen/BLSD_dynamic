# Heterogeneous disease progress of Black Leaf Streak Disease of banana at within-leaf scale

### Summary
This repository contains scripts supporting a paper aiming to investigate the within-leaf progress of Black Leaf Streak Disease (BLSD) of banana plantain. 
Disease progress was measured over time in 192 leaf fragments of 3512 mm² taken from 9 banana leaves. 
Each fragment was monitored twice a week until fully covered by the disease. 
Functional principal component analysis was used to identify the main sources of variability in the disease progress over all individual fragments. 
Four groups of fragments with similar patterns of disease progress were clustered using k-means algorithm.
The study revealed distinct temporal patterns that were significantly different from each other and linked to the predefined leaf compartments through a correspondence analysis. 

### Analyses step and corresponding scripts
<img width="4058" height="3885" alt="Git_Repo_organisation" src="https://github.com/user-attachments/assets/3adc1384-0fea-472e-8bd4-bf325b6520cf" />


### Packages required 

All analyses were undertaken using R version 4.4.2 (http://www.r-project.org) with a two-sided p<0.05 considered statistically significant. 

- Exploratory and descriptives statistics were produced using the *corrplot* (version 0.95), *arsenal* (version 3.6.3), and *multcomp* (version 1.4-28) packages.  
- Functional analyses were performed using the *refund* package (version 0.1-37) 
- Clustering was performed using the *stats* package (version 4.4.2.). Clusters' visualizations were obtained using the *cluster*, *FactoMineR* (version 2.11), and *factoextra* (version 1.0.7) packages. 
- Generalized linear mixed models were fitted using the *lme4* package (version 1.1-37) and the Estimated Marginal Means were computed using the *emmeans* package (version 1.11.1).
- Correspondance analysis was performed using the *FactoMineR* package (version 2.11). 

## Authors 
Marine Seidel ${1,4}$ *, Mathilde Chen ${2,4}$, Jacques Avelino ${2,4}$, Fabienne Ribeyre ${2,4}$, Clara Landry ${1,4}$, Catherine Abadie ${3,4}$

${1}$ CIRAD, UMR PHIM, F-97130 Capesterre-Belle-Eau, Guadeloupe, France.

${2}$ CIRAD, UMR PHIM, F-34398 Montpellier, France.

${3}$ CIRAD, UMR PHIM, 30501 Turrialba, Costa Rica.

${4}$ PHIM, CIRAD, INRAE, Institut Agro, IRD, Université de Montpellier, Montpellier, France.

*Corresponding author contact: marine.seidel@cirad.fr

Authors ORCID:
- Marine Seidel    : 0009-0002-3786-4673
- Mathilde Chen    : 0000-0002-5982-2143
- Jacques Avelino  : 0000-0003-1983-9431
- Fabienne Ribeyre : 0000-0001-9721-1485
- Catherine Abadie : 0000-0002-5075-6338

## Key words 
*Pseudocercospora fijiensis*, dynamics, spatial pattern, functional data analysis, Dominican Republic, plantain



