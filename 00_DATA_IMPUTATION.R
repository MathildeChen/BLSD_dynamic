###################

# Black Leaf Streak Disease (BLSD) dynamic on banana leaves

# Authors: M. Seidel, M. Chen, J. Avelino, C. Ababdie CIRAD, 2025

###################

# Step 1: data preparation for analyses
# - imputation of missing observations at the start of the season
# - formatting 

# ----------------------------------
# Packages
# > General data management 
library(tidyverse) # 

# > Graphics
library(cowplot) # add several panels in a plot

# ----------------------------------
# Data loading

# > Indicate the path to project
path_project <- "D:/Mes Donnees/CIRAD/COLLABORATIONS/2024_MARINE/"

# > Load data
load(paste0(path_project, "00_DONNEES/data_fragment.RData"))

# ----------------------------------
# In a *given leaf* all fragments are observed 
# as soon as one lesion appeared on at least one fragment. 

# The absence of lesion is not indicated for the other fragments
# and need to be recomputed

# 1. Identify the unique fragments
all_fragments <- data_fragment %>%
  group_by(feuille) %>%
  distinct(feuille, compartiment, fragment) ; all_fragments

# Number of fragments by compartment and by leaf
all_fragments %>% 
  group_by(feuille, compartiment) %>% 
  count() %>% 
  spread(key = compartiment, value=n)

#   feuille   `1`   `2`   `3`   `4`   `5`   `6`   `7`   `8`
#1 B1F1        2     4     2     4     2     2     2     4
#2 B1F2        2     4     2     4     2     2     2     4
#3 B1F3        2     4     2     4     2     2     2     4
#4 B2F1        2     4     2     4     2     2     2     4
#5 B2F2        2     4     2     4     2    NA     2     3
#6 B2F3        2     4     1     4     2     2     2     4
#7 B3F1        2     4     2     4     2     2     2     4
#8 B3F2        2     4     2     4     2     2     2     4
#9 B3F3        2     4     1     4     1     2     2     4
# -> some missing fragments on leaves B2F2 (3), B2F3 (1), B3F3 (2), 
# total: 6 missing fragments

# 2. Count the total number of observations per fragments
data_fragment %>% 
  group_by(feuille, compartiment, fragment) %>% 
  count() %>% 
  pivot_wider(names_from = c(compartiment, fragment), values_from=n)

#  feuille `1_z1` `1_z2` `2_z3` `2_z4` `2_z5` `2_z6` `3_z7` `3_z8` `4_z10` `4_z11` `4_z12` `4_z9` `5_z13` `5_z14` `6_z19` `6_z20` `7_z21` `7_z22` `8_z15` `8_z16` `8_z17` `8_z18`
#1 B1F1        13     13     13     13     13     13     14     14      13      13      14     13      13      12      13      12      13      14      13      13      13      13
#2 B1F2        14     13     12     11     10      9     12     11      10      11      11      9      11       8      11      11      11      10      14      10      10      14
#3 B1F3        14     14     11     11     13     14     14     13      13      13      14     13      12      12      14      14      12      13      11      12      11      12
#4 B2F1        10     10     10     10      9     10     11     10      10      10      10     10       9       9       9       9      10      10       8       8       7       9
#5 B2F2        14     14     11     13     13     13     14     14      14      11      11     14      13      11      NA      NA      14      12      13      13      11      NA
#6 B2F3        11     12     10     11      8      9     12     NA       9      10      11     11      11       7       9       9      11       8       7       9      11      11
#7 B3F1        16     16     13     13     12     12     12     13      11      14      16     13      16      13      12      13      11      13      11      12      11      12
#8 B3F2        16     15     15     15     15     14     15     17      14      15      16     13      16      13      14      15      14      14      13      13      13      13
#9 B3F3        13     15     11     12     12     13     NA     16      13      13      13     14      NA      13      14      13      12      12      13      12      13      13
# -> normally the number of observations should be the same for the 
# fragments of a same leaf

# First and last observation for each fragment by leaf
data_fragment %>% 
  group_by(feuille, compartiment, fragment) %>% 
  summarize(first_obs_frag = first(date),
            last_obs_frag = last(date)) %>% 
  ggplot(.) +
  geom_linerange(aes(xmin=first_obs_frag, xmax=last_obs_frag, y = fragment)) +
  geom_point(aes(x = first_obs_frag, y = fragment)) +
  geom_point(aes(x = last_obs_frag, y = fragment), color="red") +
  facet_grid(.~feuille) +
  theme_bw() +
  labs(x = "Première obs (noir) - dernière obs (rouge)", y = "")

# 3. Identify the different types of missing data
data_fragment %>% 
  filter(feuille=="B1F1" & compartiment %in% c(5, 8) |
         feuille=="B1F2" & compartiment %in% c(5, 8) |
         feuille=="B2F3" & compartiment %in% c(5, 8)) %>% 
  ggplot(., aes(x = date, y = surface, color = fragment_unique)) +
  geom_path() +
  geom_point(aes(y = surface)) +
  theme_bw() +
  theme(legend.position = "none") +
  facet_grid(fragment ~ feuille, scales = "free_x") 

# -> some are missing in early season (we know they are 0 if at least 1 fragment of the leaf as >0 mm² infected)
# -> some are missing during the observation period (these observations are "real" missing observations)

# 4. Build a fake table with all the observation dates for each leaf and each fragment 
full_dates_fragment <- data_fragment %>% 
  #filter(feuille=="B2F3") %>%
  split(.$feuille) %>% 
  map_dfr(., ~{
    
    # Observations dates in each leaf 
    unique_dates_feuille <- .x %>% 
      distinct(bananier, feuille, date_rang4, date)  %>% 
      arrange(feuille, date)
    
    # Build a table with all dates for all fragments in each leaf 
    all_dates_all_frag <- NULL
    for(fragment_i in unique(.x$fragment))
    {
      
      # fragment and their respective compartment
      frag_comp <- .x %>% 
        filter(fragment == fragment_i)
      
      # data to add 
      to_add <- unique_dates_feuille
      to_add$compartiment    <- unique(frag_comp$compartiment)
      to_add$fragment        <- unique(frag_comp$fragment)
      to_add$fragment_unique <- unique(frag_comp$fragment_unique)
      
      # merge 
      all_dates_all_frag <- rbind(all_dates_all_frag, 
                                  to_add)
      
    }
    
    all_dates_all_frag
    
  }, .id = "feuille")

# -> check that the number of observations is equal 
#    among the existing fragment of a leaf
full_dates_fragment %>% 
  group_by(feuille, compartiment, fragment) %>% 
  count() %>% 
  pivot_wider(names_from = c(compartiment, fragment), values_from=n)

# -> check the first observation date 
# should be the same for all fragments of a same leaf
full_dates_fragment %>% 
  group_by(feuille, compartiment, fragment) %>% 
  summarize(first_obs_frag = first(date),
            last_obs_frag = last(date)) %>% 
  ggplot(.) +
  geom_linerange(aes(xmin=first_obs_frag, xmax=last_obs_frag, y = fragment)) +
  geom_point(aes(x = first_obs_frag, y = fragment)) +
  geom_point(aes(x = last_obs_frag, y = fragment), color="red") +
  facet_grid(.~feuille) +
  theme_bw() +
  labs(x = "Première obs (noir) - dernière obs (rouge)", y = "")

# -> difference of the nb of observations before and after imputation 
# among compartments
data_fragment %>% 
  group_by(feuille, compartiment) %>% 
  summarise(n_obs_init = n()) %>%
  left_join(full_dates_fragment %>% 
              group_by(feuille, compartiment) %>% 
              summarise(n_obs_new = n())) %>% 
  ungroup() %>% 
  mutate(dif_abs = n_obs_new - n_obs_init,
         dif_rel = ((n_obs_new - n_obs_init)/n_obs_init)*100) %>% 
  mutate(id_color = if_else(dif_rel>30, 1, 0)) %>% 
  ggplot(., aes(x = feuille, y = compartiment, fill = dif_rel, label = round(dif_rel,1))) +
  geom_tile() +
  geom_text(aes(color=as.factor(id_color)), size=3) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_fill_viridis_c(direction=-1) +
  scale_color_manual(values = c("black", "white")) 

# Compartments with more than 30% of new observations! 
# - B1F2 compartments 2, 4, 5, 7
# - B2F1 compartment 8
# - B2F3 compartments 5, 6
# - B3F1 compartments 7, 8
# - B3F3 compartment 8
# - B3F3 compartment 2, 7

# -> difference of the nb of observations before and after imputation 
# among fragments
data_fragment %>% 
  group_by(feuille, compartiment, fragment) %>% 
  summarise(n_obs_init = n()) %>%
  left_join(full_dates_fragment %>% 
              group_by(feuille, compartiment, fragment) %>% 
              summarise(n_obs_new = n())) %>% 
  ungroup() %>% 
  mutate(dif_abs = n_obs_new - n_obs_init,
         dif_rel = ((n_obs_new - n_obs_init)/n_obs_init)*100) %>% 
  mutate(id_color = if_else(dif_rel>30, 1, 0)) %>% 
  ggplot(., aes(x = fragment, y = reorder(feuille, dif_rel, mean), fill = dif_rel, label = round(dif_rel,1))) +
  geom_tile() +
  geom_text(aes(color=as.factor(id_color)), size=3) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_fill_viridis_c(direction=-1) +
  scale_color_manual(values = c("black", "white")) +
  facet_grid(.~compartiment, scales = "free_x", space = "free")

# -> Some fragments with more than 50% of missing data ! 

# 5. Imputation 
# Merge the full dates with recorded dates
data_fragment_imp <- full_dates_fragment %>% 
  ungroup() %>% 
  left_join(data_fragment %>% ungroup()) %>%
  arrange(bananier, feuille, compartiment, fragment, date) %>% 
  group_by(fragment_unique) %>%
  mutate(order = 1:n()) %>% 
  #ungroup() %>% 
  # > store initial surface measure 
  mutate(surface_init = surface) %>% 
  ungroup() %>% 
  # > replace the NA as 0 when it's in early observations
  mutate(surface = case_when(
    
    is.na(surface_init) == TRUE ~ 0,
    TRUE ~ surface_init
    
  )) %>% 
  # > let the missing observations 
  mutate(surface = case_when(
    
    order > 1 & is.na(surface_init) & lag(surface) == 0 ~ 0,
    order > 1 & is.na(surface_init) & lag(surface) > 0 ~ NA,
    TRUE ~ surface
    
  )) %>% 
  # > highlight lines where 0 were added and lines with missing data 
  mutate(added_zero = if_else(is.na(surface_init) & surface == 0, 1, NA),
         is_missing = if_else(is.na(surface_init) & is.na(surface), 1, NA)) %>% 
  # > recompute the jour_depuis_rang4 variable for the added_zero and missing data 
  mutate(jour_depuis_rang4_init = jour_depuis_rang4,
         jour_depuis_rang4 = as.numeric(as.character(date-date_rang4)))

# -> Quick check of imputation on few fragments
# - Fragment with only 1 missing observation at the beginning
data_fragment_imp %>% 
  filter(feuille=="B2F1", fragment %in% c("z2")) %>% 
  dplyr::select(date, order, jour_depuis_rang4_init, jour_depuis_rang4, surface_init, surface, added_zero, is_missing)

#    date       order jour_depuis_rang4_init jour_depuis_rang4 surface_init surface added_zero is_missing
#   <date>     <int>                  <dbl>             <dbl>        <dbl>   <dbl>      <dbl>      <dbl>
# 1 2012-10-22     1                     NA                 0          NA       0           1         NA
# 2 2012-10-29     2                      7                 7         392.    392.         NA         NA
# 3 2012-11-01     3                     10                10        1528.   1528.         NA         NA
# 4 2012-11-05     4                     14                14        2161.   2161.         NA         NA
# 5 2012-11-08     5                     17                17        2494.   2494.         NA         NA
# 6 2012-11-12     6                     21                21        3512    3512          NA         NA
# 7 2012-11-15     7                     24                24        3512    3512          NA         NA
# 8 2012-11-19     8                     28                28        3512    3512          NA         NA
# 9 2012-11-22     9                     31                31        3512    3512          NA         NA
#10 2012-11-26    10                     35                35        3512    3512          NA         NA
#11 2012-11-29    11                     38                38        3512    3512          NA         NA

# - Fragment with several missing observations at the beginning
data_fragment_imp %>% 
  filter(feuille=="B2F2", fragment %in% c("z12")) %>% 
  dplyr::select(date, order, jour_depuis_rang4_init, jour_depuis_rang4, surface_init, surface, added_zero, is_missing)

#   date       order jour_depuis_rang4_init jour_depuis_rang4 surface_init surface added_zero is_missing
#   <date>     <int>                  <dbl>             <dbl>        <dbl>   <dbl>      <dbl>      <dbl>
# 1 2012-10-29     1                     NA                 0        NA       0             1         NA
# 2 2012-11-01     2                     NA                 3        NA       0             1         NA
# 3 2012-11-05     3                     NA                 7        NA       0             1         NA
# 4 2012-11-08     4                     10                10         6.19    6.19         NA         NA
# 5 2012-11-12     5                     14                14        40.1    40.1          NA         NA
# 6 2012-11-15     6                     17                17       195.    195.           NA         NA
# 7 2012-11-19     7                     21                21       387.    387.           NA         NA
# 8 2012-11-22     8                     24                24       403.    403.           NA         NA
# 9 2012-11-26     9                     28                28       982.    982.           NA         NA
#10 2012-11-29    10                     31                31      1900.   1900.           NA         NA
#11 2012-12-06    11                     38                38      2226.   2226.           NA         NA
#12 2012-12-10    12                     42                42      2877.   2877.           NA         NA
#13 2012-12-13    13                     45                45      3512    3512            NA         NA
#14 2012-12-17    14                     49                49      3512    3512            NA         NA

# - Fragment with several missing observation at the start + during the season
data_fragment_imp %>% 
  filter(feuille=="B2F3", fragment %in% c("z16")) %>% 
  dplyr::select(date, order, jour_depuis_rang4_init, jour_depuis_rang4, surface_init, surface, added_zero, is_missing)

#   date       order jour_depuis_rang4_init jour_depuis_rang4 surface_init surface added_zero is_missing
#   <date>     <int>                  <dbl>             <dbl>        <dbl>   <dbl>      <dbl>      <dbl>
# 1 2012-11-05     1                     NA                 0         NA       0            1         NA
# 2 2012-11-08     2                     NA                 3         NA       0            1         NA
# 3 2012-11-12     3                      7                 7         14.7    14.7         NA         NA
# 4 2012-11-15     4                     10                10         21.1    21.1         NA         NA
# 5 2012-11-19     5                     NA                14         NA      NA           NA          1
# 6 2012-11-26     6                     21                21         85.1    85.1         NA         NA
# 7 2012-11-29     7                     24                24        282.    282.          NA         NA
# 8 2012-12-06     8                     31                31       1824.   1824.          NA         NA
# 9 2012-12-10     9                     35                35       2826.   2826.          NA         NA
#10 2012-12-13    10                     38                38       3203.   3203.          NA         NA
#11 2012-12-17    11                     42                42       3512    3512           NA         NA
#12 2012-12-20    12                     45                45       3512    3512           NA         NA

# - Fragment with 75% more data than initial data
data_fragment_imp %>% 
  filter(feuille=="B1F2", fragment %in% c("z14")) %>% 
  dplyr::select(date, order, jour_depuis_rang4_init, jour_depuis_rang4, surface_init, surface, added_zero, is_missing)

#   date       order jour_depuis_rang4_init jour_depuis_rang4 surface_init surface added_zero is_missing
#   <date>     <int>                  <dbl>             <dbl>        <dbl>   <dbl>      <dbl>      <dbl>
# 1 2012-10-15     1                     NA                 4          NA       0           1         NA
# 2 2012-10-18     2                     NA                 7          NA       0           1         NA
# 3 2012-10-22     3                     NA                11          NA       0           1         NA
# 4 2012-10-29     4                     NA                18          NA       0           1         NA
# 5 2012-11-01     5                     NA                21          NA       0           1         NA
# 6 2012-11-05     6                     NA                25          NA       0           1         NA
# 7 2012-11-08     7                     28                28         107.    107.         NA         NA
# 8 2012-11-12     8                     32                32         376.    376.         NA         NA
# 9 2012-11-15     9                     35                35         700.    700.         NA         NA
#10 2012-11-19    10                     39                39        1453.   1453.         NA         NA
#11 2012-11-22    11                     42                42        1814.   1814.         NA         NA
#12 2012-11-26    12                     46                46        2513.   2513.         NA         NA
#13 2012-11-29    13                     49                49        2952.   2952.         NA         NA
#14 2012-12-06    14                     56                56        3512    3512          NA         NA


# - On more fragments
data_fragment_imp %>% 
  filter(feuille=="B1F1" & compartiment %in% c(5, 8) |
           feuille=="B1F2" & compartiment %in% c(5, 8) |
           feuille=="B2F3" & compartiment %in% c(5, 8)) %>% 
  ggplot(., aes(x = date, y = surface_init)) +
  geom_path(aes(y=surface), color="red") +
  geom_path() +
  geom_point(aes(y = surface), color="red") +
  geom_point(aes(y = surface_init)) +
  facet_grid(fragment ~ feuille) +
  theme_bw()

# -> Check on all the fragments 
# - cross = 0 added 
# - red points = missing data
data_fragment_imp %>% 
  ggplot(.) +
  geom_point(aes(x = date, y = fragment, 
                 color=surface, shape=as.factor(added_zero))) +
  scale_color_viridis_c(na.value = "red", direction = -1) +
  scale_shape_manual(values = 4, na.value = 19) +
  facet_grid(compartiment~feuille, scales="free_y", space ="free") +
  theme_dark() +
  theme(legend.position = "none")

# -> distribution of the surface
ggplot() + 
  geom_density(data = data_fragment, aes(x=surface), fill="blue", alpha=0.2) +
  geom_density(data = data_fragment_imp, aes(x=surface), fill="red", alpha=0.2)


# Total number of obs
dim(data_fragment)     # Initially: 2332   
dim(data_fragment_imp) # After imputation: 2730

# No NAs in the fragments id after the merge
unique(data_fragment_imp$fragment_unique) # 192
unique(data_fragment_imp$fragment)        # 22

# 6. Get the number of lesions at the first date when nb of lesions > 0
# in a fragment 
nb_lesions1_fragment <- data_fragment_imp %>% 
  filter(surface > 0) %>% 
  group_by(fragment_unique) %>% 
  filter(date == first(date)) %>% 
  ungroup() %>%
  dplyr::select(bananier, feuille, fragment_unique, fragment, "date1"="date", "nblesions1"="nblesions") 


# 7. Format data for functional analysis using refund
# - ID obs 
tab_ID <- data.frame(fragment_unique = unique(data_fragment_imp$fragment_unique),
                     ID              = 1:length(unique(data_fragment_imp$fragment_unique)))
dim(tab_ID) # should be 192 unique ID

# > Scalar variables for each of the 192 fragments
X_init <- data_fragment_imp %>% 
  dplyr::select(bananier, "feuille_unique"="feuille", compartiment, fragment, fragment_unique) %>% 
  distinct(.) %>% 
  # > ID of fragment
  left_join(., tab_ID) %>%
  # > nb of lesions at the 1rst date when more than 0 lesion
  left_join(., nb_lesions1_fragment) %>% 
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
  arrange(bananier, feuille, compartiment, fragment_unique)

str(X_init)
#tibble [192 × 9] (S3: tbl_df/tbl/data.frame)
# $ bananier       : Factor w/ 3 levels "Tree 1","Tree 2",..: 1 1 1 1 1 1 1 1 1 1 ...
# $ feuille_unique : Factor w/ 9 levels "B1F1","B1F2",..: 1 1 1 1 1 1 1 1 1 1 ...
# $ compartiment   : Factor w/ 8 levels "1","2","3","4",..: 1 1 2 2 2 2 3 3 4 4 ...
# $ fragment       : num [1:192] 1 2 3 4 5 6 7 8 10 11 ...
# $ fragment_unique: Factor w/ 192 levels "P1H1z1","P1H1z10",..: 1 12 16 17 18 19 20 21 2 3 ...
# $ ID             : int [1:192] 1 2 3 4 5 6 7 8 9 10 ...
# $ feuille        : Factor w/ 3 levels "Leaf 1","Leaf 2",..: 1 1 1 1 1 1 1 1 1 1 ...
# $ date1          : Date[1:192], format: "2012-10-08" "2012-10-08" ...
# $ nblesions1     : int [1:192] 21 60 20 9 46 36 42 79 20 12 ...

# > Functional sparse response in long format
Y_t <- data_fragment_imp %>% 
  left_join(., tab_ID) %>% 
  dplyr::select(".id"= ID, 
                #".index"=jour_depuis_app_feuille, 
                ".index"=jour_depuis_rang4, 
                ".value"=surface) 

dim(Y_t) # 2730 observations

summary(Y_t$.value)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    0.0   112.6  1245.8  1575.4  3141.8  3512.0      17 
# --> 17 NAs

# > check 
Y_t %>% 
  left_join(X_init, by=c('.id'='ID')) %>% 
  arrange(bananier, feuille, fragment) %>% 
  ggplot(data = .) + 
  geom_hline(yintercept = c(0, 3512/2, 3512), color = "lightgrey", linetype = 2) +
  geom_point(aes(x = .index, y = .value)) + 
  geom_path(aes(x = .index, y = .value, color=as.factor(.id))) + 
  theme_cowplot() +
  theme(panel.grid.major.x = element_line(color = "lightgrey", linetype = 2),
        legend.position = "none") +
  facet_wrap(. ~ feuille_unique) + 
  labs(x = "Days after rank 4", y = "Infected surface (mm²)") +
  ggtitle("Imputed values") +
  scale_color_viridis_d(na.value="red")
# 

# > initial data
data_fragment %>% 
  ggplot(data = .) + 
  geom_hline(yintercept = c(0, 3512/2, 3512), color = "lightgrey", linetype = 2) +
  geom_point(aes(x = jour_depuis_rang4, y = surface)) + 
  geom_path(aes(x = jour_depuis_rang4, y = surface, color=as.factor(fragment_unique))) + 
  theme_cowplot() +
  theme(panel.grid.major.x = element_line(color = "lightgrey", linetype = 2),
        legend.position = "none") +
  facet_wrap(. ~ feuille) + 
  labs(x = "Days after rank 4", y = "Infected surface (mm²)") +
  ggtitle("Initial values") +
  scale_color_viridis_d(na.value="red")
# 

# > check 
plot_grid(
  # > initial data
  data_fragment %>% 
    filter(feuille=="B1F1", compartiment %in% c(5,6)) %>% 
    ggplot(data = .) + 
    geom_hline(yintercept = c(0, 3512/2, 3512), color = "lightgrey", linetype = 2) +
    geom_point(aes(x = jour_depuis_rang4, y = surface)) + 
    geom_path(aes(x = jour_depuis_rang4, y = surface, color=as.factor(fragment_unique))) + 
    theme_cowplot() +
    theme(panel.grid.major.x = element_line(color = "lightgrey", linetype = 2),
          legend.position = "none") +
    facet_wrap(fragment ~ ., nrow=1) + 
    labs(x = "Days after rank 4", y = "Infected surface (mm²)") +
    ggtitle("Initial values"),
  # Imputed
  Y_t %>% 
    left_join(X_init, by=c('.id'='ID')) %>% 
    arrange(bananier, feuille, fragment) %>% 
    filter(feuille_unique=="B1F1", compartiment%in% c(5,6)) %>% 
    ggplot(data = .) + 
    geom_hline(yintercept = c(0, 3512/2, 3512), color = "lightgrey", linetype = 2) +
    geom_point(aes(x = .index, y = .value)) + 
    geom_path(aes(x = .index, y = .value, color=as.factor(fragment_unique))) + 
    theme_cowplot() +
    theme(panel.grid.major.x = element_line(color = "lightgrey", linetype = 2),
          legend.position = "none") +
    facet_wrap(fragment ~ ., nrow=1) + 
    labs(x = "Days after rank 4", y = "Infected surface (mm²)") +
    ggtitle("Imputed values"),
  # 
  ncol=1
)

# > check 
plot_grid(
  # > initial data
  data_fragment %>% 
    filter(feuille=="B1F1", compartiment %in% c(5,6)) %>% 
    ggplot(data = .) + 
    geom_hline(yintercept = c(0, 3512/2, 3512), color = "lightgrey", linetype = 2) +
    geom_point(aes(x = jour_depuis_rang4, y = surface)) + 
    geom_path(aes(x = jour_depuis_rang4, y = surface, color=as.factor(fragment_unique))) + 
    theme_cowplot() +
    theme(panel.grid.major.x = element_line(color = "lightgrey", linetype = 2),
          legend.position = "none") +
    facet_wrap(fragment ~ ., nrow=1) + 
    labs(x = "Days after rank 4", y = "Infected surface (mm²)") +
    ggtitle("Initial values"),
  # Imputed
  Y_t %>% 
    left_join(X_init, by=c('.id'='ID')) %>% 
    arrange(bananier, feuille, fragment) %>% 
    filter(feuille_unique=="B1F1", compartiment %in% c(5,6)) %>% 
    filter(is.na(.value)==F) %>%
    ggplot(data = .) + 
    geom_hline(yintercept = c(0, 3512/2, 3512), color = "lightgrey", linetype = 2) +
    geom_point(aes(x = .index, y = .value)) + 
    geom_path(aes(x = .index, y = .value, color=as.factor(fragment_unique))) + 
    theme_cowplot() +
    theme(panel.grid.major.x = element_line(color = "lightgrey", linetype = 2),
          legend.position = "none") +
    facet_wrap(fragment ~ ., nrow=1) + 
    labs(x = "Days after rank 4", y = "Infected surface (mm²)") +
    ggtitle("Imputed values without the NAs"),
  # 
  ncol=1
)


# > Functional sparse response in long format for analyses using the refund package
ydata <- Y_t
unique(ydata$.id)  
dim(ydata)         # 2730 observations

# > Temporal index on which disease is measured
#   in days after the rank 4 of the leaf
yindex <- sort(unique(ydata$.index)) ; yindex
#  -4 -3  0  3  4  7 10 11 14 17 18 21 24 25 28 31 32 35 38 39 42 45 46 49 52 56 59 60 63

# remove the 17 NAs
ydata_no_NA <- Y_t %>% filter(is.na(.value)==F)
summary(ydata_no_NA$.value)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.0   112.6  1245.8  1575.4  3141.8  3512.0 

unique(ydata_no_NA$.id)  # still 192 fragments
dim(ydata_no_NA)         # 2713 observations

testthat::expect_equal(sort(unique(ydata_no_NA$.index)), yindex)
#  ok 

# > Transformation express the surface into proportion
Y_t_trans <- data_fragment_imp %>% 
  left_join(., tab_ID) %>%
  mutate(surface_trans=surface / 3512) %>% 
  dplyr::select(".id"= ID, 
                #".index"=jour_depuis_app_feuille, 
                ".index"=jour_depuis_rang4, 
                ".value"=surface_trans) 

ydata_trans <- Y_t_trans
unique(ydata_trans$.id)  # still 192 fragments
dim(ydata_trans)         # 2730 observations

testthat::expect_equal(sort(unique(ydata_trans$.index)), yindex)
#  ok 

# 7. Save the data 
data_fragment_for_analyses <- list(
  
  # Functional response
  "ydata" = ydata,
  "ydata_no_NA" = ydata_no_NA,
  "ydata_trans" = ydata_trans,
  
  # Temporal index on which disease is measured in days after the rank 4 of the leaf
  "yindex" = yindex,
  
  # Scalar variables
  "X_init" = X_init
  
  
)

save(data_fragment_for_analyses, file = paste0(path_project, "/00_DONNEES/data_fragment_for_analyses.rda"))




