## load libraries

library(tidyverse)
library(ggbeeswarm)
library(lme4)
library(lmerTest)

## load data

# formatted gene abundance data for FMT recipients, post-intervention timepoints
load("abundances_fmt_postBL.RData")

# engraftment-dependent HGT data from 4_hgt_engrafted.R
load("engraftment_dependent_hgt.RData")

# high quality MAGs, high quality genes, high quality gene clusters
load("HQ-gene-mapping-filtered.RData")

# metadata
load("metadata_clean.RData")

# colour palette
palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "slategrey", "#CC79A7", "darkslateblue", "#D55E00")

## running

# filter HQ_cluster data for engrafted male/female HGT clusters
engrafted_female_htg_clusters <- engrafted_hgt_female %>%
  left_join(HQ_clusters)

engrafted_male_htg_clusters <- engrafted_hgt_male %>%
  left_join(HQ_clusters)

# join fmt kma abundances with engraftment-dependent HGT clusters, by centroid
hgt_cluster_abundances_female <- engrafted_female_htg_clusters %>%
  # deselect week 6 timepoint to get all timepoints for the HGT clusters
  select(Cluster_ID, centroid, Recipient_ID, COG_Functional_cat.) %>% 
  left_join(kma_abundances_fmt, by = c("centroid" = "Gene.ID", "Recipient_ID" = "Participant_ID"))

hgt_cluster_abundances_female_distinct <- hgt_cluster_abundances_female %>%
  select(-c(Sample_ID, centroid)) %>%
  filter(!is.na(CPM)) %>% # remove rows with no abundance data (TF49:Cluster 2927876, TF29:Cluster 342406)
  distinct() # remove duplicate rows - these are from multiple HGT events with same cluster (different genes, eg. TF36:Cluster 1580304)

hgt_cluster_abundances_male <- engrafted_male_htg_clusters %>%
  # deselect week 6 timepoint to get all timepoints for the HGT clusters
  select(Cluster_ID, centroid, Recipient_ID, COG_Functional_cat.) %>% 
  left_join(kma_abundances_fmt, by = c("centroid" = "Gene.ID", "Recipient_ID" = "Participant_ID"))

hgt_cluster_abundances_male_distinct <- hgt_cluster_abundances_male %>%
  select(-c(Sample_ID, centroid)) %>%
  filter(!is.na(CPM)) %>% # 0
  distinct() # remove duplicate rows

# find number of distinct engraftment-dependent HGT events at week 6
edhgt_female_wk6 <- hgt_cluster_abundances_female_distinct %>% filter(Timepoint == "Week 6") # 139 week 6 events
edhgt_male_wk6 <- hgt_cluster_abundances_male_distinct %>% filter(Timepoint == "Week 6") # 289 week 6 events

# find number of distinct week 6 engraftment-dependent HTGCs retained at week 12 and week 26
edhgt_female_wk12 <- hgt_cluster_abundances_female_distinct %>% filter(Timepoint == "Week 12") # events at week 12
edhgt_male_wk12 <- hgt_cluster_abundances_male_distinct %>% filter(Timepoint == "Week 12") # events at week 12

length(which(edhgt_female_wk12$Cluster_ID %in% edhgt_female_wk6$Cluster_ID == TRUE)) # 134 retained at week 12 in females
length(which(edhgt_male_wk12$Cluster_ID %in% edhgt_male_wk6$Cluster_ID == TRUE)) # 248 retained at week 12 in males

edhgt_female_wk26 <- hgt_cluster_abundances_female_distinct %>% filter(Timepoint == "Week 26") # events at week 26
edhgt_male_wk26 <- hgt_cluster_abundances_male_distinct %>% filter(Timepoint == "Week 26") # events at week 26

length(which(edhgt_female_wk26$Cluster_ID %in% edhgt_female_wk6$Cluster_ID == TRUE)) # 135 retained at week 26 in females
length(which(edhgt_male_wk26$Cluster_ID %in% edhgt_male_wk6$Cluster_ID == TRUE)) # 232 retained at week 26 in males


# bind male and female data
all_cluster_abundance_retention_data <- hgt_cluster_abundances_female_distinct %>%
  bind_rows(hgt_cluster_abundances_male_distinct)

# plot box plot
fig6_plot <- all_cluster_abundance_retention_data %>%
  ggplot(aes(x = factor(Timepoint, level = c("Week 6", "Week 12", "Week 26")), y = CPM, fill = Sex)) +
  geom_quasirandom(shape = 21, dodge.width = 0.75, size = 1.7) +
  geom_boxplot(alpha = 0.6, outlier.colour = NA) +
  facet_wrap(~Sex, scales = "free") +
  scale_fill_manual(values = palette[c(3, 8)]) +
  xlab(NULL) + ylab("HTGC relative abundance (CPM)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4),
        legend.position = "none")

# prepare data for linear mixed effect model
# compare back to week 6
all_cluster_abundance_retention_data$Timepoint <- factor(all_cluster_abundance_retention_data$Timepoint, 
                                                         levels = c("Week 6", "Week 12", "Week 26")) # reorder Group levels

# fit full model
full_lmer_retention <- lmer(CPM ~ Sex * Timepoint+(1|Recipient_ID),data = all_cluster_abundance_retention_data)
summary(full_lmer_retention)

# reduced model (no interactions between fixed effects)
reduced_lmer_retention <- lmer(CPM ~ Sex + Timepoint+(1|Recipient_ID),data =all_cluster_abundance_retention_data)
summary(reduced_lmer_retention)

# run anova to compare full and reduced models
anova(full_lmer_retention, reduced_lmer_retention) # p = 0.6794 
# reduced model fits better - non-significant p value and lower AIC value

# get confidence interval of reduced model fixed effect coefficients
confint(reduced_lmer_retention)