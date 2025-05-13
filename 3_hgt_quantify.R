## load libraries

library(tidyverse)
library(eulerr)
library(ggbeeswarm)

## load data

# formatted HQ gene data from 2_hgt_analysis_genes.R
load("HQ_genes_formatted.RData")

# HGT events data from 2_hgt_analysis_genes.R
load("group_hgt_events.RData")

# high quality MAGs, high quality genes, high quality gene clusters
load("HQ-gene-mapping-filtered.RData")

# metadata
load("metadata_clean.RData")

# colour palette
palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "slategrey", "#CC79A7", "darkslateblue", "#D55E00")

## running

# find relative proportions of genes involved in hgt across each fmt and placebo recipient (week 6)

# female fmt
HQ_genes_fmt_female <- HQ_genes_female_clean %>% # load all HQ female genes
  filter(Group == "FMT" & Timepoint == "Week 6") # filter for FMT week 6 
nrow(HQ_genes_fmt_female) # 616236

HQ_genes_fmt_female_individual <- HQ_genes_fmt_female %>% # load HQ female genes for FMT week 6
  group_by(Participant_ID) %>% # group by participant ID
  summarise(HQ_genes_per_indiv = n(), # get number of HQ genes for each individual
            HQ_clusters_per_indiv = n_distinct(Cluster_ID)) # get number of distinct gene clusters for each individual

hgt_genes_fmt_female_individual <- female_fmt_hgts_clean %>% # load HGT events
  mutate(Sample_ID = str_split(recipient_gene, "_", simplify = TRUE)[,1]) %>% # extract sample ID from gene ID
  mutate(Sample_ID = str_replace(Sample_ID, "wk6", "6wk")) %>% # format of sample ID from genes is reversed from meta data
  mutate(Sample_ID = str_replace(Sample_ID, "wk12", "12wk")) %>%
  mutate(Sample_ID = str_replace(Sample_ID, "wk26", "26wk")) %>%
  mutate(Sample_ID = str_replace(Sample_ID, "DF17D2a", "DF17aD2")) %>% # others have been swapped around also
  mutate(Sample_ID = str_replace(Sample_ID, "DF17D2b", "DF17bD2")) %>%
  mutate(Sample_ID = str_replace(Sample_ID, "DF16D4a", "DF16aD4")) %>%
  mutate(Sample_ID = str_replace(Sample_ID, "DF16D4b", "DF16bD4")) %>%
  left_join(reduced_meta, by = "Sample_ID") %>% # join with meta data
  group_by(Participant_ID) %>% # group by participant ID
  summarise(hgt_genes_per_indiv = n_distinct(recipient_gene), # get number of genes involved in HGT for each individual
            hgt_clusters_per_indiv = n_distinct(Cluster_ID)) %>% # get number of distinct gene clusters involved in HGT for each individual
  left_join(HQ_genes_fmt_female_individual, by = "Participant_ID") %>% # join with the "total" data from HQ genes above
  mutate(pct_hgt_genes = hgt_genes_per_indiv/HQ_genes_per_indiv*100) %>% # calculate percentage of genes at week 6 involved in HGT
  mutate(pct_hgt_clusters = hgt_clusters_per_indiv/HQ_clusters_per_indiv*100) %>% # calculate percentage of gene clusters at week 6 involved in HGT
  left_join(reduced_meta) %>% # rejoin with meta data
  filter(Sex == "Female" & Group == "FMT" & Timepoint == "Week 6") %>% # filter for relevant data
  select(-Sample_ID) # deselect irrelevant columns

# female placebo
HQ_genes_placebo_female <- HQ_genes_female_clean %>% # load all HQ female genes
  filter(Group == "Placebo" & Timepoint == "Week 6") # filter for Placebo week 6 
nrow(HQ_genes_placebo_female) # 741090

HQ_genes_placebo_female_individual <- HQ_genes_placebo_female %>% # load HQ female genes for placebo week 6
  group_by(Participant_ID) %>% # group by participant ID
  summarise(HQ_genes_per_indiv = n(), # get number of HQ genes for each individual
            HQ_clusters_per_indiv = n_distinct(Cluster_ID)) # get number of distinct gene clusters for each individual

hgt_genes_placebo_female_individual <- female_placebo_hgts_clean %>% # load HGT events
  mutate(Sample_ID = str_split(recipient_gene, "_", simplify = TRUE)[,1]) %>% # extract sample ID from gene ID
  mutate(Sample_ID = str_replace(Sample_ID, "wk6", "6wk")) %>% # format of sample ID from genes is reversed from meta data
  mutate(Sample_ID = str_replace(Sample_ID, "wk12", "12wk")) %>%
  mutate(Sample_ID = str_replace(Sample_ID, "wk26", "26wk")) %>%
  mutate(Sample_ID = str_replace(Sample_ID, "DF17D2a", "DF17aD2")) %>% # others have been swapped around also
  mutate(Sample_ID = str_replace(Sample_ID, "DF17D2b", "DF17bD2")) %>%
  mutate(Sample_ID = str_replace(Sample_ID, "DF16D4a", "DF16aD4")) %>%
  mutate(Sample_ID = str_replace(Sample_ID, "DF16D4b", "DF16bD4")) %>%
  left_join(reduced_meta, by = "Sample_ID") %>% # join with meta data
  group_by(Participant_ID) %>% # group by participant ID
  summarise(hgt_genes_per_indiv = n_distinct(recipient_gene), # get number of genes involved in HGT for each individual
            hgt_clusters_per_indiv = n_distinct(Cluster_ID)) %>% # get number of distinct gene clusters involved in HGT for each individual
  left_join(HQ_genes_placebo_female_individual, by = "Participant_ID") %>% # join with the "total" data from HQ genes above
  mutate(pct_hgt_genes = hgt_genes_per_indiv/HQ_genes_per_indiv*100) %>% # calculate percentage of genes at week 6 involved in HGT
  mutate(pct_hgt_clusters = hgt_clusters_per_indiv/HQ_clusters_per_indiv*100) %>% # calculate percentage of gene clusters at week 6 involved in HGT
  left_join(reduced_meta) %>% # rejoin with meta data
  filter(Sex == "Female" & Group == "Placebo" & Timepoint == "Week 6") %>% # filter for relevant data
  select(-Sample_ID) # deselect irrelevant columns

# male fmt
HQ_genes_fmt_male <- HQ_genes_male_clean %>% # load all HQ male genes
  filter(Group == "FMT" & Timepoint == "Week 6") # filter for FMT week 6 
nrow(HQ_genes_fmt_male) # 251198

HQ_genes_fmt_male_individual <- HQ_genes_fmt_male %>% # load HQ male genes for FMT week 6
  group_by(Participant_ID) %>% # group by participant ID
  summarise(HQ_genes_per_indiv = n(), # get number of HQ genes for each individual
            HQ_clusters_per_indiv = n_distinct(Cluster_ID)) # get number of distinct gene clusters for each individual

hgt_genes_fmt_male_individual <- male_fmt_hgts_clean %>% # load HGT events
  mutate(Sample_ID = str_split(recipient_gene, "_", simplify = TRUE)[,1]) %>% # extract sample ID from gene ID
  mutate(Sample_ID = str_replace(Sample_ID, "wk6", "6wk")) %>% # format of sample ID from genes is reversed from meta data
  mutate(Sample_ID = str_replace(Sample_ID, "wk12", "12wk")) %>%
  mutate(Sample_ID = str_replace(Sample_ID, "wk26", "26wk")) %>%
  mutate(Sample_ID = str_replace(Sample_ID, "DF17D2a", "DF17aD2")) %>% # others have been swapped around also
  mutate(Sample_ID = str_replace(Sample_ID, "DF17D2b", "DF17bD2")) %>%
  mutate(Sample_ID = str_replace(Sample_ID, "DF16D4a", "DF16aD4")) %>%
  mutate(Sample_ID = str_replace(Sample_ID, "DF16D4b", "DF16bD4")) %>%
  left_join(reduced_meta, by = "Sample_ID") %>% # join with meta data
  group_by(Participant_ID) %>% # group by participant ID
  summarise(hgt_genes_per_indiv = n_distinct(recipient_gene), # get number of genes involved in HGT for each individual
            hgt_clusters_per_indiv = n_distinct(Cluster_ID)) %>% # get number of distinct gene clusters involved in HGT for each individual
  left_join(HQ_genes_fmt_male_individual, by = "Participant_ID") %>% # join with the "total" data from HQ genes above
  mutate(pct_hgt_genes = hgt_genes_per_indiv/HQ_genes_per_indiv*100) %>% # calculate percentage of genes at week 6 involved in HGT
  mutate(pct_hgt_clusters = hgt_clusters_per_indiv/HQ_clusters_per_indiv*100) %>% # calculate percentage of gene clusters at week 6 involved in HGT
  left_join(reduced_meta) %>% # rejoin with meta data
  filter(Sex == "Male" & Group == "FMT" & Timepoint == "Week 6") %>% # filter for relevant data
  select(-Sample_ID) # deselect irrelevant columns

# male placebo
HQ_genes_placebo_male <- HQ_genes_male_clean %>% # load all HQ male genes
  filter(Group == "Placebo" & Timepoint == "Week 6") # filter for placebo week 6 
nrow(HQ_genes_placebo_male) # 353485

HQ_genes_placebo_male_individual <- HQ_genes_placebo_male %>% # load HQ male genes for placebo week 6
  group_by(Participant_ID) %>% # group by participant ID
  summarise(HQ_genes_per_indiv = n(), # get number of HQ genes for each individual
            HQ_clusters_per_indiv = n_distinct(Cluster_ID)) # get number of distinct gene clusters for each individual

hgt_genes_placebo_male_individual <- male_placebo_hgts_clean %>% # load HGT events
  mutate(Sample_ID = str_split(recipient_gene, "_", simplify = TRUE)[,1]) %>% # extract sample ID from gene ID
  mutate(Sample_ID = str_replace(Sample_ID, "wk6", "6wk")) %>% # format of sample ID from genes is reversed from meta data
  mutate(Sample_ID = str_replace(Sample_ID, "wk12", "12wk")) %>%
  mutate(Sample_ID = str_replace(Sample_ID, "wk26", "26wk")) %>%
  mutate(Sample_ID = str_replace(Sample_ID, "DF17D2a", "DF17aD2")) %>% # others have been swapped around also
  mutate(Sample_ID = str_replace(Sample_ID, "DF17D2b", "DF17bD2")) %>%
  mutate(Sample_ID = str_replace(Sample_ID, "DF16D4a", "DF16aD4")) %>%
  mutate(Sample_ID = str_replace(Sample_ID, "DF16D4b", "DF16bD4")) %>%
  left_join(reduced_meta, by = "Sample_ID") %>% # join with meta data
  group_by(Participant_ID) %>% # group by participant ID
  summarise(hgt_genes_per_indiv = n_distinct(recipient_gene), # get number of genes involved in HGT for each individual
            hgt_clusters_per_indiv = n_distinct(Cluster_ID)) %>% # get number of distinct gene clusters involved in HGT for each individual
  left_join(HQ_genes_placebo_male_individual, by = "Participant_ID") %>% # join with the "total" data from HQ genes above
  mutate(pct_hgt_genes = (hgt_genes_per_indiv/HQ_genes_per_indiv)*100) %>% # calculate percentage of genes at week 6 involved in HGT
  mutate(pct_hgt_clusters = (hgt_clusters_per_indiv/HQ_clusters_per_indiv)*100) %>% # calculate percentage of gene clusters at week 6 involved in HGT
  left_join(reduced_meta) %>% # rejoin with meta data
  filter(Sex == "Male" & Group == "Placebo" & Timepoint == "Week 6") %>% # filter for relevant data
  select(-Sample_ID) # deselect irrelevant columns

# bind rows
all_hgt_genes_individual <- bind_rows(hgt_genes_fmt_female_individual, hgt_genes_placebo_female_individual,
                                      hgt_genes_fmt_male_individual, hgt_genes_placebo_male_individual)
nrow(all_hgt_genes_individual) # 83

# plot percentage of gene clusters involved in HGT
fig3a_plot <- 
  all_hgt_genes_individual %>% # load combined data
  ggplot(aes(x = Sex, y = pct_hgt_clusters, fill = Group)) +
  geom_quasirandom(shape = 21, dodge.width = 0.75, size = 1.7) +
  geom_boxplot(alpha = 0.6, outlier.colour = NA) +
  scale_fill_manual(values = c("#b2182b", "#2166ac")) + 
  scale_colour_manual(values = c("#b2182b", "#2166ac"), name = NULL) +
  scale_y_continuous(limits = c(0, 12)) +
  xlab(NULL) + ylab("Percentage of gene clusters in HGT") + 
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4)) 

# run statistics to find impact of FMT on percentage of HTGCs (increase/ decrease?) from male/ female placebo vs FMT distributions
# check normality
all_hgt_genes_individual %>% ggplot(aes(x= pct_hgt_clusters)) + geom_histogram() # not normal
pull(all_hgt_genes_individual, pct_hgt_clusters) %>% shapiro.test() # p = 1.928e-07 (significant) == not normal
# run wilcoxon test
# female
female_box_data <- all_hgt_genes_individual %>% 
  filter(Sex == "Female") %>% # toggle between sex
  spread(Group, pct_hgt_clusters)
wilcox.test(female_box_data$FMT, female_box_data$Placebo)  # p = 0.8399
# male
male_box_data <- all_hgt_genes_individual %>% 
  filter(Sex == "Male") %>% # toggle between sex
  spread(Group, pct_hgt_clusters)
wilcox.test(male_box_data$FMT, male_box_data$Placebo)  # p = 0.63
# female and male grouped
all_box_data <- all_hgt_genes_individual %>% 
  spread(Group, pct_hgt_clusters)
wilcox.test(all_box_data$FMT, all_box_data$Placebo)  # p = 0.5581

# find overlap/ distinction between FMT and placebo HTGCs (total across fmt / placebo groups)
# female
female_fmt_hgt_clusters <- female_fmt_hgts_clean %>% distinct(Cluster_ID) %>% pull(Cluster_ID)
female_placebo_hgt_clusters <- female_placebo_hgts_clean %>% distinct(Cluster_ID) %>% pull(Cluster_ID)
female_venn_data <- euler(list("FMT" = female_fmt_hgt_clusters, "Placebo" = female_placebo_hgt_clusters))
fig3c_plot_female <- plot(female_venn_data, quantities = TRUE, 
                             fills = c("#b2182b", "#2166ac"), alpha = 0.7, legend = NULL)

# male
male_fmt_hgt_clusters <- male_fmt_hgts_clean %>% distinct(Cluster_ID) %>% pull(Cluster_ID)
male_placebo_hgt_clusters <- male_placebo_hgts_clean %>% distinct(Cluster_ID) %>% pull(Cluster_ID)
male_venn_data <- euler(list("FMT" = male_fmt_hgt_clusters, "Placebo" = male_placebo_hgt_clusters))
fig3c_plot_male <- plot(male_venn_data, quantities = TRUE, 
                           fills = c("#b2182b", "#2166ac"), alpha = 0.7, legend = NULL)

# quantify transfer events for FMT- and placebo-specific HTGCs

# subset for unique FMT clusters not in unique placebo clusters
# females
unique_female_fmt_hgt_clusters <- female_fmt_hgt_clusters[!(female_fmt_hgt_clusters %in% female_placebo_hgt_clusters)]
length(unique_female_fmt_hgt_clusters) # 4260

# males
unique_male_fmt_hgt_clusters <- male_fmt_hgt_clusters[!(male_fmt_hgt_clusters %in% male_placebo_hgt_clusters)]
length(unique_male_fmt_hgt_clusters) # 3321

# join with clean HGT data by cluster ID to get donor and recipient species
# females
length(unique_female_fmt_hgt_clusters) # 4260
unique_female_fmt_hgt_clusters_spp <- female_fmt_hgts_clean[female_fmt_hgts_clean$Cluster_ID %in% unique_female_fmt_hgt_clusters, ]
nrow(unique_female_fmt_hgt_clusters_spp) # 13595
# note: more rows due to multiple HGT events for each cluster ID - different donors, or multiple recipients with each cluster ID

# males
length(unique_male_fmt_hgt_clusters) # 3321
unique_male_fmt_hgt_clusters_spp <- male_fmt_hgts_clean[male_fmt_hgts_clean$Cluster_ID %in% unique_male_fmt_hgt_clusters, ]
nrow(unique_male_fmt_hgt_clusters_spp) # 7150

# subset for unique placebo clusters not in unique FMT clusters
unique_female_placebo_hgt_clusters <- female_placebo_hgt_clusters[!(female_placebo_hgt_clusters %in% female_fmt_hgt_clusters)] # female
length(unique_female_placebo_hgt_clusters) # 9794

unique_male_placebo_hgt_clusters <- male_placebo_hgt_clusters[!(male_placebo_hgt_clusters %in% male_fmt_hgt_clusters)] # male
length(unique_male_placebo_hgt_clusters) # 4379

# join with clean HGT data by cluster ID to get donor and recipient species
unique_female_placebo_hgt_clusters_spp <- female_placebo_hgts_clean[female_placebo_hgts_clean$Cluster_ID %in% unique_female_placebo_hgt_clusters, ]
nrow(unique_female_placebo_hgt_clusters_spp) # 54989

unique_male_placebo_hgt_clusters_spp <- male_placebo_hgts_clean[male_placebo_hgts_clean$Cluster_ID %in% unique_male_placebo_hgt_clusters, ]
nrow(unique_male_placebo_hgt_clusters_spp) # 11277


# calculate transfer events per unique FMT HTGC
unique_fmt_cluster_events_female <- unique_female_fmt_hgt_clusters_spp %>% # female
  group_by(Cluster_ID) %>%
  mutate(n_transfers_per_cluster = n()) %>%
  select(Cluster_ID, n_transfers_per_cluster) %>%
  distinct() %>%
  #pull(n_transfers_per_cluster) %>%
  mutate(Sex = "Female") %>%
  mutate(Type = "FMT-specific HTGCs")

unique_fmt_cluster_events_male <- unique_male_fmt_hgt_clusters_spp %>% # male
  group_by(Cluster_ID) %>%
  mutate(n_transfers_per_cluster = n()) %>%
  select(Cluster_ID, n_transfers_per_cluster) %>%
  distinct() %>%
  #pull(n_transfers_per_cluster) %>%
  mutate(Sex = "Male") %>%
  mutate(Type = "FMT-specific HTGCs")

# events per unique placebo HTGC
unique_placebo_cluster_events_female <- unique_female_placebo_hgt_clusters_spp %>% # female
  group_by(Cluster_ID) %>%
  mutate(n_transfers_per_cluster = n()) %>%
  select(Cluster_ID, n_transfers_per_cluster) %>%
  distinct() %>%
  mutate(Sex = "Female") %>%
  mutate(Type = "Placebo-specific HTGCs")

unique_placebo_cluster_events_male <- unique_male_placebo_hgt_clusters_spp %>% # male
  group_by(Cluster_ID) %>%
  mutate(n_transfers_per_cluster = n()) %>%
  select(Cluster_ID, n_transfers_per_cluster) %>%
  distinct() %>%
  mutate(Sex = "Male") %>%
  mutate(Type = "Placebo-specific HTGCs")

# join data
all_cluster_events_data <- unique_fmt_cluster_events_female %>% # unique female fmt
  bind_rows(unique_placebo_cluster_events_female) %>% # unique female placebo
  bind_rows(unique_fmt_cluster_events_male) %>% # unique male fmt
  bind_rows(unique_placebo_cluster_events_male) # unique male placebo

# for plot formatting, make category of n transfers > 20 == 21, to then edit as 20+
all_cluster_events_data_20plus <- all_cluster_events_data %>%
  mutate(n_transfers_20plus = ifelse(n_transfers_per_cluster <= 20, n_transfers_per_cluster, 21))

fig3b_plot <- all_cluster_events_data_20plus %>% # load data
  ggplot(aes(x = n_transfers_20plus, fill = Type)) +
  geom_histogram(alpha = 0.7, colour = "black", binwidth = 1) +
  scale_fill_manual(values = c("#b2182b", "#2166ac")) + 
  #scale_colour_manual(values = c("#b2182b", "#2166ac")) +
  facet_grid(Sex ~ Type) +
  xlab("Transfer events/HTGC") + ylab("Count") +
  #xlim(0, 20) + # arbitrary x axis limit for readability
  scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 21),
                     labels = c(0, 5, 10, 15, 20, 21)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4)) 

# regular plot
cluster_events_plot <- all_cluster_events_data %>% # load data
  ggplot(aes(x = n_transfers_per_cluster, fill = Type)) +
  geom_histogram(binwidth = 1) +
  scale_fill_manual(values = c("#b2182b", "#2166ac")) + 
  scale_colour_manual(values = c("#b2182b", "#2166ac")) +
  facet_grid(Sex ~ Type) +
  xlab("Transfer events/HTGC") + ylab("Count") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4)) 

# statistics
# combined
mean(all_cluster_events_data$n_transfers_per_cluster) # 3.99977 = 4
sd(all_cluster_events_data$n_transfers_per_cluster) # 4.864443
# FMT-specific HTGCs
FMT_filtered <- all_cluster_events_data %>%
  filter(Type == "FMT-specific HTGCs")
mean(FMT_filtered$n_transfers_per_cluster) # 2.736446 = 2.7
sd(FMT_filtered$n_transfers_per_cluster) # 3.217232 = 3.2
# placebo-specific HTGCs
placebo_filtered <- all_cluster_events_data %>%
  filter(Type == "Placebo-specific HTGCs")
mean(placebo_filtered$n_transfers_per_cluster) # 4.67551 = 4.7
sd(placebo_filtered$n_transfers_per_cluster) # 5.429017 = 5.4
