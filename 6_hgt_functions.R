## load libraries

library(tidyverse)
library(vegan)

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
palette_2 <- c("#44AA99", "#6699CC", "#88CCEE", "#DDCC77", "#CC6677", "#cc7fbf", "#888888")

## running

# define cog category classes
information_storage_and_processing <- c("J", "A", "K", "L", "B")
cellular_processes_and_signalling <- c("D", "Y", "V", "T", "M", "N", "Z", "W", "U", "O")
metabolism <- c("C", "G", "E", "F", "H", "I", "P", "Q")
poorly_characterised <- c("R", "S")

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
  distinct() # remove duplicate rows as abundance data represents the abundance of the entire cluster

# bind male and female data
all_cluster_abundance_data <- hgt_cluster_abundances_female_distinct %>%
  bind_rows(hgt_cluster_abundances_male_distinct)

# split rows with multiple COG categories
split_cluster_abundance_data <- all_cluster_abundance_data %>%
  mutate(single_category = strsplit(COG_Functional_cat., "")) %>%
  unnest(single_category)

# prepare relative abundance of functional classes for each sex / timepoint
cluster_abundance_data <- split_cluster_abundance_data %>%
  # include major COG functional category classification for colouring
  mutate(Class = case_when(grepl(paste(information_storage_and_processing, collapse="|"), single_category) ~ "Information storage and processing",
                           grepl(paste(cellular_processes_and_signalling, collapse="|"), single_category) ~ "Cellular processes and signalling",
                           grepl(paste(metabolism, collapse="|"), single_category) ~ "Metabolism",
                           grepl(paste(poorly_characterised, collapse="|"), single_category) ~ "Poorly characterised")) %>%
  select(-c(Cluster_ID, Group, COG_Functional_cat.)) %>%
  group_by(Recipient_ID, Sex, Timepoint, Class) %>%
  summarise(Total_CPM = sum(CPM)) %>% # total abundance of genes in each functional class for each sample
  ungroup()

# run statistics to find any over/ under-representation of COG functional classes for each sex/ timepoint
# format data (samples as rownames, features as columns)
permanova_data_female <- 
  cluster_abundance_data %>%
  filter(Sex == "Female") %>% # female only
  left_join(reduced_meta, by = c("Recipient_ID" = "Participant_ID", "Timepoint", "Sex")) %>% # join with meta data to get sample ID
  select(Sample_ID, Class, Total_CPM) %>% # select relevant columns
  spread(key = Class, value = Total_CPM, fill = 0) # make data wide format

permanova_data_male <- 
  cluster_abundance_data %>%
  filter(Sex == "Male") %>% # male only
  left_join(reduced_meta, by = c("Recipient_ID" = "Participant_ID", "Timepoint", "Sex")) %>% # join with meta data to get sample ID
  select(Sample_ID, Class, Total_CPM) %>% # select relevant columns
  spread(key = Class, value = Total_CPM, fill = 0) # make data wide format

# format meta data
reduced_meta_female <-
  reduced_meta %>%
  filter(Sample_ID %in% permanova_data_female$Sample_ID)

reduced_meta_male <-
  reduced_meta %>%
  filter(Sample_ID %in% permanova_data_male$Sample_ID)

# run PERMANOVA
adonis2(permanova_data_female %>% column_to_rownames(var = "Sample_ID") ~ Timepoint, 
        data = reduced_meta_female, 
        by = "margin", 
        permutations = 999, 
        method = "bray")
# timepoint p value = 0.874

adonis2(permanova_data_male %>% column_to_rownames(var = "Sample_ID") ~ Timepoint, 
        data = reduced_meta_male, 
        by = "margin", 
        permutations = 999, 
        method = "bray")
# timepoint p value = 0.916

# prepare average relative abundance of functional classes for each sex / timepoint
ave_cluster_abundance_data <- split_cluster_abundance_data %>%
  # include major COG functional category classification for colouring
  mutate(Class = case_when(grepl(paste(information_storage_and_processing, collapse="|"), single_category) ~ "Information storage and processing",
                           grepl(paste(cellular_processes_and_signalling, collapse="|"), single_category) ~ "Cellular processes and signalling",
                           grepl(paste(metabolism, collapse="|"), single_category) ~ "Metabolism",
                           grepl(paste(poorly_characterised, collapse="|"), single_category) ~ "Poorly characterised")) %>%
  select(-c(Recipient_ID, Cluster_ID, Group, COG_Functional_cat.)) %>%
  group_by(Sex, Timepoint, Class) %>%
  summarise(Average_CPM = mean(CPM))

# plot stacked bar chart
fig7_plot <- ave_cluster_abundance_data %>%
  # plot
  ggplot(aes(x = factor(Timepoint, level = c("Week 6", "Week 12", "Week 26")), y = Average_CPM, fill = Class)) +
  geom_bar(position = "fill", stat = "identity") +
  facet_grid(~ Sex) +
  xlab(NULL) + ylab("Average HTGC relative abundance (CPM)") +
  scale_fill_manual("COG functional class", values = palette_2[c(5, 2, 1, 3)]) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))


# repeating for individual cog categories
# run statistics to find any over/ under-representation of COG functional classes for each sex/ timepoint
# format data (samples as rownames, features as columns)
permanova_data_female_indiv <- 
  split_cluster_abundance_data %>%
  group_by(Recipient_ID, Sex, Timepoint, single_category) %>%
  summarise(Total_CPM = sum(CPM)) %>% # total abundance of genes in each functional category for each sample
  ungroup() %>%
  filter(Sex == "Female") %>% # female only
  left_join(reduced_meta, by = c("Recipient_ID" = "Participant_ID", "Timepoint", "Sex")) %>% # join with meta data to get sample ID
  select(Sample_ID, single_category, Total_CPM) %>% # select relevant columns
  spread(key = single_category, value = Total_CPM, fill = 0) # make data wide format

permanova_data_male_indiv <- 
  split_cluster_abundance_data %>%
  group_by(Recipient_ID, Sex, Timepoint, single_category) %>%
  summarise(Total_CPM = sum(CPM)) %>% # total abundance of genes in each functional category for each sample
  ungroup() %>%
  filter(Sex == "Male") %>% # male only
  left_join(reduced_meta, by = c("Recipient_ID" = "Participant_ID", "Timepoint", "Sex")) %>% # join with meta data to get sample ID
  select(Sample_ID, single_category, Total_CPM) %>% # select relevant columns
  spread(key = single_category, value = Total_CPM, fill = 0) # make data wide format

# run PERMANOVA
adonis2(permanova_data_female_indiv %>% column_to_rownames(var = "Sample_ID") ~ Timepoint, 
        data = reduced_meta_female, 
        by = "margin", 
        permutations = 999, 
        method = "bray")
# timepoint p value = 0.918

adonis2(permanova_data_male_indiv %>% column_to_rownames(var = "Sample_ID") ~ Timepoint, 
        data = reduced_meta_male, 
        by = "margin", 
        permutations = 999, 
        method = "bray")
# timepoint p value = 0.988

# prepare average relative abundance of individual categories
ave_cluster_abundance_data_individual <- split_cluster_abundance_data %>%
  select(-c(Recipient_ID, Cluster_ID, Group, COG_Functional_cat.)) %>%
  group_by(Sex, Timepoint, single_category) %>%
  summarise(Average_CPM = mean(CPM))

# load large colour palette
large_colour_palette <- c("hotpink", "#3B9AB2", "pink2", "cornflowerblue", "indianred3", "lightgreen",  "#666666", "orange", "turquoise", 
                          "#EBCC2A", "lightblue", "#BF5B17", "#BEAED4", "#386CB0", "khaki2", "slategray3", "#1B9E77", "#7570B3","maroon", 
                          "azure4", "goldenrod", "firebrick4", "aliceblue", "bisque4", "sienna4", "darksalmon", "steelblue4", "sandybrown", 
                          "royalblue", "thistle", "darkgreen", "yellow", "blue", "green", "red", "purple", "coral", "chartreuse4", 
                          "darkblue", "yellowgreen")

# plot stacked bar chart for individual COG functional categories
fig8_plot <- ave_cluster_abundance_data_individual %>%
  filter(!is.na(single_category)) %>% # remove rows with no COG category
  # plot
  ggplot(aes(x = factor(Timepoint, level = c("Week 6", "Week 12", "Week 26")), y = Average_CPM, fill = single_category)) +
  geom_bar(position = "fill", stat = "identity") +
  facet_grid(~ Sex) +
  xlab(NULL) + ylab("Average relative abundance (CPM)") +
  scale_fill_manual(name = "COG functional\ncategory", values = large_colour_palette) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4)) +
  guides(fill=guide_legend(ncol=2)) # split legend into two columns