## load libraries

library(tidyverse)
library(ggbeeswarm)
library(rstatix)
library(DescTools)
library(lme4)
library(lmerTest)

## load data

# lists of dataframes containing of contigs for each sample (created from WAAFLE output)
load("hgt_dfs.RData") # contigs with HGT event
load("no_hgt_dfs.RData") # contigs with no HGT event
load("unclassified_dfs.RData") # unclassified contigs

# metadata
load("metadata_clean.RData")

# metaphlan species data
load("bacterial_species.RData")

# colour palette
palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "slategrey", "#CC79A7", "darkslateblue", "#D55E00")

## running

# extract sample IDs
sample_IDs <- reduced_meta %>% pull(Sample_ID)

# find number of each type of contig for each sample
nrow_hgt_dfs <- lapply(hgt_dfs, nrow) # HGT
nrow_no_hgt_dfs <- lapply(no_hgt_dfs, nrow) # no HGT
nrow_unclassified_dfs <- lapply(unclassified_dfs, nrow) # unclassified
nrow_total <- mapply(sum, nrow_hgt_dfs, nrow_no_hgt_dfs, nrow_unclassified_dfs) # total

# find number of each type of contig for each sample
pct_hgt <- unlist(nrow_hgt_dfs)/nrow_total*100 # HGT
pct_no_hgt <- unlist(nrow_no_hgt_dfs)/nrow_total*100 # no HGT
pct_unclassified <- unlist(nrow_unclassified_dfs)/nrow_total*100 # unclassified

# make dataframes of contig counts
fig1_n_data <- data.frame(sample_IDs, unlist(nrow_hgt_dfs), unlist(nrow_no_hgt_dfs), unlist(nrow_unclassified_dfs))
colnames(fig1_n_data) <- c("Sample_ID", "n_HGT_contigs", "n_no_HGT_contigs", "n_unclassified_contigs")

# make dataframes of contig percentages
fig1_pct_data <- data.frame(sample_IDs, pct_hgt, pct_no_hgt, pct_unclassified)
colnames(fig1_pct_data) <- c("Sample_ID", "pct_HGT_contigs", "pct_no_HGT_contigs", "pct_unclassified_contigs")

# get mean, min, max percentage values for results
# HGT contigs
mean(fig1_pct_data$pct_HGT_contigs) #0.1528122
min(fig1_pct_data$pct_HGT_contigs) #0.07412074
max(fig1_pct_data$pct_HGT_contigs) #0.2556661
# no HGT contigs
mean(fig1_pct_data$pct_no_HGT_contigs) #58.6406
min(fig1_pct_data$pct_no_HGT_contigs) #38.1893
max(fig1_pct_data$pct_no_HGT_contigs) #80.33072
# unclassified contigs
mean(fig1_pct_data$pct_unclassified_contigs) #41.20659
min(fig1_pct_data$pct_unclassified_contigs) #19.54895
max(fig1_pct_data$pct_unclassified_contigs) #61.72985

# plot contig percentages for each sample
fig1_plot <- fig1_pct_data %>%
  pivot_longer(cols = `pct_HGT_contigs`:`pct_unclassified_contigs`, # reshape data to long format
               names_to = "Contig_type",
               values_to = "Percentage_value") %>%
  ggplot(aes(x = factor(Contig_type, level = c("pct_HGT_contigs", "pct_no_HGT_contigs", "pct_unclassified_contigs")), 
             y = Percentage_value, 
             fill = Contig_type)) +
  geom_quasirandom(shape = 21, dodge.width = 0.75, size = 1.7) +
  geom_boxplot(alpha = 0.6, outlier.colour = NA) +
  scale_fill_manual(values = palette, name = NULL) +
  scale_colour_manual(values = palette, name = NULL) +
  scale_y_log10() +
  scale_x_discrete(labels=c("pct_HGT_contigs" = "With HGT", "pct_no_HGT_contigs" = "Without HGT", "pct_unclassified_contigs" = "Unclassified")) +
  xlab("Contig type") + ylab("Percentage of sample contigs (log base 10)") + 
  theme_bw() +
  theme(axis.text.x = element_text(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none")

# compare FMT and placebo HGT events
metaphlan_spp <- rowSums(species != 0) # number of species in each sample (species richness)
hgt_contigs <- unlist(nrow_hgt_dfs) # number of contigs in each sample
hgt_cont_per_spp <- hgt_contigs/metaphlan_spp # number of contigs per species
richness_df <- data.frame(sample_IDs, hgt_cont_per_spp) # number of contigs per species for each sample
names(richness_df)[1] <- "Sample_ID" # name sample ID column
richness_meta <- inner_join(richness_df, reduced_meta, by = "Sample_ID") # merge with metadata

# plot HGT events per species in each recipient sample
fig2_plot <- richness_meta %>%
  filter(!grepl("Donor", Timepoint)) %>% # filter for only recipient data
  ggplot(aes(x = factor(Timepoint, level = c("Baseline", "Week 6", "Week 12", "Week 26")), y = hgt_cont_per_spp, fill = Group)) +
  geom_quasirandom(shape = 21, dodge.width = 0.75, size = 1.7) +
  geom_boxplot(alpha = 0.6, outlier.colour = NA) +
  scale_fill_manual(values = c("#b2182b", "#2166ac")) +
  scale_colour_manual(values = c("#b2182b", "#2166ac"), name = NULL) +
  scale_y_continuous(limits = c(0, 4)) +
  facet_wrap(~Sex) +
  xlab(NULL) + ylab("HGT events/species") + 
  theme_bw() +
  theme(axis.text.x = element_text(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))

# prepare data for linear mixed effect model
full_dataset <- richness_meta %>% # full dataset (include all participants)
  filter(!grepl("Donor", Timepoint)) %>% # filter for only recipient data
  select(c(Participant_ID, Timepoint, hgt_cont_per_spp)) %>% # select relevant columns
  reshape(idvar = "Participant_ID", timevar = "Timepoint", direction = "wide") %>% # make wide format
  inner_join(richness_meta %>% select(c(Participant_ID, Group, Sex))) %>% # join with selected metadata
  distinct() %>% # select only distinct rows (duplicates after inner_join)
  rename("Baseline" = "hgt_cont_per_spp.Baseline", "Week 6" = "hgt_cont_per_spp.Week 6", 
         "Week 12" = "hgt_cont_per_spp.Week 12", "Week 26" = "hgt_cont_per_spp.Week 26") %>% # rename columns
  select(c(Participant_ID, Sex, Group, Baseline, `Week 6`, `Week 12`, `Week 26`)) %>%
  gather(key = "Timepoint", value = "hgt_cont_per_spp", Baseline, `Week 6`, `Week 12`, `Week 26`) %>% # gather data back to long format
  convert_as_factor(Participant_ID, Sex, Group, Timepoint) %>%
  na.omit() # remove samples with missing data
full_dataset$Group <- factor(full_dataset$Group, levels = c("Placebo", "FMT")) # reorder Group levels

# fit full model
full_lmer <- lmer(hgt_cont_per_spp ~ Sex *Group * Timepoint+(1|Participant_ID),data = full_dataset) #full dataset ok for lmm
summary(full_lmer)

# fit reduced model (no interactions between fixed effects)
reduced_lmer <- lmer(hgt_cont_per_spp ~ Sex +Group + Timepoint+(1|Participant_ID),data =full_dataset)
summary(reduced_lmer)

# run anova to compare full and reduced models
anova(full_lmer, reduced_lmer) # p = 0.1081
# reduced model fits better - non-significant p value and lower AIC value

# get confidence interval of reduced model fixed effect coefficients
confint(reduced_lmer)

## supplementary file 1
suppl1 <- reduced_meta %>%
  inner_join(fig1_n_data) %>% # contig counts
  rename(`Number of contigs with HGT` = n_HGT_contigs,
         `Number of contigs without HGT` = n_no_HGT_contigs,
         `Number of unclassified contigs` = n_unclassified_contigs)