## load libraries

library(tidyverse)
#install_version("ggalluvial", version = "0.12.3", repos = "http://cran.us.r-project.org") # stratum labelling doesn't work with current version (0.12.4)
# press 'enter'
library(ggalluvial)
library(chisq.posthoc.test)

## load data

# formatted HQ gene data from 2_hgt_analysis_genes.R
load("HQ_genes_formatted.RData")

# HGT events data from 2_hgt_analysis_genes.R
load("group_hgt_events.RData")

# high quality MAGs, high quality genes, high quality gene clusters
load("HQ-gene-mapping-filtered.RData")

# strain engraftment data
load("strain_source_17nov2022.Rdata")

# metadata
load("metadata_clean.RData")

# colour palette
palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "slategrey", "#CC79A7", "darkslateblue", "#D55E00")

## running

# get engrafted donor hgt strains at week 6 that horizontally transferred genes to recipient species
donor_engrafted_strains <- strain_source %>%
  #filter(Treatment == "FMT") %>% # get FMT engrafted strains
  filter(Source != "absent" & Source != "other" & Source != "recipient" & Source != "multiple_donors") %>% # specific donor strain
  filter(Timepoint == "Week 6") %>% # filter for bacteria engrafted at week 6
  filter(grepl("s__", Species)) %>% # filter for species
  mutate(Species = str_split(Species, "s__", simplify = TRUE)[,2]) %>% # remove 's__' prefix
  mutate(Species = str_split(Species, ".fasta", simplify = TRUE)[,1]) %>% # remove '.fasta' suffix
  mutate(Species = str_replace_all(Species, "_", " ")) # replace '_' with ' '

engrafted_hgt_female <-
  female_fmt_hgts_clean %>% # 44789 hgt events
  # extract donor ID from donor gene
  mutate(Sample_ID = str_split(donor_gene, "_", simplify = TRUE)[,1]) %>% # extract sample ID from gene ID
  mutate(Sample_ID = str_replace(Sample_ID, "wk6", "6wk")) %>% # format of sample ID from genes is reversed from meta data
  mutate(Sample_ID = str_replace(Sample_ID, "wk12", "12wk")) %>%
  mutate(Sample_ID = str_replace(Sample_ID, "wk26", "26wk")) %>%
  mutate(Sample_ID = str_replace(Sample_ID, "DF17D2a", "DF17aD2")) %>% # others have been swapped around also
  mutate(Sample_ID = str_replace(Sample_ID, "DF17D2b", "DF17bD2")) %>%
  mutate(Sample_ID = str_replace(Sample_ID, "DF16D4a", "DF16aD4")) %>%
  mutate(Sample_ID = str_replace(Sample_ID, "DF16D4b", "DF16bD4")) %>%
  left_join(reduced_meta) %>% # join with meta data
  rename("Donor_ID" = "Participant_ID") %>%
  select(c(Cluster_ID, recipient_gene, recipient_species, donor_species, Donor_ID)) %>%
  # extract recipient ID from recipient gene
  mutate(Sample_ID = str_split(recipient_gene, "_", simplify = TRUE)[,1]) %>% # extract sample ID from gene ID
  mutate(Sample_ID = str_replace(Sample_ID, "wk6", "6wk")) %>% # format of sample ID from genes is reversed from meta data
  mutate(Sample_ID = str_replace(Sample_ID, "wk12", "12wk")) %>%
  mutate(Sample_ID = str_replace(Sample_ID, "wk26", "26wk")) %>%
  mutate(Sample_ID = str_replace(Sample_ID, "DF17D2a", "DF17aD2")) %>% # others have been swapped around also
  mutate(Sample_ID = str_replace(Sample_ID, "DF17D2b", "DF17bD2")) %>%
  mutate(Sample_ID = str_replace(Sample_ID, "DF16D4a", "DF16aD4")) %>%
  mutate(Sample_ID = str_replace(Sample_ID, "DF16D4b", "DF16bD4")) %>%
  left_join(reduced_meta) %>% # join with meta data
  rename("Recipient_ID" = "Participant_ID") %>%
  select(c(Cluster_ID, Donor_ID, donor_species, Recipient_ID, recipient_species)) %>%
  # join with engrafted strain data by Donor ID, Recipient ID and donor species
  inner_join(donor_engrafted_strains, by = c("Donor_ID" = "Source", "Recipient_ID" = "Participant_ID", "donor_species" = "Species"))

engrafted_hgt_male <-
  male_fmt_hgts_clean %>% # 12801 hgt events
  # extract donor ID from donor gene
  mutate(Sample_ID = str_split(donor_gene, "_", simplify = TRUE)[,1]) %>% # extract sample ID from gene ID
  mutate(Sample_ID = str_replace(Sample_ID, "wk6", "6wk")) %>% # format of sample ID from genes is reversed from meta data
  mutate(Sample_ID = str_replace(Sample_ID, "wk12", "12wk")) %>%
  mutate(Sample_ID = str_replace(Sample_ID, "wk26", "26wk")) %>%
  mutate(Sample_ID = str_replace(Sample_ID, "DF17D2a", "DF17aD2")) %>% # others have been swapped around also
  mutate(Sample_ID = str_replace(Sample_ID, "DF17D2b", "DF17bD2")) %>%
  mutate(Sample_ID = str_replace(Sample_ID, "DF16D4a", "DF16aD4")) %>%
  mutate(Sample_ID = str_replace(Sample_ID, "DF16D4b", "DF16bD4")) %>%
  left_join(reduced_meta) %>% # join with meta data
  rename("Donor_ID" = "Participant_ID") %>%
  select(c(Cluster_ID, recipient_gene, recipient_species, donor_species, Donor_ID)) %>%
  # extract recipient ID from recipient gene
  mutate(Sample_ID = str_split(recipient_gene, "_", simplify = TRUE)[,1]) %>% # extract sample ID from gene ID
  mutate(Sample_ID = str_replace(Sample_ID, "wk6", "6wk")) %>% # format of sample ID from genes is reversed from meta data
  mutate(Sample_ID = str_replace(Sample_ID, "wk12", "12wk")) %>%
  mutate(Sample_ID = str_replace(Sample_ID, "wk26", "26wk")) %>%
  mutate(Sample_ID = str_replace(Sample_ID, "DF17D2a", "DF17aD2")) %>% # others have been swapped around also
  mutate(Sample_ID = str_replace(Sample_ID, "DF17D2b", "DF17bD2")) %>%
  mutate(Sample_ID = str_replace(Sample_ID, "DF16D4a", "DF16aD4")) %>%
  mutate(Sample_ID = str_replace(Sample_ID, "DF16D4b", "DF16bD4")) %>%
  left_join(reduced_meta) %>% # join with meta data
  rename("Recipient_ID" = "Participant_ID") %>%
  select(c(Cluster_ID, Donor_ID, donor_species, Recipient_ID, recipient_species)) %>%
  # join with engrafted strain data by Donor ID, Recipient ID and donor species
  inner_join(donor_engrafted_strains, by = c("Donor_ID" = "Source", "Recipient_ID" = "Participant_ID", "donor_species" = "Species"))

# plot engraftment-specific HGT events
fig4_plot_female <- engrafted_hgt_female %>%
  select(Donor_ID, Recipient_ID) %>%
  group_by(Donor_ID, Recipient_ID) %>% # group by donor / recipient pairing
  mutate(n_transfers = n()) %>% # calculate number of treatment-specific transfers between fmt donor / recipient pairings
  distinct() %>%
  ggplot(aes(axis1 = factor(Donor_ID), axis2 = factor(Recipient_ID), y = n_transfers)) +
  geom_alluvium(aes(fill = Donor_ID), alpha = 0.6) + 
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_fill_manual("Donor ID", drop = FALSE, values = palette[1:4]) +
  scale_x_continuous(breaks=c(1, 2), 
                     labels=c("Donor", "FMT recipient")) +
  ylab("HGT from engrafted species") +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        legend.background = element_rect(color="black", size=0.4))

fig4_plot_male <- engrafted_hgt_male %>%
  select(Donor_ID, Recipient_ID) %>%
  group_by(Donor_ID, Recipient_ID) %>% # group by donor / recipient pairing
  mutate(n_transfers = n()) %>% # calculate number of treatment-specific transfers between fmt donor / recipient pairings
  distinct() %>%
  ggplot(aes(axis1 = factor(Donor_ID), axis2 = factor(Recipient_ID), y = n_transfers)) +
  geom_alluvium(aes(fill = Donor_ID), alpha = 0.6) + 
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_fill_manual("Donor ID", drop = FALSE, values = palette[c(5, 7:9)]) + # no DM05 in results
  scale_x_continuous(breaks=c(1, 2), 
                     labels=c("Donor", "FMT recipient")) +
  ylab("HGT from engrafted species") +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        legend.background = element_rect(color="black", size=0.4))

# investigating chunk of events between DM07 and TM04
DM07_TM04 <- engrafted_hgt_male %>%
  filter(Donor_ID == "DM07",
         Recipient_ID == "TM04") %>% # 161 events
  filter(donor_species == "Bacteroides uniformis",
         recipient_species == "Bacteroides vulgatus") # 158 events

# identify inter-/ intra-phylum HGT events between engrafted donor species and host species

# make lookup table for clean (no underscores) species and phyla names
phylum_species <- HQ_mags %>%
  select(classification) %>%
  mutate(Species = str_split(classification, "s__", simplify = TRUE)[,2]) %>% # extract species of MAG
  mutate(Species_clean = str_replace_all(Species, "_.+?", "")) %>% # '?' for non-greedy matching
  select(-Species) %>% # deselect columns
  rename(Species = Species_clean) %>%
  mutate(tmp_Phylum = str_split(classification, "p__", simplify = TRUE)[,2]) %>%
  mutate(Phylum = str_split(tmp_Phylum, ";", simplify = TRUE)[,1]) %>%
  select(c(Species, Phylum)) %>%
  mutate(across(where(is.character), ~ na_if(.,""))) %>% # replace empty cells with NA
  filter(!is.na(Species) & !is.na(Phylum)) %>%
  distinct()

# add alternative phylum names
phylum_alt_species <- phylum_species %>%
  distinct(Phylum) %>%
  mutate(Phylum_alt = c("Firmicutes", "Proteobacteria", "Actinobacteria", "Firmicutes", "Firmicutes", "Bacteroidetes", "Desulfobacterota", 
                        "Verrucomicrobia", "Fusobacteria")) %>% # entering manually based on uniprot synonyms
  right_join(phylum_species)

# plot female phyla heatmap
engrafted_female_phyla_hgt <- engrafted_hgt_female %>% # female engrafted fmt-specific hgt data
  group_by(donor_species, recipient_species) %>%
  count() %>% # count combinations of donor and recipient species
  left_join(phylum_alt_species, by = c("donor_species" = "Species")) %>% # join to get phylum names for donor species
  rename("donor_phylum" = "Phylum_alt") %>%
  select(-Phylum) %>%
  left_join(phylum_alt_species, by = c("recipient_species" = "Species")) %>% # join to get phylum names for recipient species
  rename("recipient_phylum" = "Phylum_alt") %>%
  select(-Phylum)

fig5_plot_female <- engrafted_female_phyla_hgt %>%
  ggplot(aes(x = donor_species, y = recipient_species, fill = n)) +
  geom_tile() +
  facet_grid (recipient_phylum~ donor_phylum, scales = "free", space = "free") + # rows ~ columns
  xlab("Engrafted donor species") + ylab("Recipient species") +
  scale_y_discrete(limits = rev) + # descending alphabetical order
  geom_vline(xintercept = seq(1.5, length(unique(engrafted_female_phyla_hgt$donor_species))-0.5, 1), linewidth = 1, colour = "grey") + # add vertical lines
  geom_hline(yintercept = seq(1.5, length(unique(engrafted_female_phyla_hgt$recipient_species))-0.5, 1), linewidth = 1, colour = "grey") + # add horizontal lines
  scale_fill_gradient(name = "Number of\nHGT events", low = "#f3a5ae", high = "#b2182b") +
  theme_bw() +
  theme(strip.text.x = element_text(angle = 90), # rotate x axis facel labels vertical
        strip.text.y = element_text(angle = 0), # rotate y axis facet labels horizontal 
        axis.text.x = element_text(angle = 45, hjust=1), # rotate x axis labels 45 degrees
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
engrafted_female_phyla_hgt_sum <- engrafted_female_phyla_hgt %>% group_by(donor_phylum) %>% summarise(sum(n)) # quantify HGT donors

# plot male phyla heatmap
engrafted_male_phyla_hgt <- engrafted_hgt_male %>% # male engrafted fmt-specific hgt data
  group_by(donor_species, recipient_species) %>%
  count() %>% # count combinations of donor and recipient species
  left_join(phylum_alt_species, by = c("donor_species" = "Species")) %>% # join to get phylum names for donor species
  rename("donor_phylum" = "Phylum_alt") %>%
  select(-Phylum) %>%
  left_join(phylum_alt_species, by = c("recipient_species" = "Species")) %>% # join to get phylum names for recipient species
  rename("recipient_phylum" = "Phylum_alt") %>%
  select(-Phylum)

fig5_plot_male <- engrafted_male_phyla_hgt %>%
  ggplot(aes(x = donor_species, y = recipient_species, fill = n)) +
  geom_tile() +
  facet_grid (recipient_phylum~ donor_phylum, scales = "free", space = "free") +
  xlab("Engrafted donor species") + ylab("Recipient species") +
  scale_y_discrete(limits = rev) + # descending alphabetical order
  geom_vline(xintercept = seq(1.5, length(unique(engrafted_male_phyla_hgt$donor_species))-0.5, 1), linewidth = 1, colour = "grey") + # add vertical lines
  geom_hline(yintercept = seq(1.5, length(unique(engrafted_male_phyla_hgt$recipient_species))-0.5, 1), linewidth = 1, colour = "grey") + # add horizontal lines
  scale_fill_gradient(name = "Number of\nHGT events", low = "#f3a5ae", high = "#b2182b") +
  theme_bw() +
  theme(strip.text.x = element_text(angle = 90), # rotate x axis facet labels vertical
        strip.text.y = element_text(angle = 0), # rotate y axis facet labels horizontal 
        axis.text.x = element_text(angle = 45, hjust=1), # rotate x axis labels 45 degrees
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
engrafted_male_phyla_hgt_sum <- engrafted_male_phyla_hgt %>% group_by(donor_phylum) %>% summarise(sum(n)) # quantify HGT donors

# find relative abundance of background phyla (phyla represented by HQ donor genes)
HQ_donor_genes_phyla <- HQ_genes_female_clean %>%
  bind_rows(HQ_genes_male_clean) %>% # join male and female HQ genes data
  left_join(phylum_alt_species) %>% # join to get phylum names for species
  filter(Timepoint == "Donor") %>% # donor phyla
  group_by(Sex, Phylum_alt) %>% # group data
  summarise(Freq = n()) %>% # find frequencies of phyla within each sex
  mutate(Rel_abundance = Freq/sum(Freq)) %>% # find relative abundance
  mutate(Type = "background donor phyla") # make a type column

# find relative abundance of phyla in engrafted FMT donor strains with HGT
stackedbar_engrafted_hgt_data <- engrafted_hgt_female %>% mutate(Sex = "Female") %>% # load female engrafted HGT data
  bind_rows(engrafted_hgt_male %>% mutate(Sex = "Male")) %>% # join with male engrafted HGT data
  select(donor_species, Sex) %>%
  left_join(phylum_alt_species, by = c("donor_species" = "Species")) %>% # join with species:phylum lookup table
  group_by(Sex, Phylum_alt) %>% # group by sex and phylum
  summarise(Freq = n()) %>% # find frequencies of phyla within each sex
  mutate(Rel_abundance = Freq/sum(Freq)) %>% # find relative abundance
  mutate(Type = "engrafted phyla with HGT") # make a type column

# find rel abundance of phyla in engrafted FMT donor strains
stackedbar_engrafted_all_data <- donor_engrafted_strains %>% filter(Treatment == "FMT") %>% # load week 6 engrafted strains in FMT recipients
  left_join(reduced_meta) %>% # join with meta data
  select(Species, Sex) %>%
  left_join(phylum_alt_species, by = c("Species" = "Species")) %>% # join with species:phylum lookup table
  #filter(!is.na(Species.y)) %>% # removes 56 rows - too many?
  # manually include missing phylum data, looking up on UniProt taxonomy
  mutate(new_Phylum_alt = case_when(grepl("Bacteroides", Species) ~ "Bacteroidetes",
                                    grepl("Oscillibacter", Species) ~ "Firmicutes", # Bacillota
                                    grepl("Coprococcus", Species) ~ "Firmicutes", # Bacillota
                                    grepl("Ruminococcus", Species) ~ "Firmicutes", # Bacillota
                                    grepl("Eubacterium", Species) ~ "Firmicutes", # Bacillota
                                    grepl("Alistipes", Species) ~ "Bacteroidetes",
                                    grepl("Ruminococcaceae", Species) ~ "Firmicutes", # Bacillota
                                    grepl("Lachnospiraceae", Species) ~ "Firmicutes", # Bacillota
                                    grepl("Clostridium", Species) ~ "Firmicutes", # Bacillota
  )) %>% 
  mutate(Phylum_alt_replacement = coalesce(Phylum_alt, new_Phylum_alt)) %>%
  #filter(is.na(Phylum_alt_replacement)) %>% # check remaining NAs
  # only Enterobacteria phage mEp460 (phage) and Methanobrevibacter smithii (archaea) left - ok to remove
  select(-c(Phylum_alt, new_Phylum_alt)) %>%
  rename("Phylum_alt" = "Phylum_alt_replacement") %>%
  filter(!is.na(Phylum_alt)) %>% # removed phage and archaea rows
  # format
  group_by(Sex, Phylum_alt) %>% # group by sex and phylum
  summarise(Freq = n()) %>% # find frequencies of phyla within each sex
  mutate(Rel_abundance = Freq/sum(Freq)) %>% # find relative abundance
  mutate(Type = "engrafted phyla") %>% # make a type column
  bind_rows(stackedbar_engrafted_hgt_data)%>% # bind with engrafted HGT data
  bind_rows(HQ_donor_genes_phyla) # bind with all HQ donor genes phyla

# sort colouring for stacked bar plot
phylum_hexcode <- phylum_alt_species %>%
  select(Phylum_alt) %>%
  distinct() %>%
  mutate(Hexcode = c("#44AA99", "#cc7fbf", "#DDCC77", "#88CCEE", "#CC6677", "#6699CC", "#888888"))

all_engrafted_phyla_meta_alt <- stackedbar_engrafted_all_data %>%
  left_join(phylum_hexcode)

unique_all_phyla <- unique(all_engrafted_phyla_meta_alt$Phylum_alt)
my_phylum_hexcodes <- phylum_hexcode[phylum_hexcode$Phylum_alt %in% unique_all_phyla, ]
my_phylum_hexcodes_distinct <- my_phylum_hexcodes %>% arrange(Phylum_alt) %>% distinct()

# plot stacked bar
fig5c_plot <- stackedbar_engrafted_all_data %>%
  ggplot(aes(x = Type, y = Rel_abundance, fill = Phylum_alt)) +
  geom_bar(position = "fill", stat = "identity") +
  facet_grid(~ Sex) +
  scale_fill_manual("Phylum", values = my_phylum_hexcodes_distinct$Hexcode, 
                    labels = my_phylum_hexcodes_distinct$Phylum_alt) +
  xlab(NULL) + ylab("Phylum relative abundance") + 
  scale_x_discrete(labels = c("Background", "Engrafted", "Engrafted\nwith HGT")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), #formats background
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))

# statistics to compare relative abundances of phyla in 'engrafted' and 'engrafted with HGT' distributions
female <- 
  stackedbar_engrafted_all_data %>% 
  select(-Rel_abundance) %>% 
  filter(str_detect(Type, "engrafted")) %>% 
  filter(Sex == "Female") %>% 
  pivot_wider(names_from = Type, values_from = Freq, values_fill = 0) %>% 
  ungroup() %>% 
  select(-Sex)

female_names <- 
  female %>% 
  select(Phylum_alt)%>% 
  mutate(Dimension = as.character(row_number()))

female_chisq <- 
  female %>% 
  select(-Phylum_alt) 

chisq.test(female_chisq)  #p-value = 1.247e-06

female_results <- chisq.posthoc.test(female_chisq, method = "BH", round = 17)

female_results <- 
  left_join(female_names, female_results)
female_pvalues <- 
  female_results %>% 
  filter(Value == "p values", if_any(.cols = everything(), ~ . <= 0.05))

# males
male <- 
  stackedbar_engrafted_all_data %>% 
  select(-Rel_abundance) %>% 
  filter(str_detect(Type, "engrafted")) %>% 
  filter(Sex == "Male") %>% 
  pivot_wider(names_from = Type, values_from = Freq, values_fill = 0) %>% 
  ungroup() %>% 
  select(-Sex)

male_names <- 
  male %>% 
  select(Phylum_alt)%>% 
  mutate(Dimension = as.character(row_number()))

male_chisq <- 
  male %>% 
  select(-Phylum_alt) 

chisq.test(male_chisq)  #p-value = < 2.2e-16

male_results <- chisq.posthoc.test(male_chisq, method = "BH", round = 17)

male_results <- 
  left_join(male_names, male_results)
male_pvalues <- 
  male_results %>% 
  filter(Value == "p values", if_any(.cols = everything(), ~ . <= 0.05))
