## load libraries

library(tidyverse)

## load data

# metadata
load("metadata_clean.RData")

# high quality MAGs, high quality genes, high quality gene clusters
load("HQ-gene-mapping-filtered.RData")

## load functions

# fmt hgt events
tentative_hgts_fmt <- function(HQ_genes){
  
  recipients <- 
    HQ_genes %>%
    filter(Group == "FMT") %>% 
    distinct(Participant_ID) %>% 
    pull(Participant_ID)
  
  donor_genes <- 
    HQ_genes %>%
    filter(Group == "Donor") %>%
    rename(donor_gene = Gene_ID,
           donor_species = Species) %>%
    select(Cluster_ID, donor_gene, donor_species)
  
  tentative_hgt_events_fmt <- tibble()
  for (recipient in recipients) {
    
    # make a table of recipient gene clusters
    # only include necessary columns to make operations quicker
    recipient_genes <- 
      HQ_genes %>%
      filter(Participant_ID == recipient) %>%
      rename(recipient_gene = Gene_ID,
             recipient_species = Species) %>%
      select(Cluster_ID, recipient_gene, recipient_species, Timepoint)
    
    # find gene clusters present at baseline
    baseline_genes <- 
      recipient_genes %>%
      filter(Timepoint == "Baseline") %>%
      pull(Cluster_ID)
    
    # filter out gene clusters present at baseline
    recipient_genes_filtered <-
      recipient_genes %>%
      filter(!(Cluster_ID %in% baseline_genes),
             Timepoint == "Week 6") %>%
      select(-Timepoint)
    
    # find joint gene clusters with donors
    # filter out cases where the donor and recipient species match
    tentative_hgt_events_fmt <-
      inner_join(donor_genes, 
                 recipient_genes_filtered,
                 relationship = "many-to-many") %>%
      filter(donor_species != recipient_species) %>% 
      bind_rows(tentative_hgt_events_fmt)
    
  }
  tentative_hgt_events_fmt
}

# placebo hgt events
tentative_hgts_placebo <- function(HQ_genes){
  
  recipients <- 
    HQ_genes %>%
    filter(Group == "Placebo") %>% 
    distinct(Participant_ID) %>% 
    pull(Participant_ID)
  
  donor_genes <- 
    HQ_genes %>%
    filter(Group == "Donor") %>%
    rename(donor_gene = Gene_ID,
           donor_species = Species) %>%
    select(Cluster_ID, donor_gene, donor_species)
  
  tentative_hgt_events_placebo <- tibble()
  for (recipient in recipients) {
    
    # make a table of recipient gene clusters
    # only include necessary columns to make operations quicker
    recipient_genes <- 
      HQ_genes %>%
      filter(Participant_ID == recipient) %>%
      rename(recipient_gene = Gene_ID,
             recipient_species = Species) %>%
      select(Cluster_ID, recipient_gene, recipient_species, Timepoint)
    
    # find gene clusters present at baseline
    baseline_genes <- 
      recipient_genes %>%
      filter(Timepoint == "Baseline") %>%
      pull(Cluster_ID)
    
    # filter out gene clusters present at baseline
    recipient_genes_filtered <-
      recipient_genes %>%
      filter(!(Cluster_ID %in% baseline_genes),
             Timepoint == "Week 6") %>%
      select(-Timepoint)
    
    # find joint gene clusters with donors
    # filter out cases where the donor and recipient species match
    tentative_hgt_events_placebo <-
      inner_join(donor_genes, 
                 recipient_genes_filtered,
                 relationship = "many-to-many") %>%
      filter(donor_species != recipient_species) %>% 
      bind_rows(tentative_hgt_events_placebo)
    
  }
  tentative_hgt_events_placebo
}

## running

# filter the HQ MAG data to get just MAG ID and classification
mag_spp <- 
  HQ_mags %>%
  select(MAG_ID, classification)

# join the mag species data with the genes on high quality MAGs, extract sample ID and species data
HQ_genes_spp <- 
  HQ_genes %>%
  inner_join(mag_spp) %>%
  mutate(Species = str_split(classification, "s__", simplify = TRUE)[,2]) %>% # extract species of MAG
  mutate(across(where(is.character), ~ na_if(.,""))) %>% # replace empty cells (no species) with NA
  select(-c(classification, MAG_ID)) %>%
  mutate(Sample_ID = str_split(Gene_ID, "_", simplify = TRUE)[,1]) %>% # extract sample ID from gene ID
  mutate(Sample_ID = str_replace(Sample_ID, "wk6", "6wk")) %>% # format of sample ID from genes is reversed from metadata
  mutate(Sample_ID = str_replace(Sample_ID, "wk12", "12wk")) %>%
  mutate(Sample_ID = str_replace(Sample_ID, "wk26", "26wk")) %>%
  mutate(Sample_ID = str_replace(Sample_ID, "DF17D2a", "DF17aD2")) %>% # others have been swapped around also
  mutate(Sample_ID = str_replace(Sample_ID, "DF17D2b", "DF17bD2")) %>%
  mutate(Sample_ID = str_replace(Sample_ID, "DF16D4a", "DF16aD4")) %>%
  mutate(Sample_ID = str_replace(Sample_ID, "DF16D4b", "DF16bD4"))

# join metadata with the HQ gene and MAG data
HQ_genes_meta <-
  HQ_genes_spp %>%
  left_join(reduced_meta)

# filter out rows with missing species data
HQ_genes_meta_spp <-
  HQ_genes_meta %>%
  filter(!is.na(Species))

# make separate dataframes for males and females
HQ_genes_female <-
  HQ_genes_meta_spp %>%
  filter(Sex == "Female")

HQ_genes_male <-
  HQ_genes_meta_spp %>%
  filter(Sex == "Male")

# remove underscores from genera and species
HQ_genes_female_clean <- HQ_genes_female %>%
  mutate(Species_clean = str_replace_all(Species, "_.+?", "")) %>% # '?' for non-greedy matching
  select(-c(Species, Contig_ID, Sample_ID)) %>% # deselect columns
  rename(Species = Species_clean)

HQ_genes_male_clean <- HQ_genes_male %>%
  mutate(Species_clean = str_replace_all(Species, "_.+?", "")) %>% # '?' for non-greedy matching
  select(-c(Species, Contig_ID, Sample_ID)) %>% # deselect columns
  rename(Species = Species_clean)

# get hgt events

# female fmt
female_fmt_hgts_clean <- tentative_hgts_fmt(HQ_genes_female_clean)
nrow(female_fmt_hgts_clean) # 44789 hgt events
nrow(female_fmt_hgts_clean %>% distinct(recipient_gene)) # 21376 distinct recipient genes involved

# female placebo
female_placebo_hgts_clean <- tentative_hgts_placebo(HQ_genes_female_clean)
nrow(female_placebo_hgts_clean) # 93158 hgt events
nrow(female_placebo_hgts_clean %>% distinct(recipient_gene)) # 29901 distinct recipient genes involved

# male fmt
male_fmt_hgts_clean <- tentative_hgts_fmt(HQ_genes_male_clean)
nrow(male_fmt_hgts_clean) # 12801 hgt events
nrow(male_fmt_hgts_clean %>% distinct(recipient_gene)) # 7396 distinct recipient genes involved

# male placebo
male_placebo_hgts_clean <- tentative_hgts_placebo(HQ_genes_male_clean)
nrow(male_placebo_hgts_clean) # 18115 hgt events
nrow(male_placebo_hgts_clean %>% distinct(recipient_gene)) # 8769 distinct recipient genes involved