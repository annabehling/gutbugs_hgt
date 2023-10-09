# gutbugs_hgt

## The Gut Bugs Trial

The Gut Bugs Trial was a double-blinded randomised placebo-controlled trial that assessed the efficacy of fecal microbiome transplantation (FMT) to treat adolescent obesity and improve metabolism. The Gut Bugs Trial recruited 87 adolescents (aged 14-18 years) with a BMI >30 kg/m2 and randomised them 1:1 to receive a single dose of 28 FMT or placebo capsules. FMT capsules contained concentrated fecal material derived from 4 healthy same-sex lean donors (7 capsules per donor). Placebo capsules contained saline. Recipients were clinically assessed at baseline, and at 6-, 12-, and 26-weeks post-treatment. The primary objective was a change in BMI SDS at week 6. Secondary objectives included a variety of metabolic health parameters, body composition, and gut microbiome alterations.

Donor and recipient stool samples collected at each clinical assessment underwent shotgun metagenomic sequencing. In total, 381 metagenomes were analysed.

Protocol paper: https://doi.org/10.1136/bmjopen-2018-026174

Trial paper: https://doi.org/10.1001/jamanetworkopen.2020.30415

Metagenomic data: https://doi.org/10.1186/s40168-021-01060-7

## Horizontal gene transfer analysis

This repository contains the following R scripts used to analyse metagenomic data from the Gut Bugs Trial for evidence of horizontal gene transfer (HGT).

- 1_hgt_analysis_waafle.R : analysis of [WAAFLE](https://github.com/biobakery/waafle) HGT output
- 2_hgt_analysis_genes.R : gene-based detection of HGT events
- 3_hgt_quantify.R : analysis of gene-based HGT events
- 4_hgt_engrafted.R : identification of HGT from engrafted donor strains
- 5_hgt_retention.R : quantification of transferred gene clusters
- 6_hgt_functions.R : determination of transferred gene cluster functions
