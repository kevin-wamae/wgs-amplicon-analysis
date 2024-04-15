# *****************************************************************************#
#     An R pipeline for aggregating allele frequency data from SeekDeep
# *****************************************************************************#

# clear environment
# =============================================================================#

rm(list = ls())



# load packages
# =============================================================================#

library(data.table, quietly = TRUE)
library(tidyverse, quietly = TRUE)
# source("wf-seekdeep-illumina-amplicon/scripts/functions-resistance-profile.R")



# *****************************************************************************#
# 2. import data ----
# *****************************************************************************#

# specify the study dir name
STUDY = "wf-seekdeep-illumina-amplicon/input/ssurvey_2022 - western_kenya/2024_02_23_ilri_illumina_2x300/2024_04_11-04-seekdeep/"
# STUDY = "wf-seekdeep-illumina-amplicon/input/ssurvey_2022 - western_kenya/2023_05_25_ilri_illumina_2x300/2024_04_12-01-seekdeep-dhps/"
# STUDY = "wf-seekdeep-illumina-amplicon/input/ssurvey_2022 - western_kenya/2023_05_25_ilri_illumina_2x300/2024_04_12-01-seekdeep-dhfr/"



## __a. import fastq extraction reports (by target) ----
# =============================================================================#


### ____i. import report ----
# -----------------------------------------------------------------------------#

extProfile <- read_tsv(paste0(STUDY, "reports/allExtractionProfile.tab.txt"),
                       show_col_types = FALSE) %>%
  mutate_all(list(~str_replace(., "\\(.+\\)", ""))) %>%  # remove all characters in braces
  mutate_at(.vars = 3:ncol(.), .funs = as.numeric)       # to numeric



### ____ii. qc summary - by target ----
# -----------------------------------------------------------------------------#

extProfileTarget <- extProfile %>%
  mutate(name = str_remove(name, pattern = "MID.+")) %>% # drop MID substring
  summarise(
            # generate sum across all columns by fastq and target
            across(totalMatching:last_col(), ~sum(.x, na.rm = TRUE)),
            .by = c(inputName, name)
            ) %>%
  # select relevant columns
  select(inputName, name, totalMatching, good) %>%
  # summary of totals, all data
  mutate(
         totalAll = sum(totalMatching),
         totalAllGood = sum(good),
         totalAllGoodPerc = round(totalAllGood / totalAll, 2),
         totalAllGood = base::format(totalAllGood, big.mark = ","),
         totalAllGood = paste0(totalAllGood, " [", totalAllGoodPerc, "]"),
         ) %>%
  # summary of totals by fastq
  mutate(
         totalFastq = sum(totalMatching),
         totalFastqGood = sum(good),
         totalFastqGoodPerc = round(totalFastqGood / totalFastq, 2),
         totalFastq = base::format(totalFastq, big.mark = ","),
         totalFastqGoodPerc = paste0(totalFastq, " [", totalFastqGoodPerc, "]"),
         .by = inputName
         ) %>%
  # summary of totals by fastq and target
  mutate(
         totalTargetGoodPerc = round(good / totalMatching, 2),
         totalTargetMatching = base::format(totalMatching, big.mark = ","),
         totalTargetGoodPerc = paste0(totalTargetMatching, " [", totalTargetGoodPerc, "]"),
         ) %>%
  # format: long to wide
  pivot_wider(
              id_cols = c(inputName, totalAllGood, totalFastqGoodPerc),
              names_from = "name",
              values_from = "totalTargetGoodPerc"
              )



### ____save table ----
# -----------------------------------------------------------------------------#
write_csv(extProfileTarget, paste0(STUDY, "output/qc-read-depth-target.csv"))



## __b. import analysis data ----
# =============================================================================#

selectedClustersInfo <- read_tsv(paste0(STUDY, "/selectedClustersInfo.tab.txt.gz"),
                     show_col_types = FALSE) %>%
  mutate(source = "None") %>% # define source of samples
  filter(h_AATyped != "Untranslatable")


# summarise targets available for analysis
unique(selectedClustersInfo$p_name)



## __c. extract clusters for each available ----
# =============================================================================#

###___PfAMA1 ----
# =============================================================================#

source("wf-seekdeep-illumina-amplicon/scripts/aggregate-clusters-ama1.R")


###___PfK13 ----
# =============================================================================#

source("wf-seekdeep-illumina-amplicon/scripts/aggregate-clusters-k13.R")


###___PfMDR1 ----
# =============================================================================#

source("wf-seekdeep-illumina-amplicon/scripts/aggregate-clusters-mdr1.R")


###___PfDHPS ----
# =============================================================================#

source("wf-seekdeep-illumina-amplicon/scripts/aggregate-clusters-dhps.R")
source("wf-seekdeep-illumina-amplicon/scripts/functions-resistance-profile.R")


###___PfDHFR ----
# =============================================================================#

source("wf-seekdeep-illumina-amplicon/scripts/aggregate-clusters-dhfr.R")
source("wf-seekdeep-illumina-amplicon/scripts/functions-resistance-profile.R")



# *****************************************************************************#
# 3. compute allele frequencies ----
# *****************************************************************************#


## __PfAMA1 ----
# =============================================================================#


### ____i. coi by source ----
# -----------------------------------------------------------------------------#

(
  coi_source <- clusters_AMA1 %>%
    summarise(
              min = min(s_COI),
              mean = median(s_COI),
              max = max(s_COI),
              .by = source
              )
)



### ____ii. coi by sample ----
# -----------------------------------------------------------------------------#

coi_sample <- clusters_AMA1 %>%
  distinct(s_Sample, .keep_all = TRUE) %>%
  select(source, s_Sample, s_COI)


### ____iii. save table ----
write_csv(coi_sample, paste0(STUDY, "output/coi-by-sample.csv"))



## __PfK13 ----
# =============================================================================#

##___import the wildtype sequence ----
# -----------------------------------------------------------------------------#

fasta_K13 <- read_lines("resources-genome/fasta-cds/PfK13.txt")



##___extract the wildtype alleles/haplotypes ----
# -----------------------------------------------------------------------------#

# alleles
(
  wt_alleles <- sapply(positions_K13, function(pos) substr(fasta_K13, pos, pos))
)



# haplotypes
(
  wt_haplotype <- paste(wt_alleles, collapse = "")
)



##___compute allele frequencies ----
# -----------------------------------------------------------------------------#

freqSNP_K13 <- clusters_K13 %>%
  # select relevant columns
  # ---------------------------------#
  select(s_Sample, source, starts_with("pos")) %>%
  # collapse alleles per codon per sample
  # ---------------------------------#
  reframe(source,
          across(starts_with("pos"), ~paste(unique(sort(.)), collapse = ","),
                 .names = "{.col}"),
          .by = s_Sample) %>%
  # remove duplicates
  # ---------------------------------#
  distinct(s_Sample, .keep_all = TRUE) %>%
  # transform: wide to long
  # ---------------------------------#
  pivot_longer(
               cols = starts_with("pos"),
               names_to = "codon",
               values_to = "allele"
               ) %>%
  # drop string 'pos'
  # ---------------------------------#
  mutate(codon = str_remove(codon, "pos")) %>%
  # generate allele frequencies
  # ---------------------------------#
  summarise(count=n(), .by = c(source, codon, allele)) %>%
  arrange(source, as.numeric(codon)) %>%
  mutate(
         freq = count/sum(count) * 100, .by = c(source, codon),
         freq = round(freq,1),
         allele = str_remove_all(allele, "\\d+")
         ) %>%
  mutate(codon = as.numeric(codon)) %>%
  # merge with wildtype alleles
  # ---------------------------------#
  left_join(
            data.frame(position = positions_K13, wildtype = wt_alleles),
            by = c("codon" = "position")) %>%
  # code for infection-type
  # ---------------------------------#
  mutate(codon_allele = paste0(wildtype, codon, allele),
         variant = case_when(
                             str_detect(allele, "\\,") ~ "mixed",
                             allele == wildtype ~ "wildtype",
                             allele != wildtype ~ "mutant"),
         freq = paste0(freq, " [", count, "]")
         ) %>%
  # transform: long to wide
  # ---------------------------------#
  pivot_wider(id_cols = c(source, codon, codon_allele),
              names_from = "variant",
              values_from = "freq") %>%
  # ensure potentially missing columns are present
  # ---------------------------------#
  mutate(wildtype = if ("wildtype" %in% names(.)) wildtype else NA_character_,
         mutant = if ("mutant" %in% names(.)) mutant else NA_character_,
         mixed = if ("mixed" %in% names(.)) mixed else NA_character_
         ) %>%
  # collapse rows by source and codon position
  # ---------------------------------#
  summarise(codon_allele = paste(codon_allele, collapse=", "),
            wildtype = first(na.omit(wildtype)),
            mutant = first(na.omit(mutant)),
            mixed = first(na.omit(mixed)),
            .by = c(source, codon)
            )



### ____save table ----
# -----------------------------------------------------------------------------#
write_csv(freqSNP_K13, paste0(STUDY, "output/freq-allele-k13.csv"))



## ___compute haplotype frequencies ----
# -----------------------------------------------------------------------------#

freqHap_K13 <- clusters_K13 %>%
  reframe(
          source,
          haplotype = paste(sort(unique(haplotype)), collapse = ","),
          .by = c(s_Sample)
          ) %>%
  distinct(s_Sample, .keep_all = TRUE) %>%
  summarise(count=n(), .by = c(source, haplotype)) %>%
  mutate(
         freq = count/sum(count), .by = source,
         freq = round(freq * 100, 1),
         variant = case_when(
                             haplotype == wt_haplotype ~ "wildtype",
                             str_detect(haplotype, ",") ~ "mixed",
                             TRUE ~ "mutant",
                             ),
         total = sum(count)
         ) %>%
  arrange(source, desc(freq))



### ____save table ----
# -----------------------------------------------------------------------------#
write_csv(freqHap_K13, paste0(STUDY, "output/freq-haplotype-k13.csv"))



## __PfMDR1 ----
# =============================================================================#


##___import the wildtype sequence ----
# -----------------------------------------------------------------------------#

fasta_MDR1 <- read_lines("resources-genome/fasta-cds/PfMDR1.txt")



##___extract the wildtype alleles/haplotypes ----
# -----------------------------------------------------------------------------#

# alleles
(
  wt_alleles <- sapply(positions_MDR1, function(pos) substr(fasta_MDR1, pos, pos))
)



# haplotypes
(
  wt_haplotype <- paste(wt_alleles, collapse = "")
)


##___compute allele frequencies ----
# -----------------------------------------------------------------------------#

freqSNP_MDR1 <- clusters_MDR1 %>%
  # select relevant columns
  # ---------------------------------#
  select(s_Sample, source, starts_with("pos")) %>%
          # collapse alleles per codon per sample
  # ---------------------------------#
  reframe(source,
          across(starts_with("pos"), ~paste(unique(sort(.)), collapse = ","),
                 .names = "{.col}"),
          .by = s_Sample) %>%
  # remove duplicates
  # ---------------------------------#
  distinct(s_Sample, .keep_all = TRUE) %>%
  # transform: wide to long
  # ---------------------------------#
  pivot_longer(
               cols = starts_with("pos"),
               names_to = "codon",
               values_to = "allele"
               ) %>%
  # drop string 'pos'
  # ---------------------------------#
  mutate(codon = str_remove(codon, "pos")) %>%
  # generate allele frequencies
  # ---------------------------------#
  summarise(count=n(), .by = c(source, codon, allele)) %>%
  arrange(source, as.numeric(codon)) %>%
  mutate(
         freq = count/sum(count) * 100, .by = c(source, codon),
         freq = round(freq,1),
         allele = str_remove_all(allele, "\\d+")
         ) %>%
  mutate(codon = as.numeric(codon)) %>%
  # merge with wildtype alleles
  # ---------------------------------#
  left_join(
            data.frame(position = positions_MDR1, wildtype = wt_alleles),
            by = c("codon" = "position")) %>%
  # code for infection-type
  # ---------------------------------#
  mutate(codon_allele = paste0(wildtype, codon, allele),
         variant = case_when(
                             str_detect(allele, "\\,") ~ "mixed",
                             allele == wildtype ~ "wildtype",
                             allele != wildtype ~ "mutant"),
         freq = paste0(freq, " [", count, "]")
         ) %>%
  # transform: long to wide
  # ---------------------------------#
  pivot_wider(id_cols = c(source, codon, codon_allele),
              names_from = "variant",
              values_from = "freq") %>%
  # ensure potentially missing columns are present
  # ---------------------------------#
  mutate(wildtype = if ("wildtype" %in% names(.)) wildtype else NA_character_,
         mutant = if ("mutant" %in% names(.)) mutant else NA_character_,
         mixed = if ("mixed" %in% names(.)) mixed else NA_character_
         ) %>%
  # collapse rows by source and codon position
  # ---------------------------------#
  summarise(codon_allele = paste(codon_allele, collapse=", "),
            wildtype = first(na.omit(wildtype)),
            mutant = first(na.omit(mutant)),
            mixed = first(na.omit(mixed)),
            .by = c(source, codon)
            )



### ____save table ----
# -----------------------------------------------------------------------------#
write_csv(freqSNP_MDR1, paste0(STUDY, "output/freq-allele-mdr1.csv"))



##___compute haplotype frequencies ----
# -----------------------------------------------------------------------------#

freqHap_MDR1 <- clusters_MDR1 %>%
  reframe(
          source,
          haplotype = paste(sort(unique(haplotype)), collapse = ","),
          .by = c(s_Sample)
          ) %>%
  distinct(s_Sample, .keep_all = TRUE) %>%
  summarise(count=n(), .by = c(source, haplotype)) %>%
  mutate(
         freq = count/sum(count), .by = source,
         freq = round(freq * 100, 1),
         variant = case_when(
                             haplotype == wt_haplotype ~ "wildtype",
                             str_detect(haplotype, ",") ~ "mixed",
                             TRUE ~ "mutant",
                             ),
         total = sum(count)
         ) %>%
  arrange(source, desc(freq))



### ____save table ----
# -----------------------------------------------------------------------------#
write_csv(freqHap_MDR1, paste0(STUDY, "output/freq-haplotype-mdr1.csv"))



## __PfDHPS ----
# =============================================================================#

##___import the wildtype sequence ----
# -----------------------------------------------------------------------------#
fasta_DHPS <- read_lines("resources-genome/fasta-cds/PfDHPS.txt")



##___extract the wildtype alleles/haplotypes ----
# -----------------------------------------------------------------------------#

# alleles
(
  wt_alleles <- sapply(positions_DHPS, function(pos) substr(fasta_DHPS, pos, pos))
)



# haplotypes
(
  wt_haplotype <- paste(wt_alleles, collapse = "")
)



##___compute allele frequencies ----
# -----------------------------------------------------------------------------#

freqSNP_DHPS <- clusters_DHPS %>%
  # select relevant columns
  # ---------------------------------#
  select(s_Sample, source, starts_with("pos")) %>%
  # collapse alleles per codon per sample
  # ---------------------------------#
  reframe(source,
          across(starts_with("pos"), ~paste(unique(sort(.)), collapse = ","),
                 .names = "{.col}"),
          .by = s_Sample) %>%
  # remove duplicates
  # ---------------------------------#
  distinct(s_Sample, .keep_all = TRUE) %>%
  # transform: wide to long
  # ---------------------------------#
  pivot_longer(
               cols = starts_with("pos"),
               names_to = "codon",
               values_to = "allele"
               ) %>%
  # drop string 'pos'
  # ---------------------------------#
  mutate(codon = str_remove(codon, "pos")) %>%
  # generate allele frequencies
  # ---------------------------------#
  summarise(count=n(), .by = c(source, codon, allele)) %>%
  arrange(source, as.numeric(codon)) %>%
  mutate(
         freq = count/sum(count) * 100, .by = c(source, codon),
         freq = round(freq,1),
         allele = str_remove_all(allele, "\\d+")
         ) %>%
  mutate(codon = as.numeric(codon)) %>%
  # merge with wildtype alleles
  # ---------------------------------#
  left_join(
            data.frame(position = positions_DHPS, wildtype = wt_alleles),
            by = c("codon" = "position")) %>%
  # code for infection-type
  # ---------------------------------#
  mutate(codon_allele = paste0(wildtype, codon, allele),
         variant = case_when(
                             str_detect(allele, "\\,") ~ "mixed",
                             allele == wildtype ~ "wildtype",
                             allele != wildtype ~ "mutant"),
         freq = paste0(freq, " [", count, "]")
         ) %>%
  # transform: long to wide
  # ---------------------------------#
  pivot_wider(id_cols = c(source, codon, codon_allele),
              names_from = "variant",
              values_from = "freq") %>%
  # ensure potentially missing columns are present
  # ---------------------------------#
  mutate(wildtype = if ("wildtype" %in% names(.)) wildtype else NA_character_,
         mutant = if ("mutant" %in% names(.)) mutant else NA_character_,
         mixed = if ("mixed" %in% names(.)) mixed else NA_character_
         ) %>%
  # collapse rows by source and codon position
  # ---------------------------------#
  summarise(codon_allele = paste(codon_allele, collapse=", "),
            wildtype = first(na.omit(wildtype)),
            mutant = first(na.omit(mutant)),
            mixed = first(na.omit(mixed)),
            .by = c(source, codon)
            )



### ____save table ----
# -----------------------------------------------------------------------------#
write_csv(freqSNP_DHPS, paste0(STUDY, "output/freq-allele-dhps.csv"))



##___compute haplotype frequencies and resistance profiles ----
# -----------------------------------------------------------------------------#

freqHap_DHPS <- clusters_DHPS %>%
  select(source, s_Sample, contains("pos"), haplotype) %>%
  # define level of resistance
  mutate(
         profile_437 = case_when(pos437 == "437A" ~ "wt", .default = "mut"),
         profile_540 = case_when(pos540 == "540K" ~ "wt", .default = "mut"),
         profile_581 = case_when(pos581 == "581A" ~ "wt", .default = "mut")
         ) %>%
  rowwise() %>%
  # calculate resistance profiles for dhps and total mutations for dhps
  mutate(
         mutations = sum(profile_437 == "mut", profile_540 == "mut", profile_581 == "mut", na.rm = TRUE),
         dhps_resistance = get_resistance_profile(mutations, n=sum(!is.na(c(profile_437, profile_540, profile_581)))),
         haplotype = paste0(haplotype," (", dhps_resistance, ")")
         ) %>%
  ungroup() %>%
  reframe(
          source, dhps_resistance,
          haplotype = paste(haplotype, collapse = ","),
          dhps_resistance = paste(unique(sort(dhps_resistance)), collapse = ","),
          .by = c(s_Sample)
          ) %>%
  distinct(s_Sample, .keep_all = TRUE) %>%
  summarise(count=n(), .by = c(source, haplotype)) %>%
  mutate(
         freq = count/sum(count), .by = source,
         freq = round(freq * 100, 1),
         total = sum(count)
         ) %>%
  arrange(source, desc(freq))



### ____save table ----
# -----------------------------------------------------------------------------#
write_csv(freqHap_DHPS, paste0(STUDY, "output/freq-haplotype-dhps.csv"))



## __PfDHFR ----
# =============================================================================#


##___import the wildtype sequence ----
# -----------------------------------------------------------------------------#
fasta_DHFR <- read_lines("resources-genome/fasta-cds/PfDHFR.txt")



##___extract the wildtype alleles/haplotypes ----
# -----------------------------------------------------------------------------#


# alleles
(
  wt_alleles <- sapply(positions_DHFR, function(pos) substr(fasta_DHFR, pos, pos))
)



# haplotypes
(
  wt_haplotype <- paste(wt_alleles, collapse = "")
)



##___compute allele frequencies ----
# -----------------------------------------------------------------------------#

freqSNP_DHFR <- clusters_DHFR %>%
  # select relevant columns
  # ---------------------------------#
  select(s_Sample, source, starts_with("pos")) %>%
  # collapse alleles per codon per sample
  # ---------------------------------#
  reframe(source,
          across(starts_with("pos"), ~paste(unique(sort(.)), collapse = ","),
                 .names = "{.col}"),
          .by = s_Sample) %>%
  # remove duplicates
  # ---------------------------------#
  distinct(s_Sample, .keep_all = TRUE) %>%
  # transform: wide to long
  # ---------------------------------#
  pivot_longer(
               cols = starts_with("pos"),
               names_to = "codon",
               values_to = "allele"
               ) %>%
  # drop string 'pos'
  # ---------------------------------#
  mutate(codon = str_remove(codon, "pos")) %>%
  # generate allele frequencies
  # ---------------------------------#
  summarise(count=n(), .by = c(source, codon, allele)) %>%
  arrange(source, as.numeric(codon)) %>%
  mutate(
         freq = count/sum(count) * 100, .by = c(source, codon),
         freq = round(freq,1),
         allele = str_remove_all(allele, "\\d+")
         ) %>%
  mutate(codon = as.numeric(codon)) %>%
  # merge with wildtype alleles
  # ---------------------------------#
  left_join(
            data.frame(position = positions_DHFR, wildtype = wt_alleles),
            by = c("codon" = "position")) %>%
  # code for infection-type
  # ---------------------------------#
  mutate(codon_allele = paste0(wildtype, codon, allele),
         variant = case_when(
                             str_detect(allele, "\\,") ~ "mixed",
                             allele == wildtype ~ "wildtype",
                             allele != wildtype ~ "mutant"),
         freq = paste0(freq, " [", count, "]")
         ) %>%
  # transform: long to wide
  # ---------------------------------#
  pivot_wider(id_cols = c(source, codon, codon_allele),
              names_from = "variant",
              values_from = "freq") %>%
  # ensure potentially missing columns are present
  # ---------------------------------#
  mutate(wildtype = if ("wildtype" %in% names(.)) wildtype else NA_character_,
         mutant = if ("mutant" %in% names(.)) mutant else NA_character_,
         mixed = if ("mixed" %in% names(.)) mixed else NA_character_
         ) %>%
  # collapse rows by source and codon position
  # ---------------------------------#
  summarise(codon_allele = paste(codon_allele, collapse=", "),
            wildtype = first(na.omit(wildtype)),
            mutant = first(na.omit(mutant)),
            mixed = first(na.omit(mixed)),
            .by = c(source, codon)
            )



### ____save table ----
# -----------------------------------------------------------------------------#
write_csv(freqSNP_DHFR, paste0(STUDY, "output/freq-allele-dhfr.csv"))



##___compute haplotype frequencies and resistance profiles ----
# -----------------------------------------------------------------------------#

freqHap_DHFR <- clusters_DHFR %>%
  select(source, s_Sample, contains("pos"), haplotype) %>%
  # define level of resistance
  mutate(
         profile_51 = case_when(pos51 == "51N" ~ "wt", .default = "mut"),
         profile_59 = case_when(pos59 == "59C" ~ "wt", .default = "mut"),
         profile_108 = case_when(pos108 == "108S" ~ "wt", .default = "mut")
         ) %>%
  rowwise() %>%
  # calculate resistance profiles for dhps and total mutations for dhps
  mutate(
         mutations = sum(profile_51 == "mut", profile_59 == "mut", profile_108 == "mut", na.rm = TRUE),
         dhps_resistance = get_resistance_profile(mutations, n=sum(!is.na(c(profile_51, profile_59, profile_108)))),
         haplotype = paste0(haplotype," (", dhps_resistance, ")")
         ) %>%
  ungroup() %>%
  reframe(
          source, dhps_resistance,
          haplotype = paste(haplotype, collapse = ","),
          dhps_resistance = paste(unique(sort(dhps_resistance)), collapse = ","),
          .by = c(s_Sample)
          ) %>%
  distinct(s_Sample, .keep_all = TRUE) %>%
  summarise(count=n(), .by = c(source, haplotype)) %>%
  mutate(
         freq = count/sum(count), .by = source,
         freq = round(freq * 100, 1),
         total = sum(count)
         ) %>%
  arrange(source, desc(freq))



### ____save table ----
# -----------------------------------------------------------------------------#
write_csv(freqHap_DHFR, paste0(STUDY, "output/freq-haplotype-dhfr.csv"))


