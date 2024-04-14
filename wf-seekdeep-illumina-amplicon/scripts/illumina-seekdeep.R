# This R script processes output from SeekDeep 

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
source("wf-seekdeep-illumina-amplicon/scripts/functions.R")



# *****************************************************************************#
# 2. import data ----
# *****************************************************************************#

# specify the study dir name
STUDY = "wf-seekdeep-illumina-amplicon/input/ssurvey_2022 - western_kenya/2024_02_23_ilri_illumina_2x300/2024_04_11-04-seekdeep/"
# STUDY = "wf-seekdeep-illumina-amplicon/input/ssurvey_2022 - western_kenya/2023_05_25_ilri_illumina_2x300/2024_04_12-01-seekdeep-dhps/"
# STUDY = "wf-seekdeep-illumina-amplicon/input/ssurvey_2022 - western_kenya/2023_05_25_ilri_illumina_2x300/2024_04_12-01-seekdeep-dhfr/"


## __a. import FASTQ extraction reports ----
# =============================================================================#


### ___i. extraction reports by fastq ----
# -----------------------------------------------------------------------------#

extStats <- read_tsv(paste0(STUDY, "/reports/allExtractionStats.tab.txt"),
                     show_col_types = FALSE) %>%
  mutate_all(list(~str_replace(., "\\(.+\\)", ""))) %>%  # remove all characters in braces
  mutate_at(.vars = 2:ncol(.), .funs = as.numeric)       # to numeric



### ___ii. fastq qc summary ----
# -----------------------------------------------------------------------------#

extStats %>%
  mutate(filteredOut = Total - used) %>%
  summarise(
            sumTotal = sum(Total),           # total reads-input
            sumUsed = sum(used),             # total reads-used
            sumDiscarded = sum(filteredOut)  # total reads-discarded
            ) %>%
  # apply format with thousands separator to the resulting summary
  mutate(across(everything(), ~ base::format(.x, big.mark = ",")))


# extraction reports by FASTQ and gene
# -----------------------------------------------------------------------------#

extProfile <- read_tsv(paste0(STUDY, "reports/allExtractionProfile.tab.txt"),
                       show_col_types = FALSE) %>%
  mutate_all(list(~str_replace(., "\\(.+\\)", ""))) %>%  # remove all characters in braces
  mutate_at(.vars = 3:ncol(.), .funs = as.numeric)       # to numeric



## __b. import analysis data ----
# =============================================================================#

###__extract clusters for each target ----
# # -----------------------------------------------------------------------------#
selectedClustersInfo <- read_tsv(paste0(STUDY, "/selectedClustersInfo.tab.txt.gz"),
                     show_col_types = FALSE) %>%
  mutate(source = "None") %>% # define source of samples
  filter(h_AATyped != "Untranslatable")



# summarise targets available for analysis
unique(selectedClustersInfo$p_name)



###___AMA1 ----
# =============================================================================#


# generate a vector of polymorphic codons
# -----------------------------------------------------------------------------#
(
  positions_AMA1 <- selectedClustersInfo %>%
    filter(p_name == "PFAMA1" & str_detect(h_AATyped, "^PF3D7")) %>%
    select(h_AATyped) %>%
    slice(1) %>%
    mutate(h_AATyped = str_replace(h_AATyped, "PF3D7_.+--", ""),
           numbers = str_extract_all(h_AATyped, "\\d+")) %>%
    unnest(numbers) %>%
    mutate(numbers = as.numeric(numbers)) %>%
    pull(numbers)
)



# filter variant data
# -----------------------------------------------------------------------------#
clusters_AMA1 <- selectedClustersInfo %>%
  filter(p_name == "PFAMA1") %>%
  mutate(
        purrr::map_dfc(
                      set_names(positions_AMA1, paste0("pos", positions_AMA1)),
                      ~ str_extract(string = h_AATyped, pattern = paste0(.x, "."))
                      )) %>%
  rowwise() %>%
  mutate(
         codon_pos = paste(c_across(all_of(paste0("pos", positions_AMA1))), collapse = ", "),
         codon_pos = str_remove_all(string = codon_pos, pattern = "-."),
         haplotype = str_remove_all(string = codon_pos, pattern = "\\d+|,\\s")
         ) %>%
  ungroup()



###___K13 ----
# =============================================================================#


# vector of polymorphic codons
# -----------------------------------------------------------------------------#
(
  positions_K13 <- selectedClustersInfo %>%
    filter(p_name == "PFK13" & str_detect(h_AATyped, "^PF3D7")) %>%
    select(h_AATyped) %>%
    slice(1) %>%
    mutate(h_AATyped = str_replace(h_AATyped, "PF3D7_.+--", ""),
           numbers = str_extract_all(h_AATyped, "\\d+")) %>%
    unnest(numbers) %>%
    mutate(numbers = as.numeric(numbers)) %>%
    pull(numbers)
)


# filter variant data
# -----------------------------------------------------------------------------#
clusters_K13 <- selectedClustersInfo %>%
  filter(p_name == "PFK13") %>%
  mutate(
        purrr::map_dfc(
                      set_names(positions_K13, paste0("pos", positions_K13)),
                      ~ str_extract(string = h_AATyped, pattern = paste0(.x, "."))
                      )) %>%
  rowwise() %>%
  mutate(
         codon_pos = paste(c_across(all_of(paste0("pos", positions_K13))), collapse = ", "),
         codon_pos = str_remove_all(string = codon_pos, pattern = "-."),
         haplotype = str_remove_all(string = codon_pos, pattern = "\\d+|,\\s")
         ) %>%
  ungroup()



###___MDR1 ----
# =============================================================================#


# vector of polymorphic codons
# -----------------------------------------------------------------------------#
(
  positions_MDR1 <- selectedClustersInfo %>%
    filter(p_name == "PFMDR1" & str_detect(h_AATyped, "^PF3D7")) %>%
    # mutate(h_AATyped = str_remove(h_AATyped, ":136.")) %>%  # drop codon 16
    select(h_AATyped) %>%
    slice(1) %>%
    mutate(h_AATyped = str_replace(h_AATyped, "PF3D7_.+--", ""),
           numbers = str_extract_all(h_AATyped, "\\d+")) %>%
    unnest(numbers) %>%
    mutate(numbers = as.numeric(numbers)) %>%
    pull(numbers)
)



# filter variant data
# -----------------------------------------------------------------------------#
clusters_MDR1 <- selectedClustersInfo %>%
  filter(p_name == "PFMDR1") %>%
  mutate(
        purrr::map_dfc(
                      set_names(positions_MDR1, paste0("pos", positions_MDR1)),
                      ~ str_extract(string = h_AATyped, pattern = paste0(.x, "."))
                      )) %>%
  rowwise() %>%
  mutate(
         codon_pos = paste(c_across(all_of(paste0("pos", positions_MDR1))), collapse = ", "),
         codon_pos = str_remove_all(string = codon_pos, pattern = "-."),
         haplotype = str_remove_all(string = codon_pos, pattern = "\\d+|,\\s")
         ) %>%
  ungroup()



###___DHPS ----
# =============================================================================#


# vector of polymorphic codons
(
  positions_DHPS <- selectedClustersInfo %>%
    filter(p_name == "PFDHPS" & str_detect(h_AATyped, "^PF3D7")) %>%
    select(h_AATyped) %>%
    slice(1) %>%
    mutate(h_AATyped = str_replace(h_AATyped, "PF3D7_.+--", ""),
           numbers = str_extract_all(h_AATyped, "\\d+")) %>%
    unnest(numbers) %>%
    mutate(numbers = as.numeric(numbers)) %>%
    pull(numbers)
)



# filter variant data
# -----------------------------------------------------------------------------#
clusters_DHPS <- selectedClustersInfo %>%
  filter(p_name == "PFDHPS") %>%
  mutate(
        purrr::map_dfc(
                      set_names(positions_DHPS, paste0("pos", positions_DHPS)),
                      ~ str_extract(string = h_AATyped, pattern = paste0(.x, "."))
                      )) %>%
  rowwise() %>%
  mutate(
         codon_pos = paste(c_across(all_of(paste0("pos", positions_DHPS))), collapse = ", "),
         codon_pos = str_remove_all(string = codon_pos, pattern = "-."),
         haplotype = str_remove_all(string = codon_pos, pattern = "\\d+|,\\s")
         ) %>%
  ungroup()



###___DHFR ----
# =============================================================================#


# vector of polymorphic codons
(
  positions_DHFR <- selectedClustersInfo %>%
    filter(p_name == "PFDHFR" & str_detect(h_AATyped, "^PF3D7")) %>%
    mutate(h_AATyped = str_remove(h_AATyped, "16A:")) %>%  # drop codon 16
    select(h_AATyped) %>%
    slice(1) %>%
    mutate(h_AATyped = str_replace(h_AATyped, "PF3D7_.+--", ""),
           numbers = str_extract_all(h_AATyped, "\\d+")) %>%
    unnest(numbers) %>%
    mutate(numbers = as.numeric(numbers)) %>%
    pull(numbers)
)



# filter variant data
# -----------------------------------------------------------------------------#
clusters_DHFR <- selectedClustersInfo %>%
  filter(p_name == "PFDHFR") %>%
  mutate(
        purrr::map_dfc(
                      set_names(positions_DHFR, paste0("pos", positions_DHFR)),
                      ~ str_extract(string = h_AATyped, pattern = paste0(.x, "."))
                      )) %>%
  rowwise() %>%
  mutate(
         codon_pos = paste(c_across(all_of(paste0("pos", positions_DHFR))), collapse = ", "),
         codon_pos = str_remove_all(string = codon_pos, pattern = "-."),
         haplotype = str_remove_all(string = codon_pos, pattern = "\\d+|,\\s"),
         ) %>%
  ungroup()



# *****************************************************************************#
# 3. compute allele frequencies ----
# *****************************************************************************#


## __AMA1 ----
# =============================================================================#


# coi by source
# -----------------------------------------------------------------------------#
(
  coi_AMA1 <- clusters_AMA1 %>%
    summarise(
              min = min(s_COI),
              mean = median(s_COI),
              max = max(s_COI),
              .by = source
              )
)



## __K13 ----
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
  select(s_Sample, source, starts_with("pos")) %>%
  reframe(
          # collapse alleles per codon per sample
          source,
          across(starts_with("pos"), ~paste(unique(sort(.)), collapse = ","),
                 .names = "{.col}"),
          .by = s_Sample) %>%
  distinct(s_Sample, .keep_all = TRUE) %>%
  pivot_longer(
               cols = starts_with("pos"), # dynamically select columns that start with "pos"
               names_to = "codon",
               values_to = "allele"
               ) %>%
  mutate(codon = str_remove(codon, "pos")) %>%
  summarise(count=n(), .by = c(source, codon, allele)) %>%
  arrange(source, as.numeric(codon)) %>%
  mutate(
         freq = count/sum(count) * 100, .by = c(source, codon),
         freq = round(freq,1),
         allele = str_remove_all(allele, "\\d+")
         ) %>%
  mutate(codon = as.numeric(codon)) %>%
  left_join(data.frame(position = positions_K13, wildtype = wt_alleles),
            by = c("codon" = "position")) %>%
  mutate(
         codon_allele = paste0(wildtype, codon, allele),
         variant = case_when(
                             str_detect(allele, "\\,") ~ "mixed",
                             allele == wildtype ~ "wildtype",
                             allele != wildtype ~ "mutant"),
         freq = paste0(freq, " [", count, "]")
         ) %>%
  pivot_wider(
              id_cols = c(source, codon, codon_allele),
              names_from = "variant",
              values_from = "freq") %>%
  summarise(
            codon_allele = paste(codon_allele, collapse=", "), 
            wildtype = first(na.omit(wildtype)),
            mutant = first(na.omit(mutant)),
            mixed = first(na.omit(mixed)),
            .by = c(source, codon)
            )



##___compute haplotype frequencies ----
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



## __MDR1 ----
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
  select(s_Sample, source, starts_with("pos")) %>%
  reframe(
          # collapse alleles per codon per sample
          source,
          across(starts_with("pos"), ~paste(unique(sort(.)), collapse = ","),
                 .names = "{.col}"),
          .by = s_Sample) %>%
  distinct(s_Sample, .keep_all = TRUE) %>%
  pivot_longer(
               cols = starts_with("pos"), # dynamically select columns that start with "pos"
               names_to = "codon",
               values_to = "allele"
               ) %>%
  mutate(codon = str_remove(codon, "pos")) %>%
  summarise(count=n(), .by = c(source, codon, allele)) %>%
  arrange(source, as.numeric(codon)) %>%
  mutate(
         freq = count/sum(count) * 100, .by = c(source, codon),
         freq = round(freq,1),
         allele = str_remove_all(allele, "\\d+")
         ) %>%
  mutate(codon = as.numeric(codon)) %>%
  left_join(data.frame(position = positions_MDR1, wildtype = wt_alleles),
            by = c("codon" = "position")) %>%
  mutate(
         codon_allele = paste0(wildtype, codon, allele),
         variant = case_when(
                             str_detect(allele, "\\,") ~ "mixed",
                             allele == wildtype ~ "wildtype",
                             allele != wildtype ~ "mutant"),
         freq = paste0(freq, " [", count, "]")
         ) %>%
  pivot_wider(
              id_cols = c(source, codon, codon_allele),
              names_from = "variant",
              values_from = "freq") %>%
  summarise(
            codon_allele = paste(codon_allele, collapse=", "), 
            wildtype = first(na.omit(wildtype)),
            mutant = first(na.omit(mutant)),
            mixed = first(na.omit(mixed)),
            .by = c(source, codon)
            )



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



## __DHPS ----
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

###____table ----
# --------------------------------------------#

freqSNP_DHPS <- clusters_DHPS %>%
  select(s_Sample, source, starts_with("pos")) %>%
  reframe( # collapse alleles per codon per sample
          source,
          across(starts_with("pos"), ~paste(unique(sort(.)), collapse = ","),
                 .names = "{.col}"),
          .by = s_Sample) %>%
  distinct(s_Sample, .keep_all = TRUE) %>%
  pivot_longer(
               cols = starts_with("pos"), # dynamically select columns that start with "pos"
               names_to = "codon",
               values_to = "allele"
               ) %>%
  mutate(codon = str_remove(codon, "pos")) %>%
  summarise(count=n(), .by = c(source, codon, allele)) %>%
  arrange(source, as.numeric(codon)) %>%
  mutate(
         freq = count/sum(count) * 100, .by = c(source, codon),
         freq = round(freq,1),
         allele = str_remove_all(allele, "\\d+")
         ) %>%
  mutate(codon = as.numeric(codon)) %>%
  left_join(data.frame(position = positions_DHPS, wildtype = wt_alleles),
            by = c("codon" = "position")) %>%
  mutate(
         codon_allele = paste0(wildtype, codon, allele),
         variant = case_when(
                             str_detect(allele, "\\,") ~ "mixed",
                             allele == wildtype ~ "wildtype",
                             allele != wildtype ~ "mutant"),
         freq = paste0(freq, " [", count, "]")
         ) %>%
  pivot_wider(
              id_cols = c(source, codon, codon_allele),
              names_from = "variant",
              values_from = "freq") %>%
  summarise(
            codon_allele = paste(codon_allele, collapse=", "), 
            wildtype = first(na.omit(wildtype)),
            mutant = first(na.omit(mutant)),
            mixed = first(na.omit(mixed)),
            .by = c(source, codon)
            )



###____plot ----
# --------------------------------------------#

(
  my_plot <- clusters_DHPS %>%
    select(s_Sample, source, starts_with("pos")) %>%
    reframe( # collapse alleles per codon per sample
            source,
            across(starts_with("pos"), ~paste(unique(sort(.)), collapse = ","),
                   .names = "{.col}"),
            .by = s_Sample) %>%
    distinct(s_Sample, .keep_all = TRUE) %>%
    pivot_longer(
                 cols = starts_with("pos"), # dynamically select columns that start with "pos"
                 names_to = "codon",
                 values_to = "allele"
                 ) %>%
    mutate(codon = str_remove(codon, "pos")) %>%
    summarise(count=n(), .by = c(source, codon, allele)) %>%
    arrange(source, as.numeric(codon)) %>%
    mutate(
           freq = count/sum(count) * 100, .by = c(source, codon),
           freq = round(freq,1),
           allele = str_remove_all(allele, "\\d+")
           ) %>%
    mutate(codon = as.numeric(codon)) %>%
    left_join(data.frame(position = positions_DHPS, wildtype = wt_alleles),
              by = c("codon" = "position")) %>%
    mutate(variant = case_when(
                               str_detect(allele, "\\,") ~ "mixed",
                               allele == wildtype ~ "wildtype",
                               allele != wildtype ~ "mutant"),
           codon_allele = paste(wildtype, codon, allele, sep = "")
           ) %>%
    # drop positions not associated with resistance
    filter(! codon %in% c(431, 436)) %>%
    ggplot(aes(x=as.factor(codon), y=freq, fill=variant)) +
      geom_bar(stat = "identity", position = "stack", width = 0.7) +
      ggprism::theme_prism() +
      scale_fill_manual(values=c("grey","maroon", "black")) +
      facet_wrap(. ~ source, scales = "free_x") +
      labs(x = "Codon Position", y = "Relative Frequency", fill = "Genotype") +
      scale_y_continuous(limits=c(0, 101), breaks=seq(0, 100, by=10))
)



###____save plot ----
# -----------------------------------------------------------------------------#

jpeg(filename = "plot/dhps-snps.jpeg",
     width = 8, height = 5, units = 'in', res = 600,
     # type = "cairo",       #options::"cairo", "Xlib", "quartz"
     # antialias = "none" #options::"default", "none", "gray", "subpixel"
); print(my_plot); dev.off()



##___compute haplotype frequencies and resistance profiles ----
# -----------------------------------------------------------------------------#

###____table ----
# --------------------------------------------#

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
          dhps_resistance = paste(dhps_resistance, collapse = ","),
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



###____plot ----
# --------------------------------------------#

(
  my_plot <- clusters_DHPS %>%
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
           dhfr_resistance = get_resistance_profile(mutations, n=sum(!is.na(c(profile_437, profile_540, profile_581)))),
           haplotype = paste0(haplotype," (", dhfr_resistance, ")")
           ) %>%
    ungroup() %>%
    reframe(
            source, dhfr_resistance,
            haplotype = paste(haplotype, collapse = ","),
            dhfr_resistance = paste(dhfr_resistance, collapse = ","),
            .by = c(s_Sample)
            ) %>%
    distinct(s_Sample, .keep_all = TRUE) %>%
    summarise(count=n(), .by = c(source, dhfr_resistance)) %>%
    mutate(
           profile = case_when(str_detect(dhfr_resistance, ",") ~ "mixed",
                               .default = dhfr_resistance),
           freq = count/sum(count) * 100, .by = source
           ) %>%
    ggplot(aes(x=source, y=freq, group=profile, fill=profile)) +
      geom_bar(stat = "identity", position = "stack", width = 0.7) +
      ggprism::theme_prism() +
      scale_fill_manual(values=c("navy blue", "grey","orange", "maroon", "black")) +
      facet_wrap(. ~ source, scales = "free_x") +
      labs(x = "Haplotype", y = "Relative Frequency", fill = "Haplotype") +
      scale_y_continuous(limits=c(0, 101), breaks=seq(0, 100, by=10))
)



###____save plot ----
# -----------------------------------------------------------------------------#

jpeg(filename = "plot/dhps-haplotypes.jpeg",
     width = 5, height = 5, units = 'in', res = 600,
     # type = "cairo",       #options::"cairo", "Xlib", "quartz"
     # antialias = "none" #options::"default", "none", "gray", "subpixel"
); print(my_plot); dev.off()



## __DHFR ----
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

###____table ----
# --------------------------------------------#

freqSNP_DHFR <- clusters_DHFR %>%
  select(s_Sample, source, starts_with("pos")) %>%
  reframe(
          # collapse alleles per codon per sample
          source,
          across(starts_with("pos"), ~paste(unique(sort(.)), collapse = ","),
                 .names = "{.col}"),
          .by = s_Sample) %>%
  distinct(s_Sample, .keep_all = TRUE) %>%
  pivot_longer(
               cols = starts_with("pos"), # dynamically select columns that start with "pos"
               names_to = "codon",
               values_to = "allele"
               ) %>%
  mutate(codon = str_remove(codon, "pos")) %>%
  summarise(count=n(), .by = c(source, codon, allele)) %>%
  arrange(source, as.numeric(codon)) %>%
  mutate(
         freq = count/sum(count) * 100, .by = c(source, codon),
         freq = round(freq,1),
         allele = str_remove_all(allele, "\\d+")
         ) %>%
  mutate(codon = as.numeric(codon)) %>%
  left_join(data.frame(position = positions_DHFR, wildtype = wt_alleles),
            by = c("codon" = "position")) %>%
  mutate(
         codon_allele = paste0(wildtype, codon, allele),
         variant = case_when(
                             str_detect(allele, "\\,") ~ "mixed",
                             allele == wildtype ~ "wildtype",
                             allele != wildtype ~ "mutant"),
         freq = paste0(freq, " [", count, "]")
         ) %>%
  pivot_wider(
              id_cols = c(source, codon, codon_allele),
              names_from = "variant",
              values_from = "freq") %>%
  summarise(
            codon_allele = paste(codon_allele, collapse=", "), 
            wildtype = first(na.omit(wildtype)),
            mutant = first(na.omit(mutant)),
            mixed = first(na.omit(mixed)),
            .by = c(source, codon)
            )



###____plot ----
# --------------------------------------------#

(
  my_plot <- clusters_DHFR %>%
    select(s_Sample, source, starts_with("pos")) %>%
    reframe( # collapse alleles per codon per sample
            source,
            across(starts_with("pos"), ~paste(unique(sort(.)), collapse = ","),
                   .names = "{.col}"),
            .by = s_Sample) %>%
    distinct(s_Sample, .keep_all = TRUE) %>%
    pivot_longer(
                 cols = starts_with("pos"), # dynamically select columns that start with "pos"
                 names_to = "codon",
                 values_to = "allele"
                 ) %>%
    mutate(codon = str_remove(codon, "pos")) %>%
    summarise(count=n(), .by = c(source, codon, allele)) %>%
    arrange(source, as.numeric(codon)) %>%
    mutate(
           freq = count/sum(count) * 100, .by = c(source, codon),
           freq = round(freq,1),
           allele = str_remove_all(allele, "\\d+")
           ) %>%
    mutate(codon = as.numeric(codon)) %>%
    left_join(data.frame(position = positions_DHFR, wildtype = wt_alleles),
              by = c("codon" = "position")) %>%
    mutate(variant = case_when(
                               str_detect(allele, "\\,") ~ "mixed",
                               allele == wildtype ~ "wildtype",
                               allele != wildtype ~ "mutant"),
           codon_allele = paste(wildtype, codon, allele, sep = "")
           ) %>%
    ggplot(aes(x=as.factor(codon), y=freq, fill=variant)) +
      geom_bar(stat = "identity", position = "stack", width = 0.7) +
      ggprism::theme_prism() +
      scale_fill_manual(values=c("grey","maroon", "black")) +
      facet_wrap(. ~ source, scales = "free_x") +
      labs(x = "Codon Position", y = "Relative Frequency", fill = "Genotype") +
      scale_y_continuous(limits=c(0, 101), breaks=seq(0, 100, by=10))
)



###____save plot ----
# -----------------------------------------------------------------------------#

jpeg(filename = "plot/dhfr-snps.jpeg",
     width = 8, height = 5, units = 'in', res = 600,
     # type = "cairo",       #options::"cairo", "Xlib", "quartz"
     # antialias = "none" #options::"default", "none", "gray", "subpixel"
); print(my_plot); dev.off()



##___compute haplotype frequencies and resistance profiles ----
# -----------------------------------------------------------------------------#

###____table ----
# --------------------------------------------#

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
          dhps_resistance = paste(dhps_resistance, collapse = ","),
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



###____plot ----
# --------------------------------------------#

(
  my_plot <- clusters_DHFR %>%
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
            dhps_resistance = paste(dhps_resistance, collapse = ","),
            .by = c(s_Sample)
            ) %>%
    distinct(s_Sample, .keep_all = TRUE) %>%
    summarise(count=n(), .by = c(source, dhps_resistance)) %>%
    mutate(
           profile = case_when(str_detect(dhps_resistance, ",") ~ "mixed",
                               .default = dhps_resistance),
           freq = count/sum(count) * 100, .by = source
           ) %>%
    ggplot(aes(x=source, y=freq, group=profile, fill=profile)) +
      geom_bar(stat = "identity", position = "stack", width = 0.7) +
      ggprism::theme_prism() +
      scale_fill_manual(values=c("navy blue","grey", "maroon", "black")) +
      facet_wrap(. ~ source, scales = "free_x") +
      labs(x = "Haplotype", y = "Relative Frequency", fill = "Haplotype") +
      scale_y_continuous(limits=c(0, 101), breaks=seq(0, 100, by=10))
)



###____save plot ----
# -----------------------------------------------------------------------------#

jpeg(filename = "plot/dhfr-haplotypes.jpeg",
     width = 5, height = 5, units = 'in', res = 600,
     # type = "cairo",       #options::"cairo", "Xlib", "quartz"
     # antialias = "none" #options::"default", "none", "gray", "subpixel"
); print(my_plot); dev.off()



