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



# *****************************************************************************#
# 1. define directory paths ----
# *****************************************************************************#


# make sure file-paths terminate with `/`
# -----------------------------------------------------------------------------#
PATH_STUDY = "input/ssurvey_2022 - western_kenya/"
PATH_RUN = "2023_05_25_ilri_illumina_2x300/"
PATH_DATE = "2024_04_12-01-seekdeep-dhfr/"



# create the output directory for generated reports/tables
# -----------------------------------------------------------------------------#
dir.create(paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/"),
           recursive = TRUE,
           showWarnings = FALSE)



# *****************************************************************************#
# 2. import read-extraction qc reports ----
# *****************************************************************************#


## ___i. illumina ----
# -----------------------------------------------------------------------------#

source("../resources-src/functions/aggregate-illumina-extraction-qc.R")



## ___save table
# -----------------------------------------------------------------------------#

write_csv(raw_extProfileTarget, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/qc-target-read-depth.csv"))



## __c. import analysis data ----
# =============================================================================#

### ____i. one file ----
# -----------------------------------------------------------------------------#

raw_selectedClustersInfo <- read_tsv(paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "selectedClustersInfo.tab.txt.gz"),
                     show_col_types = FALSE) %>%
  mutate(source = "None") %>% # define source of samples
  # filter untranslatable and show targets available for analysis
  filter(h_AATyped != "Untranslatable")



# check what markers are available
(
  available_markers <- unique(raw_selectedClustersInfo$p_name)
)



### ____ii. multiple files ----
# -----------------------------------------------------------------------------#

# start by listing files with the specific prefix
file_list <- list.files(path = paste0(PATH_STUDY, PATH_RUN, PATH_DATE),
                    pattern = "^selectedClustersInfo_.*\\.gz$",
                    full.names = TRUE)



# read and bind the data frames
raw_selectedClustersInfo <- file_list %>%
  map_dfr(~ read_tsv(.x, col_types = cols())) %>%  # adjust col_types as needed
  mutate(source = "None") %>% # define source of samples
  # filter untranslatable and show targets available for analysis
  filter(h_AATyped != "Untranslatable"); rm(file_list)



# check what markers are available
(
  available_markers <- unique(raw_selectedClustersInfo$p_name)
)



## __d. add sample origin, if avaialable ----
# =============================================================================#

source(paste0(PATH_STUDY, "scripts/add_sample_source.R"))



## __e. identify samples without data ----
# =============================================================================#

# import sample names
raw_sampleNames <- read_tsv(paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "/info/sampNames.tab.txt"),
                        show_col_types = FALSE,
                        col_names = c("inputName", "s_Sample", "mid"))



# dynamically create columns based on available_markers
for (marker in available_markers) {
  raw_sampleNames[[paste0("gene_", marker)]] <- paste0(toupper(marker), raw_sampleNames$mid)
}



# merge sample names with read-extraction profile
raw_sampleNamesProfile <- raw_sampleNames %>%
  pivot_longer(                                    # transform: wide to long
               cols = starts_with("gene_"),
               names_to = "target",
               values_to = "name"
               ) %>%
  left_join(
            raw_extProfile,                        # merge with raw_extProfile 
            by = c("inputName", "name")
            ) %>% 
  pivot_wider(                                     # transform: long to wide
              id_cols = c(inputName, s_Sample),
              names_from = "target",
              values_from = "good"
              ) %>%
  mutate(across(.cols = starts_with("gene_"),      # add thousands separator and replace_na
                ~ifelse(is.na(.x), "0", format(.x,
                                               big.mark = ",",
                                               decimal.mark = ".",
                                               nsmall = 0)))); rm(marker)



### save table
write_csv(raw_sampleNamesProfile, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/qc-read-ext-profile.csv"))



# =============================================================================#
## ___Clusters - PfAMA1 ----
# =============================================================================#

### ____aggregate clusters ----
# -----------------------------------------------------------------------------#

STRING_TARGET = "^PFAMA1"
STRING_GENOME = "^PF3D7"
source("../resources-src/functions/aggregate-clusters-target.R")



### ____coi by source ----
# -----------------------------------------------------------------------------#

(
  df_coi_source <- df_clusters_Target %>%
    distinct(s_Sample, .keep_all = TRUE) %>%  # de-duplicate sample entries
    summarise(
      sample_size=n(),                # determine sample size
      min = min(s_COI),
      mean = median(s_COI),
      max = max(s_COI),
      .by = source                    # group by sample-origin
    )
)



### save table
write_csv(df_coi_source, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/coi-by-souce.csv"))



### ____coi by sample ----
# -----------------------------------------------------------------------------#

df_coi_sample <- df_clusters_Target %>%
  distinct(s_Sample, .keep_all = TRUE) %>%
  select(source, s_Sample, s_COI)



### save table
write_csv(df_coi_sample, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/coi-by-sample.csv"))



# =============================================================================#
## ___Clusters - PfK13 ----
# =============================================================================#

### ____aggregate clusters ----
# -----------------------------------------------------------------------------#

STRING_TARGET = "^PFK13-469"
STRING_GENOME = "^PF3D7"
source("../resources-src/functions/aggregate-clusters-target.R")



### ____import the wildtype sequence ----
# -----------------------------------------------------------------------------#

fasta_file <- read_lines("../resources-genome/fasta-cds/PfK13.txt")



### ____extract the wildtype alleles/haplotypes ----
# -----------------------------------------------------------------------------#

# alleles
(
  wt_alleles <- sapply(positions_Target, function(pos) substr(fasta_file, pos, pos))
)



# haplotypes
(
  wt_haplotype <- paste(wt_alleles, collapse = "")
)



### ____compute allele frequencies ----
# -----------------------------------------------------------------------------#

source("../resources-src/functions/compute-snpfreq-target.R")



### ____save table
# -----------------------------------------------------------------------------#

write_csv(df_freqSNP_All, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-allele-k13-all.csv"))
write_csv(df_freqSNP_Source, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-allele-k13-source.csv"))
write_csv(df_freqSNP_Sample, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-allele-k13-sample.csv"))




# =============================================================================#
## ___Clusters - PfMDR1 ----
# =============================================================================#

### ____aggregate clusters ----
# -----------------------------------------------------------------------------#

STRING_TARGET = "^PFMDR1"
STRING_GENOME = "^PF3D7"
source("../resources-src/functions/aggregate-clusters-target.R")



### ____import the wildtype sequence ----
# -----------------------------------------------------------------------------#

fasta_file <- read_lines("../resources-genome/fasta-cds/PfMDR1.txt")



### ____extract the wildtype alleles/haplotypes ----
# -----------------------------------------------------------------------------#

# alleles
(
  wt_alleles <- sapply(positions_Target, function(pos) substr(fasta_file, pos, pos))
)



# haplotypes
(
  wt_haplotype <- paste(wt_alleles, collapse = "")
)



### ____compute allele frequencies ----
# -----------------------------------------------------------------------------#

source("../resources-src/functions/compute-snpfreq-target.R")



### ____save table
# -----------------------------------------------------------------------------#

write_csv(df_freqSNP_All, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-allele-mdr1-all.csv"))
write_csv(df_freqSNP_Source, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-allele-mdr1-source.csv"))
write_csv(df_freqSNP_Sample, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-allele-mdr1-sample.csv"))



### ____compute haplotype frequencies ----
# -----------------------------------------------------------------------------#

source("../resources-src/functions/compute-hapfreq-target-clonal.R")



### ____save table ----
# -----------------------------------------------------------------------------#
write_csv(df_freqHap_MDR1_All, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-haplotype-mdr1-all.csv"))
write_csv(df_freqHap_MDR1_Source, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-haplotype-mdr1-source.csv"))


# =============================================================================#
## ___Clusters - PfDHPS ----
# =============================================================================#

### ____aggregate clusters ----
# -----------------------------------------------------------------------------#

STRING_TARGET = "^PFDHPS"
STRING_GENOME = "^PF3D7"
source("../resources-src/functions/aggregate-clusters-target.R")
source("../resources-src/functions/functions-resistance-profile.R")



### ____import the wildtype sequence ----
# -----------------------------------------------------------------------------#

fasta_file <- read_lines("../resources-genome/fasta-cds/PfDHPS.txt")



### ____extract the wildtype alleles/haplotypes ----
# -----------------------------------------------------------------------------#

# alleles
(
  wt_alleles <- sapply(positions_Target, function(pos) substr(fasta_file, pos, pos))
)



# haplotypes
(
  wt_haplotype <- paste(wt_alleles, collapse = "")
)



### ____compute allele frequencies ----
# -----------------------------------------------------------------------------#

source("../resources-src/functions/compute-snpfreq-target.R")



### ____save table
# -----------------------------------------------------------------------------#

write_csv(df_freqSNP_All, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-allele-dhps-all.csv"))
write_csv(df_freqSNP_Source, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-allele-dhps-source.csv"))
write_csv(df_freqSNP_Sample, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-allele-dhps-sample.csv"))



### ____compute haplotype frequencies ----
# -----------------------------------------------------------------------------#

source("../resources-src/functions/compute-hapfreq-dhps-clonal.R")



### ____save table
# -----------------------------------------------------------------------------#

write_csv(df_freqHap_All, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-haplotype-dhps-all.csv"))
write_csv(df_freqHap_Source, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-haplotype-dhps-source.csv"))



### ____save table
# -----------------------------------------------------------------------------#

write_csv(df_freqHap_DHPS_All, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-haplotype-dhps-all.csv"))
write_csv(df_freqHap_DHPS_Source, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-haplotype-dhps-source.csv"))



# =============================================================================#
## ___Clusters - PfDHFR ----
# =============================================================================#

### ____aggregate clusters ----
# -----------------------------------------------------------------------------#

STRING_TARGET = "^PFDHFR"
STRING_GENOME = "^PF3D7"
source("../resources-src/functions/aggregate-clusters-target.R")
source("../resources-src/functions/functions-resistance-profile.R")



### ____import the wildtype sequence ----
# -----------------------------------------------------------------------------#

fasta_file <- read_lines("../resources-genome/fasta-cds/PfDHFR.txt")



### ____extract the wildtype alleles/haplotypes ----
# -----------------------------------------------------------------------------#

# alleles
(
  wt_alleles <- sapply(positions_Target, function(pos) substr(fasta_file, pos, pos))
)



# haplotypes
(
  wt_haplotype <- paste(wt_alleles, collapse = "")
)



### ____compute allele frequencies ----
# -----------------------------------------------------------------------------#

source("../resources-src/functions/compute-snpfreq-target.R")



### ____save table
# -----------------------------------------------------------------------------#

write_csv(df_freqSNP_All, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-allele-dhfr-all.csv"))
write_csv(df_freqSNP_Source, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-allele-dhfr-source.csv"))
write_csv(df_freqSNP_Sample, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-allele-dhfr-sample.csv"))



### ____compute haplotype frequencies ----
# -----------------------------------------------------------------------------#

source("../resources-src/functions/compute-hapfreq-dhfr-clonal.R")



### ____save table
# -----------------------------------------------------------------------------#

write_csv(df_freqHap_All, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-haplotype-dhfr-all.csv"))
write_csv(df_freqHap_Source, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-haplotype-dhfr-source.csv"))
