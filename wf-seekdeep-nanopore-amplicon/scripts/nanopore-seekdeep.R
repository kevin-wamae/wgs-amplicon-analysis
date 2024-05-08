# *****************************************************************************#
#     An R pipeline for aggregating allele frequency data from SeekDeep
# *****************************************************************************#

# clear environment
# =============================================================================#

rm(list = ls())



# =============================================================================#
# load packages
library(data.table, quietly = TRUE)
library(tidyverse, quietly = TRUE)



# *****************************************************************************#
# 1. import data ----
# *****************************************************************************#

## __a. specify input-dir path ----

# make sure to end file paths with `/`
PATH_STUDY = "input/tes_busia-2018-2019/"
PATH_RUN = "2024_04_12_nanopore_r10.4.1/"
PATH_DATE = "2024_04_24-01-seekdeep/"



# create the output directory for generated reports
dir.create(paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/"),
           recursive = TRUE,
           showWarnings = FALSE)



## __a. import FASTQ extraction reports ----
# =============================================================================#


# extraction reports by FASTQ
# -----------------------------------------------------------------------------#
raw_extProfile <- read_tsv(paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "reports/allExtractionStats.tab.txt"),
                     show_col_types = FALSE) %>%
  mutate_at(.vars = 2:ncol(.), .funs = as.numeric) # to numeric, from "totalReadsProcessed" to "passedFrac"



# extraction reports by FASTQ and gene
# -----------------------------------------------------------------------------#
raw_extProfileTarget <- read_tsv(paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "reports/allExtractionProfile.tab.txt"),
                       show_col_types = FALSE) %>%
  mutate(
         total = sum(count), .by = sample,
         freq = as.character(round(count / total, 2)),
         freq = paste0(count, " [", freq, "]")
         ) %>%
  pivot_wider(id_cols = sample,
              names_from = target,
              values_from = freq,
              names_sep = "_"
              )



### ____save table ----
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
  filter(h_AATyped != "Untranslatable")



# check what markers are available
(
  available_markers <- unique(raw_selectedClustersInfo$p_name)
)



## __d. add sample origin, if avaialable ----
# =============================================================================#

source(paste0(PATH_STUDY, "scripts/add_sample_source.R"))



## __f. extract clusters for each available gene ----
# =============================================================================#

###___PfAMA1 ----
# =============================================================================#

source("../resources-src/functions/aggregate-clusters-ama1.R")


###___PfK13 ----
# =============================================================================#

source("../resources-src/functions/aggregate-clusters-k13.R")


###___PfMDR1 ----
# =============================================================================#

source("../resources-src/functions/aggregate-clusters-mdr1.R")


###___PfDHPS ----
# =============================================================================#

source("../resources-src/functions/aggregate-clusters-dhps.R")
source("../resources-src/functions/functions-resistance-profile.R")


###___PfDHFR ----
# =============================================================================#

source("../resources-src/functions/aggregate-clusters-dhfr.R")
source("../resources-src/functions/functions-resistance-profile.R")



# *****************************************************************************#
# 2. compute allele frequencies ----
# *****************************************************************************#


## __PfAMA1 ----
# =============================================================================#


### ____i. coi by source ----
# -----------------------------------------------------------------------------#

(
  df_coi_source <- df_clusters_AMA1 %>%
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



### ____ii. coi by sample ----
# -----------------------------------------------------------------------------#

df_coi_sample <- df_clusters_AMA1 %>%
  distinct(s_Sample, .keep_all = TRUE) %>%
  select(source, s_Sample, s_COI)



### save table
write_csv(df_coi_sample, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/coi-by-sample.csv"))



## __PfK13 ----
# =============================================================================#

###___import the wildtype sequence ----
# -----------------------------------------------------------------------------#

fasta_file <- read_lines("../resources-genome/fasta-cds/PfK13.txt")



###___extract the wildtype alleles/haplotypes ----
# -----------------------------------------------------------------------------#

# alleles
(
  wt_alleles <- sapply(positions_K13, function(pos) substr(fasta_file, pos, pos))
)



# haplotypes
(
  wt_haplotype <- paste(wt_alleles, collapse = "")
)



###___compute allele frequencies ----
# -----------------------------------------------------------------------------#

source("../resources-src/functions/compute-snpfreq-k13.R")



### ____save table ----
# -----------------------------------------------------------------------------#

write_csv(df_freqSNP_K13_All, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-allele-k13-all.csv"))
write_csv(df_freqSNP_K13_Source, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-allele-k13-source.csv"))
write_csv(df_freqSNP_K13_Sample, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-allele-k13-sample.csv"))



## __PfMDR1 ----
# =============================================================================#


##___import the wildtype sequence ----
# -----------------------------------------------------------------------------#

fasta_file <- read_lines("../resources-genome/fasta-cds/PfMDR1.txt")



##___extract the wildtype alleles/haplotypes ----
# -----------------------------------------------------------------------------#

# alleles
(
  wt_alleles <- sapply(positions_MDR1, function(pos) substr(fasta_file, pos, pos))
)



# haplotypes
(
  wt_haplotype <- paste(wt_alleles, collapse = "")
)


##___compute allele frequencies ----
# -----------------------------------------------------------------------------#

source("../resources-src/functions/compute-snpfreq-mdr1.R")



### ____save table ----
# -----------------------------------------------------------------------------#

write_csv(df_freqSNP_MDR1_All, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-allele-mdr1-all.csv"))
write_csv(df_freqSNP_MDR1_Source, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-allele-mdr1-source.csv"))



##___compute haplotype frequencies ----
# -----------------------------------------------------------------------------#

source("../resources-src/functions/compute-hapfreq-mdr1-clonal.R")



### ____save table ----
# -----------------------------------------------------------------------------#
write_csv(df_freqHap_MDR1_All, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-haplotype-mdr1-all.csv"))
write_csv(df_freqHap_MDR1_Source, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-haplotype-mdr1-source.csv"))



## __PfDHPS ----
# =============================================================================#

##___import the wildtype sequence ----
# -----------------------------------------------------------------------------#

fasta_file <- read_lines("../resources-genome/fasta-cds/PfDHPS.txt")



##___extract the wildtype alleles/haplotypes ----
# -----------------------------------------------------------------------------#

# alleles
(
  wt_alleles <- sapply(positions_DHPS, function(pos) substr(fasta_file, pos, pos))
)



# haplotypes
(
  wt_haplotype <- paste(wt_alleles, collapse = "")
)



##___compute allele frequencies ----
# -----------------------------------------------------------------------------#

source("../resources-src/functions/compute-snpfreq-dhps.R")



### ____save table ----
# -----------------------------------------------------------------------------#

write_csv(df_freqSNP_DHPS_All, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-allele-dhps-all.csv"))
write_csv(df_freqSNP_DHPS_Source, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-allele-dhps-source.csv"))



##___compute haplotype frequencies and resistance profiles ----
# -----------------------------------------------------------------------------#

source("../resources-src/functions/compute-hapfreq-dhps.R")


### ____save table ----
# -----------------------------------------------------------------------------#

write_csv(df_freqHap_DHPS_All, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-haplotype-dhps-all.csv"))
write_csv(df_freqHap_DHPS_Source, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-haplotype-dhps-source.csv"))



## __PfDHFR ----
# =============================================================================#


##___import the wildtype sequence ----
# -----------------------------------------------------------------------------#

fasta_file <- read_lines("../resources-genome/fasta-cds/PfDHFR.txt")



##___extract the wildtype alleles/haplotypes ----
# -----------------------------------------------------------------------------#


# alleles
(
  wt_alleles <- sapply(positions_DHFR, function(pos) substr(fasta_file, pos, pos))
)



# haplotypes
(
  wt_haplotype <- paste(wt_alleles, collapse = "")
)



##___compute allele frequencies ----
# -----------------------------------------------------------------------------#

source("../resources-src/functions/compute-snpfreq-dhfr.R")



### ____save table ----
# -----------------------------------------------------------------------------#

write_csv(df_freqSNP_DHFR_All, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-allele-dhfr-all.csv"))
write_csv(df_freqSNP_DHFR_Source, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-allele-dhfr-source.csv"))



##___compute haplotype frequencies and resistance profiles ----
# -----------------------------------------------------------------------------#

source("../resources-src/functions/compute-hapfreq-dhfr.R")



### ____save table ----
# -----------------------------------------------------------------------------#

write_csv(df_freqHap_DHFR_All, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-haplotype-dhfr-all.csv"))
write_csv(df_freqHap_DHFR_Source, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-haplotype-dhfr-source.csv"))



