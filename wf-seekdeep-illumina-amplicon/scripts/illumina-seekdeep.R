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
# 2. import data ----
# *****************************************************************************#

## __a. specify input-dir path ----

# make sure to end file paths with `/`
PATH_STUDY = "wf-seekdeep-illumina-amplicon/input/ssurvey_2022 - western_kenya/"
PATH_RUN = "2024_02_23_ilri_illumina_2x300/"
PATH_DATE = "2024_04_24-01-seekdeep/"



# create the output directory for generated reports
dir.create(paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/"),
           recursive = TRUE,
           showWarnings = FALSE)



## __b. import fastq extraction reports (by target) ----
# =============================================================================#

### ____i. import report ----
# -----------------------------------------------------------------------------#

raw_extProfile <- read_tsv(paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "reports/allExtractionProfile.tab.txt"),
                       show_col_types = FALSE) %>%
  mutate_all(list(~str_replace(., "\\(.+\\)", ""))) %>%       # remove all characters in braces
  mutate_at(.vars = 3:ncol(.), .funs = as.numeric) %>%        # to numeric
  mutate(inputName = str_remove(inputName, "_R1|_S.+_L.+"))   # drop read or miseq-lane suffix



### ____ii. qc summary - by target ----
# -----------------------------------------------------------------------------#

raw_extProfileTarget <- raw_extProfile %>%
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
         totalAllGoodPerc = paste0(totalAllGood, " [", totalAllGoodPerc, "]"),
         ) %>%
  # summary of totals by fastq
  mutate(
         totalByFastq = sum(totalMatching),
         totalByFastqGood = sum(good),
         totalByFastqGoodPerc = round(totalByFastqGood / totalByFastq, 2),
         totalByFastq = base::format(totalByFastq, big.mark = ","),
         totalByFastqGoodPerc = paste0(totalByFastq, " [", totalByFastqGoodPerc, "]"),
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
              id_cols = c(inputName, totalAllGoodPerc, totalByFastqGoodPerc),
              names_from = "name",
              values_from = "totalTargetGoodPerc"
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



### ____iii. add sample origin, if avaialable ----
# -----------------------------------------------------------------------------#

source(paste0(PATH_STUDY, "scripts/add_sample_source.R"))



## __d. identify samples without data ----
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
                                               nsmall = 0))))



### save table
write_csv(raw_sampleNamesProfile, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/read-extraction-profile.csv"))



## __d. extract clusters for each available gene ----
# =============================================================================#

###___PfAMA1 ----
# =============================================================================#

source("wf-seekdeep-illumina-amplicon/scripts/functions/aggregate-clusters-ama1.R")


###___PfK13 ----
# =============================================================================#

source("wf-seekdeep-illumina-amplicon/scripts/functions/aggregate-clusters-k13.R")


###___PfMDR1 ----
# =============================================================================#

source("wf-seekdeep-illumina-amplicon/scripts/functions/aggregate-clusters-mdr1.R")


###___PfDHPS ----
# =============================================================================#

source("wf-seekdeep-illumina-amplicon/scripts/functions/aggregate-clusters-dhps.R")
source("wf-seekdeep-illumina-amplicon/scripts/functions/functions-resistance-profile.R")


###___PfDHFR ----
# =============================================================================#

source("wf-seekdeep-illumina-amplicon/scripts/functions/aggregate-clusters-dhfr.R")
source("wf-seekdeep-illumina-amplicon/scripts/functions/functions-resistance-profile.R")



# *****************************************************************************#
# 3. compute allele frequencies ----
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

##___import the wildtype sequence ----
# -----------------------------------------------------------------------------#

fasta_file <- read_lines("resources-genome/fasta-cds/PfK13.txt")



##___extract the wildtype alleles/haplotypes ----
# -----------------------------------------------------------------------------#

# alleles
(
  wt_alleles <- sapply(positions_K13, function(pos) substr(fasta_file, pos, pos))
)



# haplotypes
(
  wt_haplotype <- paste(wt_alleles, collapse = "")
)



##___compute allele frequencies ----
# -----------------------------------------------------------------------------#

source("wf-seekdeep-illumina-amplicon/scripts/functions/compute-snpfreq-k13.R")



### ____save table ----
# -----------------------------------------------------------------------------#

write_csv(df_freqSNP_K13_All, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-allele-k13-all.csv"))
write_csv(df_freqSNP_K13_Source, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-allele-k13-source.csv"))



## __PfMDR1 ----
# =============================================================================#


##___import the wildtype sequence ----
# -----------------------------------------------------------------------------#

fasta_file <- read_lines("resources-genome/fasta-cds/PfMDR1.txt")



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

source("wf-seekdeep-illumina-amplicon/scripts/functions/compute-snpfreq-mdr1.R")



### ____save table ----
# -----------------------------------------------------------------------------#

write_csv(df_freqSNP_MDR1_All, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-allele-mdr1-all.csv"))
write_csv(df_freqSNP_MDR1_Source, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-allele-mdr1-source.csv"))



##___compute haplotype frequencies ----
# -----------------------------------------------------------------------------#

source("wf-seekdeep-illumina-amplicon/scripts/functions/compute-hapfreq-mdr1.R")



### ____save table ----
# -----------------------------------------------------------------------------#
write_csv(df_freqHap_MDR1_All, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-haplotype-mdr1-all.csv"))
write_csv(df_freqHap_MDR1_Source, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-haplotype-mdr1-source.csv"))



## __PfDHPS ----
# =============================================================================#

##___import the wildtype sequence ----
# -----------------------------------------------------------------------------#

fasta_file <- read_lines("resources-genome/fasta-cds/PfDHPS.txt")



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

source("wf-seekdeep-illumina-amplicon/scripts/functions/compute-snpfreq-dhps.R")



### ____save table ----
# -----------------------------------------------------------------------------#

write_csv(df_freqSNP_DHPS_All, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-allele-dhps-all.csv"))
write_csv(df_freqSNP_DHPS_Source, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-allele-dhps-source.csv"))



##___compute haplotype frequencies and resistance profiles ----
# -----------------------------------------------------------------------------#

source("wf-seekdeep-illumina-amplicon/scripts/functions/compute-hapfreq-dhps.R")


### ____save table ----
# -----------------------------------------------------------------------------#

write_csv(df_freqHap_DHPS_All, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-haplotype-dhps-all.csv"))
write_csv(df_freqHap_DHPS_Source, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-haplotype-dhps-source.csv"))



## __PfDHFR ----
# =============================================================================#


##___import the wildtype sequence ----
# -----------------------------------------------------------------------------#

fasta_file <- read_lines("resources-genome/fasta-cds/PfDHFR.txt")



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

source("wf-seekdeep-illumina-amplicon/scripts/functions/compute-snpfreq-dhfr.R")



### ____save table ----
# -----------------------------------------------------------------------------#

write_csv(df_freqSNP_DHFR_All, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-allele-dhfr-all.csv"))
write_csv(df_freqSNP_DHFR_Source, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-allele-dhfr-source.csv"))



##___compute haplotype frequencies and resistance profiles ----
# -----------------------------------------------------------------------------#

source("wf-seekdeep-illumina-amplicon/scripts/functions/compute-hapfreq-dhfr.R")



### ____save table ----
# -----------------------------------------------------------------------------#

write_csv(df_freqHap_DHFR_All, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-haplotype-dhfr-all.csv"))
write_csv(df_freqHap_DHFR_Source, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-haplotype-dhfr-source.csv"))


