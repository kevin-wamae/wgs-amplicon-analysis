# *****************************************************************************#
#     An R Pipeline for Aggregating Allele Frequency Data from SeekDeep
# *****************************************************************************#
# This script is part of a robust pipeline designed to aggregate and analyze allele frequency data
# from SeekDeep outputs, related to the genomic study of malaria parasites. It systematically processes
# multiple types of genomic data across several malaria genes. The pipeline includes data cleaning,
# importing, processing, and saving allele and haplotype frequencies, enabling detailed genetic analysis.

# Workflow Overview:

# 1. Environment Setup:
#    - Clears the current R environment to ensure a clean slate for data processing.
#    - Loads necessary R packages from the tidyverse and data.table for efficient data manipulation.

# 2. Directory and File Management:
#    - Establishes directory paths for input and output data, ensuring all file paths terminate with `/`.
#    - Creates an output directory to store generated reports and tables.

# 3. Data Importation and Processing:
#    - Imports quality control reports for Illumina and Nanopore data to assess read extraction efficiency.
#    - Processes SeekDeep analysis data for targeted genetic clusters, adding sample source information
#      and identifying samples without data.

# 4. Gene-Specific Analysis:
#    - Performs targeted aggregation of cluster data for specific malaria genes (e.g., PfAMA1, PfK13, PfMDR1,
#      PfDHPS, PfDHFR) including:
#         * Importing raw cluster data and computing allele frequencies.
#         * Aggregating clusters to extract wildtype sequences and compute haplotype frequencies.
#    - Saves the processed data into CSV files for each gene, facilitating easy access and further analysis.

# 5. Final Data Output:
#    - Outputs include allele frequency tables and haplotype frequency tables for each gene studied,
#      categorized by overall data, source-specific, and sample-specific breakdowns.

# This pipeline is essential for researchers studying genetic variations and drug resistance in malaria parasites,
# providing a structured approach to handle and analyze high-throughput sequencing data.

# Dependencies: R (>= 4.0), tidyverse, data.table
# *****************************************************************************#



# clear environment
# =============================================================================#

rm(list = ls())



# load packages
# =============================================================================#

library(tidyverse, quietly = TRUE)



# *****************************************************************************#
# 1. define directory paths ----
# *****************************************************************************#


# make sure file-paths terminate with `/`
# -----------------------------------------------------------------------------#
PATH_STUDY = "input/ssurvey_2022_western_kenya/"
PATH_RUN = "2024_02_23_ilri_illumina_2x300/"
PATH_DATE = "2024_04_24-01-seekdeep/"



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

source("../resources-src/aggregate-illumina-extraction-qc.R")


## ___ii. nanopore ----
# -----------------------------------------------------------------------------#

source("../resources-src/aggregate-nanopore-extraction-qc.R")



## ___save table
# -----------------------------------------------------------------------------#

write_csv(raw_extProfileTarget, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/qc-target-read-depth.csv"))



# *****************************************************************************#
# 3. import SeekDeep analysis data ----
# *****************************************************************************#

# import

source("../resources-src/import-target-clusters.R")


# add sample origin, if available

source(paste0(PATH_STUDY, "scripts/add_sample_source.R"))



# *****************************************************************************#
# 4. identify samples without data ----
# *****************************************************************************#

# aggregate read-extraction data
source("../resources-src/aggregate-illumina-extraction-qc-target.R")



# save table
write_csv(raw_sampleNamesProfile, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/qc-read-ext-profile.csv"))



# =============================================================================#
## ___Clusters - PfAMA1 ----
# =============================================================================#

### ____aggregate clusters ----
# -----------------------------------------------------------------------------#

STRING_TARGET = "^PFAMA1"
STRING_GENOME = "^PF3D7"
source("../resources-src/aggregate-clusters-target.R")



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
source("../resources-src/aggregate-clusters-target.R")



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

source("../resources-src/compute-snpfreq-target.R")



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
source("../resources-src/aggregate-clusters-target.R")



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

source("../resources-src/compute-snpfreq-target.R")



### ____save table
# -----------------------------------------------------------------------------#

write_csv(df_freqSNP_All, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-allele-mdr1-all.csv"))
write_csv(df_freqSNP_Source, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-allele-mdr1-source.csv"))
write_csv(df_freqSNP_Sample, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-allele-mdr1-sample.csv"))



### ____compute haplotype frequencies ----
# -----------------------------------------------------------------------------#

source("../resources-src/compute-hapfreq-target-clonal.R")



### ____save table ----
# -----------------------------------------------------------------------------#
write_csv(df_freqHap_All, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-haplotype-mdr1-all.csv"))
write_csv(df_freqHap_Source, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-haplotype-mdr1-source.csv"))


# =============================================================================#
## ___Clusters - PfDHPS ----
# =============================================================================#

### ____aggregate clusters ----
# -----------------------------------------------------------------------------#

STRING_TARGET = "^PFDHPS"
STRING_GENOME = "^PF3D7"
source("../resources-src/aggregate-clusters-target.R")
source("../resources-src/functions-resistance-profile.R")



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

source("../resources-src/compute-snpfreq-target.R")



### ____save table
# -----------------------------------------------------------------------------#

write_csv(df_freqSNP_All, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-allele-dhps-all.csv"))
write_csv(df_freqSNP_Source, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-allele-dhps-source.csv"))
write_csv(df_freqSNP_Sample, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-allele-dhps-sample.csv"))



### ____compute haplotype frequencies ----
# -----------------------------------------------------------------------------#

source("../resources-src/compute-hapfreq-dhps-clonal.R")



### ____save table
# -----------------------------------------------------------------------------#

write_csv(df_freqHap_All, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-haplotype-dhps-all.csv"))
write_csv(df_freqHap_Source, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-haplotype-dhps-source.csv"))



# =============================================================================#
## ___Clusters - PfDHFR ----
# =============================================================================#

### ____aggregate clusters ----
# -----------------------------------------------------------------------------#

STRING_TARGET = "^PFDHFR"
STRING_GENOME = "^PF3D7"
source("../resources-src/aggregate-clusters-target.R")
source("../resources-src/functions-resistance-profile.R")



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

source("../resources-src/compute-snpfreq-target.R")



### ____save table
# -----------------------------------------------------------------------------#

write_csv(df_freqSNP_All, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-allele-dhfr-all.csv"))
write_csv(df_freqSNP_Source, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-allele-dhfr-source.csv"))
write_csv(df_freqSNP_Sample, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-allele-dhfr-sample.csv"))



### ____compute haplotype frequencies ----
# -----------------------------------------------------------------------------#

source("../resources-src/compute-hapfreq-dhfr-clonal.R")



### ____save table
# -----------------------------------------------------------------------------#

write_csv(df_freqHap_All, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-haplotype-dhfr-all.csv"))
write_csv(df_freqHap_Source, paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/freq-haplotype-dhfr-source.csv"))
