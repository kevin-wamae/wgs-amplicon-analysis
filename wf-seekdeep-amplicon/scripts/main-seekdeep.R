# *****************************************************************************#
#     An R Pipeline for Aggregating Allele Frequency Data from SeekDeep
# *****************************************************************************#
# main-seekdeep.R

# Description:
# This script is part of a robust pipeline designed to aggregate and analyze allele
# frequency data from SeekDeep outputs related to the genomic study of malaria 
# parasites. It systematically processes multiple types of genomic data across
# several malaria genes. The pipeline includes data cleaning, importing, 
# processing, and saving allele and haplotype frequencies, enabling detailed 
# genetic analysis.

# Workflow Overview:

# 1. Environment Setup:
#    - Clear the current R environment to ensure a clean slate for data processing
#    - Load necessary R packages from the tidyverse for efficient data manipulation

# 2. Directory and File Management:
#    - Establish directory paths for input and output data
#    - Create an output directory to store generated reports and tables

# 3. Data Importation and Processing:
#    - Import quality control reports for Illumina and Nanopore data 
#    - Process SeekDeep analysis data for targeted genetic clusters
#    - Add sample source information and identify samples without data

# 4. Gene-Specific Analysis:
#    - Perform targeted aggregation for specific malaria genes:
#        * Import cluster data and compute allele frequencies
#        * Extract wild-type sequences and compute haplotype frequencies
#    - Save processed data into CSV files for each gene

# 5. Final Data Output:
#    - Generate allele and haplotype frequency tables for each gene
#    - Categorize by overall data, source-specific, and sample-specific metrics


# Dependencies:
#    - R (>= 4.0)
#    - tidyverse


# clear environment
# =============================================================================#

rm(list = ls())



# load packages
# =============================================================================#

library(tidyverse, quietly = TRUE)



# *****************************************************************************#
# 1. DEFINE DIRECTORY PATHS ----
# *****************************************************************************#


# make sure file-paths terminate with `/`
# -----------------------------------------------------------------------------#
PATH_STUDY = "input/dir1/"  # e.g. "input/ssurvey_western_kenya/"
PATH_RUN = "dir2/"          # e.g. "2025_01_01_illumina_2x300/"
PATH_ANALYSIS = "dir3/"     # e.g. "2024_04_12-01-seekdeep-dhps/"



# create the output directory for generated reports/tables
# -----------------------------------------------------------------------------#
dir.create(paste0(PATH_STUDY, PATH_RUN, PATH_ANALYSIS, "output/"),
           recursive = TRUE,
           showWarnings = FALSE)



# *****************************************************************************#
# 2. IMPORT READ-EXTRACTION QC REPORTS ----

# NOTE: Before proceeding, ensure you run either the Illumina or Nanopore script
# below depending on your sequencing data type.
# *****************************************************************************#

## __i. illumina ----
# -----------------------------------------------------------------------------#

source("scripts/functions/quality-control/qc-extraction-illumina.R")


## __ii. nanopore ----
# -----------------------------------------------------------------------------#

source("scripts/functions/quality-control/qc-extraction-nanopore.R")



## ___save table
# -----------------------------------------------------------------------------#

write_csv(raw_extractionFastq,
          paste0(PATH_STUDY, PATH_RUN, PATH_ANALYSIS, "output/qc-target-read-depth.csv"))



# *****************************************************************************#
# 3. IMPORT SEEKDEEP CLUSTERS DATA ----
# *****************************************************************************#

## __i. import ----
# -----------------------------------------------------------------------------#
source("scripts/functions/sequence-clusters/clusters-aggregate.R")



## __ii. add samples' geographical origin (IF AVAILABLE) ----
# -----------------------------------------------------------------------------#
# NOTE: Before proceeding, please:
# 1. Create a samples/ directory in your PATH_STUDY location
# 2. samples_source.csv with the following columns:
#     - s_Sample: Unique sample identifier matching s_Sample in sequence data
#     - source: Geographical origin or sampling location
#
# Place this file at: PATH_STUDY/samples/samples_source.csv

# If the file does not exist, downstream analysis will proceed without sample source info.
# -----------------------------------------------------------------------------#
source("scripts/functions/study-related/add_samples_source.R")



# *****************************************************************************#
# 4. IDENTIFY SAMPLES WITHOUT DATA ----
# *****************************************************************************#

# aggregate read-extraction data
source("scripts/functions/quality-control/qc-samples-missing-data.R"); head(df_missingDataSamples)



# save table
write_csv(df_missingDataSamples,
          paste0(PATH_STUDY, PATH_RUN, PATH_ANALYSIS, "output/qc-read-ext-profile.csv"))



# remove temporary objects
rm(df_missingDataSamples)



# =============================================================================#
## __A) CLUSTERS - PfAMA1 ----
# =============================================================================#

### ___aggregate/filter clusters ----
# -----------------------------------------------------------------------------#

STRING_TARGET = "^PFAMA1"
STRING_GENOME = "^PF3D7"
source("scripts/functions/sequence-clusters/clusters-filter.R")



### ___coi by source ----
# -----------------------------------------------------------------------------#

(
  df_coi_source <- df_clusters_Target %>%
    distinct(s_Sample, .keep_all = TRUE) %>%  # de-duplicate sample entries
    mutate(s_COI = as.numeric(s_COI)) %>%
    summarise(
      sample_size=n(),       # determine sample size
      min = min(s_COI),      # compute mean, median and max COI
      mean = median(s_COI),
      max = max(s_COI),
      .by = source           # group by sample-origin
    )
)



### save table
write_csv(df_coi_source,
          paste0(PATH_STUDY, PATH_RUN, PATH_ANALYSIS, "output/coi-by-souce.csv"))



### ___coi by sample ----
# -----------------------------------------------------------------------------#

df_coi_sample <- df_clusters_Target %>%
  distinct(s_Sample, .keep_all = TRUE) %>%
  select(source, s_Sample, s_COI); head(df_coi_sample)



### save table
write_csv(df_coi_sample,
          paste0(PATH_STUDY, PATH_RUN, PATH_ANALYSIS, "output/coi-by-sample.csv"))



# remove temporary objects
rm(df_clusters_Target, df_clusters_Segregating, df_coi_source, df_coi_sample)



# =============================================================================#
## __B) CLUSTERS - PfK13 ----
# =============================================================================#

### ___aggregate/filter clusters ----
# -----------------------------------------------------------------------------#

STRING_TARGET = "^PFK13-675"
STRING_GENOME = "^PF3D7"
source("scripts/functions/sequence-clusters/clusters-filter.R")



### ___import the wildtype sequence ----
# -----------------------------------------------------------------------------#

fasta_file <- read_lines("../resources-genome/fasta-protein/PfK13.txt")



### ___define alleles/haplotypes by extracting alleles at specified positions ----
# -----------------------------------------------------------------------------#

# alleles
(
  wt_alleles <- sapply(
    # get unique numeric positions for polymorphic sites
    unique(as.numeric(positions_df$position)),
    # extract a single character (allele) from fasta_file at each position
    function(pos) substr(fasta_file, pos, pos)
  )
)



# haplotypes
(
  wt_haplotype <- paste(wt_alleles, collapse = "")
)



#### ____compute allele frequencies (no weighting) ----
# -----------------------------------------------------------------------------#

source("scripts/functions/var-snps/snpfreq-target.R")



#### ____save table
# -----------------------------------------------------------------------------#

write_csv(df_freqSNP_All,
          paste0(PATH_STUDY, PATH_RUN, PATH_ANALYSIS, "output/freq-allele-k13-all.csv"))
write_csv(df_freqSNP_Source,
          paste0(PATH_STUDY, PATH_RUN, PATH_ANALYSIS, "output/freq-allele-k13-source.csv"))
write_csv(df_freqSNP_Sample,
          paste0(PATH_STUDY, PATH_RUN, PATH_ANALYSIS, "output/freq-allele-k13-sample.csv"))



#### ____compute allele frequencies (with weighting) ----
# -----------------------------------------------------------------------------#

source("scripts/functions/var-snps/snpfreq-target-weighted.R")



#### ____save table
# -----------------------------------------------------------------------------#

write_csv(df_freqSNP_All,
          paste0(PATH_STUDY, PATH_RUN, PATH_ANALYSIS, "output/freq-allele-k13-all-weighted.csv"))
write_csv(df_freqSNP_Source,
          paste0(PATH_STUDY, PATH_RUN, PATH_ANALYSIS, "output/freq-allele-k13-source-weighted.csv"))



#### ____compute haplotype frequencies (no weighting) ----
# -----------------------------------------------------------------------------#

source("scripts/functions/var-haplotype/hapfreq-target-clonal.R")



#### ____save table
# -----------------------------------------------------------------------------#
write_csv(df_freqHap_All,
          paste0(PATH_STUDY, PATH_RUN, PATH_ANALYSIS, "output/freq-haplotype-k13-all.csv"))
write_csv(df_freqHap_Source,
          paste0(PATH_STUDY, PATH_RUN, PATH_ANALYSIS, "output/freq-haplotype-k13-source.csv"))



#### ____compute haplotype frequencies (with weighting) ----
# -----------------------------------------------------------------------------#

source("scripts/functions/var-haplotype/hapfreq-target-clonal-weighted.R")



#### ____save table
# -----------------------------------------------------------------------------#
write_csv(df_freqHap_All,
          paste0(PATH_STUDY, PATH_RUN, PATH_ANALYSIS, "output/freq-haplotype-k13-all-weighted.csv"))
write_csv(df_freqHap_Source,
          paste0(PATH_STUDY, PATH_RUN, PATH_ANALYSIS, "output/freq-haplotype-k13-source-weighted.csv"))



# =============================================================================#
## __C) CLUSTERS - PfMDR1 ----
# =============================================================================#

### ___aggregate/filter clusters ----
# -----------------------------------------------------------------------------#

STRING_TARGET = "^PFMDR1"
STRING_GENOME = "^PF3D7"
source("scripts/functions/sequence-clusters/clusters-filter.R")



### ___import the wildtype sequence ----
# -----------------------------------------------------------------------------#

fasta_file <- read_lines("../resources-genome/fasta-protein/PfMDR1.txt")



### ___define alleles/haplotypes by extracting alleles at specified positions ----
# -----------------------------------------------------------------------------#

# alleles
(
  wt_alleles <- sapply(
    # get unique numeric positions for polymorphic sites
    unique(as.numeric(positions_df$position)),
    # extract a single character (allele) from fasta_file at each position
    function(pos) substr(fasta_file, pos, pos)
  )
)



# haplotypes
(
  wt_haplotype <- paste(wt_alleles, collapse = "")
)



#### ____compute allele frequencies (no weighting) ----
# -----------------------------------------------------------------------------#

source("scripts/functions/var-snps/snpfreq-target.R")



#### ____save table
# -----------------------------------------------------------------------------#

write_csv(df_freqSNP_All,
          paste0(PATH_STUDY, PATH_RUN, PATH_ANALYSIS, "output/freq-allele-mdr1-all.csv"))
write_csv(df_freqSNP_Source,
          paste0(PATH_STUDY, PATH_RUN, PATH_ANALYSIS, "output/freq-allele-mdr1-source.csv"))
write_csv(df_freqSNP_Sample,
          paste0(PATH_STUDY, PATH_RUN, PATH_ANALYSIS, "output/freq-allele-mdr1-sample.csv"))



#### ____compute allele frequencies (with weighting) ----
# -----------------------------------------------------------------------------#

source("scripts/functions/var-snps/snpfreq-target-weighted.R")



#### ____save table
# -----------------------------------------------------------------------------#

write_csv(df_freqSNP_All,
          paste0(PATH_STUDY, PATH_RUN, PATH_ANALYSIS, "output/freq-allele-mdr1-all-weighted.csv"))
write_csv(df_freqSNP_Source,
          paste0(PATH_STUDY, PATH_RUN, PATH_ANALYSIS, "output/freq-allele-mdr1-source-weighted.csv"))



#### ____compute haplotype frequencies (no weighting) ----
# -----------------------------------------------------------------------------#

source("scripts/functions/var-haplotype/hapfreq-target-clonal.R")



#### ____save table
# -----------------------------------------------------------------------------#
write_csv(df_freqHap_All,
          paste0(PATH_STUDY, PATH_RUN, PATH_ANALYSIS, "output/freq-haplotype-mdr1-all.csv"))
write_csv(df_freqHap_Source,
          paste0(PATH_STUDY, PATH_RUN, PATH_ANALYSIS, "output/freq-haplotype-mdr1-source.csv"))



#### ____compute haplotype frequencies (with weighting) ----
# -----------------------------------------------------------------------------#

source("scripts/functions/var-haplotype/hapfreq-target-clonal-weighted.R")



#### ____save table
# -----------------------------------------------------------------------------#
write_csv(df_freqHap_All,
          paste0(PATH_STUDY, PATH_RUN, PATH_ANALYSIS, "output/freq-haplotype-mdr1-all-weighted.csv"))
write_csv(df_freqHap_Source,
          paste0(PATH_STUDY, PATH_RUN, PATH_ANALYSIS, "output/freq-haplotype-mdr1-source-weighted.csv"))



# =============================================================================#
## __D) CLUSTERS - PfDHPS ----
# =============================================================================#

### ___aggregate/filter clusters ----
# -----------------------------------------------------------------------------#

STRING_TARGET = "^PFDHPS"
STRING_GENOME = "^PF3D7"
source("scripts/functions/sequence-clusters/clusters-filter.R")
source("scripts/functions/var-snps/sp-resistance-profile.R")



### ___import the wildtype sequence ----
# -----------------------------------------------------------------------------#

fasta_file <- read_lines("../resources-genome/fasta-protein/PfDHPS.txt")



### ___define alleles/haplotypes by extracting alleles at specified positions ----
# -----------------------------------------------------------------------------#

# alleles
(
  wt_alleles <- sapply(
    # get unique numeric positions for polymorphic sites
    unique(as.numeric(positions_df$position)),
    # extract a single character (allele) from fasta_file at each position
    function(pos) substr(fasta_file, pos, pos)
  )
)



# haplotypes
(
  wt_haplotype <- paste(wt_alleles, collapse = "")
)



#### ____compute allele frequencies (no weighting) ----
# -----------------------------------------------------------------------------#

source("scripts/functions/var-snps/snpfreq-target.R")



#### ____save table
# -----------------------------------------------------------------------------#

write_csv(df_freqSNP_All,
          paste0(PATH_STUDY, PATH_RUN, PATH_ANALYSIS, "output/freq-allele-dhps-all.csv"))
write_csv(df_freqSNP_Source,
          paste0(PATH_STUDY, PATH_RUN, PATH_ANALYSIS, "output/freq-allele-dhps-source.csv"))
write_csv(df_freqSNP_Sample,
          paste0(PATH_STUDY, PATH_RUN, PATH_ANALYSIS, "output/freq-allele-dhps-sample.csv"))



#### ____compute allele frequencies (with weighting) ----
# -----------------------------------------------------------------------------#

source("scripts/functions/var-snps/snpfreq-target-weighted.R")



#### ____save table
# -----------------------------------------------------------------------------#

write_csv(df_freqSNP_All,
          paste0(PATH_STUDY, PATH_RUN, PATH_ANALYSIS, "output/freq-allele-dhps-all-weighted.csv"))
write_csv(df_freqSNP_Source,
          paste0(PATH_STUDY, PATH_RUN, PATH_ANALYSIS, "output/freq-allele-dhps-source-weighted.csv"))



#### ____compute haplotype frequencies (no weighting) ----
# -----------------------------------------------------------------------------#

source("scripts/functions/var-haplotype/hapfreq-dhps-clonal.R")



#### ____save table
# -----------------------------------------------------------------------------#

write_csv(df_freqHap_All,
          paste0(PATH_STUDY, PATH_RUN, PATH_ANALYSIS, "output/freq-haplotype-dhps-all.csv"))
write_csv(df_freqHap_Source,
          paste0(PATH_STUDY, PATH_RUN, PATH_ANALYSIS, "output/freq-haplotype-dhps-source.csv"))



#### ____compute haplotype frequencies (with weighting) ----
# -----------------------------------------------------------------------------#

source("scripts/functions/var-haplotype/hapfreq-target-clonal-weighted.R")



#### ____save table
# -----------------------------------------------------------------------------#
write_csv(df_freqHap_All,
          paste0(PATH_STUDY, PATH_RUN, PATH_ANALYSIS, "output/freq-haplotype-dhps-all-weighted.csv"))
write_csv(df_freqHap_Source,
          paste0(PATH_STUDY, PATH_RUN, PATH_ANALYSIS, "output/freq-haplotype-dhps-source-weighted.csv"))



# =============================================================================#
## __E) CLUSTERS - PfDHFR ----
# =============================================================================#

### ___aggregate/filter clusters ----
# -----------------------------------------------------------------------------#

STRING_TARGET = "^PFDHFR"
STRING_GENOME = "^PF3D7"
source("scripts/functions/sequence-clusters/clusters-filter.R")
source("scripts/functions/var-snps/sp-resistance-profile.R")



### ___import the wildtype sequence ----
# -----------------------------------------------------------------------------#

fasta_file <- read_lines("../resources-genome/fasta-protein/PfDHFR.txt")



### ___define alleles/haplotypes by extracting alleles at specified positions ----
# -----------------------------------------------------------------------------#

# alleles
(
  wt_alleles <- sapply(
    # get unique numeric positions for polymorphic sites
    unique(as.numeric(positions_df$position)),
    # extract a single character (allele) from fasta_file at each position
    function(pos) substr(fasta_file, pos, pos)
  )
)



# haplotypes
(
  wt_haplotype <- paste(wt_alleles, collapse = "")
)



#### ____compute allele frequencies (no weighting) ----
# -----------------------------------------------------------------------------#

source("scripts/functions/var-snps/snpfreq-target.R")



#### ____save table
# -----------------------------------------------------------------------------#

write_csv(df_freqSNP_All,
          paste0(PATH_STUDY, PATH_RUN, PATH_ANALYSIS, "output/freq-allele-dhfr-all.csv"))
write_csv(df_freqSNP_Source,
          paste0(PATH_STUDY, PATH_RUN, PATH_ANALYSIS, "output/freq-allele-dhfr-source.csv"))
write_csv(df_freqSNP_Sample,
          paste0(PATH_STUDY, PATH_RUN, PATH_ANALYSIS, "output/freq-allele-dhfr-sample.csv"))



#### ____compute allele frequencies (with weighting) ----
# -----------------------------------------------------------------------------#

source("scripts/functions/var-snps/snpfreq-target-weighted.R")



#### ____save table
# -----------------------------------------------------------------------------#

write_csv(df_freqSNP_All,
          paste0(PATH_STUDY, PATH_RUN, PATH_ANALYSIS, "output/freq-allele-dhfr-all-weighted.csv"))
write_csv(df_freqSNP_Source,
          paste0(PATH_STUDY, PATH_RUN, PATH_ANALYSIS, "output/freq-allele-dhfr-source-weighted.csv"))



#### ____compute haplotype frequencies (no weighting) ----
# -----------------------------------------------------------------------------#

source("scripts/functions/var-haplotype/hapfreq-dhfr-clonal.R")



#### ____save table
# -----------------------------------------------------------------------------#

write_csv(df_freqHap_All,
          paste0(PATH_STUDY, PATH_RUN, PATH_ANALYSIS, "output/freq-haplotype-dhfr-all.csv"))
write_csv(df_freqHap_Source,
          paste0(PATH_STUDY, PATH_RUN, PATH_ANALYSIS, "output/freq-haplotype-dhfr-source.csv"))



#### ____compute haplotype frequencies (with weighting) ----
# -----------------------------------------------------------------------------#

source("scripts/functions/var-haplotype/hapfreq-target-clonal-weighted.R")



#### ____save table
# -----------------------------------------------------------------------------#
write_csv(df_freqHap_All,
          paste0(PATH_STUDY, PATH_RUN, PATH_ANALYSIS, "output/freq-haplotype-dhfr-all-weighted.csv"))
write_csv(df_freqHap_Source,
          paste0(PATH_STUDY, PATH_RUN, PATH_ANALYSIS, "output/freq-haplotype-dhfr-source-weighted.csv"))


