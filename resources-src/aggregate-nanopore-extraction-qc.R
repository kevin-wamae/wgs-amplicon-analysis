# DESCRIPTION OF R SCRIPT
# =============================================================================#
# This R script is specifically tailored to import and process FASTQ extraction reports
# for SeekDeep data within a defined study, run, and date. The script undertakes several
# crucial operations to ensure the data is prepared for downstream analysis:
#
# 1. Import and Process General Extraction Reports:
#    Reads a tabulated file containing overall extraction statistics for FASTQ files.
#    This includes converting all numeric data from string format to numeric to ensure
#    accurate computations and data handling. The operation focuses on fields from
#    "totalReadsProcessed" to the last column of the dataset, ensuring all are treated as
#    numeric values.
#
# 2. Import and Enhance Target-Specific Extraction Reports:
#    Processes another set of extraction reports that detail the extraction statistics
#    by FASTQ files and specific genetic targets. The script:
#      - Calculates the total counts per sample to use as a denominator in proportion calculations.
#      - Computes the frequency of each target as a proportion of the total counts, formatted
#        as a percentage alongside the absolute count.
#      - Reshapes the data from long to wide format, where each target becomes a separate column
#        with combined frequency and count information for easier comparison and analysis.
#
# This script utilizes the 'readr', 'dplyr', and 'tidyr' packages from the tidyverse to handle
# data importation, transformation, and reshaping tasks efficiently. It is essential for researchers
# or analysts looking to preprocess extraction data for quality checks or preliminary analysis before
# deeper genetic evaluations.
#
# Dependencies: readr, dplyr, tidyr
# =============================================================================#



# i. extraction reports by FASTQ
# -----------------------------------------------------------------------------#
raw_extProfile <- read_tsv(paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "reports/allExtractionStats.tab.txt"),
                     show_col_types = FALSE) %>%
  mutate_at(.vars = 2:ncol(.), .funs = as.numeric) # to numeric, from "totalReadsProcessed" to "passedFrac"



# ii. extraction reports by FASTQ and gene
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
