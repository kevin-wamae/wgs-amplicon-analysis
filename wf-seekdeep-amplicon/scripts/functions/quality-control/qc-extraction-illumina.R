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

raw_extractionTarget <- read_tsv(paste0(PATH_STUDY, PATH_RUN, PATH_ANALYSIS, "reports/allExtractionProfile.tab.txt"),
                       show_col_types = FALSE) %>%
  mutate_all(list(~str_replace(., "\\(.+\\)", ""))) %>%       # remove all characters in braces
  mutate_at(.vars = 3:ncol(.), .funs = as.numeric) %>%        # to numeric
  mutate(inputName = str_remove(inputName, "_R1|_S.+_L.+"))   # drop read or miseq-lane suffix



# ii. extraction reports by FASTQ and gene
# -----------------------------------------------------------------------------#

raw_extractionFastq <- raw_extractionTarget %>%
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



##___print a message in the console ----
# -----------------------------------------------------------------------------#

# Using yellow for the border
cat("\033[1m\033[33m", "\n##############################################################", "\033[0m")

# Using magenta for the "NOTE TO USER" heading
cat("\033[1m\033[35m", "\nNOTE TO USER:", "\033[0m")

# Using green for the table descriptions
cat("\033[1m\033[32m", "\n1. raw_extractionFastq  - This table reports stats on sequence-reads extraction per FASTQ file", "\033[0m")
cat("\033[1m\033[32m", "\n2. raw_extractionTarget - This table reports stats on sequence-reads extraction per target/marker", "\033[0m")

# Yellow for the border again
cat("\033[1m\033[33m", "\n##############################################################\n", "\033[0m")
