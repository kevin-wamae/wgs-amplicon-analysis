# DESCRIPTION OF R SCRIPT
# =============================================================================#
# This R script is designed for identify samples without read-data by target.

# It performs several key functions:

# 1. Importing Sample Names:
#    Reads sample names and related identifiers from a tab-separated values file,
#    setting custom column names to ensure data is appropriately labeled for
#    downstream processing.
#
# 2. Dynamic Column Creation:
#    For each available genetic marker, a new column is dynamically created. This
#    column combines the marker name with a unique identifier from the sample data
#    to create a new identifier that incorporates genetic information.
#
# 3. Data Merging and Transformation:
#    The script merges sample names with their corresponding read-extraction profiles,
#    transforming the data structure from wide to long format and back to wide again.
#    This process facilitates the analysis by aligning sample information with
#    extraction data.
#
# 4. Data Formatting:
#    Applies formatting to numerical data to improve readability, including adding
#    thousands separators and handling missing values by replacing them with a
#    default value.
#
# The script uses functions from the 'tidyverse' package, particularly 'dplyr' and 
# 'tidyr', to manipulate and format the data efficiently. This script is crucial for
# researchers looking to integrate and prepare genomic data for detailed analysis.
#
# Dependencies: readr, dplyr, tidyr
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