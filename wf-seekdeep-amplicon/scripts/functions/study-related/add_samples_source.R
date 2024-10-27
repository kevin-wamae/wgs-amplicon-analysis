# An R script for determining the source of samples based on sample ids
# =============================================================================#


# Import the dataset containing sample geographical origins
study_county <- read_csv(paste0(PATH_STUDY, "samples/", "sample_sources.csv"),
                         col_types = cols(.default = col_character()),  # specify all columns as character type
                         show_col_types = FALSE                         # suppress column type display in console
                         ) %>%
  # select only the relevant columns: sample origin ('source') and sample identifier ('s_Sample')
  select(source, s_Sample)



# replace source in raw_selectedClustersInfo with actual sample-origins
raw_selectedClustersInfo <- raw_selectedClustersInfo %>%
  select(-source) %>%
  left_join(study_county, by = "s_Sample", relationship = "many-to-many")



# remove temporary objects
rm(study_county)

