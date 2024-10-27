# An R script for determining the source of samples based on sample ids
# =============================================================================#

# Define file path for sample sources
file_sample_sources <- paste0(PATH_STUDY, "samples/", "samples_source.csv")

# Check if the samples_source.csv file exists
if (file.exists(file_sample_sources)) {
  
  # import the dataset containing sample geographical origins
  samples_source <- read_csv(
    file_sample_sources,
    col_types = cols(.default = col_character()),  # specify all columns as character type
    show_col_types = FALSE                         # suppress column type display in console
  ) %>%
    # Select only the relevant columns: sample origin ('source') and sample identifier ('s_Sample')
    select(source, s_Sample)
  
  # replace source in raw_selectedClustersInfo with actual sample origins
  raw_selectedClustersInfo <- raw_selectedClustersInfo %>%
    select(-source) %>%
    left_join(samples_source, by = "s_Sample", relationship = "many-to-many")
  
  # remove temporary objects
  rm(samples_source)
  
} else {
  # display a warning message in the console if the file is missing
  warning(
    "The file 'sample_sources.csv' is missing. Sample origins were not added. ",
    "Ensure the file is located in the directory: '", PATH_STUDY, "/samples/'. ",
    "The file should contain two columns: 'source' for sample source and 's_Sample' for sample names."
  )
}
