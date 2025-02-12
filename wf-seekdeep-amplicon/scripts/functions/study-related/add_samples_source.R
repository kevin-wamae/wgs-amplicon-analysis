# An R script for mapping sample identifiers to their geographical origins
# =============================================================================#
# This script enhances the dataset `raw_selectedClustersInfo` by adding a new
# column that indicates the geographical origin ('source') of each sample based
# on their sample IDs ('s_Sample'). It reads the sample source information from
# a CSV file and merges it with the main dataset.

# Workflow Overview:

# 1. Define the file path to the `samples_source.csv` file, which should be
#    located in the `samples` directory within the study path (`PATH_STUDY`).

# 2. Check if the `samples_source.csv` file exists:
#    - If it exists:
#        a. Import the dataset containing sample geographical origins.
#        b. Read all columns as character type to maintain data consistency.
#        c. Select only the relevant columns: 'source' (geographical origin)
#           and 's_Sample' (sample ID).
#        d. Merge the sample source data with `raw_selectedClustersInfo` using
#           a left join to add the 'source' information to each sample.
#        e. Remove temporary objects to clean up the environment.
#    - If it does not exist:
#        a. Display a warning message indicating that the file is missing.
#        b. Inform the user that sample origins were not added and provide the
#           expected file location and required columns.

# 3. Remove the temporary `file_sample_sources` object to finalize the script.

# Notes:
# - Ensure that the `samples_source.csv` file contains the following columns:
#     * 'source': the geographical origin of each sample.
#     * 's_Sample': the sample identifiers matching those in `raw_selectedClustersInfo`.
# - Place the `samples_source.csv` file in the directory:
#   `PATH_STUDY/samples/` for the script to locate it correctly.
# - Missing the sample source information may lead to incomplete downstream
#   analysis, as geographical origin is a key variable.

# Dependencies:
# - Packages: `dplyr`, `readr`
# =============================================================================#

# define file path for sample sources
file_sample_sources <- paste0(PATH_STUDY, "samples/", "samples_source.csv")

# check if the samples_source.csv file exists
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
    "Ensure the file is located in the directory: '", PATH_STUDY, "samples/'. ",
    "The file should contain two columns: 'source' for sample source and 's_Sample' for sample names."
  )
}

# remove temporary objects
rm(file_sample_sources)
