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
  filter(h_AATyped != "Untranslatable"); rm(file_list)



# check what markers are available
(
  available_markers <- unique(raw_selectedClustersInfo$p_name)
)

