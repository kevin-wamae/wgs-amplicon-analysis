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



# =============================================================================#
# load packages
library(data.table, quietly = TRUE)
library(tidyverse, quietly = TRUE)



# *****************************************************************************#
# 1. import data ----
# *****************************************************************************#

## __a. specify input-dir path ----

# make sure to end file paths with `/`
PATH_STUDY = "wf-seekdeep-nanopore-amplicon/input/turkana_embatalk/"
PATH_RUN = "2024_04_12_nanopore_r10.4.1/"
PATH_DATE = "2024_04_24-01-seekdeep/"



# create the output directory for generated reports
dir.create(paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "output/"),
           recursive = TRUE,
           showWarnings = FALSE)



## __a. import FASTQ extraction reports ----
# =============================================================================#


# extraction reports by FASTQ
# -----------------------------------------------------------------------------#
raw_extProfile <- read_tsv(paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "reports/allExtractionStats.tab.txt"),
                     show_col_types = FALSE) %>%
  mutate_at(.vars = 2:ncol(.), .funs = as.numeric) # to numeric, from "totalReadsProcessed" to "passedFrac"



# extraction reports by FASTQ and gene
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



## __d. add sample origin, if avaialable ----
# =============================================================================#

source(paste0(PATH_STUDY, "scripts/add_sample_source.R"))



## __f. extract clusters for each available gene ----
# =============================================================================#

###___PfAMA1 ----
# =============================================================================#

source("src/functions/aggregate-clusters-ama1.R")


###___PfK13-469 ----
# =============================================================================#

source("src/functions/aggregate-clusters-k13-469.R")


###___PfK13-675 ----
# =============================================================================#

source("src/functions/aggregate-clusters-k13-675.R")


###___PfMDR1 ----
# =============================================================================#

source("src/functions/aggregate-clusters-mdr1.R")


###___PfDHPS ----
# =============================================================================#

source("src/functions/aggregate-clusters-dhps.R")
source("src/functions/functions-resistance-profile.R")


###___PfDHFR ----
# =============================================================================#

source("src/functions/aggregate-clusters-dhfr.R")
source("src/functions/functions-resistance-profile.R")
