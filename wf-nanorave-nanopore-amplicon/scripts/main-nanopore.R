# This R script processes TSV files for genomic data analysis. 

# It starts by defining a list of directory paths and corresponding algorithms.
# Each directory contains TSV files associated with a specific genomic sequencing algorithm.

# The script includes a custom function 'read_and_label' which:
#   1. Lists all TSV files in a given directory path.
#   2. Reads each TSV file into a data frame.
#   3. Adds the basename of each file (the filename without path) to its corresponding data frame.
#   4. Extracts the barcode from the filename (e.g., 'barcode17') and adds it as a column.
#   5. Adds a column for the algorithm based on the directory of the file.
#   6. Removes the now redundant filename column.

# The 'read_and_label' function is then applied to each directory and algorithm pair.
# As a result, it creates a combined data frame that consolidates all files,
# categorizing them by algorithm and labeling each with its appropriate barcode.

# *****************************************************************************#
# 1. load packages ----
# *****************************************************************************#

# -----------------------------------------------------------------------------#
# clear environment
rm(list = ls())


# -----------------------------------------------------------------------------#
# load packages
library(janitor, quietly = TRUE)
library(data.table, quietly = TRUE)
library(tidyverse, quietly = TRUE)


# *****************************************************************************#
# 2. import data ----
# *****************************************************************************#

# specify the study dir name
STUDY = "wf-nanorave-nanopore-amplicon/input/project-turkana-embatalk-k13/"


## __a. import genome bed file ----
# =============================================================================#

# P. falciparum bed file
genome_bed <- read_delim(
                         "resources-genome/bed/Pf3D7_v51.bed",
                         delim = "\t",
                         col_names = c("chrom", "start", "end", "name", "score",
                                       "strand", "source", "feature_type",
                                       "attributes", "description"),
                        show_col_types = FALSE) %>%
  mutate(description = str_remove(description, ".+description=")) %>%
  select(name, description)


## __b. import sample names ----
# =============================================================================#
sample_ids <- read_csv(paste0(STUDY, "/study/sample_names.csv"),
                      show_col_types = FALSE)


## __c. import TSVs (clair3 and medaka) ----
# =============================================================================#


# -----------------------------------------------------------------------------#
# define the folder paths and corresponding algorithm names
paths_and_algorithms <- tibble(
      path = c(
              paste0(STUDY, "/vcf_to_tsv/clair3/"),
              paste0(STUDY, "/vcf_to_tsv/medaka_haploid_variant/"),
              paste0(STUDY, "/vcf_to_tsv/medaka_variant/")
      ),
      algorithm = c(
                    "clair3",
                    "medaka_haploid_variant",
                    "medaka_variant"
                    )
      )


# -----------------------------------------------------------------------------#
# function to read files from a given path and add algorithm and barcode columns
read_and_label <- function(path, algorithm) {
  # -------------------------------------------#
  # list all TSV files in the given path
  files <- dir(path, pattern = "*.tsv$", full.names = TRUE)
  # -------------------------------------------#
  # read each TSV file and create a list of data frames
  data_list <- lapply(files, function(file) {
    data_frame <- read_tsv(file, col_types = cols(.default = col_character()))
    data_frame$filename <- basename(file)  # Add the filename to each data frame
    return(data_frame)
  })
  # -------------------------------------------#
  # combine all data frames into one
  data <- bind_rows(data_list)
  # -------------------------------------------#
  # add algorithm and barcode columns, then remove the filename column
  data <- data %>%
    mutate(
      algorithm = algorithm,  # Add algorithm column
      barcode = str_extract(filename, "barcode\\d+")  # Extract barcode (e.g., barcode17) from filename
    ) %>%
    select(-filename)  # Remove the filename column if not needed
  return(data)
}


# -----------------------------------------------------------------------------#
# apply the function to each path-algorithm pair and combine the results into one data frame
variants_v1 <- paths_and_algorithms %>%
  pmap_dfr(~ read_and_label(..1, ..2)) %>%
  janitor::clean_names()


# -----------------------------------------------------------------------------#
# optionally, print the combined data to view
glimpse(variants_v1, width = 3)


## __d. import TSVs (freebayes) ----
# =============================================================================#


# -----------------------------------------------------------------------------#
# Set the data_path variable to the folder containing the TSV files
data_path <- paste0(STUDY, "/vcf_to_tsv/freebayes/")



# -----------------------------------------------------------------------------#
# get a list of all TSV files in data_path
files <- dir(data_path, pattern = "*.tsv$")


# -----------------------------------------------------------------------------#
# function to read TSV and rename columns
read_and_rename <- function(filename, data_path) {
  # -------------------------------------------#
  # read the file
  data <- fread(file.path(data_path, filename), colClasses = "character")
  # -------------------------------------------#
  # extract barcode from the filename
  barcode <- str_extract(filename, "barcode\\d+")
  # -------------------------------------------#
  # rename columns based on barcode
  colnames(data)[which(colnames(data) == paste0(barcode, ".GT"))] <- "sample.GT"
  colnames(data)[which(colnames(data) == paste0(barcode, ".AD"))] <- "sample.AD"
  colnames(data)[which(colnames(data) == paste0(barcode, ".DP"))] <- "sample.DP"
  colnames(data)[which(colnames(data) == paste0(barcode, ".GQ"))] <- "sample.GQ"
  colnames(data)[which(colnames(data) == paste0(barcode, ".PGT"))] <- "sample.PGT"
  colnames(data)[which(colnames(data) == paste0(barcode, ".PID"))] <- "sample.PID"
  colnames(data)[which(colnames(data) == paste0(barcode, ".PL"))] <- "sample.PL"
  colnames(data)[which(colnames(data) == paste0(barcode, ".PS"))] <- "sample.PS"
  
  return(data)
}


# -----------------------------------------------------------------------------#
# create a tibble with a column "filename" containing the file names
variants_v2 <- tibble(barcode = files,
                      algorithm = "freebayes") %>%
  mutate(data = map(barcode, ~ read_and_rename(., data_path))) %>%  # read and rename columns for each file
  unnest(cols = c(data)) %>%                                        # unnest the list column "data" into separate rows
  mutate(barcode = str_remove(barcode, "-freebayes.tsv")) %>%
  janitor::clean_names()


# -----------------------------------------------------------------------------#
# optionally, print the combined data to view
glimpse(variants_v2)



## __e. merge datasets ----
# =============================================================================#
variants_v3 <- bind_rows(variants_v1, variants_v2) %>%
  inner_join(., genome_bed,               # merge with genome bed file
             by = c("gene_id" = "name")) %>%
  mutate(description = case_when(description == "apical membrane antigen 1" ~ "ama1",
                                 description == "chloroquine resistance transporter" ~ "crt",
                                 description == "bifunctional dihydrofolate reductase-thymidylate synthase" ~ "dhfr",
                                 description == "hydroxymethyldihydropterin pyrophosphokinase-dihydropteroate synthase" ~ "dhps",
                                 description == "multidrug resistance protein 1" ~ "mdr1",
                                 description == "kelch protein K13" ~ "k13",
                                 TRUE ~ description)) %>%
  filter(description %in% c("crt", "dhps", "dhfr", "mdr1", "ama1")) %>%
  inner_join(., sample_ids, by = "barcode")


# *****************************************************************************#
# 3 . reshape data ----
# *****************************************************************************#


# -----------------------------------------------------------------------------#
# variables
READ_DEPTH = 250 # mininum read-depth to retain variant
AGREEMENT = 1  # minimum number of algorithms to retain variant



# -----------------------------------------------------------------------------#
# reshape from wide to long and filter based on algorithm-concordance
variants_v4 <- variants_v3 %>%
  select(barcode, sample_id, description, type, sample_dp, sample_ad, hgvs_c, hgvs_p, algorithm) %>%
  separate(sample_ad,
           into = c("coverage_wt", "coverage_alt"),
           sep = ",",
           remove = FALSE,
           extra = "merge") %>%
  mutate_at(.vars = c("sample_dp", "coverage_wt", "coverage_alt"), .funs = as.numeric) %>%
  mutate(
         freq_alt = coverage_alt/(coverage_wt + coverage_alt) * 100,
         freq_alt = round(freq_alt, 2),
         codon = as.numeric(str_extract(hgvs_p, "\\d+")),
         sample_dp = coverage_wt + coverage_alt,
         freq_depth = paste0(coverage_alt, " [", freq_alt, "]")
         ) %>%
  filter(sample_dp >= READ_DEPTH & type == "SNP" & freq_alt >= 5) %>%
  pivot_wider(
              # removed sample_ad from id columns
              id_cols = c(barcode, sample_id, description, codon, hgvs_p),
              names_from = "algorithm",
              values_from = "freq_depth",
              # values_fill = list(sample_ad = list(NA_character_))  # fill missing values with NA
              ) %>%
  arrange(barcode, description, codon) %>%
  rowwise() %>%             # apply the following operations to each row
  mutate(concordance = sum(!is.na(c(
                                    clair3,
                                    medaka_variant,
                                    medaka_haploid_variant
                                    # freebayes
                                    )))) %>%
  ungroup() %>%             # remove the row-wise grouping
  filter(concordance >= AGREEMENT)  # filter alleles detected by "x" programs


# -----------------------------------------------------------------------------#
# cleanup: remove variables to free memory
rm(STUDY, data_path, files, read_and_label, read_and_rename, paths_and_algorithms)

# *****************************************************************************#
# 4 . allele frequencies ----
# *****************************************************************************#

df_allele_freq <- variants_v4 %>%
  reframe(n=n(), .by = c(description, hgvs_p))
