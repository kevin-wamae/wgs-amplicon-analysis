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

## __a. specify input-dir path ----

# make sure to end file paths with `/`
PATH_STUDY = "input/turkana_embatalk/"
PATH_RUN = "2024_02_23_illumina_2x300_fail/"
PATH_DATE = "2024_04_21-bcftools/"



## __b. import genome bed file ----
# =============================================================================#

# P. falciparum bed file
genome_bed <- read_delim(
                         "../resources-genome/bed/Pf3D7_v51.bed",
                         delim = "\t",
                         col_names = c("chrom", "start", "end", "name", "score",
                                       "strand", "source", "feature_type",
                                       "attributes", "description"),
                        show_col_types = FALSE) %>%
  mutate(description = str_remove(description, ".+description=")) %>%
  select(name, description)



## __c. import sample names ----
# =============================================================================#
# sample_ids <- read_csv(paste0(STUDY, "/study/sample_names.csv"),
#                       show_col_types = FALSE)


## __d. import TSVs (clair3 and medaka) ----
# =============================================================================#


# define the folder paths
# -----------------------------------------------------------------------------#
file_list <- dir(path = paste0(PATH_STUDY, PATH_RUN, PATH_DATE,
                               "8_extracted_variants"),
                 pattern = "*.tsv$",
                 full.names = TRUE)



# import data
# -----------------------------------------------------------------------------#
variants_v1 = tibble(source = file_list) %>%
  tidyr::extract(source, "file", "(8_extracted_variants.+)", remove = FALSE) %>%
  mutate(data = lapply(source, function(x) {
    dt <- fread(x, skip=0, blank.lines.skip=TRUE, fill=TRUE,sep="\t")
    dt[] <- lapply(dt, as.character)  # Convert all columns to character
    return(dt)
  })) %>%
  unnest(data) %>%
  select(-source) %>%
  mutate(file = str_remove_all(file, "8_extracted_variants/|.tsv")) %>%
  janitor::clean_names()



## __e. append target names ----
# =============================================================================#
variants_v2 <- variants_v1 %>%
  left_join(., genome_bed,               # merge with genome bed file
             by = c("geneid" = "name")) %>%
  mutate(description = case_when(description == "apical membrane antigen 1" ~ "ama1",
                                 description == "chloroquine resistance transporter" ~ "crt",
                                 description == "bifunctional dihydrofolate reductase-thymidylate synthase" ~ "dhfr",
                                 description == "hydroxymethyldihydropterin pyrophosphokinase-dihydropteroate synthase" ~ "dhps",
                                 description == "multidrug resistance protein 1" ~ "mdr1",
                                 description == "kelch protein K13" ~ "k13",
                                 TRUE ~ description)) %>%
  filter(description %in% c("crt", "dhps", "dhfr", "mdr1", "ama1"))




# *****************************************************************************#
# 3 . reshape data ----
# *****************************************************************************#


# -----------------------------------------------------------------------------#
# variables
READ_DEPTH = 250 # mininum read-depth to retain variant
AGREEMENT = 1  # minimum number of algorithms to retain variant



# -----------------------------------------------------------------------------#
# reshape from wide to long and filter based on algorithm-concordance
variants_v3 <- variants_v2 %>%
  select(file, chrom, pos, geneid, description, hgvs_c, hgvs_p, cds_pos, aa_pos, gt, ad, af_ref, af_alt) %>%
  separate(ad,
           into = c("coverage_wt", "coverage_alt"),
           sep = ",",
           remove = FALSE,
           extra = "merge") %>%
  mutate_at(.vars = c("ad", "coverage_wt", "coverage_alt"), .funs = as.numeric) %>%
  mutate(
         sample_dp = coverage_wt + coverage_alt,
         #freq_depth = paste0(coverage_alt, " [", freq_alt, "]")
         ) #%>%
  # filter(sample_dp >= READ_DEPTH & type == "SNP" & freq_alt >= 5) %>%
  # pivot_wider(
  #             # removed sample_ad from id columns
  #             id_cols = c(barcode, sample_id, description, codon, hgvs_p),
  #             names_from = "algorithm",
  #             values_from = "freq_depth",
  #             # values_fill = list(sample_ad = list(NA_character_))  # fill missing values with NA
  #             ) %>%
  # arrange(barcode, description, codon) %>%
  # rowwise() %>%             # apply the following operations to each row
  # mutate(concordance = sum(!is.na(c(
  #                                   clair3,
  #                                   medaka_variant,
  #                                   medaka_haploid_variant
  #                                   # freebayes
  #                                   )))) %>%
  # ungroup() %>%             # remove the row-wise grouping
  # filter(concordance >= AGREEMENT)  # filter alleles detected by "x" programs

# cleanup: remove variables to free memory
rm(STUDY, data_path, files, read_and_label, read_and_rename, paths_and_algorithms)


# *****************************************************************************#
# 4. allele frequencies ----
# *****************************************************************************#

## ____check list of available markers ----
# -----------------------------------------------------------------------------#
(
  unique(variants_v3$description)
)



## __a. PFAMA1 ----
# -----------------------------------------------------------------------------#

clusters_AMA1 <- variants_v3 %>%
  filter(description == "ama1") %>%
  pivot_wider(id_cols = file, names_from = aa_pos, values_from = af_alt)



## __b. PFCRT ----
# -----------------------------------------------------------------------------#

clusters_CRT <- variants_v3 %>%
  filter(description == "crt") %>%
  pivot_wider(id_cols = file, names_from = aa_pos, values_from = af_alt)



## __c. PFMDR1 ----
# -----------------------------------------------------------------------------#

clusters_MDR1 <- variants_v3 %>%
  filter(description == "crt") %>%
  pivot_wider(id_cols = file, names_from = aa_pos, values_from = af_alt)



## __d. PFDHPS ----
# -----------------------------------------------------------------------------#

clusters_DHPS <- variants_v3 %>%
  filter(description == "dhps") %>%
  pivot_wider(id_cols = file, names_from = aa_pos, values_from = af_alt)



## __e. PFDHFR ----
# -----------------------------------------------------------------------------#

clusters_DHFR <- variants_v3 %>%
  filter(description == "dhfr") %>%
  pivot_wider(id_cols = file, names_from = aa_pos, values_from = af_alt)



## __f. PFK13 ----
# -----------------------------------------------------------------------------#

clusters_K13 <- variants_v3 %>%
  filter(description == "k13") %>%
  pivot_wider(id_cols = file, names_from = aa_pos, values_from = af_alt)

