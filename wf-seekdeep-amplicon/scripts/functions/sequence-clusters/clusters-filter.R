# DESCRIPTION OF R SCRIPT
# =============================================================================#
# This R script is designed for extracting data from raw_selectedClustersInto matching
# the target define in the object: STRING_TARGET
#
# 1. Generation of a Vector of Polymorphic Codons:
#    This section of the script identifies specific polymorphic positions from
#    genomic sequences that match predefined target strings in gene names and
#    genome annotations. It extracts these positions to create a numeric vector
#    representing the codons of interest. The process involves:
#      - Filtering sequences by specific gene name and genome annotation patterns.
#      - Extracting numeric codon positions from the sequences.
#      - Cleansing and converting extracted data to a numeric format.
#
# 2. Filtering Variant Data and Extracting Codon Information:
#    The second part of the script uses the vector of codon positions to filter
#    and extract specific variant data from the sequences. It constructs new columns
#    for each codon position, extracts the corresponding variant, and assembles a
#    haplotype string for each sequence. Key steps include:
#      - Mapping each codon position to extract specific variants using a dynamic
#        pattern matching.
#      - Collapsing extracted codon positions into a single string to form the haplotype.
#      - Cleaning and formatting the final codon position and haplotype strings.
#
# The script utilizes the `dplyr`, `stringr`, and `purrr` libraries to handle data
# manipulation, string operations, and functional programming techniques respectively,
# ensuring efficient processing of genomic data.
# =============================================================================#



# Note that SeekDeep may also report user-supplied alleles for positions that
# are not truly segregating, which could result in some positions appearing
# without variation, so we'll generate df_clusters_Target to include all and
# df_clusters_Segregating to include only segregating positions.
# -----------------------------------------------------------------------------#


# Define a function to parse each entry
# -----------------------------------------#
parse_entry <- function(entry) {
  if (str_detect(entry, "^(\\d+)([A-Z]+)$")) {
      # Case 1: Number(s) followed by Letter(s)
      m <- str_match(entry, "^(\\d+)([A-Z]+)$")
      position <- as.numeric(m[2])
      ref_allele <- m[3]
      mut_allele <- NA
      mutation <- FALSE
  } else if (str_detect(entry, "^([A-Z]+)(\\d+)$")) {
      # Case 2: Letter(s) followed by Number(s)
      m <- str_match(entry, "^([A-Z]+)(\\d+)$")
      position <- as.numeric(m[3])
      ref_allele <- m[2]
      mut_allele <- NA
      mutation <- FALSE
  } else if (str_detect(entry, "^([A-Z]+)(\\d+)([A-Z]+)$")) {
      # Case 3: Letter(s), Number(s), then Letter(s)
      m <- str_match(entry, "^([A-Z]+)(\\d+)([A-Z]+)$")
      position <- as.numeric(m[3])
      ref_allele <- m[2]
      mut_allele <- m[4]
      mutation <- TRUE
  } else {
      position <- NA
      ref_allele <- NA
      mut_allele <- NA
      mutation <- NA
  }
  tibble(position = position, ref_allele = ref_allele, mut_allele = mut_allele, mutation = mutation)
}



# process the h_AATyped data
# -----------------------------------------#
positions_df <- raw_selectedClustersInfo %>%
  mutate(id = row_number(),
         h_AATyped = str_replace(h_AATyped, "PF3D7_.+--", "")) %>%
  select(id, h_AATyped) %>%
  mutate(entries = str_split(h_AATyped, ":", simplify = FALSE)) %>%
  unnest(entries) %>%
  mutate(parsed = map(entries, parse_entry)) %>%
  unnest_wider(parsed) %>%
  mutate(allele_to_use = ifelse(mutation, mut_allele, ref_allele))



# pivot the data to wide format
# -----------------------------------------#
positions_df_wide <- positions_df %>%
  select(id, position, allele_to_use) %>%
  pivot_wider(names_from = position, names_prefix = "pos", values_from = allele_to_use)



# join back to the original data and create haplotype
# -----------------------------------------#
df_clusters_Target <- raw_selectedClustersInfo %>%
  mutate(id = row_number()) %>%
  left_join(positions_df_wide, by = "id") %>%
  rowwise() %>%
  mutate(
    codon_pos = paste(c_across(starts_with("pos")), collapse = ", "),
    haplotype = paste(c_across(starts_with("pos")), collapse = "")
  ) %>%
  ungroup()



# generate a vector of polymorphic codons that are truly segregating
# -----------------------------------------------------------------------------#
# extract columns related to codon positions from the dataframe
df_filtered <- df_clusters_Target %>% select(starts_with("pos"))



# identify columns that have more than one unique value (indicating polymorphism)
positions_Segregating <- names(df_filtered)[sapply(df_filtered, function(col) length(unique(col)) > 1)]



# Remove the 'pos' prefix and convert the result to a numeric vector representing
# polymorphic codon positions
positions_Segregating <- as.numeric(gsub("pos", "", positions_Segregating))



# filter variant data at segregating alleles, only
# -----------------------------------------------------------------------------#
df_clusters_Segregating <- raw_selectedClustersInfo %>%
  filter(str_detect(p_name, STRING_TARGET)) %>%
  mutate(
    purrr::map_dfc(
      set_names(positions_Segregating, paste0("pos", positions_Segregating)),
      ~ str_extract(string = h_AATyped, pattern = paste0(.x, "."))
    )) %>%
  rowwise() %>%
  mutate(
    codon_pos = paste(c_across(all_of(paste0("pos", positions_Segregating))), collapse = ", "),
    codon_pos = str_remove_all(string = codon_pos, pattern = "-."),
    haplotype = str_remove_all(string = codon_pos, pattern = "\\d+|,\\s")
  ) %>%
  ungroup()



# remove temporary objects
rm(df_filtered)



##___print a message in the console ----
# -----------------------------------------------------------------------------#

# Using yellow for the border
cat("\033[1m\033[33m", "\n##############################################################", "\033[0m")

# Using magenta for the table descriptions
cat("\033[1m\033[35m", "\nData Summary:", "\033[0m")

# Using cyan for the first table description
cat("\033[1m\033[36m", "\n1. df_clusters_Segregating  - This table contains only the segregating codons, i.e., positions with true variation.", "\033[0m")

# Using cyan for the second table description
cat("\033[1m\033[36m", "\n2. df_clusters_Target - This table includes all initial variants identified. It may include both segragating and
                        non-segregating positions based on what the user requested for.", "\033[0m")

# Yellow for the border again
cat("\033[1m\033[33m", "\n##############################################################\n", "\033[0m")

