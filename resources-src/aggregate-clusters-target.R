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



# generate a vector of polymorphic codons reported by SeekDeep
#  - note that SeekDeep may also report user-supplied alleles for positions that
#    are not truly segregating, which could result in some positions appearing
#    without variation
# -----------------------------------------------------------------------------#
(
  positions_Target <- raw_selectedClustersInfo %>%
    filter(str_detect(p_name, STRING_TARGET), str_detect(h_AATyped, STRING_GENOME)) %>%
    select(h_AATyped) %>%
    slice(1) %>%
    mutate(h_AATyped = str_replace(h_AATyped, "PF3D7_.+--", ""),
           numbers = str_extract_all(h_AATyped, "\\d+")) %>%
    unnest(numbers) %>%
    mutate(numbers = as.numeric(numbers)) %>%
    pull(numbers)
)



# filter variant data
# -----------------------------------------------------------------------------#
df_clusters_Target <- raw_selectedClustersInfo %>%
  filter(str_detect(p_name, STRING_TARGET)) %>%
  mutate(
        purrr::map_dfc(
                      set_names(positions_Target, paste0("pos", positions_Target)),
                      ~ str_extract(string = h_AATyped, pattern = paste0(.x, "."))
                      )) %>%
  rowwise() %>%
  mutate(
         codon_pos = paste(c_across(all_of(paste0("pos", positions_Target))), collapse = ", "),
         codon_pos = str_remove_all(string = codon_pos, pattern = "-."),
         haplotype = str_remove_all(string = codon_pos, pattern = "\\d+|,\\s")
         ) %>%
  ungroup()



# generate a vector of polymorphic codons that are truly segregating
# -----------------------------------------------------------------------------#
# extract columns related to codon positions from the dataframe
df_filtered <- df_clusters_Target %>% select(starts_with("pos"))

# identify columns that have more than one unique value (indicating polymorphism)
positions_Segrating <- names(df_filtered)[sapply(df_filtered, function(col) length(unique(col)) > 1)]

# Remove the 'pos' prefix and convert the result to a numeric vector representing
# polymorphic codon positions
positions_Segrating <- as.numeric(gsub("pos", "", positions_Segrating))
