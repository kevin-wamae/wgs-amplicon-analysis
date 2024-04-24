# vector of polymorphic codons
# -----------------------------------------------------------------------------#
(
  positions_K13 <- raw_selectedClustersInfo %>%
    filter(str_detect(p_name, "^PFK13"), str_detect(h_AATyped, "^PF3D7")) %>%
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
df_clusters_K13 <- raw_selectedClustersInfo %>%
  filter(str_detect(p_name, "^PFK13")) %>%
  mutate(
         purrr::map_dfc(
                        set_names(positions_K13, paste0("pos", positions_K13)),
                        ~ str_extract(string = h_AATyped, pattern = paste0(.x, "."))
                        )) %>%
  rowwise() %>%
  mutate(
    codon_pos = paste(c_across(all_of(paste0("pos", positions_K13))), collapse = ", "),
    codon_pos = str_remove_all(string = codon_pos, pattern = "-."),
    haplotype = str_remove_all(string = codon_pos, pattern = "\\d+|,\\s")
  ) %>%
  ungroup()
