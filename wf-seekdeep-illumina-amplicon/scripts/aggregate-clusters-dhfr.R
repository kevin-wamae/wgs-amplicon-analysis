# vector of polymorphic codons
(
  positions_DHFR <- selectedClustersInfo %>%
    filter(p_name == "PFDHFR" & str_detect(h_AATyped, "^PF3D7")) %>%
    mutate(h_AATyped = str_remove(h_AATyped, "16A:")) %>%  # drop codon 16
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
clusters_DHFR <- selectedClustersInfo %>%
  filter(p_name == "PFDHFR") %>%
  mutate(
        purrr::map_dfc(
                      set_names(positions_DHFR, paste0("pos", positions_DHFR)),
                      ~ str_extract(string = h_AATyped, pattern = paste0(.x, "."))
                      )) %>%
  rowwise() %>%
  mutate(
         codon_pos = paste(c_across(all_of(paste0("pos", positions_DHFR))), collapse = ", "),
         codon_pos = str_remove_all(string = codon_pos, pattern = "-."),
         haplotype = str_remove_all(string = codon_pos, pattern = "\\d+|,\\s"),
         ) %>%
  ungroup()