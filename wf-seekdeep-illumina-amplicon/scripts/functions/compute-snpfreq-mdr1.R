##___compute allele frequencies, regardless source ----
# -----------------------------------------------------------------------------#

df_freqSNP_MDR1_All <- df_clusters_MDR1 %>%
  # select relevant columns
  # ---------------------------------#
  select(s_Sample, starts_with("pos")) %>%
  # collapse alleles per codon per sample
  # ---------------------------------#
  reframe(
    across(starts_with("pos"), ~paste(unique(sort(.)), collapse = ","),
           .names = "{.col}"),
    .by = s_Sample) %>%
  # remove duplicates
  # ---------------------------------#
  distinct(s_Sample, .keep_all = TRUE) %>%
  # transform: wide to long
  # ---------------------------------#
  pivot_longer(
    cols = starts_with("pos"),
    names_to = "codon",
    values_to = "allele"
  ) %>%
  # drop string 'pos'
  # ---------------------------------#
  mutate(codon = str_remove(codon, "pos")) %>%
  # generate allele frequencies
  # ---------------------------------#
  summarise(count=n(), .by = c(codon, allele)) %>%
  arrange(as.numeric(codon)) %>%
  mutate(
    freq = count/sum(count) * 100, .by = c(codon),
    freq = round(freq,1),
    allele = str_remove_all(allele, "\\d+")
  ) %>%
  mutate(codon = as.numeric(codon)) %>%
  # merge with wildtype alleles
  # ---------------------------------#
  left_join(
    data.frame(position = positions_MDR1, wildtype = wt_alleles),
    by = c("codon" = "position")) %>%
  # code for infection-type
  # ---------------------------------#
  mutate(codon_allele = paste0(wildtype, codon, allele),
         variant = case_when(
           str_detect(allele, "\\,") ~ "mixed",
           allele == wildtype ~ "wildtype",
           allele != wildtype ~ "mutant"),
         freq = paste0(freq, " [", count, "]")
  ) %>%
  # transform: long to wide
  # ---------------------------------#
  pivot_wider(id_cols = c(codon, codon_allele),
              names_from = "variant",
              values_from = "freq") %>%
  # ensure potentially missing columns are present
  # ---------------------------------#
  mutate(wildtype = if ("wildtype" %in% names(.)) wildtype else NA_character_,
         mutant = if ("mutant" %in% names(.)) mutant else NA_character_,
         mixed = if ("mixed" %in% names(.)) mixed else NA_character_
  ) %>%
  # collapse rows by source and codon position
  # ---------------------------------#
  summarise(codon_allele = paste(codon_allele, collapse=", "),
            wildtype = first(na.omit(wildtype)),
            mutant = first(na.omit(mutant)),
            mixed = first(na.omit(mixed)),
            .by = c(codon)) %>%
  # replace na in allele frequencies with 0
  # ---------------------------------#
  mutate_at(.vars = 3:5, replace_na, "0 [0]")



##___compute allele frequencies, by source ----
# -----------------------------------------------------------------------------#

df_freqSNP_MDR1_Source <- df_clusters_MDR1 %>%
  # select relevant columns
  # ---------------------------------#
  select(s_Sample, source, starts_with("pos")) %>%
  # collapse alleles per codon per sample
  # ---------------------------------#
  reframe(source,
          across(starts_with("pos"), ~paste(unique(sort(.)), collapse = ","),
                 .names = "{.col}"),
          .by = s_Sample) %>%
  # remove duplicates
  # ---------------------------------#
  distinct(s_Sample, .keep_all = TRUE) %>%
  # transform: wide to long
  # ---------------------------------#
  pivot_longer(
    cols = starts_with("pos"),
    names_to = "codon",
    values_to = "allele"
  ) %>%
  # drop string 'pos'
  # ---------------------------------#
  mutate(codon = str_remove(codon, "pos")) %>%
  # generate allele frequencies
  # ---------------------------------#
  summarise(count=n(), .by = c(source, codon, allele)) %>%
  arrange(source, as.numeric(codon)) %>%
  mutate(
    freq = count/sum(count) * 100, .by = c(source, codon),
    freq = round(freq,1),
    allele = str_remove_all(allele, "\\d+")
  ) %>%
  mutate(codon = as.numeric(codon)) %>%
  # merge with wildtype alleles
  # ---------------------------------#
  left_join(
    data.frame(position = positions_MDR1, wildtype = wt_alleles),
    by = c("codon" = "position")) %>%
  # code for infection-type
  # ---------------------------------#
  mutate(codon_allele = paste0(wildtype, codon, allele),
         variant = case_when(
           str_detect(allele, "\\,") ~ "mixed",
           allele == wildtype ~ "wildtype",
           allele != wildtype ~ "mutant"),
         freq = paste0(freq, " [", count, "]")
  ) %>%
  # transform: long to wide
  # ---------------------------------#
  pivot_wider(id_cols = c(source, codon, codon_allele),
              names_from = "variant",
              values_from = "freq") %>%
  # ensure potentially missing columns are present
  # ---------------------------------#
  mutate(wildtype = if ("wildtype" %in% names(.)) wildtype else NA_character_,
         mutant = if ("mutant" %in% names(.)) mutant else NA_character_,
         mixed = if ("mixed" %in% names(.)) mixed else NA_character_
  ) %>%
  # collapse rows by source and codon position
  # ---------------------------------#
  summarise(codon_allele = paste(codon_allele, collapse=", "),
            wildtype = first(na.omit(wildtype)),
            mutant = first(na.omit(mutant)),
            mixed = first(na.omit(mixed)),
            .by = c(source, codon)) %>%
  # replace na in allele frequencies with 0
  # ---------------------------------#
  mutate_at(.vars = 3:5, replace_na, "0 [0]")
