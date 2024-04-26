##___determine allele changes by codon-position ----
# -----------------------------------------------------------------------------#
df_Alleles <- df_clusters_MDR1 %>%
  select(starts_with("pos")) %>%
  # transform: wide to long format
  # ---------------------------------#
  pivot_longer(
               cols = starts_with("pos"),
               names_to = "codon",
               values_to = "allele"
               ) %>%
  # remove duplicates to reduce redundancy
  # ---------------------------------#
  distinct(codon, allele, .keep_all = TRUE) %>%
  # extract numerical part and letter part of allele into separate columns
  # ---------------------------------#
  mutate(
         codon = as.numeric(str_extract(allele, "\\d+")),
         allele = str_extract(allele, "[A-Z]"),
         ) %>%
  # merge with reference data to add wildtype allele information
  # ---------------------------------#
  left_join(
            data.frame(position = positions_MDR1, wildtype = wt_alleles),
            by = c("codon" = "position")
            ) %>%
  # apply filter based on group count: if more than one, filter by allele difference and retain relevant codons
  # ---------------------------------#
  filter(if (n() > 1) {allele != wildtype | is.na(wildtype)} else {TRUE}, .by = codon) %>%
  mutate(codon_residue = paste0(wildtype, codon, allele)) %>%
  select(codon, codon_residue)



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
  mutate_at(.vars = 3:5, replace_na, "0 [0]") %>%
  # replace codon allele
  # ---------------------------------#
  select(-codon_allele) %>%
  left_join(df_Alleles, by = "codon") %>%
  # drop alleles with 100% wildtype frequency or missing allele information (stop codons)
  # ---------------------------------#
  filter(! str_detect(wildtype, "100 \\[") & ! is.na(codon_residue) | c(86, 184, 1246)) %>%
  mutate(codon = codon_residue) %>%
  select(-codon_residue)



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
  mutate_at(.vars = 4:6, replace_na, "0 [0]") %>%
  # replace codon allele
  # ---------------------------------#
  select(-codon_allele) %>%
  left_join(df_Alleles, by = "codon") %>%
  # drop alleles with 100% wildtype frequency or missing allele information (stop codons)
  # ---------------------------------#
  filter(! str_detect(wildtype, "100 \\[") & ! is.na(codon_residue) | codon %in% c(86, 184, 1246)) %>%
  mutate(codon = codon_residue) %>%
  select(-codon_residue)

