##___determine allele changes by codon-position ----
# -----------------------------------------------------------------------------#
df_WT_Alleles <- df_clusters_Target %>%
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
            data.frame(position = positions_Target, wildtype = wt_alleles),
            by = c("codon" = "position")
            ) %>%
  # apply filter based on group count: if more than one, filter by allele difference
  # ---------------------------------#
  filter(if (n() > 1) {allele != wildtype | is.na(wildtype)} else {TRUE}, .by = codon) %>%
  mutate(aa_change = paste0(wildtype, codon, allele)) %>%
  select(codon, aa_change)



##___compute allele frequencies, regardless source ----
# -----------------------------------------------------------------------------#

df_freqSNP_All <- df_clusters_Target %>%
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
            data.frame(position = positions_Target, wildtype = wt_alleles),
            by = c("codon" = "position")) %>%
  # code for infection-type
  # ---------------------------------#
  mutate(
         aa_change = paste0(wildtype, codon, allele),
         variant = case_when(
                             str_detect(allele, "\\,") ~ "mixed",
                             allele == wildtype ~ "wildtype",
                             allele != wildtype ~ "mutant"),
         freq = paste0(freq, " [", count, "]")
         ) %>%
  # transform: long to wide
  # ---------------------------------#
  pivot_wider(
              id_cols = c(codon, aa_change),
              names_from = "variant",
              values_from = "freq") %>%
  # ensure potentially missing columns are present
  # ---------------------------------#
  mutate(
         wildtype = if ("wildtype" %in% names(.)) wildtype else NA_character_,
         mutant = if ("mutant" %in% names(.)) mutant else NA_character_,
         mixed = if ("mixed" %in% names(.)) mixed else NA_character_
         ) %>%
  arrange(codon, str_length(aa_change)) %>%
  # collapse rows by source and codon position
  # ---------------------------------#
  summarise(
            aa_change = paste(aa_change, collapse=", "),
            wildtype = ifelse(all(is.na(wildtype)), NA_character_,
                              paste(wildtype[!is.na(wildtype)], collapse=", ")),
            mutant = ifelse(all(is.na(mutant)), NA_character_,
                            paste(mutant[!is.na(mutant)], collapse=", ")),
            mixed = ifelse(all(is.na(mixed)), NA_character_,
                           paste(mixed[!is.na(mixed)], collapse=", ")),
            .by = c(codon)
            ) %>%
  # replace na in allele frequencies with 0
  # ---------------------------------#
  mutate_at(.vars = 3:5, replace_na, "0 [0]") %>%
  # drop alleles with 100% wildtype frequency or missing allele information (stop codons)
  # ---------------------------------#
  filter(! str_detect(wildtype, "100 \\["))



##___compute allele frequencies, by source ----
# -----------------------------------------------------------------------------#

df_freqSNP_Source <- df_clusters_Target %>%
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
            data.frame(position = positions_Target, wildtype = wt_alleles),
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
            .by = c(source, codon)
            ) %>%
  # replace na in allele frequencies with 0
  # ---------------------------------#
  mutate_at(.vars = 4:6, replace_na, "0 [0]") %>%
  # replace codon allele
  # ---------------------------------#
  select(-codon_allele) %>%
  left_join(df_WT_Alleles, by = "codon") %>%
  # drop alleles with 100% wildtype frequency or missing allele information (stop codons)
  # ---------------------------------#
  filter(! str_detect(wildtype, "100 \\[") & ! is.na(aa_change)) %>%
  mutate(codon = aa_change) %>%
  select(-aa_change)



##___compute allele frequencies, by sample ----
# -----------------------------------------------------------------------------#

df_freqSNP_Sample <- df_clusters_Target %>%
  # select relevant columns
  # ---------------------------------#
  select(R1_name, s_Sample, source, starts_with("pos")) %>%
  # collapse alleles per codon per sample
  # ---------------------------------#
  reframe(source, R1_name,
          across(starts_with("pos"), ~paste(unique(sort(.)), collapse = ","),
                 .names = "{.col}"),
          .by = s_Sample) %>%
  # remove duplicates
  # ---------------------------------#
  distinct(s_Sample, .keep_all = TRUE)



##___remove temporary objects ----
# -----------------------------------------------------------------------------#

rm(df_WT_Alleles)
