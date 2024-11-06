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
            data.frame(
                       position = unique(as.numeric(positions_df$position)),
                       wildtype = wt_alleles),
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
            aa_change_wildtype = paste(aa_change[!is.na(wildtype)], collapse=", "),
            aa_change_mutant = paste(aa_change[!is.na(mutant)], collapse=", "),
            aa_change_mixed = paste(aa_change[!is.na(mixed)], collapse=", "),
            wildtype = ifelse(all(is.na(wildtype)), NA_character_,
                              paste(wildtype[!is.na(wildtype)], collapse=", ")),
            mutant = ifelse(all(is.na(mutant)), NA_character_,
                            paste(mutant[!is.na(mutant)], collapse=", ")),
            mixed = ifelse(all(is.na(mixed)), NA_character_,
                           paste(mixed[!is.na(mixed)], collapse=", ")),
            .by = c(codon)
            ) %>%
    mutate(
      aa_change = paste(
                        ifelse(aa_change_wildtype != "", aa_change_wildtype, NA_character_),
                        ifelse(aa_change_mutant != "", aa_change_mutant, NA_character_),
                        ifelse(aa_change_mixed != "", aa_change_mixed, NA_character_),
                        sep = " / "
                        ) %>%
        str_replace_all("NA / ", "") %>%
        str_replace_all(" / NA", "") %>%
        str_replace_all("NA", "")
      ) %>%
  select(-c(aa_change_wildtype, aa_change_mutant, aa_change_mixed)) %>%
  relocate(aa_change, .after = codon) %>%
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
            data.frame(
                       position = unique(as.numeric(positions_df$position)),
                       wildtype = wt_alleles),
            by = c("codon" = "position")) %>%
  # code for infection-type
  # ---------------------------------#
  mutate(aa_change = paste0(wildtype, codon, allele),
         variant = case_when(
                             str_detect(allele, "\\,") ~ "mixed",
                             allele == wildtype ~ "wildtype",
                             allele != wildtype ~ "mutant"),
         freq = paste0(freq, " [", count, "]")
         ) %>%
  # transform: long to wide
  # ---------------------------------#
  pivot_wider(id_cols = c(source, codon, aa_change),
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
  summarise(
            aa_change = paste(aa_change, collapse=", "),
            wildtype = ifelse(all(is.na(wildtype)), NA_character_,
                              paste(wildtype[!is.na(wildtype)], collapse=", ")),
            mutant = ifelse(all(is.na(mutant)), NA_character_,
                            paste(mutant[!is.na(mutant)], collapse=", ")),
            mixed = ifelse(all(is.na(mixed)), NA_character_,
                           paste(mixed[!is.na(mixed)], collapse=", ")),
            .by = c(source, codon)
            ) %>%
  # replace na in allele frequencies with 0
  # ---------------------------------#
  mutate_at(.vars = 4:6, replace_na, "0 [0]") %>%
  # drop alleles with 100% wildtype frequency or missing allele information (stop codons)
  # ---------------------------------#
  filter(! str_detect(wildtype, "100 \\[")) %>%
  # collapse aa_change to remove redundancy
  # ---------------------------------#
  mutate(
         wt_allele = sapply(str_extract_all(aa_change, "(?<=, |^)[A-Z](?=[0-9])"), function(x) paste(unique(x), collapse = "/")),
         mut_allele = sapply(seq_along(aa_change), function(i) {
           wt <- unique(str_extract_all(aa_change[i], "(?<=, |^)[A-Z](?=[0-9])")[[1]])
           all_alleles <- unique(unlist(str_extract_all(aa_change[i], "[A-Z]")))
           mut <- setdiff(all_alleles, wt)
           paste(mut, collapse = "/")}
           ),
         position = str_extract(aa_change, "[0-9]+"),
         aa_change = paste0(wt_allele, position, mut_allele)
         ) %>%
  select(-c(wt_allele, position, mut_allele))



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



##___print a message in the console ----
# -----------------------------------------------------------------------------#

# Using yellow for the border
cat("\033[1m\033[33m", "\n##############################################################", "\033[0m")

# Using magenta for the table descriptions
cat("\033[1m\033[35m", "\nData Summary:", "\033[0m")

# Using cyan for the first table description
cat("\033[1m\033[36m", "\n1. df_freqSNP_All    - This table shows the aggregated SNP frequencies across all geographical regions", "\033[0m")

# Using cyan for the second table description
cat("\033[1m\033[36m", "\n2. df_freqSNP_Source - This table shows the SNP frequencies by geographical region", "\033[0m")

# Using cyan for the third table description
cat("\033[1m\033[36m", "\n3. df_freqSNP_Sample - This table shows the SNP frequencies per sample", "\033[0m")

# Yellow for the border again
cat("\033[1m\033[33m", "\n##############################################################\n", "\033[0m")


