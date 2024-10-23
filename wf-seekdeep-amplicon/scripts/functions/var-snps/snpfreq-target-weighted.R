##___compute allele frequencies, regardless of source ----
# -----------------------------------------------------------------------------#

df_freqSNP_All <- df_clusters_Target %>%
  # # transform: wide to long
  # ---------------------------------#
  pivot_longer(
               cols = starts_with("pos"),
               names_to = "codon",
               values_to = "allele"
               ) %>%
  # nest data by location and codon
  tidyr::nest(allele_data = -c(codon)) %>%
  # compute sample size
  nplyr::nest_mutate(
                     allele_data,
                     pop_size = n_distinct(s_Sample),
                     across(.cols = everything())
                     ) %>%
  # compute allele frequencies
  nplyr::nest_group_by( # grouping
                       .nest_data = allele_data,
                       allele) %>%
  nplyr::nest_mutate(.nest_data = allele_data, # allele frequency calculation
                     freq = round(sum(as.numeric(c_AveragedFrac)/pop_size),2)
                     ) %>%
  tidyr::unnest(cols = c(allele_data)) %>%
  # remove duplicates
  distinct(codon, allele, pop_size, freq) %>%
  # drop string 'pos' and digits in allele
  # ---------------------------------#
  mutate(
         codon = str_remove(codon, "pos"),
         allele = str_remove(allele, "\\d+"),
         count = round(pop_size * freq, 0)
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
  pivot_wider(id_cols = c(codon, aa_change),
              names_from = "variant",
              values_from = "freq") %>%
  summarise(
            aa_change = paste(aa_change, collapse=", "),
            wildtype = ifelse(all(is.na(wildtype)), NA_character_,
                              paste(wildtype[!is.na(wildtype)], collapse=", ")),
            mutant = ifelse(all(is.na(mutant)), NA_character_,
                            paste(mutant[!is.na(mutant)], collapse=", ")),
            .by = c(codon)
            ) %>%
  # replace na in allele frequencies with 0
  # ---------------------------------#
  mutate_at(.vars = 3:4, replace_na, "0 [0]") %>%
  # drop alleles with 100% wildtype frequency or missing allele information (stop codons)
  # ---------------------------------#
  filter( !grepl("^0 \\[0\\]$", wildtype) & !grepl("^0 \\[0\\]$", mutant) ) %>%
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



##___compute allele frequencies, by source ----
# -----------------------------------------------------------------------------#

df_freqSNP_Source <- df_clusters_Target %>%
  # # transform: wide to long
  # ---------------------------------#
  pivot_longer(
               cols = starts_with("pos"),
               names_to = "codon",
               values_to = "allele"
               ) %>%
  # nest data by location and codon
  tidyr::nest(allele_data = -c(source, codon)) %>%
  # compute sample size
  nplyr::nest_mutate(
                     allele_data,
                     pop_size = n_distinct(s_Sample),
                     across(.cols = everything())
                     ) %>%
  # compute allele frequencies
  nplyr::nest_group_by( # grouping
                       .nest_data = allele_data,
                       allele) %>%
  nplyr::nest_mutate(.nest_data = allele_data, # allele frequency calculation
                     freq = round(sum(as.numeric(c_AveragedFrac)/pop_size),2)
                     ) %>%
  tidyr::unnest(cols = c(allele_data)) %>%
  # remove duplicates
  distinct(source, codon, allele, pop_size, freq) %>%
  # drop string 'pos' and digits in allele
  # ---------------------------------#
  mutate(
         codon = str_remove(codon, "pos"),
         allele = str_remove(allele, "\\d+"),
         count = round(pop_size * freq, 0)
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
              id_cols = c(source, codon, aa_change),
              names_from = "variant",
              values_from = "freq") %>%
  summarise(
            aa_change = paste(aa_change, collapse=", "),
            wildtype = ifelse(all(is.na(wildtype)), NA_character_,
                              paste(wildtype[!is.na(wildtype)], collapse=", ")),
            mutant = ifelse(all(is.na(mutant)), NA_character_,
                            paste(mutant[!is.na(mutant)], collapse=", ")),
            .by = c(source, codon)
            ) %>%
  # replace na in allele frequencies with 0
  # ---------------------------------#
  mutate_at(.vars = 4:5, replace_na, "0 [0]") %>%
  # drop alleles with 100% wildtype frequency or missing allele information (stop codons)
  # ---------------------------------#
  filter( !grepl("^0 \\[0\\]$", wildtype) & !grepl("^0 \\[0\\]$", mutant) ) %>%
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

# Yellow for the border again
cat("\033[1m\033[33m", "\n##############################################################\n", "\033[0m")
