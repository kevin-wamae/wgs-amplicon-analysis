##___compute allele frequencies, regardless source ----
# -----------------------------------------------------------------------------#

df_freqHap_DHFR_All <- df_clusters_DHFR %>%
  select(s_Sample, contains("pos"), haplotype) %>%
  # define level of resistance
  mutate(
         profile_51 = case_when(pos51 == "51N" ~ "wt", .default = "mut"),
         profile_59 = case_when(pos59 == "59C" ~ "wt", .default = "mut"),
         profile_108 = case_when(pos108 == "108S" ~ "wt", .default = "mut")
         ) %>%
  rowwise() %>%
  # calculate resistance profiles for dhfr and total mutations for dhfr
  mutate(
         mutations = sum(profile_51 == "mut", profile_59 == "mut", profile_108 == "mut", na.rm = TRUE),
         dhfr_resistance = get_resistance_profile(mutations, n=sum(!is.na(c(profile_51, profile_59, profile_108)))),
         haplotype = paste0(haplotype," (", dhfr_resistance, ")")
         ) %>%
  ungroup() %>%
  reframe(
          dhfr_resistance,
          haplotype = paste(sort(unique(haplotype)), collapse = ","),
          dhfr_resistance = paste(unique(sort(dhfr_resistance)), collapse = ","),
          .by = c(s_Sample)
          ) %>%
  distinct(s_Sample, .keep_all = TRUE) %>%
  summarise(count=n(), .by = c(haplotype)) %>%
  mutate(
         freq = count/sum(count),
         freq = round(freq * 100, 1),
         freq = paste0(freq, " [", count, "]"),
         variant = case_when(
                             haplotype == wt_haplotype ~ "wildtype",
                             str_detect(haplotype, ",") ~ "mixed",
                             TRUE ~ "mutant",
                             ),
         total = sum(count)
         ) %>%
  select(-count) %>%
  arrange(desc(freq))



##___compute allele frequencies, by source ----
# -----------------------------------------------------------------------------#

df_freqHap_DHFR_Source <- df_clusters_DHFR %>%
  select(source, s_Sample, contains("pos"), haplotype) %>%
  # define level of resistance
  mutate(
         profile_51 = case_when(pos51 == "51N" ~ "wt", .default = "mut"),
         profile_59 = case_when(pos59 == "59C" ~ "wt", .default = "mut"),
         profile_108 = case_when(pos108 == "108S" ~ "wt", .default = "mut")
         ) %>%
  rowwise() %>%
  # calculate resistance profiles for dhfr and total mutations for dhfr
  mutate(
         mutations = sum(profile_51 == "mut", profile_59 == "mut", profile_108 == "mut", na.rm = TRUE),
         dhfr_resistance = get_resistance_profile(mutations, n=sum(!is.na(c(profile_51, profile_59, profile_108)))),
         haplotype = paste0(haplotype," (", dhfr_resistance, ")")
         ) %>%
  ungroup() %>%
  reframe(
          source, dhfr_resistance,
          haplotype = paste(sort(unique(haplotype)), collapse = ","),
          dhfr_resistance = paste(unique(sort(dhfr_resistance)), collapse = ","),
          .by = c(s_Sample)
          ) %>%
  distinct(s_Sample, .keep_all = TRUE) %>%
  summarise(count=n(), .by = c(source, haplotype)) %>%
  mutate(
         freq = count/sum(count), .by = source,
         freq = round(freq * 100, 1),
         freq = paste0(freq, " [", count, "]"),
         variant = case_when(
                             haplotype == wt_haplotype ~ "wildtype",
                             str_detect(haplotype, ",") ~ "mixed",
                             TRUE ~ "mutant",
                             ),
         total = sum(count)
         ) %>%
  select(-count) %>%
  arrange(source, desc(freq))


