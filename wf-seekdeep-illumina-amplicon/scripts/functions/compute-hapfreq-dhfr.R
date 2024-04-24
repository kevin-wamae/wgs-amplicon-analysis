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
  # calculate resistance profiles for dhps and total mutations for dhps
  mutate(
         mutations = sum(profile_51 == "mut", profile_59 == "mut", profile_108 == "mut", na.rm = TRUE),
         dhps_resistance = get_resistance_profile(mutations, n=sum(!is.na(c(profile_51, profile_59, profile_108)))),
         haplotype = paste0(haplotype," (", dhps_resistance, ")")
         ) %>%
  ungroup() %>%
  reframe(
          dhps_resistance,
          haplotype = paste(sort(unique(haplotype)), collapse = ","),
          dhps_resistance = paste(unique(sort(dhps_resistance)), collapse = ","),
          .by = c(s_Sample)
          ) %>%
  distinct(s_Sample, .keep_all = TRUE) %>%
  summarise(count=n(), .by = c(haplotype)) %>%
  mutate(
         freq = count/sum(count),
         freq = round(freq * 100, 1),
         total = sum(count)
         ) %>%
  arrange(desc(freq))


##___compute allele frequencies, by source ----
# -----------------------------------------------------------------------------#

df_freqHap_DHPS_Source <- df_clusters_DHPS %>%
  select(source, s_Sample, contains("pos"), haplotype) %>%
  # define level of resistance
  mutate(
         profile_437 = case_when(pos437 == "437A" ~ "wt", .default = "mut"),
         profile_540 = case_when(pos540 == "540K" ~ "wt", .default = "mut"),
         profile_581 = case_when(pos581 == "581A" ~ "wt", .default = "mut")
         ) %>%
  rowwise() %>%
  # calculate resistance profiles for dhps and total mutations for dhps
  mutate(
         mutations = sum(profile_437 == "mut", profile_540 == "mut", profile_581 == "mut", na.rm = TRUE),
         dhps_resistance = get_resistance_profile(mutations, n=sum(!is.na(c(profile_437, profile_540, profile_581)))),
         haplotype = paste0(haplotype," (", dhps_resistance, ")")
         ) %>%
  ungroup() %>%
  reframe(
          source, dhps_resistance,
          haplotype = paste(sort(unique(haplotype)), collapse = ","),
          dhps_resistance = paste(unique(sort(dhps_resistance)), collapse = ","),
          .by = c(s_Sample)
          ) %>%
  distinct(s_Sample, .keep_all = TRUE) %>%
  summarise(count=n(), .by = c(source, haplotype)) %>%
  mutate(
         freq = count/sum(count), .by = source,
         freq = round(freq * 100, 1),
         total = sum(count)
         ) %>%
  arrange(source, desc(freq))
