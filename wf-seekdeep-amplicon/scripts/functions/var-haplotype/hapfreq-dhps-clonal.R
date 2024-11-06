##___compute allele frequencies, regardless source ----
# -----------------------------------------------------------------------------#

df_freqHap_All <- df_clusters_Target %>%
  select(s_Sample, contains("pos"), haplotype) %>%
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
          dhps_resistance,
          haplotype,
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
                             )
         ) %>%
  select(-count) %>%
  arrange(desc(freq))


##___compute allele frequencies, by source ----
# -----------------------------------------------------------------------------#

df_freqHap_Source <- df_clusters_Target %>%
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
          haplotype,
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
                             )
         ) %>%
  select(-count) %>%
  arrange(source, desc(freq))



##___print a message in the console ----
# -----------------------------------------------------------------------------#

# Using yellow for the border
cat("\033[1m\033[33m", "\n##############################################################", "\033[0m")

# Using magenta for the table descriptions
cat("\033[1m\033[35m", "\nData Summary:", "\033[0m")

# Using cyan for the first table description
cat("\033[1m\033[36m", "\n1. df_freqHap_All - This table shows the aggregated haplotype-frequencies across all geographical regions", "\033[0m")

# Using cyan for the second table description
cat("\033[1m\033[36m", "\n2. df_freqHap_Source - This table shows the haplotype-frequencies by geographical region", "\033[0m")

# Yellow for the border again
cat("\033[1m\033[33m", "\n##############################################################\n", "\033[0m")

