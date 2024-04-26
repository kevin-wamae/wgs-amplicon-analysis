##___compute allele frequencies, regardless source ----
# -----------------------------------------------------------------------------#

df_freqHap_MDR1_All <- df_clusters_MDR1 %>%
  reframe(
          haplotype = paste(sort(unique(haplotype)), collapse = ","),
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
         total = sum(count),
         ) %>%
  select(-count) %>%
  arrange(desc(freq))


##___compute allele frequencies, by source ----
# -----------------------------------------------------------------------------#

df_freqHap_MDR1_Source <- df_clusters_MDR1 %>%
  reframe(
          source,
          haplotype = paste(sort(unique(haplotype)), collapse = ","),
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
