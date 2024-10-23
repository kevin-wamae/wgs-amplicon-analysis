# compute haplotype frequencies - with weighting
# ------------------------------------------------------------------------------
df_freqHap_All <- df_clusters_Target %>%
  select(s_Sample, codon_pos, haplotype, c_AveragedFrac, s_COI, source) %>%
  # sum c_AveragedFrac for identical haplotypes within each sample
  mutate(across(c(c_AveragedFrac, s_COI), as.numeric)) %>%
  group_by(s_Sample, haplotype) %>%
  mutate(c_AveragedFrac = sum(c_AveragedFrac)) %>%
  distinct(haplotype, .keep_all = TRUE) %>%
  ungroup() %>%
  # compute the sample size
  mutate(pop_size = n_distinct(s_Sample)) %>%
  # calculate haplotype frequencies across all samples
  group_by(haplotype) %>%
  summarize(freq = round(sum(c_AveragedFrac) / first(pop_size) * 100, 2), .groups = 'drop')



# haplotype frequencies - with weighting
# ------------------------------------------------------------------------------
df_freqHap_Source <- df_clusters_Target %>%
  select(s_Sample, codon_pos, haplotype, c_AveragedFrac, s_COI, source) %>%
  # sum c_AveragedFrac for identical haplotypes
  mutate(across(c(c_AveragedFrac, s_COI), as.numeric)) %>%
  mutate(c_AveragedFrac = sum(c_AveragedFrac), .by = c(s_Sample, haplotype)) %>%
  # drop duplicate haplotypes, within samples
  group_by(s_Sample, haplotype) %>% distinct(haplotype, .keep_all = TRUE) %>% ungroup() %>%
  # nest data by source and codon
  tidyr::nest(country_data = -source) %>%
  # compute sample size
  nplyr::nest_mutate(country_data,
                     pop_size = n_distinct(s_Sample)) %>%
  # compute allele frequencies
  nplyr::nest_group_by(.nest_data = country_data,
                       haplotype) %>%
  nplyr::nest_mutate(.nest_data = country_data,
                     freq = round(sum(c_AveragedFrac/pop_size),3)) %>%
  tidyr::unnest(cols = c(country_data)) %>%
  # retain allele frequencies
  distinct(source, haplotype, freq)



##___print a message in the console ----
# -----------------------------------------------------------------------------#

cat("\033[1m", "\n##############################################################", "\033[0m")

cat("\033[1m", "\n1. df_freqSNP_All    - This table shows the aggregated SNP frequencies across all geographical regions", "\033[0m")

cat("\033[1m", "\n2. df_freqSNP_Source - This table shows the SNP frequencies by geographical region", "\033[0m")  # No newline before last line

cat("\033[1m", "\n##############################################################\n", "\033[0m")

