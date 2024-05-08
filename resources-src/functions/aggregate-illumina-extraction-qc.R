## __b. import fastq extraction reports (by target) ----
# =============================================================================#

### ____i. import report ----
# -----------------------------------------------------------------------------#

raw_extProfile <- read_tsv(paste0(PATH_STUDY, PATH_RUN, PATH_DATE, "reports/allExtractionProfile.tab.txt"),
                       show_col_types = FALSE) %>%
  mutate_all(list(~str_replace(., "\\(.+\\)", ""))) %>%       # remove all characters in braces
  mutate_at(.vars = 3:ncol(.), .funs = as.numeric) %>%        # to numeric
  mutate(inputName = str_remove(inputName, "_R1|_S.+_L.+"))   # drop read or miseq-lane suffix



### ____ii. qc summary - by target ----
# -----------------------------------------------------------------------------#

raw_extProfileTarget <- raw_extProfile %>%
  mutate(name = str_remove(name, pattern = "MID.+")) %>% # drop MID substring
  summarise(
            # generate sum across all columns by fastq and target
            across(totalMatching:last_col(), ~sum(.x, na.rm = TRUE)),
            .by = c(inputName, name)
            ) %>%
  # select relevant columns
  select(inputName, name, totalMatching, good) %>%
  # summary of totals, all data
  mutate(
         totalAll = sum(totalMatching),
         totalAllGood = sum(good),
         totalAllGoodPerc = round(totalAllGood / totalAll, 2),
         totalAllGood = base::format(totalAllGood, big.mark = ","),
         totalAllGoodPerc = paste0(totalAllGood, " [", totalAllGoodPerc, "]"),
         ) %>%
  # summary of totals by fastq
  mutate(
         totalByFastq = sum(totalMatching),
         totalByFastqGood = sum(good),
         totalByFastqGoodPerc = round(totalByFastqGood / totalByFastq, 2),
         totalByFastq = base::format(totalByFastq, big.mark = ","),
         totalByFastqGoodPerc = paste0(totalByFastq, " [", totalByFastqGoodPerc, "]"),
         .by = inputName
         ) %>%
  # summary of totals by fastq and target
  mutate(
         totalTargetGoodPerc = round(good / totalMatching, 2),
         totalTargetMatching = base::format(totalMatching, big.mark = ","),
         totalTargetGoodPerc = paste0(totalTargetMatching, " [", totalTargetGoodPerc, "]"),
         ) %>%
  # format: long to wide
  pivot_wider(
              id_cols = c(inputName, totalAllGoodPerc, totalByFastqGoodPerc),
              names_from = "name",
              values_from = "totalTargetGoodPerc"
              )
