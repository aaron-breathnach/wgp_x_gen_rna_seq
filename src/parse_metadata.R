library(tidyverse)

xlsx <- "Meta data file for CHJ RNA seq.xlsx"

cols <- c("mouse", "sample_id", "training_status", "genotype")

df_1 <- readxl::read_excel(xlsx, range = "B2:E18") %>%
  setNames(cols) %>%
  fill(all_of(cols)) %>%
  mutate(stim = "Unstimulated") %>%
  select(2, 1, 4, 3, 5)

df_2 <- readxl::read_excel(xlsx, range = "B23:E35") %>%
  setNames(cols) %>%
  fill(all_of(cols)) %>%
  mutate(stim = "Poly I:C Stimulated") %>%
  select(2, 1, 4, 3, 5)

metadata <- rbind(df_1, df_2) %>%
  mutate(genotype = ifelse(genotype == "Wild-Type", "WT", "KO")) %>%
  mutate(training_status = ifelse(training_status == "Untrained", "Untrained", "Trained")) %>%
  mutate(genotype_x_training_status = paste0(genotype, "_", training_status))

write_tsv(metadata, "metadata.tsv")
