library(tidyverse)

make_counts_table <- function(out_dir = ".") {
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  edb <- EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79
  k <- AnnotationDbi::keys(edb, keytype = "TXNAME")
  tx2gene <- AnnotationDbi::select(edb, k, "GENEID", "TXNAME")
  
  files <- list.files(
    "quants",
    pattern = ".sf",
    recursive = TRUE,
    full.names = TRUE
  )
  
  txi <- tximport::tximport(
    files,
    type = "salmon",
    tx2gene = tx2gene,
    ignoreTxVersion = TRUE
  )
  
  sample_ids <- files %>%
    purrr::map(function(x) {
      str_split(x, "\\/") %>%
        unlist() %>%
        nth(2) %>%
        str_replace("_quant", "")}) %>%
    unlist()
  
  col_ord <- read_delim("metadata.tsv")$sample_id
  
  counts <- txi$counts
  colnames(counts) <- sample_ids
  counts <- counts[,col_ord]
  
  ensgene_to_symbol <- annotables::grcm38[,c(1, 3)] %>%
    filter(nchar(symbol) > 0)
  
  counts <- counts %>%
    as.data.frame() %>%
    rownames_to_column("ensgene") %>%
    inner_join(ensgene_to_symbol, by = "ensgene") %>%
    select(-ensgene) %>%
    pivot_longer(!symbol, names_to = "sample_id", values_to = "expr") %>%
    group_by(sample_id, symbol) %>%
    summarise(expr = sum(expr)) %>%
    pivot_wider(names_from = "sample_id", values_from = "expr")
  
  output <- sprintf("%s/counts.tsv", out_dir)
  write_tsv(counts, output)
  
}

make_counts_table()
