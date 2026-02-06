library(tidyverse)

main <- function(path2metadata, path2counts, out_dir = "plots") {
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  metadata <- read_delim(path2metadata) %>%
    filter(training_status == "Trained" & mouse != "M4") %>%
    column_to_rownames("sample_id") %>%
    select(genotype) %>%
    mutate_if(is.character, as.factor)
  
  cols <- c("symbol", rownames(metadata))
  
  counts <- read_delim(path2counts) %>%
    select(all_of(cols)) %>%
    column_to_rownames("symbol") %>%
    round() %>%
    mutate_if(is.numeric, as.integer)
  
  ## run DESeq2
  de <- run_deseq2(counts, metadata)
  write_tsv(de, sprintf("%s/deseq2.tsv", out_dir))
  
  ## GO enrichment analysis
  enriched_go <- run_enrichment_analysis(de, db = "GO")
  write_tsv(enriched_go, sprintf("%s/enriched_go.tsv", out_dir))
  p1 <- make_enrichment_plot(enriched_go, "GO", 25)
  
  ## KEGG enrichment analysis
  enriched_ko <- run_enrichment_analysis(de, db = "KEGG")
  write_tsv(enriched_ko, sprintf("%s/enriched_ko.tsv", out_dir))
  p2 <- make_enrichment_plot(enriched_ko, "KEGG", 25)
  
  ## PCA
  p3 <- make_pca_plot(counts, metadata)
  
  ## volcano plot
  p4 <- make_volcano_plot(de)
  
  ## save the plots
  filenames <- paste0(
    sprintf("%s/", out_dir), 
    c(
      "enrichment_plot_go.png",
      "enrichment_plot_kegg.png",
      "pca_plot.png",
      "volcano_plot.png"
    )
  )
  
  dimensions <- list(
    c(7.5, 5),
    c(7.5, 5),
    c(5, 5),
    c(5, 5)
  )
  
  plot_list <- list(p1, p2, p3, p4)
  
  purrr::map(1:4, function(x) ggsave(
    filenames[[x]],
    plot_list[[x]],
    width = dimensions[[x]][1],
    height = dimensions[[x]][2]
  ))
  
}
