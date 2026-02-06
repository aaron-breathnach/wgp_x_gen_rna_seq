library(tidyverse)

get_top_table <- function(n, fit) {
  
  cols <- c("intercept",
            "genotype",
            "training_status",
            "genotype_x_training_status")
  
  coefficient <- cols[n]
  
  limma::topTable(fit, coef = n, sort.by = "P", n = Inf) %>%
    rownames_to_column("gene") %>%
    as_tibble() %>%
    mutate(coefficient = coefficient) %>%
    select(ncol(.), 1:ncol(.)-1)
  
}

make_training_status_volcano_plot <- function(res) {
  
  dat <- res %>%
    filter(coefficient == "training_status") %>%
    select(gene, logFC, P.Value, adj.P.Val) %>%
    mutate(colour = case_when(
      sign(logFC) == -1 & adj.P.Val <= 0.05 ~ "\U02193 WGP",
      adj.P.Val > 0.05 ~ "n.s.",
      sign(logFC) == +1 & adj.P.Val <= 0.05 ~ "\U02191 WGP"
    ))
  
  y_intercept <- dat %>%
    filter(adj.P.Val <= 0.05) %>%
    filter(P.Value == max(P.Value)) %>%
    pull(P.Value) %>%
    unique()
  
  pal <- c("salmon", "steelblue", "darkgrey")
  
  ann <- dat %>%
    top_n(-20, P.Value)
  
  ggplot(dat, aes(x = logFC, y = -log10(P.Value))) +
    geom_hline(yintercept = -log10(y_intercept), colour = "lightgrey", linetype = "dashed") +
    geom_vline(xintercept = c(-2, -1, 1, 2), colour = "lightgrey", linetype = "dashed") +
    geom_point(aes(colour = colour, fill = colour), pch = 21, alpha = 0.75) +
    ggrepel::geom_text_repel(data = ann,
                             aes(x = logFC, y = -log10(P.Value), label = gene),
                             min.segment.length = 1,
                             size = 3,
                             fontface = "italic") +
    theme_classic() +
    theme(axis.title.x = ggtext::element_markdown(),
          axis.title.y = ggtext::element_markdown(),
          plot.title = element_text(face = "bold", hjust = 0.5),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          legend.position = "top",
          legend.title = element_blank(),
          legend.text = ggtext::element_markdown()) +
    scale_colour_manual(values = pal) +
    scale_fill_manual(values = pal) +
    labs(x = "**log<sub>2</sub> fold change**",
         y = "**-log<sub>10</sub>(*p*)**")
  
}

make_interaction_boxplot <- function(counts, res, metadata) {
  
  if (!dir.exists("plots")) dir.create("plots")
  
  sig <- res %>%
    filter(adj.P.Val <= 0.05 & coefficient == "genotype_x_training_status") %>%
    pull(gene)
  
  p_inp <- counts[sig,] %>%
    rownames_to_column("gene") %>%
    filter(gene %in% res$gene) %>%
    pivot_longer(!gene, names_to = "sample_id", values_to = "expr") %>%
    inner_join(metadata, by = "sample_id")
  
  ggplot(p_inp, aes(x = genotype, y = expr)) +
    facet_wrap(~ gene, scales = "free_y") +
    geom_boxplot(aes(colour = training_status), outlier.shape = NA) +
    geom_point(aes(colour = training_status, fill = training_status),
               position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.125),
               alpha = 0.5,
               pch = 21) +
    theme_bw(base_size = 12.5) +
    theme(axis.title = element_text(face = "bold"),
          legend.title = element_text(face = "bold"),
          panel.grid = element_blank(),
          strip.text = element_text(face = "bold.italic")) +
    labs(x = "Genotype",
         y = "Abundance",
         colour = "Training status",
         fill = "Training status")
  
}

run_limma_voom <- function(counts, metadata, stimulation) {
  
  dge <- counts %>%
    edgeR::DGEList() %>%
    edgeR::calcNormFactors()
  
  dge <- dge[which(apply(edgeR::cpm(dge), 1, max) < 1),]
  
  genotype <- metadata$genotype
  training_status <- metadata$training_status
  design <- model.matrix(~ genotype + training_status + genotype:training_status)
  
  obj <- limma::voom(dge, design, plot = FALSE)
  
  dup_cor <- limma::duplicateCorrelation(obj, design, block = metadata$mouse)
  
  fit <- limma::lmFit(
    obj,
    design,
    block = metadata$mouse,
    correlation = dup_cor$consensus
  ) %>%
    limma::eBayes()
  
  purrr::map(2:4, function(x) get_top_table(x, fit)) %>%
    bind_rows()
  
}

get_gsea_inputs <- function(de, gs_cat = "H", gs_sub = NULL) {
  
  dat <- de %>%
    select(gene, logFC) %>%
    setNames(c("gene_id", "logFC"))
  
  gene_list <- dat$logFC
  names(gene_list) <- dat$gene_id
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  gene_set <- msigdbr::msigdbr(species = "mouse",
                               category = gs_cat,
                               subcategory = gs_sub) %>%
    dplyr::select(3:6)
  
  output <- list("gene_list" = gene_list, "gene_set" = gene_set)
  
  return(output)
  
}

.run_gsea <- function(de, gs_cat = "H", gs_sub = NULL) {
  
  gsea_inputs <- get_gsea_inputs(de, gs_cat, gs_sub)
  
  gene_list   <- gsea_inputs$gene_list
  gene_set    <- gsea_inputs$gene_set
  
  clusterProfiler::GSEA(gene_list,
                        TERM2GENE = gene_set[, c(1, 2)],
                        pvalueCutoff = 1) %>%
    slot("result") %>%
    as_tibble()
  
}

run_gsea <- function(coef, de, gs_cat = "H", gs_sub = NULL) {
  
  de <- filter(de, coefficient == coef)
  
  .run_gsea(de, gs_cat, gs_sub) %>%
    mutate(coefficient = coef) %>%
    select(ncol(.), 1:ncol(.) - 1)
  
}

write_sig_res <- function() {
  
  read_delim("tables/limma_voom.unstimulated.tsv") %>%
    filter(adj.P.Val <= 0.05) %>%
    write_tsv("tables/limma_voom.sig_res.unstimulated.tsv")
  
  read_delim("tables/limma_voom.poly_i_c_stimulated.tsv") %>%
    filter(adj.P.Val <= 0.05) %>%
    write_tsv("tables/limma_voom.sig_res.poly_i_c_stimulated.tsv")
  
  read_delim("tables/limma_voom.poly_i_c_stimulated.tsv") %>%
    filter(adj.P.Val <= 0.05 & coefficient == "genotype_x_training_status")
  
}

wrapper <- function(stimulation) {
  
  purrr::map(c("plots", "tables"), function(x) dir.create(x, FALSE, FALSE))
  
  suffix <- str_to_lower(str_replace_all(stimulation, "[^[:alnum:]]", "_"))
  
  filenames <- c(
    sprintf("tables/limma_voom.%s.tsv", suffix),
    sprintf("tables/gsea.%s.tsv", suffix),
    sprintf("plots/interaction_effects.%s.png", suffix),
    sprintf("plots/training_status.%s.png", suffix)
  )
  
  metadata <- read_delim("data/metadata.tsv") %>%
    filter(stim == stimulation & mouse != "M4")
  
  sample_ids <- metadata$sample_id
  
  counts <- read.delim("data/counts.tsv", row.names = 1)
  counts <- counts[, sample_ids]
  
  de <- run_limma_voom(counts, metadata, stimulation)
  write_tsv(de, filenames[1])
  
  coefficients <- c("genotype", "training_status", "genotype_x_training_status")
  
  gsea <- purrr::map(coefficients, function(x) run_gsea(x, de)) %>%
    bind_rows() %>%
    suppressWarnings()
  write_tsv(gsea, filenames[2])
  
  p1 <- make_interaction_boxplot(counts, de, metadata)
  ggsave(filenames[3], p1, width = 12.5, height = 7.5)
  
  p2 <- make_training_status_volcano_plot(de)
  ggsave(filenames[4], p2, width = 5, height = 5)
  
}

wrapper("Unstimulated")
wrapper("Poly I:C Stimulated")
write_sig_res()
