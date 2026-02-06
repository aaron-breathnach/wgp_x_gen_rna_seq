#########
## PCA ##
#########

make_pca_plot <- function(counts, metadata) {
  
  pca <- prcomp(log1p(t(counts)), center = TRUE)
  
  perc <- 100 * summary(pca)$importance[2,][1:2]
  x_lab <- paste0("PC1 [", round(perc[[1]], 2), "%]")
  y_lab <- paste0("PC2 [", round(perc[[2]], 2), "%]")
  
  pca <- prcomp(log1p(t(counts)), center = TRUE)$x[, 1:2] %>%
    as.data.frame() %>%
    rownames_to_column("sample_id") %>%
    inner_join(rownames_to_column(metadata, "sample_id"), by = "sample_id")
  
  ggplot(pca, aes(x = PC1, y = PC2)) +
    geom_point(aes(fill = genotype),
               colour = "black",
               pch = 21,
               size = 3) +
    theme_classic(base_size = 12.5) +
    theme(axis.title = element_text(face = "bold"),
          legend.title = element_text(face = "bold")) +
    labs(x = x_lab, y = y_lab, fill = "Genotype")
  
}

################
## run DESeq2 ##
################

run_deseq2 <- function(counts, metadata) {
  
  obj <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData = metadata,
    design = ~ genotype
  )
  
  keep <- rowSums(BiocGenerics::counts(obj)) >= 10
  
  obj <- DESeq2::DESeq(obj[keep,])
  
  DESeq2::results(obj) %>%
    as.data.frame() %>%
    rownames_to_column("symbol") %>%
    as_tibble() %>%
    inner_join(annotables::grcm38[, 2:3], by = "symbol") %>%
    select(ncol(.), 1:ncol(.) - 1)
  
}

############################
## GO enrichment analysis ##
############################

run_enrichment_analysis <- function(de, db) {
  
  n <- ifelse(db == "KEGG", 1, 2)
  
  gene <- de %>%
    filter(padj <= 0.05) %>%
    pull(n)
  
  if (db == "KEGG") {
  
    enriched <- clusterProfiler::enrichKEGG(
      gene = gene,
      organism = "mmu"
    )
    
  } else {
    
    enriched <- clusterProfiler::enrichGO(
      gene = gene,
      OrgDb = org.Mm.eg.db::org.Mm.eg.db,
      keyType = "SYMBOL",
      ont = "BP",
      readable = TRUE
    )
    
  }
  
  enriched@result %>%
    as_tibble() %>%
    separate(GeneRatio, c("x", "y"), sep = "/") %>%
    mutate(GeneRatio = as.numeric(x) / as.numeric(y)) %>%
    dplyr::select(-c(x, y))
  
}

make_enrichment_plot <- function(enriched, db, top = 10) {
  
  p_inp <- enriched %>%
    top_n(-top, qvalue) %>%
    mutate(GO = str_trunc(paste0(ID, ": ", Description), 50))
  
  y_lab <- ifelse(db == "KEGG", "**KEGG pathway**", "**GO term**")
  
  ggplot(p_inp, aes(x = GeneRatio, y = reorder(GO, GeneRatio))) +
    geom_point(aes(fill = -log10(pvalue), size = Count), pch = 21, colour = "black") +
    scale_fill_distiller(palette = "Reds", direction = 1) +
    theme_classic() +
    theme(legend.title = ggtext::element_markdown(),
          axis.title.x = ggtext::element_markdown(),
          axis.title.y = ggtext::element_markdown()) +
    labs(x = "**Gene ratio**",
         y = y_lab,
         fill = "**-log<sub>10</sub>(*p*)**",
         size = "**Count**")
  
}

#######################
## make volcano plot ##
#######################

make_volcano_plot <- function(de, q = 0.05, pal = c("red", "blue", "grey")) {
  
  groups <- c(
    paste("\U2191", "in KO"),
    paste("\U2193", "in KO"),
    "non-significant"
  )
  
  tab <- de %>%
    drop_na() %>%
    mutate(colour = case_when(
      padj <  0.05 & log2FoldChange > +0001 ~ groups[1],
      padj <  0.05 & log2FoldChange < -0001 ~ groups[2],
      padj > 0.05 | abs(log2FoldChange) < 1 ~ groups[3]
    )) %>%
    mutate(colour = factor(colour, levels = groups))
  
  x_min <- min(tab$log2FoldChange)
  x_max <- max(tab$log2FoldChange)
  x_val <- max(abs(c(x_min, x_max)))
  x_lim <- c(-x_val, x_val)
  
  ann <- tab %>%
    drop_na() %>%
    filter(abs(log2FoldChange) > 2) %>%
    mutate(direction = ifelse(log2FoldChange > 0, "up", "down")) %>%
    top_n(-10, padj)
  
  names(pal) <- groups
  
  tmp_max_sig <- tab %>%
    filter(padj <= q)
  
  tmp_volcano <- ggplot(tab, aes(x = log2FoldChange, y = -log10(pvalue)))
  
  if (nrow(tmp_max_sig) > q) {
    
    max_sig <- filter(tmp_max_sig, pvalue == max(pvalue)) %>%
      .[[1, "pvalue"]]
    
    tmp_volcano <- tmp_volcano +
      geom_hline(yintercept = -log10(max_sig), linetype = "dashed", colour = "darkgrey")
    
  }
  
  tmp_volcano +
    geom_vline(xintercept = c(-2, -1, 1, 2), linetype = "dashed", colour = "darkgrey") +
    geom_point(aes(colour = colour)) +
    theme_bw(base_size = 12.5) +
    theme(axis.title.x = ggtext::element_markdown(),
          axis.title.y = ggtext::element_markdown(),
          axis.line = element_line(colour = "black"),
          plot.title = element_text(face = "bold", hjust = 0.5),
          panel.border = element_blank(),
          legend.position = "top") +
    scale_colour_manual(values = pal, breaks = groups[1:2]) +
    labs(x = "**log<sub>2</sub> fold change**",
         y = "**-log<sub>10</sub>(*p*)**",
         colour = NULL) +
    xlim(x_lim) +
    ggrepel::geom_text_repel(data = ann,
                             aes(
                               x = log2FoldChange,
                               y = -log10(pvalue),
                               label = symbol
                             ),
                             fontface = "italic",
                             min.segment.length = 0) +
    guides(colour = guide_legend(override.aes = list(size = 3)))
  
}
