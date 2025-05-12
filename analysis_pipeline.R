
getwd()
setwd("/Users/sophiewilliams/fyp2")

# Load required libraries
library(DESeq2)
library(EnhancedVolcano)
library(pheatmap)
library(clusterProfiler)
library(org.At.tair.db)
library(pathview)
library(fgsea)
library(tibble)
library(dplyr)

# === Load and prepare data ===
counts <- read.delim("counts_table_2.txt", header = FALSE, skip = 1,
                     col.names = c("Gene_id", "Col0_day_1", "Col0_day_2", "Col0_night_1", "Col0_night_2",
                                   "jmj13_day_1", "jmj13_day_2", "jmj13_night_1", "jmj13_night_2"), row.names = 1)

sample_info <- data.frame(
  row.names = colnames(counts),
  genotype = rep(c("Col0", "jmj13"), each = 4),
  condition = rep(c("day", "day", "night", "night"), times = 2)
)

sample_info$genotype <- factor(sample_info$genotype)
sample_info$condition <- factor(sample_info$condition)

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = sample_info,
                              design = ~ genotype + condition + genotype:condition)

dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)

rld <- rlog(dds, blind = FALSE)

# === Retrieve comparisons ===
res_day <- results(dds, name = "genotype_jmj13_vs_Col0")
res_interaction <- results(dds, name = "genotypejmj13.conditionnight")
res_condition <- results(dds, name = "condition_night_vs_day")

# 1. jmj13 vs Col0 during day
res_1 <- res_day

# 2. jmj13 vs Col0 during night = main + interaction
res_2 <- res_day
res_2$log2FoldChange <- res_day$log2FoldChange + res_interaction$log2FoldChange
res_2$pvalue <- pmax(res_day$pvalue, res_interaction$pvalue)

# 3. Col0 night vs day = main effect
res_3 <- res_condition

# 4. jmj13 night vs day = main + interaction
res_4 <- res_condition
res_4$log2FoldChange <- res_condition$log2FoldChange + res_interaction$log2FoldChange
res_4$pvalue <- pmax(res_condition$pvalue, res_interaction$pvalue)

results_list <- list(res_1, res_2, res_3, res_4)
names(results_list) <- c("jmj13_vs_Col0_day", "jmj13_vs_Col0_night", "Col0_night_vs_day", "jmj13_night_vs_day")

# === Define plotting & enrichment functions ===
volcano_plot <- function(res, title) {
  EnhancedVolcano(res,
                  lab = rownames(res),
                  x = 'log2FoldChange',
                  y = 'pvalue',
                  title = title)
}

heatmap_top <- function(res, rld, title) {
  top_genes <- head(order(res$pvalue), 50)
  mat <- assay(rld)[top_genes, ]
  mat <- mat - rowMeans(mat)
  pheatmap(mat, annotation_col = as.data.frame(colData(dds)), main = title)
}

go_enrich <- function(res, title) {
  sig_genes <- rownames(res)[which(res$padj < 0.05)]
  ego <- enrichGO(gene = sig_genes,
                  OrgDb = org.At.tair.db,
                  keyType = "TAIR",
                  ont = "BP",
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.05,
                  readable = TRUE)
  dotplot(ego, showCategory = 20, title = paste("GO Enrichment -", title))
}

kegg_enrich <- function(res, title) {
  sig_genes <- rownames(res)[which(res$padj < 0.05)]
  genes_entrez <- mapIds(org.At.tair.db, keys = sig_genes, keytype = "TAIR", column = "ENTREZID")
  kegg <- enrichKEGG(gene = na.omit(genes_entrez), organism = 'ath')
  dotplot(kegg, showCategory = 20, title = paste("KEGG Pathways -", title))
}

gsea_run <- function(res, title) {
  gene_ranks <- res$log2FoldChange
  names(gene_ranks) <- mapIds(org.At.tair.db, keys = rownames(res), keytype = "TAIR", column = "ENTREZID")
  gene_ranks <- sort(na.omit(gene_ranks), decreasing = TRUE)
  gsea <- gseGO(geneList = gene_ranks,
                OrgDb = org.At.tair.db,
                ont = "BP",
                keyType = "ENTREZID",
                pAdjustMethod = "BH",
                verbose = FALSE)
  gseaplot2(gsea, geneSetID = 1, title = paste("GSEA -", title))
}

# === Run all plots for each comparison ===
for (i in seq_along(results_list)) {
  name <- names(results_list)[i]
  res <- results_list[[i]]
  cat("\nRunning analyses for:", name, "\n")
  volcano_plot(res, title = name)
  heatmap_top(res, rld, title = name)
  go_enrich(res, title = name)
  kegg_enrich(res, title = name)
  gsea_run(res, title = name)
}

#fix KEGG error 
kegg_enrich <- function(res, title) {
  sig_genes <- rownames(res)[which(res$padj < 0.05)]
  genes_entrez <- mapIds(org.At.tair.db,
                         keys = sig_genes,
                         keytype = "TAIR",
                         column = "ENTREZID",
                         multiVals = "first")
  
  genes_entrez <- na.omit(genes_entrez)
  
  if (length(genes_entrez) == 0) {
    cat("No significant genes mapped to Entrez IDs for KEGG enrichment -", title, "\n")
    return(NULL)
  }
  
  kegg <- enrichKEGG(gene = genes_entrez, organism = 'ath')
  
  if (is.null(kegg) || nrow(kegg) == 0) {
    cat("No KEGG pathways enriched for -", title, "\n")
    return(NULL)
  }
  
  dotplot(kegg, showCategory = 20, title = paste("KEGG Pathways -", title))
}

#fix KEGG take 2 
kegg_enrich <- function(res, title) {
  sig_genes <- rownames(res)[which(res$padj < 0.05)]
  genes_entrez <- mapIds(org.At.tair.db,
                         keys = sig_genes,
                         keytype = "TAIR",
                         column = "ENTREZID",
                         multiVals = "first")
  
  genes_entrez <- na.omit(genes_entrez)
  genes_entrez <- as.character(genes_entrez)  # Ensure proper format
  
  if (length(genes_entrez) == 0) {
    cat("No significant genes mapped to Entrez IDs for KEGG enrichment -", title, "\n")
    return(NULL)
  }
  
  kegg <- enrichKEGG(gene = genes_entrez, organism = 'ath')
  
  if (is.null(kegg) || nrow(kegg) == 0) {
    cat("No KEGG pathways enriched for -", title, "\n")
    return(NULL)
  }
  
  print(dotplot(kegg, showCategory = 20, title = paste("KEGG Pathways -", title)))
}

#fix GSEA error 
gsea_run <- function(res, title) {
  gene_ids <- mapIds(org.At.tair.db,
                     keys = rownames(res),
                     keytype = "TAIR",
                     column = "ENTREZID",
                     multiVals = "first")
  
  gene_ranks <- res$log2FoldChange
  names(gene_ranks) <- gene_ids
  
  # Remove genes with NA names (Entrez IDs)
  gene_ranks <- gene_ranks[!is.na(names(gene_ranks))]
  
  gene_ranks <- sort(gene_ranks, decreasing = TRUE)
  
  if (length(gene_ranks) == 0) {
    cat("No valid gene ranks for GSEA -", title, "\n")
    return(NULL)
  }
  
  gsea <- gseGO(geneList = gene_ranks,
                OrgDb = org.At.tair.db,
                ont = "BP",
                keyType = "ENTREZID",
                pAdjustMethod = "BH",
                verbose = FALSE)
  
  print(gseaplot2(gsea, geneSetID = 1, title = paste("GSEA -", title)))
}

#save all plots 
#volcano plot 
volcano_plot <- function(res, title) {
  p <- EnhancedVolcano(res,
                       lab = rownames(res),
                       x = 'log2FoldChange',
                       y = 'pvalue',
                       title = title)
  ggsave(paste0(title, "_volcano.png"), p, width = 8, height = 6)
  print(p)
}

#heatmap of top genes 
heatmap_top <- function(res, rld, title) {
  top_genes <- head(order(res$pvalue), 50)
  mat <- assay(rld)[top_genes, ]
  mat <- mat - rowMeans(mat)
  
  pdf(paste0(title, "_heatmap.pdf"), width = 8, height = 8)
  pheatmap(mat, annotation_col = as.data.frame(colData(dds)), main = title)
  dev.off()
}

#GO enrichment 
go_enrich <- function(res, title) {
  sig_genes <- rownames(res)[which(res$padj < 0.05)]
  ego <- enrichGO(gene = sig_genes,
                  OrgDb = org.At.tair.db,
                  keyType = "TAIR",
                  ont = "BP",
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.05,
                  readable = TRUE)
  
  if (is.null(ego) || nrow(ego) == 0) {
    cat("No GO terms enriched for -", title, "\n")
    return(NULL)
  }
  
  p <- dotplot(ego, showCategory = 20, title = paste("GO Enrichment -", title))
  ggsave(paste0(title, "_GO_enrichment.png"), p, width = 8, height = 6)
  print(p)
}

#KEGG enrichment 
kegg_enrich <- function(res, title) {
  sig_genes <- rownames(res)[which(res$padj < 0.05)]
  genes_entrez <- mapIds(org.At.tair.db,
                         keys = sig_genes,
                         keytype = "TAIR",
                         column = "ENTREZID",
                         multiVals = "first")
  
  genes_entrez <- na.omit(genes_entrez)
  genes_entrez <- as.character(genes_entrez)
  
  if (length(genes_entrez) == 0) {
    cat("No significant genes mapped to Entrez IDs for KEGG enrichment -", title, "\n")
    return(NULL)
  }
  
  kegg <- enrichKEGG(gene = genes_entrez, organism = 'ath')
  
  if (is.null(kegg) || nrow(kegg) == 0) {
    cat("No KEGG pathways enriched for -", title, "\n")
    return(NULL)
  }
  
  p <- dotplot(kegg, showCategory = 20, title = paste("KEGG Pathways -", title))
  ggsave(paste0(title, "_KEGG_enrichment.png"), p, width = 8, height = 6)
  print(p)
}

#GSEA Plot 
gsea_run <- function(res, title) {
  gene_ids <- mapIds(org.At.tair.db,
                     keys = rownames(res),
                     keytype = "TAIR",
                     column = "ENTREZID",
                     multiVals = "first")
  
  gene_ranks <- res$log2FoldChange
  names(gene_ranks) <- gene_ids
  gene_ranks <- gene_ranks[!is.na(names(gene_ranks))]
  gene_ranks <- sort(gene_ranks, decreasing = TRUE)
  
  if (length(gene_ranks) == 0) {
    cat("No valid gene ranks for GSEA -", title, "\n")
    return(NULL)
  }
  
  gsea <- gseGO(geneList = gene_ranks,
                OrgDb = org.At.tair.db,
                ont = "BP",
                keyType = "ENTREZID",
                pAdjustMethod = "BH",
                verbose = FALSE)
  
  if (is.null(gsea) || nrow(gsea) == 0) {
    cat("No GSEA results -", title, "\n")
    return(NULL)
  }
  
  p <- gseaplot2(gsea, geneSetID = 1, title = paste("GSEA -", title))
  ggsave(paste0(title, "_GSEA.png"), p, width = 8, height = 6)
  print(p)
}


getwd()


print(file.path(getwd(), paste0(title, "_volcano.png")))
p <- EnhancedVolcano(...)
ggsave(paste0(title, "_volcano.png"), plot = p)

#trying again as error came up 
volcano_plot <- function(res, title) {
  p <- EnhancedVolcano(
    res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'pvalue',
    title = title
  )
  file_name <- paste0(title, "_volcano.png")
  file_path <- file.path(getwd(), file_name)
  ggsave(file_path, plot = p, width = 8, height = 6)
  cat("Saved volcano plot to:", file_path, "\n")
}
heatmap_top <- function(res, rld, title) {
  top_genes <- head(order(res$pvalue), 50)
  mat <- assay(rld)[top_genes, ]
  mat <- mat - rowMeans(mat)
  file_path <- file.path(getwd(), paste0(title, "_heatmap.pdf"))
  pdf(file_path, width = 8, height = 10)
  pheatmap(mat, annotation_col = as.data.frame(colData(rld)), main = title)
  dev.off()
  cat("Saved heatmap to:", file_path, "\n")
}
go_enrich <- function(res, title) {
  sig_genes <- rownames(res)[which(res$padj < 0.05)]
  ego <- enrichGO(
    gene = sig_genes,
    OrgDb = org.At.tair.db,
    keyType = "TAIR",
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  
  if (is.null(ego) || nrow(ego) == 0) {
    cat("No GO enrichment found for -", title, "\n")
    return(NULL)
  }
  
  p <- dotplot(ego, showCategory = 20, title = paste("GO Enrichment -", title))
  file_path <- file.path(getwd(), paste0(title, "_GO_enrichment.png"))
  ggsave(file_path, plot = p, width = 8, height = 6)
  cat("Saved GO enrichment plot to:", file_path, "\n")
}
kegg_enrich <- function(res, title) {
  sig_genes <- rownames(res)[which(res$padj < 0.05)]
  genes_entrez <- mapIds(
    org.At.tair.db,
    keys = sig_genes,
    keytype = "TAIR",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  genes_entrez <- na.omit(genes_entrez)
  if (length(genes_entrez) == 0) {
    cat("No significant genes mapped to Entrez IDs for KEGG enrichment -", title, "\n")
    return(NULL)
  }
  
  kegg <- enrichKEGG(gene = genes_entrez, organism = 'ath')
  
  if (is.null(kegg) || nrow(kegg) == 0) {
    cat("No KEGG pathways enriched for -", title, "\n")
    return(NULL)
  }
  
  p <- dotplot(kegg, showCategory = 20, title = paste("KEGG Pathways -", title))
  file_path <- file.path(getwd(), paste0(title, "_KEGG.png"))
  ggsave(file_path, plot = p, width = 8, height = 6)
  cat("Saved KEGG enrichment plot to:", file_path, "\n")
}
gsea_run <- function(res, title) {
  gene_ranks <- res$log2FoldChange
  names(gene_ranks) <- mapIds(org.At.tair.db, keys = rownames(res), keytype = "TAIR", column = "ENTREZID")
  gene_ranks <- sort(na.omit(gene_ranks), decreasing = TRUE)
  
  if (any(is.na(names(gene_ranks)))) {
    cat("GSEA: Some genes could not be mapped to Entrez IDs for -", title, "\n")
    gene_ranks <- gene_ranks[!is.na(names(gene_ranks))]
  }
  
  if (length(gene_ranks) == 0) {
    cat("GSEA: No valid gene ranks for -", title, "\n")
    return(NULL)
  }
  
  gsea <- gseGO(
    geneList = gene_ranks,
    OrgDb = org.At.tair.db,
    ont = "BP",
    keyType = "ENTREZID",
    pAdjustMethod = "BH",
    verbose = FALSE
  )
  
  if (is.null(gsea) || nrow(gsea) == 0) {
    cat("No GSEA enrichment found for -", title, "\n")
    return(NULL)
  }
  
  p <- gseaplot2(gsea, geneSetID = 1, title = paste("GSEA -", title))
  file_path <- file.path(getwd(), paste0(title, "_GSEA.png"))
  ggsave(file_path, plot = p, width = 8, height = 6)
  cat("Saved GSEA plot to:", file_path, "\n")
}
for (i in seq_along(results_list)) {
  name <- names(results_list)[i]
  res <- results_list[[i]]
  cat("\nRunning analyses for:", name, "\n")
  volcano_plot(res, title = name)
  heatmap_top(res, rld, title = name)
  go_enrich(res, title = name)
  kegg_enrich(res, title = name)
  gsea_run(res, title = name)
}

#fix gsea 
gsea_run <- function(res, title) {
  gene_ranks <- res$log2FoldChange
  names(gene_ranks) <- mapIds(org.At.tair.db, keys = rownames(res), keytype = "TAIR", column = "ENTREZID")
  gene_ranks <- sort(na.omit(gene_ranks), decreasing = TRUE)
  
  if (any(is.na(names(gene_ranks)))) {
    cat("GSEA: Some genes could not be mapped to Entrez IDs for -", title, "\n")
    gene_ranks <- gene_ranks[!is.na(names(gene_ranks))]
  }
  
  if (length(gene_ranks) == 0) {
    cat("GSEA: No valid gene ranks for -", title, "\n")
    return(NULL)
  }
  
  # Perform GSEA using clusterProfiler
  gsea <- gseGO(
    geneList = gene_ranks,
    OrgDb = org.At.tair.db,
    ont = "BP",
    keyType = "ENTREZID",
    pAdjustMethod = "BH",
    verbose = FALSE
  )
  
  if (is.null(gsea) || nrow(gsea) == 0) {
    cat("No GSEA enrichment found for -", title, "\n")
    return(NULL)
  }
  
  # Use gseaplot2 from clusterProfiler
  p <- gseaplot(gsea, geneSetID = 1, title = paste("GSEA -", title))  # gseaplot instead of gseaplot2
  file_path <- file.path(getwd(), paste0(title, "_GSEA.png"))
  ggsave(file_path, plot = p, width = 8, height = 6)
  cat("Saved GSEA plot to:", file_path, "\n")
}

file_path <- file.path(getwd(), paste0(title, "_GSEA.png"))
cat("File will be saved to:", file_path, "\n")
ggsave(file_path, plot = p, width = 8, height = 6)

cat("Title is:", title, "\n")
file_path <- file.path(getwd(), paste0(title, "_GSEA.png"))

title <- "jmj13_vs_Col0_day"  # Example title, ensure this is being passed as a string
file_path <- file.path(getwd(), paste0(title, "_GSEA.png"))
cat("File will be saved to:", file_path, "\n")

# Generate the GSEA plot
p <- gseaplot(gsea, geneSetID = 1, title = paste("GSEA -", title))

# Save the plot to the specified file path
ggsave(file_path, plot = p, width = 8, height = 6)

# Confirm the plot has been saved
cat("Saved GSEA plot to:", file_path, "\n")

#fix gsea again
gsea_run <- function(res, title) {
  # Create a ranked list of genes based on log2FoldChange
  gene_ranks <- res$log2FoldChange
  names(gene_ranks) <- mapIds(org.At.tair.db, keys = rownames(res), keytype = "TAIR", column = "ENTREZID")
  
  # Sort by log2FoldChange
  gene_ranks <- sort(na.omit(gene_ranks), decreasing = TRUE)
  
  if (any(is.na(names(gene_ranks)))) {
    cat("GSEA: Some genes could not be mapped to Entrez IDs for -", title, "\n")
    gene_ranks <- gene_ranks[!is.na(names(gene_ranks))]
  }
  
  if (length(gene_ranks) == 0) {
    cat("GSEA: No valid gene ranks for -", title, "\n")
    return(NULL)
  }
  
  # Perform GSEA using clusterProfiler
  gsea <- gseGO(
    geneList = gene_ranks,
    OrgDb = org.At.tair.db,
    ont = "BP",
    keyType = "ENTREZID",
    pAdjustMethod = "BH",
    verbose = FALSE
  )
  
  # Check if gsea object is valid
  if (is.null(gsea) || nrow(gsea) == 0) {
    cat("No GSEA enrichment found for -", title, "\n")
    return(NULL)
  }
  
  # If successful, create the plot
  p <- gseaplot(gsea, geneSetID = 1, title = paste("GSEA -", title))
  
  # Save the plot to the file path
  file_path <- file.path(getwd(), paste0(title, "_GSEA.png"))
  ggsave(file_path, plot = p, width = 8, height = 6)
  
  cat("Saved GSEA plot to:", file_path, "\n")
}
#code for gsea_run 
# Ensure that required libraries are loaded
library(clusterProfiler)  # For GSEA
library(org.At.tair.db)   # For Arabidopsis annotation
library(AnnotationDbi)    # For mapping gene IDs

# Define the GSEA run function
gsea_run <- function(res, title) {
  # Create a ranked list of genes based on log2FoldChange
  gene_ranks <- res$log2FoldChange
  names(gene_ranks) <- mapIds(org.At.tair.db, keys = rownames(res), keytype = "TAIR", column = "ENTREZID")
  
  # Sort by log2FoldChange
  gene_ranks <- sort(na.omit(gene_ranks), decreasing = TRUE)
  
  if (any(is.na(names(gene_ranks)))) {
    cat("GSEA: Some genes could not be mapped to Entrez IDs for -", title, "\n")
    gene_ranks <- gene_ranks[!is.na(names(gene_ranks))]
  }
  
  if (length(gene_ranks) == 0) {
    cat("GSEA: No valid gene ranks for -", title, "\n")
    return(NULL)
  }
  
  # Perform GSEA using clusterProfiler
  gsea <- gseGO(
    geneList = gene_ranks,
    OrgDb = org.At.tair.db,
    ont = "BP",
    keyType = "ENTREZID",
    pAdjustMethod = "BH",
    verbose = FALSE
  )
  
  # Check if gsea object is valid
  if (is.null(gsea) || nrow(gsea) == 0) {
    cat("No GSEA enrichment found for -", title, "\n")
    return(NULL)
  }
  
  # If successful, create the plot
  p <- gseaplot(gsea, geneSetID = 1, title = paste("GSEA -", title))
  
  # Save the plot to the file path
  file_path <- file.path(getwd(), paste0(title, "_GSEA.png"))
  ggsave(file_path, plot = p, width = 8, height = 6)
  
  cat("Saved GSEA plot to:", file_path, "\n")
}

# Example: Run GSEA for the comparison 'jmj13_vs_Col0_day'
title <- "jmj13_vs_Col0_day"
res <- results_list[["jmj13_vs_Col0_day"]]  # Assuming 'results_list' contains the DESeq2 results

# Call the gsea_run function for the selected comparison
gsea_run(res, title)





########code for all required plots 
# Ensure necessary libraries are loaded
library(clusterProfiler)  # For GO enrichment, GSEA
library(org.At.tair.db)   # For Arabidopsis annotation
library(AnnotationDbi)    # For mapping gene IDs
library(EnhancedVolcano)  # For volcano plots
library(DESeq2)           # For DESeq2 results
library(pheatmap)         # For heatmaps

# Define the functions for GO enrichment, Volcano, Heatmap, KEGG enrichment, and GSEA
# (Include all functions you have already defined: volcano_plot, go_enrich, heatmap_top, kegg_enrich, gsea_run)

# Run analyses for each comparison and save plots

# jmj13 vs Col0 day
title <- "jmj13_vs_Col0_day"
res_day <- results_list[["jmj13_vs_Col0_day"]]  # DESeq2 results for jmj13 vs Col0 day

# Volcano plot for jmj13 vs Col0 day
volcano_plot(res_day, title = title)
# GO enrichment for jmj13 vs Col0 day
go_enrich(res_day, title = title)
# Heatmap for jmj13 vs Col0 day
heatmap_top(res_day, rld, title = title)
# KEGG pathway for jmj13 vs Col0 day
kegg_enrich(res_day, title = title)
# GSEA for jmj13 vs Col0 day
gsea_run(res_day, title = title)

# jmj13 vs Col0 night
title <- "jmj13_vs_Col0_night"
res_night <- results_list[["jmj13_vs_Col0_night"]]  # DESeq2 results for jmj13 vs Col0 night

# Volcano plot for jmj13 vs Col0 night
volcano_plot(res_night, title = title)
# GO enrichment for jmj13 vs Col0 night
go_enrich(res_night, title = title)
# KEGG pathway for jmj13 vs Col0 night
kegg_enrich(res_night, title = title)
# GSEA for jmj13 vs Col0 night
gsea_run(res_night, title = title)

# jmj13 night vs jmj13 day
title <- "jmj13_night_vs_jmj13_day"
res_night_day <- results_list[["jmj13_night_vs_jmj13_day"]]  # DESeq2 results for jmj13 night vs jmj13 day

# Volcano plot for jmj13 night vs jmj13 day
volcano_plot(res_night_day, title = title)
# GO enrichment for jmj13 night vs jmj13 day
go_enrich(res_night_day, title = title)
# KEGG pathway for jmj13 night vs jmj13 day
kegg_enrich(res_night_day, title = title)
# GSEA for jmj13 night vs jmj13 day
gsea_run(res_night_day, title = title)

# col0 night vs col0 day
title <- "col0_night_vs_col0_day"
res_col0_night_day <- results_list[["col0_night_vs_col0_day"]]  # DESeq2 results for col0 night vs col0 day

# Volcano plot for col0 night vs col0 day
volcano_plot(res_col0_night_day, title = title)
# GO enrichment for col0 night vs col0 day
go_enrich(res_col0_night_day, title = title)
# KEGG pathway for col0 night vs col0 day
kegg_enrich(res_col0_night_day, title = title)
# GSEA for col0 night vs col0 day
gsea_run(res_col0_night_day, title = title)

# KEGG pathway for all 4 comparisons

# jmj13 vs Col0 day
kegg_enrich(res_day, title = "jmj13_vs_Col0_day")
# jmj13 vs Col0 night
kegg_enrich(res_night, title = "jmj13_vs_Col0_night")
# jmj13 night vs jmj13 day
kegg_enrich(res_night_day, title = "jmj13_night_vs_jmj13_day")
# col0 night vs col0 day
kegg_enrich(res_col0_night_day, title = "col0_night_vs_col0_day")

# Print message to confirm completion
cat("\nAll plots and enrichment analyses have been saved!\n")


#fix volcano plot for jmj13 night vs day 
str(res_night_day)

#generate res_night_day properly 
res_night_day <- results(dds, contrast = c("condition", "jmj13_night", "jmj13_day"))
resultsNames(dds)
res_night_day <- results(dds, list(
  c("condition_night_vs_day", "genotypejmj13.conditionnight")
))

#generate res_col0_night_day properly & run 
res_col0_night_day <- results(dds, name = "condition_night_vs_day")
volcano_plot(res_col0_night_day, title = "Col0_night_vs_day")
heatmap_top(res_col0_night_day, rld, title = "Col0_night_vs_day")
go_enrich(res_col0_night_day, title = "Col0_night_vs_day")
kegg_enrich(res_col0_night_day, title = "Col0_night_vs_day")
gsea_run(res_col0_night_day, title = "Col0_night_vs_day")



