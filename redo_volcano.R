getwd()
###try again 
library(EnhancedVolcano)
library(readr)
library(dplyr)
library(ggplot2)

# Optional: setwd("/Users/sophiewilliams/fyp2")

files <- list(
  "res_night.csv" = "jmj13_vs_Col0_night",
  "res_day.csv" = "jmj13_vs_Col0_day",
  "res_jmj13_night_vs_day.csv" = "jmj13_night_vs_day",
  "res_col0_night_vs_day.csv" = "Col0_night_vs_day"
)

generate_volcano_plot <- function(filename, title) {
  res <- read_csv(filename)
  
  sig_res <- res %>% filter(padj < 0.05)
  top_up <- sig_res %>% filter(log2FoldChange > 0) %>% slice_max(log2FoldChange, n = 10)
  top_down <- sig_res %>% filter(log2FoldChange < 0) %>% slice_min(log2FoldChange, n = 10)
  top_genes <- unique(c(top_up$geneName, top_down$geneName))
  
  volcano <- EnhancedVolcano(res,
                             lab = ifelse(res$geneName %in% top_genes, res$geneName, ""),
                             x = 'log2FoldChange',
                             y = 'padj',
                             title = title,
                             subtitle = 'EnhancedVolcano',
                             xlab = bquote(~Log[2]~ 'fold change'),
                             ylab = bquote(~-Log[10]~ 'adjusted p-value'),
                             pCutoff = 0.05,
                             FCcutoff = 1,
                             pointSize = 2.0,
                             labSize = 4.0,
                             col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                             colAlpha = 0.9,
                             legendLabels = c('NS', 'Log2 FC', 'p-value', 'p-value & Log2 FC'),
                             legendPosition = 'top',
                             legendLabSize = 12,
                             legendIconSize = 4.0,
                             drawConnectors = TRUE,
                             widthConnectors = 0.8,
                             colConnectors = 'grey30',
                             max.overlaps = Inf
  )
  
  # Save using ggsave
  outfile <- paste0("volcano_", title, ".png")
  ggsave(outfile, plot = volcano, width = 8, height = 6, dpi = 300)
  message("✅ Saved volcano plot to: ", normalizePath(outfile))
}

for (file in names(files)) {
  generate_volcano_plot(file, files[[file]])
}

##redo GO dotplots ---> DIDNT WORK 


# Named list of files and titles


#try again 
library(ggplot2)
library(readr)
library(dplyr)
library(forcats)

# Named list of files and titles
go_files <- list(
  "DEGs_GO_col0_vs_jmj13_day.csv" = "jmj13_vs_Col0_day",
  "DEGs_GO_col0_vs_jmj13_night.csv" = "jmj13_vs_Col0_night",
  "DEGs_GO_jmj13_day_vs_night.csv" = "jmj13_night_vs_day",
  "DEGs_GO_col0_day_vs_night.csv" = "Col0_night_vs_day"
)

# Plotting function
plot_go_dotplot <- function(filename, title) {
  df <- read_csv(filename)
  
  # Ensure padj is numeric for arranging
  df$padj <- as.numeric(df$padj)
  
  # Calculate GeneRatio if not available
  # Assuming the total number of genes is the total rows in df or a fixed number
  total_genes <- nrow(df)  # Update if you have a specific total number
  
  # Calculate GeneRatio: (Count of genes in the GO term) / (Total genes)
  # Assuming each row corresponds to a gene associated with a GO term
  df$GeneRatio <- 1 / total_genes  # If each gene counts as 1 for its associated GO term
  
  # Simplify term names for plotting
  df$GO_Names <- as.factor(df$GO_Names)
  df <- df %>%
    arrange(padj) %>%
    mutate(GO_Names = fct_reorder(GO_Names, GeneRatio))
  
  p <- ggplot(df, aes(x = GeneRatio, y = GO_Names)) +
    geom_point(aes(size = GeneRatio * total_genes, color = padj)) +  # Count can be GeneRatio * total_genes
    scale_color_gradient(low = "red", high = "blue") +
    theme_bw() +
    labs(
      title = paste("GO Enrichment -", title),
      x = "GeneRatio",
      y = NULL,
      color = "padj",
      size = "Count"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.text = element_text(size = 12)
    )
  
  outfile <- paste0("GO_", title, ".png")
  ggsave(outfile, plot = p, width = 10, height = 8, dpi = 300)
  message("✅ Saved GO plot to: ", normalizePath(outfile))
}

# Run for all GO files
for (file in names(go_files)) {
  plot_go_dotplot(file, go_files[[file]])
}



