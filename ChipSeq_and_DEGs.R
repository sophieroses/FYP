#to intersect methylated genes (from ChIP-seq) & DEGs (from RNA-seq)
#DEGs that are only methylated in jmj13 (exclusively), suggests jmj13 loss = gain of h3k27me3 --> repression
getwd()
#To get: table summarising total DEGs, total JMJ13-exclusive methylated genes, how many overlap 
# Load libraries
library(dplyr)
library(readr)
library(VennDiagram)

# Load DEGs: Col0 vs JMJ13 during the day
degs_col0_vs_jmj13_day <- read_csv("DEGs_GO_col0_vs_jmj13_day.csv")  # Must contain 'gene_id'

# Load methylation data
jmj13_exclusive_day <- read_csv("jmj13exclusive_methylated_genes_day.csv")    # Must contain 'gene_id'
col0_exclusive_day  <- read_csv("Col0exclusive_methylated_genes_day.csv")     # Must contain 'gene_id'

# Intersections
gained_methylation_overlap <- inner_join(
  degs_col0_vs_jmj13_day,
  jmj13_exclusive_day,
  by = c("gene_id" = "ensembl_gene_id")
)

lost_methylation_overlap <- inner_join(
  degs_col0_vs_jmj13_day,
  col0_exclusive_day,
  by = c("gene_id" = "ensembl_gene_id")
)


# Summary table
summary_overlap <- tibble(
  Comparison = "Col0 vs JMJ13 (Day)",
  Total_DEGs = nrow(degs_col0_vs_jmj13_day),
  JMJ13_Gain_Methylated = nrow(jmj13_exclusive_day),
  Col0_Gain_Methylated = nrow(col0_exclusive_day),
  DEGs_with_Gained_Methylation = nrow(gained_methylation_overlap),
  DEGs_with_Lost_Methylation = nrow(lost_methylation_overlap)
)

print(summary_overlap)

# Optional: Save overlaps
write_csv(gained_methylation_overlap, "DEGs_with_JMJ13exclusive_methylation_day.csv")
write_csv(lost_methylation_overlap, "DEGs_with_Col0exclusive_methylation_day.csv")

# Venn diagram (example for gain of methylation in JMJ13)
venn.plot <- draw.pairwise.venn(
  area1 = nrow(degs_col0_vs_jmj13_day),
  area2 = nrow(jmj13_exclusive_day),
  cross.area = nrow(gained_methylation_overlap),
  category = c("DEGs: Col0 vs JMJ13 (Day)", "JMJ13-Exclusive H3K27me3"),
  fill = c("skyblue", "salmon"),
  alpha = 0.5,
  cex = 1,
  cat.cex = 0.7,
  main = "Overlap: DEGs & JMJ13-Only Methylated Genes (Day)"
)

#filter & view the 5 genes to those downregulated in mutant 
# Genes with gained H3K27me3 AND downregulated in jmj13 mutant
# (log2FC < 0 means higher in Col-0 → down in jmj13)
direct_targets <- gained_methylation_overlap %>%
  filter(log2FoldChange < 0) %>%
  arrange(log2FoldChange)

# View the result
print(direct_targets)

# Optional: save to file
write_csv(direct_targets, "direct_targets_jmj13_H3K27me3_downregulated.csv")

#simple table of 5 overlapping genes showing their regulation in jmj13 mutant 
# Create a column for direction of regulation based on log2FoldChange
gained_methylation_overlap <- gained_methylation_overlap %>%
  mutate(regulation = ifelse(log2FoldChange < 0, "Down in jmj13", "Up in jmj13"))

# Select and display relevant columns only
gene_overlap_table <- gained_methylation_overlap %>%
  select(gene_id, log2FoldChange, padj, regulation)

# View in R
print(gene_overlap_table)

# Optional: save to CSV
write_csv(gene_overlap_table, "overlapping_genes_jmj13_H3K27me3_and_DEGs.csv")




###same for night --> load & rename
jmj13_exclusive_night <- read_csv("jmj13exclusive_methylated_genes_night.csv") %>%
  rename(gene_id = ensembl_gene_id)

degs_col0_vs_jmj13_night <- read_csv("DEGs_GO_col0_vs_jmj13_night.csv")

# Overlap
gained_methylation_night_overlap <- inner_join(degs_col0_vs_jmj13_night, jmj13_exclusive_night, by = "gene_id")
# Add regulation direction
gained_methylation_night_overlap <- gained_methylation_night_overlap %>%
  mutate(regulation = ifelse(log2FoldChange < 0, "Down in jmj13", "Up in jmj13"))

#venn diagram 
venn.plot <- draw.pairwise.venn(
  area1 = nrow(degs_col0_vs_jmj13_night),
  area2 = nrow(jmj13_exclusive_night),
  cross.area = nrow(gained_methylation_overlap),
  category = c("DEGs: Col0 vs JMJ13 (Night)", "JMJ13-Exclusive H3K27me3"),
  fill = c("skyblue", "salmon"),
  alpha = 0.5,
  cex = 1,
  cat.cex = 0.7,
  main = "Overlap: DEGs & JMJ13-Only Methylated Genes (Night)"
)






#create table
gene_overlap_night_table <- gained_methylation_night_overlap %>%
  select(gene_id, log2FoldChange, padj, regulation)

print(gene_overlap_night_table)

write_csv(gene_overlap_night_table, "gene_overlap_night.csv")

#compare day and night overlaps
intersect(gained_methylation_overlap$gene_id, gained_methylation_night_overlap$gene_id)
