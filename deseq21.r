install.packages("BiocManager")
BiocManager::install("DESeq2", force = TRUE)

# Load required libraries
library(DESeq2)
library(ggplot2)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)


# count matrix
countdata <- read.csv("/media/shubham/NV21/PRJ2025JOB/counts/output.csv", row.names = 1)

# col data
sampledata <- read.csv("/media/shubham/NV21/PRJ2025JOB/meta.csv", row.names = 1)

sampledata$condition <- factor(sampledata$condition, levels = c("tmz", "combo"))

# match sampledata and countdata
stopifnot(all(rownames(sampledata) == colnames(countdata)))

# DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = sampledata,
                              design = ~ condition) 

# Run DESeq2
dds <- DESeq(dds)
res <- results(dds)

# Remove NA values
res <- na.omit(res)

# Order by adjusted p-value
res_ordered <- res[order(res$padj), ]

# Filter significant DEGs (adjusted p-value < 0.05)
res_sig <- subset(res_ordered, padj < 0.05)

# Save significant DEGs to CSV
write.csv(as.data.frame(res_sig), file = "DEG_results_annotated.csv")

# Prepare data for volcano plot
res_df <- as.data.frame(res)
res_df <- res_df[order(res_df$padj), ]

res_df$log2FoldChange[is.na(res_df$log2FoldChange)] <- 0
res_df$padj[is.na(res_df$padj)] <- 1

# Define threshold groups
res_df$threshold <- ifelse(res_df$padj < 0.05 & res_df$log2FoldChange > 1, "Up",
                    ifelse(res_df$padj < 0.05 & res_df$log2FoldChange < -1, "Down", "NS"))
res_df$threshold <- as.factor(res_df$threshold)


# Volcano plot without labels
volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = threshold)) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value",
       color = "Expression") +
  theme(legend.position = "top")

# Create volcano plot with top 10 gene labels
top_genes <- res_df[order(res_df$padj), ][1:20, ]
top_genes$gene <- rownames(top_genes)

#all genes 
# All genes ordered by adjusted p-value
all_genes <- res_df[order(res_df$padj), ]
all_genes$gene <- rownames(all_genes)
allg <- write.csv(all_genes, "allgenes.csv")

#get up/down genes only exclude ns
res_df_up_down <- res_df %>%
  filter(threshold %in% c("Up", "Down"))



library(org.Hs.eg.db)
library(AnnotationDbi)
# Add gene symbols column if not present
res_df_up_down$symbol <- rownames(res_df_up_down)

# Map gene symbols to Entrez IDs
res_df_up_down$entrez <- mapIds(org.Hs.eg.db,
                           keys = res_df_up_down$symbol,
                           column = "ENTREZID",
                           keytype = "SYMBOL",
                           multiVals = "first")

forfunctional <- write.csv(res_df_up_down, "genesupdownforfunctionalen.csv")

#save top 20 gene
top20 <- write.csv(top_genes, "top20.csv")

volcano_labeled <- volcano_plot +
  geom_text(data = top_genes, aes(label = gene), vjust = 1, hjust = 0.5, size = 3, check_overlap = TRUE)


# Display labeled plot
print(volcano_labeled)

# Save labeled volcano plot
ggsave("volcano_plot_labeled.jpg", plot = volcano_labeled, width = 8, height = 6)

