# Load libraries
library(dplyr)
library(DESeq2)
library(tidyverse)
library(GEOquery)

#Read count RNA - Seq data
count_data <- read_tsv('F:/Bversity/Mini Project/GSE244282_raw_counts_GRCh38.p13_NCBI.tsv/GSE244282_raw_counts_GRCh38.p13_NCBI.tsv')
counts_data <- as.data.frame(count_data)
rownames(counts_data) <- counts_data$GeneID  # Set 'gene_id' as row names
counts_data <- counts_data[ , -1] 

#metadata creation
SampleID <- colnames(counts_data)  # Get column names as sample IDs
Type <- c("treatment_1", "treatment_1", "treatment_1", "treatment_2", "treatment_2", "treatment_2")  # Define conditions

# Combine into a metadata dataframe
meta_data <- data.frame(SampleID = SampleID, Type = factor(Type), row.names = 1)

# View metadata
print(metadata)


colnames(counts_data)
rownames(meta_data)
#making sure the column names in counts_data is same as sample _data
all(colnames(counts_data) %in% rownames(meta_data))

# making sure it is present in same order
all(colnames(counts_data) == rownames(metadata))
counts_data <- as.matrix(counts_data)
str(counts_data)
row_sums <- rowSums(counts_data)
# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = counts_data, colData = meta_data, design = ~ Type)

#prefiltering for DESeq dataset [removing rows less than 10 reads]
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds
# Run DESeq normalization and differential expression analysis
dds <- DESeq(dds)

# View results
res <- results(dds)
res
res_df <- as.data.frame(res)
head(rownames(res_df))

write.csv(res_df, "F:/Bversity/Mini Project/GSE244282_raw_counts_GRCh38.p13_NCBI.tsv/DESeq2_results.csv", row.names = TRUE)
cat("DESeq2 results saved to 'DESeq2_results.csv'\n")

# Exploring the results of DESeq analysis
summary(res)

#MA plot
plotMA(res)

library(tidyverse)
library(ggplot2)
library(EnhancedVolcano)

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = 'Volcano Plot of DEGs')

# Order results by adjusted p-value (padj)
res_sorted <- res[order(res$padj), ]

# View top 10 most significant genes
head(res_sorted, 10)
sig_genes <- subset(res_sorted, padj < 0.05)

# Filter upregulated genes (logFC > 0, padj < 0.05)
upregulated_genes <- res[which(res$log2FoldChange > 0 & res$padj < 0.05), ]

# Filter downregulated genes (logFC < 0, padj < 0.05)
downregulated_genes <- res[which(res$log2FoldChange < 0 & res$padj < 0.05), ]
# View upregulated genes
upregulated_genes

# View downregulated genes
downregulated_genes



# Extract normalized counts
normalized_counts <- assay(rlog(dds))

# Select significant genes for heatmap
sig_genes_data <- normalized_counts[rownames(sig_genes), ]

# Plot heatmap
library(pheatmap)
library(RColorBrewer)

# Transpose the data to place significant genes in columns
sig_genes_data_transposed <- t(sig_genes_data)

# Define sample annotations (samples in rows now)
sample_annotations <- data.frame(
  Treatment = rep(c("Treatment A", "Treatment B"), each = 3)
)
rownames(sample_annotations) <- rownames(sig_genes_data_transposed)  # Match rows with sample names

# Generate the heatmap
pheatmap(sig_genes_data_transposed,
         cluster_rows = TRUE,          # Clustering for samples (rows)
         cluster_cols = TRUE,          # Clustering for genes (columns)
         annotation_row = sample_annotations, # Add treatment annotations
         show_colnames = TRUE,         # Show gene names as column names
         color = custom_colors,        # Custom color palette
         main = "Heatmap of Significant Genes Across Treatments") # Add a title

sample_annotations <- data.frame(
  Treatment = rep(c("Treatment A", "Treatment B"), each = 6)
)
rownames(sample_annotations) <- colnames(sig_genes_data)

pheatmap(sig_genes_data,
         cluster_rows = TRUE,          # Clustering for genes (rows)
         cluster_cols = TRUE,          # Clustering for samples (columns)
         annotation_col = sample_annotations, # Add treatment annotations
         show_rownames = TRUE,         # Show gene names
         color = custom_colors,        # Custom color palette
         main = "Heatmap of Significant Genes Across Treatments") # Add a title
#Perform variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)

# Perform PCA using prcomp
pca <- prcomp(t(assay(vsd)))

# Get percentage of variance explained by each principal component
percent_variance <- (pca$sdev^2 / sum(pca$sdev^2)) * 100

# Prepare PCA data for plotting
pca_data <- as.data.frame(pca$x)
pca_data$condition <- colData(vsd)$condition  # Add condition information

# Plot PCA with % variance in labels
ggplot(pca_data, aes(x = PC1, y = PC2, color = Type)) +
  geom_point(size = 3) +
  labs(
    title = "PCA of Samples",
    x = paste0("PC1 (", round(percent_variance[1], 1), "%)"),
    y = paste0("PC2 (", round(percent_variance[2], 1), "%)")
  ) +
  theme_minimal()

# Create a volcano plot
library(ggplot2)
ggplot(res, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05 & abs(log2FoldChange) > 1)) +
  scale_color_manual(values = c("black", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10(padj)") +
  theme(legend.position = "none")

# Apply regularized log transformation
rld <- rlog(dds, blind = TRUE)  # blind = TRUE ignores experimental design
rlog_data <- assay(rld)  # Extract transformed data as a matrix
head(rlog_data)

# generation of distance matrix
sampleDists <- dist(t(assay(rld)))
library(RColorBrewer)
library(pheatmap)
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues"))) (255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col=colors)

#Dispersion plot
plotDispEsts(dds, main = 'Dispersion Plot')

library(AnnotationDbi)
library(org.Hs.eg.db)
res_df$ensembl_id <- mapIds(org.Hs.eg.db,
                            keys = rownames(res_df),
                            column = "ENSEMBL",      # Map to Ensembl IDs
                            keytype = "ENTREZID",      # Replace 'SYMBOL' with your gene ID format
                            multiVals = "first")
res_df$gene_name <- mapIds(org.Hs.eg.db,
                           keys = rownames(res_df),
                           column = "GENENAME",     # Map to Gene Names (descriptions)
                           keytype = "ENTREZID",      # Replace 'SYMBOL' with your input ID type
                           multiVals = "first")
res_df$gene_symbol <- mapIds(org.Hs.eg.db,
                             keys = rownames(res_df),
                             column = "SYMBOL",      # Map to Gene Symbols
                             keytype = "ENTREZID",     # Replace 'SYMBOL' with the current gene ID type
                             multiVals = "first")

head(res_df)
write.csv(res_df, "F:/Bversity/Mini Project/GSE244282_raw_counts_GRCh38.p13_NCBI.tsv/DESeq2_results_with_annotations.csv", row.names = TRUE)
