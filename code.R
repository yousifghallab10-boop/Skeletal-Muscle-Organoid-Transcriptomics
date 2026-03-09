# Read the file without row.names first
data_raw <- read.table(path, 
                       header = TRUE, 
                       sep = "\t", 
                       check.names = FALSE)

# Check the first few rows and columns
head(data_raw[, 1:5])

# 1. Aggregate duplicates by taking the mean of FPKM values
data_clean <- aggregate(. ~ gene_name, data = data_raw, FUN = mean)

# 2. Set gene_name as row names and then remove the column
rownames(data_clean) <- data_clean$gene_name
data_clean$gene_name <- NULL

# 3. Check the final clean data
head(data_clean[, 1:5])
dim(data_clean)

# Transform the data
data_log <- log2(data_clean + 1)

# Check the distribution with a Boxplot
boxplot(data_log, 
        las = 2, 
        main = "Log2 FPKM Distribution", 
        col = "steelblue",
        ylab = "Log2(FPKM + 1)")

# 1. Define groups based on sample names
group <- factor(ifelse(grepl("4wk", colnames(data_log)), "Week4",
                       ifelse(grepl("16wk", colnames(data_log)), "Week16", "Other")))

# 2. Filter data to keep only Week 4 and Week 16 for comparison
keep <- group %in% c("Week4", "Week16")
data_sub <- data_log[, keep]
group_sub <- factor(group[keep])

# 3. Create Design Matrix
design <- model.matrix(~0 + group_sub)
colnames(design) <- levels(group_sub)

# Install BiocManager if you don't have it
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install limma
BiocManager::install("limma")

# Now try to load it again
library(limma)

# 1. Create groups (Week 4 vs Week 16)
samples <- colnames(data_log)
group <- factor(ifelse(grepl("4wk", samples), "Week4",
                       ifelse(grepl("16wk", samples), "Week16", "Other")))

# 2. Keep only the samples we want to compare
keep <- group %in% c("Week4", "Week16")
data_sub <- data_log[, keep]
group_sub <- factor(group[keep])

# 3. Define the experimental design
design <- model.matrix(~0 + group_sub)
colnames(design) <- levels(group_sub)

# 4. Run the linear model
fit <- lmFit(data_sub, design)
cont_matrix <- makeContrasts(Week16vsWeek4 = Week16 - Week4, levels = design)
fit2 <- contrasts.fit(fit, cont_matrix)
fit2 <- eBayes(fit2)

# 5. Extract the results table
results <- topTable(fit2, coef = 1, number = Inf)

# Check how many genes are up or down
# Threshold: Log Fold Change > 1 and Adjusted P-value < 0.05
table(results$adj.P.Val < 0.05 & abs(results$logFC) > 1)

# Define significance thresholds
logFC_cutoff <- 1
p_cutoff <- 0.05

# Create classification column
results$diffexpressed <- "Not Significant"
results$diffexpressed[results$logFC > logFC_cutoff & results$adj.P.Val < p_cutoff] <- "Upregulated"
results$diffexpressed[results$logFC < -logFC_cutoff & results$adj.P.Val < p_cutoff] <- "Downregulated"

# Check counts for each category
table(results$diffexpressed)

library(ggplot2)

# Create the plot
volcano_p <- ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), color = diffexpressed)) +
  geom_point(alpha = 0.4, size = 1.5) + # alpha to handle overlapping points
  theme_minimal() +
  scale_color_manual(values = c("Downregulated" = "blue", 
                                "Not Significant" = "grey", 
                                "Upregulated" = "red")) +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), col = "black", linetype = "dashed") +
  geom_hline(yintercept = -log10(p_cutoff), col = "black", linetype = "dashed") +
  labs(title = "Volcano Plot: Differentiation Stages (16wk vs 4wk)",
       subtitle = paste0("Total DEGs: ", sum(results$diffexpressed != "Not Significant")),
       x = "log2(Fold Change)",
       y = "-log10(Adjusted P-value)",
       color = "Status") +
  theme(legend.position = "bottom")

# Show the plot
print(volcano_p)

# Save the plot in high quality for LinkedIn
# ggsave("Volcano_Plot.png", plot = volcano_p, width = 8, height = 6, dpi = 300)

# Filter significant genes only
sig_genes <- results[results$diffexpressed != "Not Significant", ]

# Save to CSV to keep it for the next steps
write.csv(sig_genes, "Significant_DEGs_Week16vs4.csv")

# Show the top 10 most significant genes
head(sig_genes[order(sig_genes$adj.P.Val), ])


if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"))

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# 1. Get the list of Upregulated genes (you can do the same for Downregulated later)
up_genes <- rownames(sig_genes[sig_genes$diffexpressed == "Upregulated", ])

# 2. Run GO Enrichment for Biological Process (BP)
go_results <- enrichGO(gene          = up_genes,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'SYMBOL',
                       ont           = "BP", # Biological Process
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05)

# 3. Show top results
head(go_results)

# A. Barplot of top enriched pathways
barplot(go_results, showCategory = 15) + 
  ggtitle("Top Biological Processes - Upregulated Genes")

# B. Dotplot (very popular in papers)
dotplot(go_results, showCategory = 15) + 
  ggtitle("Enriched Pathways: Week 16 vs Week 4")

# Select top 200 upregulated genes
top_up_genes <- sig_genes[sig_genes$diffexpressed == "Upregulated", ]
top_up_genes <- top_up_genes[order(top_up_genes$adj.P.Val), ][1:200, ]

# Write to a simple text file
write.table(rownames(top_up_genes), "Top_200_Genes_for_STRING.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
