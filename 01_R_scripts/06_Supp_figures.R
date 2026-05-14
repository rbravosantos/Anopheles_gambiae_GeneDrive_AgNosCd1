
library(tximport)
library(readr)
library(rhdf5)
library(DESeq2)
library(readxl)
library(patchwork)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

rm(list = ls())

# Define the main directory
main_dir <- "02_kallisto_output"

# Get a list of all subdirectories
subdirs <- list.dirs(main_dir, recursive = FALSE)

# Create a named vector of file paths to the abundance files
files <- file.path(subdirs, "abundance.h5")
names(files) <- basename(subdirs)

# Import the abundance files using tximport
txi <- tximport(files, type = "kallisto", txOut = TRUE)

# Print the imported data
#print(txi)

metadata <- read_excel("metadata.xlsx")

metadata$sex <- as.factor(metadata$sex)
metadata$diet <- as.factor((metadata$diet))

# Label genotypes
metadata$genotype <- factor(metadata$genotype, 
                            levels = c("WT", "GD", "KO"),
                            labels = c("AngG3", "AgNosCd-1", "Agcd-deletion11"))


# Exclude sample D4 from metadata - Did not pass QC before sequencing
metadata <- metadata[metadata$sample != "D4", ]

# Update rownames of metadata to match colnames of txi$counts
rownames(metadata) <- metadata$sample

# Filter txi to match the samples in metadata (excluding D4)
txi$counts <- txi$counts[, colnames(txi$counts) %in% rownames(metadata)]
txi$abundance <- txi$abundance[, colnames(txi$abundance) %in% rownames(metadata)]
txi$length <- txi$length[, colnames(txi$length) %in% rownames(metadata)]



# Create DESeq2 object for normalization
dds <- DESeqDataSetFromTximport(txi, 
                                colData = metadata, 
                                design = ~ genotype + sex + diet)

# Normalize the data
dds <- DESeq(dds)


# Get normalized counts with variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)

# Calculate sample-to-sample correlation matrix
cor_matrix <- cor(assay(vsd), method = "pearson")

# Create annotation dataframe for heatmap
annotation_df <- data.frame(
  Genotype = metadata$genotype,
  Sex = metadata$sex,
  Diet = metadata$diet,
  row.names = rownames(metadata)
)


# Define colors for annotations based on actual factor levels
genotype_colors <- c("AngG3" = "#619CFF", 
                     "AgNosCd-1" = "#00BA38", 
                     "Agcd-deletion11" = "#F8766D")
sex_colors <- c("Male" = "dodgerblue2", "Female" = "hotpink1")
diet_colors <- c("Sugarfed" = "tan2", "Bloodfed" = "red4")

# Create annotation object
ha_top = HeatmapAnnotation(
  Genotype = annotation_df$Genotype,
  Sex = annotation_df$Sex,
  Diet = annotation_df$Diet,
  col = list(
    Genotype = genotype_colors,
    Sex = sex_colors,
    Diet = diet_colors
  ),
  annotation_name_side = "left"
)

# Create the sample correlation heatmap (without coefficient labels)
sample_correlation_heatmap <- Heatmap(
  cor_matrix,
  name = "Correlation",
  top_annotation = ha_top,
  
  # Color scheme for correlation values
  col = colorRamp2(c(0.7, 0.85, 1), c("blue", "white", "red")),
  
  # Clustering
  clustering_distance_rows = "euclidean",
  clustering_distance_columns = "euclidean",
  clustering_method_rows = "complete",
  clustering_method_columns = "complete",
  
  # Heatmap title and labels
  column_title = "Sample Correlation Heatmap",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  
  # Sample labels
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  
  # Heatmap dimensions
  width = unit(12, "cm"),
  height = unit(12, "cm")
)

# Create gene expression heatmap with top variable genes
# Select top variable genes
top_genes <- 1000                  # Number of top variable genes to display
rv <- rowVars(assay(vsd))
select_genes <- order(rv, decreasing = TRUE)[seq_len(min(top_genes, length(rv)))]

# Get the expression matrix for selected genes
gene_matrix <- assay(vsd)[select_genes, ]

# Scale the gene expression data (z-score normalization)
gene_matrix_scaled <- t(scale(t(gene_matrix)))

# Create gene expression heatmap
gene_expression_heatmap <- Heatmap(
  gene_matrix_scaled,
  name = "Z-score",
  top_annotation = ha_top,
  
  # Color scheme for z-scores
  col = colorRamp2(c(-3, 0, 3), c("blue", "white", "red")),
  
  # Clustering
  clustering_distance_rows = "euclidean",
  clustering_distance_columns = "euclidean",
  clustering_method_rows = "complete",
  clustering_method_columns = "complete",
  
  # Hide gene names (too many to display)
  show_row_names = FALSE,
  
  # Heatmap title and labels
  column_title = paste("Top", top_genes, "Variable Genes Expression Heatmap"),
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  
  # Sample labels
  column_names_gp = gpar(fontsize = 10),
  
  # Heatmap dimensions
  width = unit(14, "cm"),
  height = unit(16, "cm")
)

print(sample_correlation_heatmap)
print(gene_expression_heatmap)




