# Load necessary libraries
library(tximport)
library(readr)
library(rhdf5)
library(DESeq2)
library(readxl)
library(patchwork)
library(ggplot2)
library(ggpubr)
library(ggrepel)

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

# Exclude females
metadata <- metadata[metadata$sex != "Female", ]

# Update rownames of metadata to match colnames of txi$counts
rownames(metadata) <- metadata$sample

# Filter txi to exclude the counts for the sample "D4"
txi$counts <- txi$counts[, colnames(txi$counts) %in% rownames(metadata)]
txi$abundance <- txi$abundance[, colnames(txi$abundance) %in% rownames(metadata)]
txi$length <- txi$length[, colnames(txi$length) %in% rownames(metadata)]

rownames(metadata) <- colnames(txi$counts) # to set up every row for every sample

dds <- DESeqDataSetFromTximport(txi, 
                                colData = metadata, 
                                design = ~ condition)


# normalize the read counts,estimate dispersions, and fit the linear model
dds <- DESeq(dds,betaPrior = T)

# Normalization factors
#sizeFactors(dds)
normalizationFactors(dds) # for tximport

# Plot the dispersions
#plotDispEsts(dds)

# Plot PCA
vsd <- vst(dds, blind = FALSE)

pca_condition <- plotPCA(vsd, intgroup = "condition") + 
  scale_color_manual("Genotype", labels = c("AngG3", "AgNosCd-1", "Agcd-deletion11"),
                     values = c("AngG3" = "#619CFF", 
                                "AgNosCd-1"  = "#00BA38",
                                "Agcd-deletion11" = "#F8766D")) +
  coord_fixed() +
  xlim(-20, 20) +
  ylim(-20, 20)

#plot(pca_condition)

##########################################
###   PCA males: AngG3 vs. AgNosCd-1   ###
##########################################

# Define data for second analysis

txi2 <- txi

metadata2 <- metadata[metadata$genotype != "Agcd-deletion11",]

# Update rownames of metadata to match colnames of txi$counts
rownames(metadata2) <- metadata2$sample

# Filter txi to exclude the counts for the sample "D4"
txi2$counts <- txi2$counts[, colnames(txi2$counts) %in% rownames(metadata2)]
txi2$abundance <- txi2$abundance[, colnames(txi2$abundance) %in% rownames(metadata2)]
txi2$length <- txi2$length[, colnames(txi2$length) %in% rownames(metadata2)]

rownames(metadata2) <- colnames(txi2$counts) # to set up every row for every sample


dds2 <- DESeqDataSetFromTximport(txi2, 
                                colData = metadata2, 
                                design = ~ condition)


# normalize the read counts,estimate dispersions, and fit the linear model
dds2 <- DESeq(dds2,betaPrior = T)

# Normalization factors
#sizeFactors(dds)
normalizationFactors(dds2) # for tximport

# Plot the dispersions
#plotDispEsts(dds2)




################################################
###   PCA males: AngG3 vs. Agcd-deletion11   ###
################################################

metadata3 <- metadata

metadata3$sex <- as.factor(metadata3$sex)
metadata3$diet <- as.factor((metadata3$diet))

# Exclude females
metadata3 <- metadata3[metadata3$sex != "Female", ]
metadata3 <- metadata3[metadata3$genotype != "AgNosCd-1",]

# Update rownames of metadata to match colnames of txi$counts
rownames(metadata3) <- metadata3$sample

# Import the abundance files using tximport
txi3 <- txi

# Filter txi to exclude the counts for the sample "D4"
txi3$counts <- txi3$counts[, colnames(txi3$counts) %in% rownames(metadata3)]
txi3$abundance <- txi3$abundance[, colnames(txi3$abundance) %in% rownames(metadata3)]
txi3$length <- txi3$length[, colnames(txi3$length) %in% rownames(metadata3)]

# Create a new group column in metadata

rownames(metadata3) <- colnames(txi3$counts) # to set up every row for every sample

dds3 <- DESeqDataSetFromTximport(txi3, 
                                 colData = metadata3, 
                                 design = ~ condition)


# normalize the read counts,estimate dispersions, and fit the linear model
dds3 <- DESeq(dds3,betaPrior = T)

# Normalization factors
#sizeFactors(dds)
normalizationFactors(dds3) # for tximport

# Plot the dispersions
#plotDispEsts(dds3)

####################################################
###   PCA males: AgNosCd-1 vs. Agcd-deletion11   ###
####################################################


metadata4 <- metadata

metadata4$sex <- as.factor(metadata4$sex)
metadata4$diet <- as.factor((metadata4$diet))

# Exclude females
metadata4 <- metadata4[metadata4$sex != "Female", ]
metadata4 <- metadata4[metadata4$genotype != "AngG3",]

# Update rownames of metadata to match colnames of txi$counts
rownames(metadata4) <- metadata4$sample

# Import the abundance files using tximport
txi4 <- txi

# Filter txi to exclude the counts for the sample "D4"
txi4$counts <- txi4$counts[, colnames(txi4$counts) %in% rownames(metadata4)]
txi4$abundance <- txi4$abundance[, colnames(txi4$abundance) %in% rownames(metadata4)]
txi4$length <- txi4$length[, colnames(txi4$length) %in% rownames(metadata4)]

# Create a new group column in metadata
#metadata$sex_diet <- factor(paste0(metadata$sex, metadata$diet))
#metadata$group <- factor(paste0(metadata$sex, metadata$diet, metadata$condition))

rownames(metadata4) <- colnames(txi4$counts) # to set up every row for every sample

#metadata$group <- factor(paste0(metadata$sex, metadata$diet, metadata$condition))

dds4 <- DESeqDataSetFromTximport(txi4, 
                                 colData = metadata4, 
                                 design = ~ condition)


# normalize the read counts,estimate dispersions, and fit the linear model
# betaPrior = T opton tells DESeq2 to squeeze the log Fold Changes of lowly expressed genes toward zero.
dds4 <- DESeq(dds4,betaPrior = T)

# Normalization factors
#sizeFactors(dds)
normalizationFactors(dds4) # for tximport

# Plot the dispersions
#plotDispEsts(dds3)


#################################
### Volcano plots: WT vs GD   ###
#################################


# GD vs WT (instead of WT vs GD)
res1 <- results(dds,
                contrast = c("condition","AgNosCd-1","AngG3"),
                lfcThreshold = 0.58,
                alpha = 0.05)

# KO vs WT
res2 <- results(dds,
                contrast = c("condition","Agcd-deletion11","AngG3"),
                lfcThreshold = 0.58,
                alpha = 0.05)

# KO vs GD
res3 <- results(dds,
                contrast = c("condition","Agcd-deletion11","AgNosCd-1"),
                lfcThreshold = 0.58,
                alpha = 0.05)


res1 <- as.data.frame(res1)
res2 <- as.data.frame(res2)
res3 <- as.data.frame(res3)


# Remove the specified substring from row names

rownames(res1) <- sub(".*(?=XP)", "", rownames(res1), perl = TRUE)
rownames(res1) <- sub(".*(?=YP)", "", rownames(res1), perl = TRUE)
rownames(res1) <- sub("_[^_]*$", "", rownames(res1), perl = TRUE)

rownames(res2) <- sub(".*(?=XP)", "", rownames(res2), perl = TRUE)
rownames(res2) <- sub(".*(?=YP)", "", rownames(res2), perl = TRUE)
rownames(res2) <- sub("_[^_]*$", "", rownames(res2), perl = TRUE)

rownames(res3) <- sub(".*(?=XP)", "", rownames(res3), perl = TRUE)
rownames(res3) <- sub(".*(?=YP)", "", rownames(res3), perl = TRUE)
rownames(res3) <- sub("_[^_]*$", "", rownames(res3), perl = TRUE)


library(clusterProfiler)
library(AnnotationDbi)
library(org.Ag.eg.db)

res1$cleaned_keys <- as.character(sub("\\..*", "", row.names(res1)))
res2$cleaned_keys <- as.character(sub("\\..*", "", row.names(res2)))
res3$cleaned_keys <- as.character(sub("\\..*", "", row.names(res3)))


res1$symbol <- mapIds(org.Ag.eg.db,
                      keys=res1$cleaned_keys, 
                      column="SYMBOL",
                      keytype="REFSEQ",
                      multiVals="first")
res1$entrez = mapIds(org.Ag.eg.db,
                     keys=res1$cleaned_keys, 
                     column="ENTREZID",
                     keytype="REFSEQ",
                     multiVals="first")
res1$name =   mapIds(org.Ag.eg.db,
                     keys=res1$cleaned_keys, 
                     column="GENENAME",
                     keytype="REFSEQ",
                     multiVals="first")
res1$refseq = mapIds(org.Ag.eg.db,
                     keys=res1$cleaned_keys, 
                     column="REFSEQ",
                     keytype="REFSEQ",
                     multiVals="first")


res2$symbol <- mapIds(org.Ag.eg.db,
                      keys=res2$cleaned_keys, 
                      column="SYMBOL",
                      keytype="REFSEQ",
                      multiVals="first")
res2$entrez = mapIds(org.Ag.eg.db,
                     keys=res2$cleaned_keys, 
                     column="ENTREZID",
                     keytype="REFSEQ",
                     multiVals="first")
res2$name =   mapIds(org.Ag.eg.db,
                     keys=res2$cleaned_keys, 
                     column="GENENAME",
                     keytype="REFSEQ",
                     multiVals="first")
res2$refseq = mapIds(org.Ag.eg.db,
                     keys=res2$cleaned_keys, 
                     column="REFSEQ",
                     keytype="REFSEQ",
                     multiVals="first")


res3$symbol <- mapIds(org.Ag.eg.db,
                      keys=res3$cleaned_keys, 
                      column="SYMBOL",
                      keytype="REFSEQ",
                      multiVals="first")
res3$entrez = mapIds(org.Ag.eg.db,
                     keys=res3$cleaned_keys, 
                     column="ENTREZID",
                     keytype="REFSEQ",
                     multiVals="first")
res3$name =   mapIds(org.Ag.eg.db,
                     keys=res3$cleaned_keys, 
                     column="GENENAME",
                     keytype="REFSEQ",
                     multiVals="first")
res3$refseq = mapIds(org.Ag.eg.db,
                     keys=res3$cleaned_keys, 
                     column="REFSEQ",
                     keytype="REFSEQ",
                     multiVals="first")



# Define differentially expressed genes with broader thresholds for coloring
res1$diffexpressed <- "NO"
res1$diffexpressed[res1$log2FoldChange > 0.58 & res1$padj < 0.05] <- "UP"
res1$diffexpressed[res1$log2FoldChange < -0.58 & res1$padj < 0.05] <- "DOWN"

# Create a new column "delabel" to contain the name of differentially expressed genes with stricter thresholds for labeling
res1$delabel <- NA
res1$delabel <- ifelse(res1$log2FoldChange > 2.98 & res1$padj < 2.93E-06, res1$name, res1$delabel)
res1$delabel <- ifelse(res1$log2FoldChange < -2.40 & res1$padj < 1.12E-09, res1$name, res1$delabel)


# Create the volcano plot
p1a <- ggplot(data=res1, aes(x=log2FoldChange, 
                            y=-log10(padj), 
                            col=diffexpressed, 
                            label=delabel)) + 
  geom_point() +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())  +
  geom_vline(xintercept=c(-0.58, 0.58), col="red", linetype = "dashed") +
  geom_hline(yintercept=-log10(0.05), col="red", linetype = "dashed") + 
  geom_label_repel(fill = "white", 
                   max.overlaps = Inf, 
                   aes(hjust=0, segment.size = 0.75),
                   min.segment.length = unit(0,"lines"),
                   color = "black",
                   #nudge_x = -5,
                   nudge_y = 60
                   ) +
  scale_color_manual(values=c("#619CFF", "black", "#00BA38")) +
  labs(color = "Differentially expressed") +
  ylim(0,150) +
  xlim(-15,15) +
  ggtitle("AngG3 vs AgNosCd1") +
  xlab("Log2 Fold Change") + ylab("-Log10 p-value") +
  theme(legend.position = "bottom", legend.box = "vertical")

p1b <- ggplot(data=res1, aes(x=log2FoldChange, 
                            y=-log10(padj), 
                            col=diffexpressed, 
                            #label=delabel
                            )) + 
  geom_point() +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())  +
  geom_vline(xintercept=c(-0.58, 0.58), col="red", linetype = "dashed") +
  geom_hline(yintercept=-log10(0.05), col="red", linetype = "dashed") + 
#  geom_label_repel(fill = "white", 
#                   max.overlaps = Inf, 
#                   aes(hjust=0, segment.size = 0.75),
#                   min.segment.length = unit(0,"lines"),
#                   color = "black",
#                   #nudge_x = -5,
#                   nudge_y = 60) +
  scale_color_manual(values=c("#619CFF", "black", "#00BA38")) +
  labs(color = "Differentially expressed") +
  ylim(0,150) +
  xlim(-15,15) +
  ggtitle("AngG3 vs AgNosCd1") +
  xlab("Log2 Fold Change") + ylab("-Log10 p-value") +
  theme(legend.position = "bottom", legend.box = "vertical")


# List of genes to highlight
highlight_genes1 <- c("XP_311486", "XP_061500553", "XP_061502559", "XP_001231036", "XP_061513319",   
                      "XP_320163", "XP_309021", "XP_001238133", "XP_001237381", "XP_001238571")      

# Create a new column to indicate if the gene should be highlighted
res1$delabel2 <- ifelse(res1$refseq %in% highlight_genes1, res1$name, NA)

# Create a label column just for the highlighted genes
#res1$delabel_highlight <- ifelse(rownames(res1) %in% highlight_genes1, rownames(res1), NA)

# Create the volcano plot with highlights
p1c <- ggplot(data=res1, aes(x=log2FoldChange, 
                             y=-log10(padj), 
                             color=diffexpressed,
                             label=delabel2)) + 
  geom_point() +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())  +
  geom_vline(xintercept=c(-0.58, 0.58), col="red", linetype = "dashed") +
  geom_hline(yintercept=-log10(0.05), col="red", linetype = "dashed") + 
  geom_label_repel(fill = "white", 
                   force = 10,
                   max.overlaps = Inf, 
                   aes(hjust=0, segment.size = 0.75),
                   min.segment.length = unit(0,"lines"),
                   color = "black",
                   nudge_y = 25,
                   direction = "both") +
  scale_color_manual(values=c("#619CFF", "black", "#00BA38")) +
  labs(color = "Differentially expressed") +
  ylim(0,200) +
  xlim(-20,20) +
  ggtitle("AngG3 vs AgNosCd1") +
  xlab("Log2 Fold Change") + ylab("-Log10 p-value") +
  theme(legend.position = "bottom", legend.box = "vertical")

# Plot the volcano plot
#plot(p1a)
#plot(p1b)
#plot(p1c)



# Define differentially expressed genes with broader thresholds for coloring
res2$diffexpressed <- "NO"
res2$diffexpressed[res2$log2FoldChange > 0.58 & res2$padj < 0.05] <- "UP"
res2$diffexpressed[res2$log2FoldChange < -0.58 & res2$padj < 0.05] <- "DOWN"

# Create a new column "delabel" to contain the name of differentially expressed genes with stricter thresholds for labeling
res2$delabel <- NA
res2$delabel <- ifelse(res2$log2FoldChange > 2.98 & res2$padj < 2.93E-06, res2$name, res2$delabel)
res2$delabel <- ifelse(res2$log2FoldChange < -2.66 & res2$padj < 5.61E-08, res2$name, res2$delabel)


# Create the volcano plot
p2a <- ggplot(data=res2, aes(x=log2FoldChange, 
                            y=-log10(padj), 
                            col=diffexpressed, 
                            label=delabel)) + 
  geom_point() +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())  +
  geom_vline(xintercept=c(-0.58, 0.58), col="red", linetype = "dashed") +
  geom_hline(yintercept=-log10(0.05), col="red", linetype = "dashed") + 
  geom_label_repel(fill = "white", 
                   max.overlaps = Inf, 
                   aes(hjust=0, segment.size = 0.75),
                   min.segment.length = unit(0,"lines"),
                   color = "black",
                   #nudge_x = -5,
                   nudge_y = 10) +
  scale_color_manual(values=c("#619CFF", "black", "#00BA38")) +
  labs(color = "Differentially expressed") +
  ylim(0,150) +
  xlim(-15,15) +
  ggtitle("AngG3 vs Agcd-deletion11") +
  xlab("Log2 Fold Change") + ylab("-Log10 p-value") +
  theme(legend.position = "bottom", legend.box = "vertical")

p2b <- ggplot(data=res2, aes(x=log2FoldChange, 
                             y=-log10(padj), 
                             col=diffexpressed, 
                             #label=delabel
                             )) + 
  geom_point() +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())  +
  geom_vline(xintercept=c(-0.58, 0.58), col="red", linetype = "dashed") +
  geom_hline(yintercept=-log10(0.05), col="red", linetype = "dashed") + 
 # geom_label_repel(fill = "white", 
#                   max.overlaps = Inf, 
#                   aes(hjust=0, segment.size = 0.75),
#                   min.segment.length = unit(0,"lines"),
#                   color = "black",
                   #nudge_x = -5,
#                   nudge_y = 60) +
  scale_color_manual(values=c("#619CFF", "black", "#00BA38")) +
  labs(color = "Differentially expressed") +
  ylim(0,150) +
  xlim(-15,15) +
  ggtitle("AngG3 vs Agcd-deletion11") +
  xlab("Log2 Fold Change") + ylab("-Log10 p-value") +
  theme(legend.position = "bottom", legend.box = "vertical")


# List of genes to highlight
highlight_genes2 <- c("XP_061501637", "XP_061505062", "XP_061505062", "XP_061497888", "XP_003435766",    
                      "XP_320163","XP_061497007","XP_061519826","XP_061506202", "XP_061498453")       

# Create a new column to indicate if the gene should be highlighted
res2$delabel2 <- ifelse(res3$refseq %in% highlight_genes2, res2$name, NA)

# Create the volcano plot with highlights
p2c <- ggplot(data=res2, aes(x=log2FoldChange, 
                             y=-log10(padj), 
                             color=diffexpressed,
                             label=delabel2)) + 
  geom_point() +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())  +
  geom_vline(xintercept=c(-0.58, 0.58), col="red", linetype = "dashed") +
  geom_hline(yintercept=-log10(0.05), col="red", linetype = "dashed") + 
  geom_label_repel(fill = "white", 
                   force = 15,
                   max.overlaps = Inf, 
                   aes(hjust=0, segment.size = 0.75),
                   min.segment.length = unit(0,"lines"),
                   color = "black",
                   nudge_y = 30,
                   direction = "both") +
  scale_color_manual(values=c("#619CFF", "black", "#00BA38")) +
  labs(color = "Differentially expressed") +
  ylim(0,200) +
  xlim(-20,20) +
  ggtitle("AngG3 vs Agcd-deletion11") +
  xlab("Log2 Fold Change") + ylab("-Log10 p-value") +
  theme(legend.position = "bottom", legend.box = "vertical")


# Plot the volcano plot
#plot(p2a)
#plot(p2b)
plot(p2c)



# Define differentially expressed genes with broader thresholds for coloring
res3$diffexpressed <- "NO"
res3$diffexpressed[res3$log2FoldChange > 0.58 & res3$padj < 0.05] <- "UP"
res3$diffexpressed[res3$log2FoldChange < -0.58 & res3$padj < 0.05] <- "DOWN"

# Create a new column "delabel" to contain the name of differentially expressed genes with stricter thresholds for labeling
res3$delabel <- NA
res3$delabel <- ifelse(res3$log2FoldChange > 2.98 & res3$padj < 2.93E-06, res3$name, res3$delabel)
res3$delabel <- ifelse(res3$log2FoldChange < -3.86 & res3$padj < 9.91E-20, res3$name, res3$delabel)


# Create the volcano plot
p3a <- ggplot(data=res3, aes(x=log2FoldChange, 
                            y=-log10(padj), 
                            col=diffexpressed, 
                            label=delabel)) + 
  geom_point() +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())  +
  geom_vline(xintercept=c(-0.58, 0.58), col="red", linetype = "dashed") +
  geom_hline(yintercept=-log10(0.05), col="red", linetype = "dashed") + 
  geom_label_repel(fill = "white", 
                   max.overlaps = Inf, 
                   aes(hjust=0, segment.size = 0.75),
                   min.segment.length = unit(0,"lines"),
                   color = "black",
                   #nudge_x = -5,
                   nudge_y = 60) +
  scale_color_manual(values=c("#619CFF", "black", "#00BA38")) +
  labs(color = "Differentially expressed") +
  ylim(0,150) +
  xlim(-15,15) +
  ggtitle("AngG3 vs Agcd-deletion11") +
  xlab("Log2 Fold Change") + ylab("-Log10 p-value") +
  theme(legend.position = "bottom", legend.box = "vertical")

p3b <- ggplot(data=res3, aes(x=log2FoldChange, 
                             y=-log10(padj), 
                             col=diffexpressed, 
                             #label=delabel
                             )) + 
  geom_point() +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())  +
  geom_vline(xintercept=c(-0.58, 0.58), col="red", linetype = "dashed") +
  geom_hline(yintercept=-log10(0.05), col="red", linetype = "dashed") + 
#  geom_label_repel(fill = "white", 
#                   max.overlaps = Inf, 
#                   aes(hjust=0, segment.size = 0.75),
#                   min.segment.length = unit(0,"lines"),
#                   color = "black",
                   #nudge_x = -5,
#                   nudge_y = 60) +
  scale_color_manual(values=c("#619CFF", "black", "#00BA38")) +
  labs(color = "Differentially expressed") +
  ylim(0,200) +
  xlim(-20,20) +
  ggtitle("AngG3 vs Agcd-deletion11") +
  xlab("Log2 Fold Change") + ylab("-Log10 p-value") +
  theme(legend.position = "bottom", legend.box = "vertical")

# List of genes to highlight
highlight_genes3 <- c("XP_061505062", "XP_061511274", "XP_309133","XP_319873", "XP_061501637",       
                      "XP_311486","XP_061498453", "XP_061500553", "XP_061509182", "XP_061515653")             

# Create a new column to indicate if the gene should be highlighted
res3$delabel2 <- ifelse(res3$refseq %in% highlight_genes3, res2$name, NA)

p3c <- ggplot(data=res3, aes(x=log2FoldChange, 
                             y=-log10(padj), 
                             color=diffexpressed,
                             label=delabel2)) + 
  geom_point() +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())  +
  geom_vline(xintercept=c(-0.58, 0.58), col="red", linetype = "dashed") +
  geom_hline(yintercept=-log10(0.05), col="red", linetype = "dashed") + 
  geom_label_repel(fill = "white", 
                   force = 5,
                   max.overlaps = Inf, 
                   aes(hjust=0, segment.size = 0.75),
                   min.segment.length = unit(0,"lines"),
                   color = "black",
                   nudge_y = 15) +
  scale_color_manual(values=c("#619CFF", "black", "#00BA38")) +
  labs(color = "Differentially expressed") +
  ylim(0,200) +
  xlim(-20,20) +
  ggtitle("AngG3 vs Agcd-deletion11") +
  xlab("Log2 Fold Change") + ylab("-Log10 p-value") +
  theme(legend.position = "bottom", legend.box = "vertical")


# Plot the volcano plot
plot(p3a)
plot(p3b)




#--------------------
###    HEATMAP   ### 
#-------------------

res4 <- results(dds, contrast = list("condition_A_vs_B", 
                                     "condition_A_vs_C", 
                                     "condition_B_vs_C"))


# Load necessary libraries
library(DESeq2)
library(pheatmap)

# Check if normalized_counts has row names
if (is.null(rownames(normalized_counts))) {
  stop("normalized_counts does not have row names. Ensure it contains gene identifiers.")
}

# Filter significant genes (adjust the threshold as needed)
sig_genes1 <- res1[res1$padj < 0.05 & abs(res1$log2FoldChange) > 0.58, ]
sig_genes2 <- res1[res2$padj < 0.05 & abs(res2$log2FoldChange) > 0.58, ]
sig_genes3 <- res1[res3$padj < 0.05 & abs(res3$log2FoldChange) > 0.58, ]

# Extract normalized counts for significant genes
normalized_counts <- counts(dds, normalized = TRUE)

# Remove the specified substring from row names
rownames(normalized_counts) <- sub(".*(?=XP)", "", rownames(normalized_counts), perl = TRUE)
rownames(normalized_counts) <- sub(".*(?=YP)", "", rownames(normalized_counts), perl = TRUE)
rownames(normalized_counts) <- sub("_[^_]*$", "", rownames(normalized_counts), perl = TRUE)

# Extract normalized counts for significant genes
sig_gene_counts <- normalized_counts[rownames(normalized_counts) %in% rownames(sig_genes1) |
                                       rownames(normalized_counts) %in% rownames(sig_genes2) |
                                       rownames(normalized_counts) %in% rownames(sig_genes3), ]

# Remove rows with NA or Inf values
sig_gene_counts <- sig_gene_counts[complete.cases(sig_gene_counts), ]
sig_gene_counts <- sig_gene_counts[apply(sig_gene_counts, 1, function(row) all(is.finite(row))), ]

# Check if sig_gene_counts is empty after filtering
if (nrow(sig_gene_counts) == 0) {
  stop("No significant genes found after filtering for NA/Inf values.")
}

# Define colors for each genotype
genotype_colors <- c("AngG3" = "#619CFF", 
                     "AgNosCd-1"  = "#00BA38",
                     "Agcd-deletion11" = "#F8766D")

# Create a list of annotation colors
ann_colors <- list(genotype = genotype_colors)

# Create a heatmap
heatmap_plot <- pheatmap(sig_gene_counts,
                         cluster_rows = TRUE,
                         cluster_cols = TRUE,
                         show_rownames = FALSE,
                         show_colnames = TRUE,
                         #annotation_col = as.data.frame(colData(dds)[, c("condition", "genotype")]),
                         annotation_col = as.data.frame(colData(dds)[, "genotype", drop = FALSE]),
                         annotation_colors = ann_colors,
                         scale = "row",
                         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                         main = "Heatmap of Differentially Expressed Genes")

