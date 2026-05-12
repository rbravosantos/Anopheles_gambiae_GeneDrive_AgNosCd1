rm(list = ls())

# Load necessary libraries
library(tximport)
library(readr)
library(rhdf5)
library(ggplot2)
library(DESeq2)
library(readxl)

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


metadata$genotype <- factor(metadata$genotype, 
                            levels = c("WT", "GD", "KO"),
                            labels = c("AngG3", "AgNosCd-1", "Agcd-deletion11"))

# Update rownames of metadata to match colnames of txi$counts
rownames(metadata) <- metadata$sample

# Filter txi to exclude the counts for the sample "D4"
txi$counts <- txi$counts[, colnames(txi$counts) %in% rownames(metadata)]
txi$abundance <- txi$abundance[, colnames(txi$abundance) %in% rownames(metadata)]
txi$length <- txi$length[, colnames(txi$length) %in% rownames(metadata)]

# Create a new group column in metadata
metadata$sex_diet <- factor(paste0(metadata$sex, metadata$diet))
metadata$group <- factor(paste0(metadata$sex, metadata$diet, metadata$condition))

rownames(metadata) <- colnames(txi$counts) # to set up every row for every sample

metadata$group <- factor(paste0(metadata$sex, metadata$diet, metadata$condition))

dds <- DESeqDataSetFromTximport(txi, 
                                colData = metadata, 
                                design = ~ group)


# normalize the read counts,estimate dispersions, and fit the linear model
dds <- DESeq(dds,betaPrior = T)

# Normalization factors
sizeFactors(dds)

# Plot the dispersions
plotDispEsts(dds)

# Plot PCA ---------------------------------------------------------------------
#library(ggpubr)
vsd <- vst(dds, blind = FALSE)

unique(colData(vsd)$genotype)
levels(colData(vsd)$diet)

colData(vsd)$group_shape <- with(colData(vsd),
                                 ifelse(sex == "Male", "Male",
                                        ifelse(sex == "Female" & diet == "Sugarfed", "Sugarfed Female",
                                               ifelse(sex == "Female" & diet == "Bloodfed", "Bloodfed Female", NA))))

colData(vsd)$group_shape <- factor(colData(vsd)$group_shape,
                                   levels = c("Male", "Sugarfed Female", "Bloodfed Female")
)

pca_final <- plotPCA(vsd, intgroup = c("genotype", "group_shape")) + 
  aes(color = genotype, shape = group_shape) +
  scale_color_manual("Genotype",
                     values = c("Agcd-deletion11" = "#619CFF", 
                                "AgNosCd-1"  = "#00BA38",
                                "AngG3" = "#F8766D")) +
  scale_shape_manual("Sex-Diet",
                     values = c("Male" = 16,
                                "Sugarfed Female" = 17,
                                "Bloodfed Female" = 15)) +
  coord_fixed() +
  xlim(-40, 40) +
  ylim(-30, 30)


