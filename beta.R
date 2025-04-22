# Install necessary packages (only run if not already installed)
install.packages("tidyverse")
install.packages("vegan")
install.packages("devtools")
library(devtools)
devtools::install_github("jbisanz/qiime2R")

# Load the libraries
library(tidyverse)
library(vegan)
library(qiime2R)

# Set the working directory to the data location
setwd("~/Desktop/ANSC516-repro/data/moving-pictures")

# List all files to check the directory
list.files()

# Create output directory if it doesn't exist
if (!dir.exists("output")) dir.create("output")

# Load metadata
metadata2 <- read.delim("Assignment7/sample-metadata.tsv", sep = "\t", header = TRUE, quote = "", stringsAsFactors = FALSE)
metadata2 <- metadata2[-1, ]  # Remove the first row
metadata <- read_q2metadata("Assignment7/sample-metadata.tsv")

# Rename columns in the metadata
colnames(metadata)[3] <- "body.site"
colnames(metadata)[8] <- "reported.antibiotic.usage"
colnames(metadata)[9] <- "days.since.experiment.start"

# Set row names for metadata
row.names(metadata) <- metadata$SampleID

# Load Bray-Curtis, Unweighted UniFrac, and Weighted UniFrac PCoA results
bc_PCoA <- read_qza("Assignment7/bray_curtis_pcoa_results.qza")
wu_PCoA <- read_qza("Assignment7/weighted_unifrac_pcoa_results.qza")
uwu_PCoA <- read_qza("Assignment7/unweighted_unifrac_pcoa_results.qza")
jaccard_PCoA <- read_qza("Assignment7/jaccard_pcoa_results.qza")

# Concatenate Bray-Curtis data with metadata
# Remove leading and trailing whitespace from all columns in the metadata
metadata <- metadata %>% 
  mutate(across(everything(), ~str_trim(., side = "both")))
bc_meta <- bc_PCoA$data$Vectors %>% 
  select(SampleID, PC1, PC2) %>% 
  inner_join(metadata, by = c("SampleID" = "SampleID"))

# Create Bray-Curtis PCoA plot
ggplot(bc_meta, aes(x = PC1, y = PC2, color = treatment)) +
  geom_point() + 
  theme_q2r() +
  xlab(paste0("PC1 (", round(100 * bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) + 
  ylab(paste0("PC2 (", round(100 * bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values = c("Blue", "Black", "Green", "Gray"), name = "Treatment") +
  ggsave("Assignment7/output/Bray_Curtis_PCoA.pdf", height = 3, width = 4.5)

# Concatenate Weighted UniFrac data with metadata
wu_meta <- wu_PCoA$data$Vectors %>% 
  select(SampleID, PC1, PC2) %>% 
  inner_join(metadata, by = c("SampleID" = "SampleID"))

# Create Weighted UniFrac PCoA plot
ggplot(wu_meta, aes(x = PC1, y = PC2, color = treatment)) +
  geom_point() + 
  theme_q2r() +
  xlab(paste0("PC1 (", round(100 * wu_PCoA$data$ProportionExplained[1], digits = 2), "%)")) + 
  ylab(paste0("PC2 (", round(100 * wu_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values = c("Blue", "Black", "Green", "Gray"), name = "treatment") +
  ggsave("Assignment7/output/Weighted_Unifrac_PCoA.pdf", height = 3, width = 4.5)

# Concatenate Unweighted UniFrac data with metadata
uwu_meta <- uwu_PCoA$data$Vectors %>% 
  select(SampleID, PC1, PC2) %>% 
  inner_join(metadata, by = c("SampleID" = "SampleID"))

# Create Unweighted UniFrac PCoA plot
ggplot(uwu_meta, aes(x = PC1, y = PC2, color = treatment)) +
  geom_point() + 
  theme_q2r() +
  xlab(paste0("PC1 (", round(100 * uwu_PCoA$data$ProportionExplained[1], digits = 2), "%)")) + 
  ylab(paste0("PC2 (", round(100 * uwu_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values = c("Blue", "Black", "Green", "Gray"), name = "Treatment") +
  ggsave("Assignment7/output/Unweighted_Unifrac_PCoA.pdf", height = 3, width = 4.5)

# Concatenate Jaccard data with metadata
jaccard_meta <- jaccard_PCoA$data$Vectors %>% 
  select(SampleID, PC1, PC2) %>% 
  inner_join(metadata, by = c("SampleID" = "SampleID"))

# Create Jaccard PCoA plot
ggplot(jaccard_meta, aes(x = PC1, y = PC2, color = treatment)) +
  geom_point() + 
  theme_q2r() +
  xlab(paste0("PC1 (", round(100 * jaccard_PCoA$data$ProportionExplained[1], digits = 2), "%)")) + 
  ylab(paste0("PC2 (", round(100 * jaccard_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values = c("Blue", "Black", "Green", "Gray"), name = "Treatment") +
  ggsave("Assignment7/output/Jaccard_PCoA.pdf", height = 3, width = 4.5)

# Run PERMANOVA for each distance matrix (Bray-Curtis, Jaccard, Weighted UniFrac, and Unweighted UniFrac)
bc_dist_mat <- read_qza("Assignment7/bray_curtis_distance_matrix.qza")
wu_dist_mat <- read_qza("Assignment7/weighted_unifrac_distance_matrix.qza")
uwu_dist_mat <- read_qza("Assignment7/unweighted_unifrac_distance_matrix.qza")
jaccard_dist_mat <- read_qza("Assignment7/jaccard_distance_matrix.qza")

# Convert distance matrices to matrices for PERMANOVA
bc_dm <- as.matrix(bc_dist_mat$data)
wu_dm <- as.matrix(wu_dist_mat$data)
uwu_dm <- as.matrix(uwu_dist_mat$data)
jaccard_dm <- as.matrix(jaccard_dist_mat$data)

# Match sample names in metadata
metadata_sub <- metadata[match(rownames(bc_dm), metadata$SampleID), ]

# Perform PERMANOVA on each distance matrix
PERMANOVA_bc <- adonis2(bc_dm ~ treatment, data = metadata_sub)
PERMANOVA_wu <- adonis2(wu_dm ~ treatment, data = metadata_sub)
PERMANOVA_uwu <- adonis2(uwu_dm ~ treatment, data = metadata_sub)
PERMANOVA_jaccard <- adonis2(jaccard_dm ~ treatment, data = metadata_sub)

# Save PERMANOVA results
write.table(PERMANOVA_bc, "Assignment7/output/PERMANOVA_Bray_Curtis.csv", sep = ",", row.names = TRUE)
write.table(PERMANOVA_wu, "Assignment7/output/PERMANOVA_Weighted_Unifrac.csv", sep = ",", row.names = TRUE)
write.table(PERMANOVA_uwu, "Assignment7/output/PERMANOVA_Unweighted_Unifrac.csv", sep = ",", row.names = TRUE)
write.table(PERMANOVA_jaccard, "Assignment7/output/PERMANOVA_Jaccard.csv", sep = ",", row.names = TRUE)

# Perform pairwise comparisons (optional, similar to above)
pairwise_bc <- pairwise.adonis2(bc_dm ~ treatment, data = metadata_sub)
pairwise_wu <- pairwise.adonis2(wu_dm ~ treatment, data = metadata_sub)
pairwise_uwu <- pairwise.adonis2(uwu_dm ~ treatment, data = metadata_sub)
pairwise_jaccard <- pairwise.adonis2(jaccard_dm ~ treatment, data = metadata_sub)

# Save pairwise results
write.table(pairwise_bc, "Assignment7/output/Pairwise_Bray_Curtis.csv", sep = ",", row.names = TRUE)
write.table(pairwise_wu, "Assignment7/output/Pairwise_Weighted_Unifrac.csv", sep = ",", row.names = TRUE)
write.table(pairwise_uwu, "Assignment7/output/Pairwise_Unweighted_Unifrac.csv", sep = ",", row.names = TRUE)
write.table(pairwise_jaccard, "Assignment7/output/Pairwise_Jaccard.csv", sep = ",", row.names = TRUE)



# Install and load required packages
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("vegan")) install.packages("vegan")
if (!require("devtools")) install.packages("devtools")
if (!require("qiime2R")) devtools::install_github("jbisanz/qiime2R")

# Load libraries
library(tidyverse)
library(vegan)
library(qiime2R)

# Set working directory
setwd("~/Desktop/ANSC516-repro/data/moving-pictures")
if (!dir.exists("output")) dir.create("output")
getwd()

# Load metadata and clean whitespace
metadata <- read_q2metadata("Assignment7/sample-metadata.tsv") %>%
  mutate(across(everything(), trimws))
colnames(metadata)[3] <- "body.site"
colnames(metadata)[8] <- "reported.antibiotic.usage"
colnames(metadata)[9] <- "days.since.experiment.start"
row.names(metadata) <- metadata$SampleID

# Define colors and plotting column
body_colors <- c("Black", "Blue", "Green", "Gray")
my_column <- "treatment"

### Bray-Curtis Plot ###
bc_PCoA <- read_qza("Assignment7/bray_curtis_pcoa_results.qza")
bc_meta <- bc_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = "SampleID")

ggplot(bc_meta, aes(x = PC1, y = PC2, color = get(my_column))) +
  geom_point(aes(shape = treatment), size = 3) +
  stat_ellipse(level = 0.95, type = "t") +
  theme_q2r() +
  xlab(paste0("PC1 (", round(100 * bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100 * bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values = body_colors, name = my_column) +
  ggtitle("Bray-Curtis PCoA")
ggsave("Assignment7/output/Bray_Curtis_ellipse.pdf", height = 3, width = 4.5)

### Jaccard Plot ###
jac_PCoA <- read_qza("Assignment7/jaccard_pcoa_results.qza")
jac_meta <- jac_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = "SampleID")

ggplot(jac_meta, aes(x = PC1, y = PC2, color = get(my_column))) +
  geom_point(aes(shape = treatment), size = 3) +
  stat_ellipse(level = 0.95, type = "t") +
  theme_q2r() +
  xlab(paste0("PC1 (", round(100 * jac_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100 * jac_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values = body_colors, name = my_column) +
  ggtitle("Jaccard PCoA")
ggsave("Assignment7/output/Jaccard_ellipse.pdf", height = 3, width = 4.5)

### Weighted UniFrac Plot ###
wuni_PCoA <- read_qza("Assignment7/weighted_unifrac_pcoa_results.qza")
wuni_meta <- wuni_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = "SampleID")

ggplot(wuni_meta, aes(x = PC1, y = PC2, color = get(my_column))) +
  geom_point(aes(shape = treatment), size = 3) +
  stat_ellipse(level = 0.95, type = "t") +
  theme_q2r() +
  xlab(paste0("PC1 (", round(100 * wuni_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100 * wuni_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values = body_colors, name = my_column) +
  ggtitle("Weighted UniFrac PCoA")
ggsave("Assignment7/output/Weighted_Unifrac_ellipse.pdf", height = 3, width = 4.5)

### Unweighted UniFrac Plot ###
uwuni_PCoA <- read_qza("Assignment7/unweighted_unifrac_pcoa_results.qza")
uwuni_meta <- uwuni_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = "SampleID")

ggplot(uwuni_meta, aes(x = PC1, y = PC2, color = get(my_column))) +
  geom_point(aes(shape = treatment), size = 3) +
  stat_ellipse(level = 0.95, type = "t") +
  theme_q2r() +
  xlab(paste0("PC1 (", round(100 * uwuni_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100 * uwuni_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values = body_colors, name = my_column) +
  ggtitle("Unweighted UniFrac PCoA")
ggsave("Assignment7/output/Unweighted_Unifrac_ellipse.pdf", height = 3, width = 4.5)

cat("All plots saved to the 'output' folder.\n")

