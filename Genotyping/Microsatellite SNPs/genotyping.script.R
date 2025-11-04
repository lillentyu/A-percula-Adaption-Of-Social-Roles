####This R script contains the code to organise, and filter the haplotype matrix based on read depth to construct genotype matrix for sibship analysis.

#by: Lili Vizer - 2025 April

# Libraries
library(tidyverse)
library(tidyr)
library(dplyr)
library(stringr)
library(ape)


# Read in dataset
haplotype_matrix <- read.csv("Aperc_hap_genotype_filt_26Mar24.csv", header = T, check.names = F)
head(haplotype_matrix)

#Locus                                                                        Haplotypes     12PA__1     12PB__2      12S__3     16P1__4      16P2__5
#1 Aperc_c1_16                  1(0.22);2(0.19);3(0.17);4(0.14);5(0.11);6(0.08);7(0.06);8(0.03); 3/1:373,158 1/4:146,105 2/5:230,126 2/5:256,139  3/6:923,332
#2  Aperc_c1_4                                          1(0.29);2(0.24);3(0.22);4(0.19);5(0.06);   2/1:42,11   3/1:42,14    4/1:34,7   4/1:38,12    2/1:90,29
#3  Aperc_c1_6                                  1(0.39);2(0.21);4(0.14);3(0.12);5(0.08);6(0.06); 1/6:317,302     1/1:512 2/3:320,196 2/3:265,172      1/1:415
#4 Aperc_c1_12          1(0.33);2(0.14);3(0.14);7(0.11);4(0.10);5(0.07);6(0.07);8(0.03);9(0.01);      4/4:46 5/2:334,224     1/1:298     1/1:294  1/3:683,291
#5  Aperc_c2_7                                          1(0.49);2(0.18);3(0.15);4(0.10);5(0.08); 1/4:156,126  3/5:240,67  3/1:213,84 2/1:179,115  1/4:413,299
#6 Aperc_c2_44 1(0.28);2(0.18);3(0.10);4(0.08);7(0.08);6(0.07);8(0.07);9(0.06);10(0.04);5(0.04);  3/10:73,51  6/4:282,82     7/7:110  1/2:219,71 3/10:167,134

#In this dataset:
#locus column is the given locus name 
#haplotypes are the alleles and their frequencies across all samples
#and the subsequent columns are the sampled individuals respective haplotypes and their counts for each allele

# No filtering of the haplotype matrix was done, as sibship construction without filtering resulted to most accurate results

#information to match samples to clutches
clutch_data <- read.csv("./Sample2clutch.csv", header = T)

# carry out parantage analysis using Sequoia in R as it allows for multi-allelic data
# need to sort and filter data

# Step 1: Reorganise the data
haplotype_matrix_long <- haplotype_matrix %>%
      pivot_longer(
             cols = -(1:2),
             names_to = "Sample",
             values_to = "Genotype")
head(haplotype_matrix_long)

# Reshape: wide format, sample as rows and locus as columns
haplotype_matrix_wide <- haplotype_matrix_long %>%
       select(Locus, Sample, Genotype) %>%
       pivot_wider(
             names_from = Locus,
             values_from = Genotype)
head(haplotype_matrix_wide)
#write.csv(haplotype_matrix_wide, "Aperc_hap_genotype_matrix_transposed_2025_04_24.csv", row.names = F)

# Extract the allele counts and frequencies
# Generate 2 columns for each locus: one for allele counts and one for frequencies
# Get all locus column names (everything except "Sample")
# Get all columns to split (excluding Sample)
data_cols <- setdiff(names(haplotype_matrix_wide), "Sample")

# Split each data column into two new columns: _allele and _freq
haplotype_matrix_split <- haplotype_matrix_wide %>%
  pivot_longer(-Sample, names_to = "Locus", values_to = "Value") %>%
  separate(Value, into = c("Allele", "Freq"), sep = ":", fill = "right") %>%
  pivot_wider(names_from = Locus, values_from = c(Allele, Freq), names_sep = "_")
head(haplotype_matrix_split)
str(haplotype_matrix_split)

# Creating 2 separate datasets, one for alleles and one for allele frequencies
# Alleles only (keeping Sample column)
allele_matrix <- haplotype_matrix_split %>%
  select(Sample, starts_with("Allele_")) %>%
  rename_with(~ gsub("^Allele_", "", .), starts_with("Allele_"))
#write.csv(allele_matrix, "Aperc_hap_genotype_matrix_alleles2025_04_24.csv", row.names = F)

# Frequencies only (keeping Sample column)
freq_matrix <- haplotype_matrix_split %>%
  select(Sample, starts_with("Freq_")) %>%
  rename_with(~ gsub("^Freq_", "", .), starts_with("Freq_"))
#write.csv(freq_matrix, "Aperc_hap_genotype_matrix_freqs2025_04_24.csv", row.names = F)


# Modify allele matrix by replace / with a space
geno_matrix <- allele_matrix %>%
  mutate(across(-Sample, ~ gsub("/", " ", .)))
nrow(geno_matrix)  # Check number of rows = 36
ncol(geno_matrix)  # Check number of columns = 61



# Remove everything after the colon (if present)
geno_clean <- geno_matrix %>%
  mutate(across(-Sample, ~ sub(":.*", "", .)))

# Add in clutch information to each sample as _<clutch_ID>
head(clutch_data)
head(geno_clean)
# Merge the datasets of genotype and clutch data by Sample
geno_clean <- geno_clean %>%
  left_join(clutch_data, by = "Sample") %>%
  mutate(Sample = paste0(Sample, "__", Clutch_ID)) %>%
  select(-Clutch_ID)  # Remove Clutch_ID column

geno_long <- geno_clean %>%
  pivot_longer(-Sample, names_to = "Locus", values_to = "Alleles") %>%
  separate(Alleles, into = c("A1", "A2"), sep = " ") %>%
  mutate(
    Genotype = paste0(A1, A2)
  ) %>%
  select(Sample, Locus, Genotype) %>%
  pivot_wider(names_from = Locus, values_from = Genotype)

# Set rownames for clustering
geno_long <- as.data.frame(geno_long)
rownames(geno_long) <- geno_long$Sample
geno_long$Sample <- NULL

geno_long[is.na(geno_long)] <- "00"  # or "NA" if you prefer

dist_mat <- dist.gene(as.matrix(geno_long), method = "pairwise", pairwise.deletion = TRUE)
hc <- hclust(dist_mat, method = "ward.D2")
plot(hc, main = "Relatedness Dendrogram of Fish Clutches", xlab = "Sample", sub = "")



## Running the dendogram with corrected clutch IDs
clutch_data_corrected <- read.csv("Sample2clutch_corrected.csv", header = T)

# Combine the geno matrix data with the corrected clutch data
# Remove everything after the colon (if present)
geno_clean <- geno_matrix %>%
  mutate(across(-Sample, ~ sub(":.*", "", .)))

# Add in clutch information to each sample as _<clutch_ID>
head(clutch_data_corrected)
head(geno_clean)
# Merge the datasets of genotype and clutch data by Sample
geno_clean_corrected <- geno_clean %>%
  left_join(clutch_data_corrected, by = "Sample") %>%
  mutate(Sample = paste0(Sample, "__", Clutch_ID)) %>%
  select(-Clutch_ID)  # Remove Clutch_ID column

geno_long_corrected <- geno_clean_corrected %>%
  pivot_longer(-Sample, names_to = "Locus", values_to = "Alleles") %>%
  separate(Alleles, into = c("A1", "A2"), sep = " ") %>%
  mutate(
    Genotype = paste0(A1, A2)
  ) %>%
  select(Sample, Locus, Genotype) %>%
  pivot_wider(names_from = Locus, values_from = Genotype)

# Set rownames for clustering
geno_long_corrected <- as.data.frame(geno_long_corrected)
rownames(geno_long_corrected) <- geno_long_corrected$Sample
geno_long_corrected$Sample <- NULL
geno_long_corrected[is.na(geno_long_corrected)] <- "00"  # or "NA" if you prefer
dist_mat_corrected <- dist.gene(as.matrix(geno_long_corrected), method = "pairwise", pairwise.deletion = TRUE)
hc_corrected <- hclust(dist_mat_corrected, method = "ward.D2")
plot(hc_corrected, main = "Relatedness Dendrogram of Fish Clutches Corrected", xlab = "Sample", sub = "")


library(gridGraphics)
library(grid)
library(cowplot)

# Function to capture base plot as a grob
grab_plot <- function(expr) {
  grid.newpage()
  grid.echo(expr)
  grid.grab()
}

# Capture the dendrogram plots
plot1 <- grab_plot(function() plot(hc, main = "", xlab = "Sample", sub = ""))
plot2 <- grab_plot(function() plot(hc_corrected, main = "", xlab = "Sample", sub = ""))

# Combine using cowplot
combined_plot <- plot_grid(plot1, plot2, labels = c("A", "B"), ncol = 1, align = "v", rel_heights = c(1, 1))
print(combined_plot)

# Save the combined plot
ggsave("combined_dendrogram_plot.png", plot = combined_plot, width = 8, height = 8, dpi = 1200, bg = "white", units = "in", limitsize = FALSE)
