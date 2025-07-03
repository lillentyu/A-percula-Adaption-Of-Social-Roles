## This script carries out targetted exploration of gene sets with known roles ##
## it uses A percula geneID and the gene`s name to generate figures ##

######################################Set up######################################
### Load libraries
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")
library(DESeq2)

#BiocManager::install("genefilter")
library(genefilter)

#BiocManager::install("DEGreport")
library(DEGreport)

#BiocManager::install("vsn")
library(vsn)

#BiocManager::install("apeglm")
library(apeglm)

#BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)

#install.packages("devtools")
library(devtools)

#install_github("stephens999/ashr")
library(ashr)

#install.packages("factoextra")
library(factoextra)

#install.packages("tidyverse")
library(tidyverse)

library(stringr)
library(dplyr)
library(RColorBrewer)
library(pheatmap)
library(dendextend)
library(ggplot2)
library(gplots)
library(gridExtra)
library(ggrepel)
library(vegan) #for adonis


### Clean up the workspace
rm(list = ls()) 

### Set directory
setwd(wd)

### Read file
data <- read.table("../Data/nemo_cleaned_count_data.txt", header = TRUE, sep = "\t", dec = ".")
colnames(data) <- gsub("^X", "", colnames(data))
head(data)
# File looks like this: (first column is your gene, following columns are counts for each sample)
#target_id   1S  10S  11S 13S  14S  15S   2S  20S  21S  23S  26S   3S   5S   6S   7S  1P1 10P1 11P1 13P1 14P1 15P1
#1 ENSAPEG00000000002   43   46   57  49   34   53   49  109  105   87   56   17   35   41   39   31   40   80   58   57   37
#2 ENSAPEG00000000003    0    0    0   0    0    0    0    0    0    0    0    1    0    0    0    0    0    0    0    0    0
#3 ENSAPEG00000000004    4    6   10   1    3    5    1    8    4    4    9    4    7    2    4    7    1    5    4    7    3
#4 ENSAPEG00000000005    1    2    7   5    4    7    6    3    8    6    9   11    0    3    3    4    1    1    1    7    6


### Format data frame
rownames(data) <- data$target_id
data <- data[, -1]

### Get counts
countData <- data # %>% round() 

### Read metadata
sampleInfo <- read.csv("../Data/metadata.csv", header = TRUE)
rownames(sampleInfo) <- sampleInfo$Sample_ID
sampleInfo$Sample_ID <- NULL

### Verify metadata
all(rownames(sampleInfo) == colnames(countData))

### Make data factors
sampleInfo$fish_type2 <- factor(sampleInfo$fish_type2)
sampleInfo$fish_type <- factor(sampleInfo$fish_type)
sampleInfo$clutch_ID <- factor(sampleInfo$clutch_ID)


### Create DeSeq object
dds <- DESeqDataSetFromMatrix(countData = as.matrix(countData), colData = sampleInfo, design = ~ clutch_ID + fish_type2)

### Pre-filtering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

### Get normalized counts
dds <- estimateSizeFactors(dds)  #normalization using median-of-ratios

normCounts <- counts(dds, normalized = TRUE) %>% data.frame()
colnames(normCounts) <- gsub("^X", "", colnames(normCounts))
normCounts$target_id <- rownames(normCounts)
normCounts <- normCounts[, c("target_id", names(normCounts)[-ncol(normCounts)])]
rownames(normCounts) <- NULL

### Transform data
vsd <- vst(dds, blind = TRUE) #variance stabilization across mean counts

### Input is matrix of transformed values
vsdMatrix <- assay(vsd)  #extracts transformed values in the matrix to be used downstream


######################################Plot PCAs for each gene set######################################
### Read files with gene sets

dataFiles <- c("percula.appetite.regulation.txt",
               "percula.digestion.bile.acids.txt",
               "percula.digestion.digestive.enzyme.carbohydrase.txt",
               "percula.digestion.digestive.enzyme.lipase.txt",
               "percula.digestion.digestive.enzyme.protease.txt",
               "percula.gastrointestinal.function.txt",
               "percula.hormone.signaling.corticoids.txt",
               "percula.hormone.signaling.thyroid.hormones.txt",
               "percula.metabolism.beta.oxidation.txt",
               "percula.metabolism.cholesterol.biosynthesis.txt",
               "percula.metabolism.fatty.acid.synthesis.txt",
               "percula.metabolism.glycolysis.txt",
               "percula.metabolism.lactic.fermentation.txt",
               "percula.metabolism.pdh.complex.txt",
               "percula.metabolism.tca.cycle.txt",
               "percula.osmoregulation.txt",
               "percula.osmoregulation.coupled.to.acid.base.regulation.txt",
               "percula.osmoregulation.coupled.to.ammonium.transport.txt",
               "percula.osmoregulation.and.permeability.tight.junctions.txt",
               "percula.ossification.txt",
               "percula.pigmentation.components.of.melanosomes.txt",
               "percula.pigmentation.iridophores.txt",
               "percula.pigmentation.leucophores.txt",
               "percula.pigmentation.melanocyte.development.txt",
               "percula.pigmentation.melanogenesis.regulation.txt",
               "percula.pigmentation.melanophore.development.txt",
               "percula.pigmentation.melanosome.biogenesis.txt",
               "percula.pigmentation.melanosome.transport.txt",
               "percula.pigmentation.pigment.cell.differentiation.txt",
               "percula.pigmentation.pteridine.synthesis.txt",
               "percula.pigmentation.xanthophore.development.txt",
               "percula.visual.perception.opsins.txt",
               "percula.visual.perception.phototransduction.txt")


### Loop through the files and add significance using adonis PERMANOVA 
for (file in dataFiles) {
  ### Read data
  data <- read.table(file.path("../Data/gene_lists", file), header = TRUE, sep = "\t")
 
  ### Filter vsd_matrix
  vsdDf <- as.data.frame(vsdMatrix)
  vsdDf <- rownames_to_column(vsdDf, var = "target_id") 
  geneCounts <- vsdDf %>% filter(target_id %in% data$target_id)
  row.names(geneCounts) <- geneCounts[, 1]
  geneCounts <- geneCounts[, -1]
  
  ### Perform PCA clustering based on specific gene set
  pca <- prcomp(t(geneCounts), scale. = T)
  percentVar <- pca$sdev^2 / sum(pca$sdev^2)
  df <- cbind(sampleInfo, pca$x)
  df$fish_type2 <- factor(df$fish_type2, levels = c("S", "P1", "P2"))
  
  ### Make PERMANOVA using adonis
  pca_adonis <- prcomp(t(assay(vsd)), center = TRUE, scale. = FALSE)
  treat_pca <- adonis2(pca$x ~ 
                         vsd$fish_type2 +
                         vsd$clutch_ID,
                       method = 'eu')
  
  ### Extract significant p-value from the adonis results
  p_value <- treat_pca$`Pr(>F)`[1]
  
  ### Plot PC1 vs PC2
  my_colors <- c("S" = "grey40", "P1" = "magenta1", "P2" = "#b700ff")
  
  pca_plot <- ggplot(df, aes(x = PC1, y = PC2, color = fish_type2, group = fish_type2)) + 
    geom_point(size = 4) + 
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.5, color = "lightgrey") + 
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.5, color = "lightgrey")+
    scale_color_manual(values = my_colors) + 
    stat_ellipse(lwd = 0.2) +
    xlab(sprintf("PC1 (%.2f%%)", percentVar[1] * 100)) +
    ylab(sprintf("PC2 (%.2f%%)", percentVar[2] * 100)) +
    theme_bw() +
    theme(axis.title = element_text(colour = "black", size = 13),
          axis.text = element_text(colour = "black", size = 12),
          axis.ticks = element_line(colour = "black", linewidth = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_rect(fill = NA, colour = "black", size = 0.8))+
    #panel.border = element_rect(fill = NA, colour = "black", size = 0.8)) +
    guides(color = guide_legend(title = "Fish Category"))+
    ggtitle(paste("PERMANOVA fish_category p-value:", p_value))
  
  ### Save each PCA plot"
  ggsave(pca_plot, filename = paste0("../Figures/2025-03-29/PCA_", gsub("percula\\.|\\.txt", "", file), "_scale_sig", ".pdf"), width = 8, height = 6) ## Change plot dimensions here
} 


######################################Plot Expression Levels######################################
### Run ANOVAs to see differences in expression levels (use normalized counts)

### Loop through all gene sets
for (file in dataFiles) {
  ### Read data
  data <- read.table(file.path("../Data/gene_lists", file), header = TRUE, sep = "\t")
  
  ### Format data frames
  geneCounts <- normCounts %>% inner_join(data, by = "target_id") %>%
    dplyr::select(gene_id, everything(), -target_id) %>%
    arrange(gene_id)
 
  ### Make gene IDs unique
  geneCounts <- geneCounts %>%
    mutate(gene_id = ifelse(duplicated(gene_id), paste0(gene_id, "_", row_number()), gene_id))
  
  geneCounts <- geneCounts %>%
    pivot_longer(-gene_id, names_to = "sample", values_to = "Value") %>%
    mutate(treatment = case_when(
      str_detect(sample, "P1$") ~ "P1",
      str_detect(sample, "P2$") ~ "P2",
      str_detect(sample, "S$") ~ "S",
      TRUE ~ "Other"
    )) %>%
    dplyr::select(sample, treatment, gene_id, Value) %>%
    pivot_wider(names_from = gene_id, values_from = Value)
  
  geneCounts$treatment <- factor(geneCounts$treatment, levels = c("P1", "P2", "S"))
  
  ### Run ANOVA and save results to a file
  genes <- colnames(geneCounts)[!(colnames(geneCounts) %in% c("sample", "treatment"))]
  outputFile <- paste0("../Data/anova_expression_levels/2025_07_01/anova_", gsub("percula\\.", "", file))
  
  sink(outputFile)
  
  for (gene in genes) {
    formula <- as.formula(paste(gene, "~ treatment"))
    cat("Gene:", gene, "\n")
    anovaResult <- aov(formula, data = geneCounts)
    anovaSummary <- summary(anovaResult)
    cat("ANOVA Summary:\n", file = outputFile, append = TRUE)
    print(anovaSummary, quote = FALSE, file = outputFile, append = TRUE)
    tukeyResults <- TukeyHSD(anovaResult)
    cat("Tukey HSD Results:\n", file = outputFile, append = TRUE)
    print(tukeyResults, quote = FALSE, file = outputFile, append = TRUE)
    cat("\n", file = outputFile, append = TRUE)
  }
  sink()
}

### Make plots
## Create empty list
geneCountsList <- list()

## Loop through data files
for (file in dataFiles) {
  ### Read the data
  data <- read.table(file.path("../Data/gene_lists", file), header = TRUE, sep = "\t")
  
  ### Format data frames
  geneCounts <- normCounts %>% inner_join(data, by = "target_id") %>%
    dplyr::select(gene_id, everything(), -target_id) %>%
    arrange(gene_id) %>%
    mutate_at(vars(-gene_id), ~ log2(. + 1))
  
  ### Make gene IDs unique
  geneCounts <- geneCounts %>%
    mutate(gene_id = ifelse(duplicated(gene_id), paste0(gene_id, "_", row_number()), gene_id))
  
  geneCounts <- geneCounts %>%
    pivot_longer(-gene_id, names_to = "sample", values_to = "Value") %>%
    mutate(treatment = case_when(
      str_detect(sample, "P1$") ~ "P1",
      str_detect(sample, "P2$") ~ "P2",
      str_detect(sample, "S$") ~ "S",
      TRUE ~ "Other"
    )) %>%
    dplyr::select(sample, treatment, gene_id, Value) %>%
    pivot_wider(names_from = gene_id, values_from = Value)
  
  geneCounts$treatment <- factor(geneCounts$treatment, levels = c("P1", "P2", "S"))
  
  
  ### Save results in the list
  geneCountsList[[file]] <- geneCounts
}

### Generate and save plots for each gene set
for (file in dataFiles) {
  geneCounts <- geneCountsList[[file]]
  
  ### Extract relevant columns for plotting
  cols_to_plot <- colnames(geneCounts)[!(colnames(geneCounts) %in% c("sample", "treatment"))]
  
  ### Split the columns into groups of 12 (or depending in how many plots per page you want)
  col_groups <- split(cols_to_plot, rep(1:(ceiling(length(cols_to_plot) / 12)), each = 12, length.out = length(cols_to_plot)))
  
  ### Generate and save plots for each group 
  for (i in seq_along(col_groups)) {
    cols <- col_groups[[i]]
    
    plots <- list()
    for (col in cols) {
      p <- ggplot(data = geneCounts, aes(x = treatment, y = !!sym(col), color = treatment, fill = treatment)) +
        geom_jitter(position = position_jitter(width = 0.2, height = 0), size = 1, show.legend = FALSE, alpha = 0.5) + # Use geom_jitter for raw data points
        geom_boxplot(alpha = 0.5, outlier.shape = NA, alpha=0.3) + # Use geom_boxplot for box plots, hide outliers to avoid overlap with jitter
        #group_by(treatment) %>%
        #summarize(mean_value = mean(!!sym(col)), aes(x = treatment, y = mean_value), size = 0.5, color = "black", group = 1) +
        #stat_summary(fun = mean, geom = "crossbar", width = 0.3, color = "black", fatten = 0.5)
        # Error bars = mean ± SEM
        #stat_summary( fun.data = mean_se, geom = "errorbar",aes(group = treatment), color = "black", width = 0.2  ) +
        # stat_summary(fun = mean, geom = "point", aes(group = treatment), color = "black", shape = 18, size =2) +
        labs(title = paste(col, ""), x = "", y = "Expression level") +
        scale_color_manual(values = c("S" = "grey40", "P1" = "magenta", "P2" = "#b700ff")) +
        scale_fill_manual(values = c("S" = "grey40", "P1" = "magenta", "P2" = "#b700ff")) +
        theme_classic() + theme(
          plot.title = element_text(size = 9),
          axis.title = element_text(colour = "black", size = 8),
          axis.text = element_text(colour = "black", size = 7),
          axis.ticks = element_line(colour = "black", linewidth = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.position = "none"
        )
      
      plots[[col]] <- p
    }
    
    ### Save plots
    filename <- gsub("percula\\.|\\.txt", "", file)
    pdf(paste0("../Figures/2025_03_01_expression_levels_", filename, "_", i, ".pdf"))
    grid.arrange(grobs = plots, ncol = 3, nrow = 4)
    dev.off()
  }
}


###running DESeq2
dds <- DESeq(dds) # differential expression analysis on gamma-poisson distribution

### the three pairwise contrasts:
P1_con_P2 <- results(dds, contrast = c("fish_type2", "P2", "P1"), alpha = 0.05)
P1_con_S <- results(dds, contrast = c("fish_type2", "S", "P1"), alpha = 0.05)
P2_con_S <- results(dds, contrast = c("fish_type2", "S", "P2"), alpha = 0.05)

######-------getting p-vals & p-adjusted from pairwise contrasts

#P1 vs P2   #133
valP1vsP2 <- data.frame(row.names = row.names(P1_con_P2), pval.P1vsP2 = P1_con_P2$pvalue, padj.P1vsP2 = P1_con_P2$padj)
head(valP1vsP2)
table(complete.cases(valP1vsP2))

#P1 vs S
valP1vsS <- data.frame(row.names = row.names(P1_con_S), pval.P1vsS = P1_con_S$pvalue, padj.P1vsS = P1_con_S$padj)
head(valP1vsS)
table(complete.cases(valP1vsS))

#P2 vs S
valP2vsS <- data.frame(row.names = row.names(P2_con_S), pval.P2vsS = P2_con_S$pvalue, padj.P2vsS = P2_con_S$padj)
head(valP2vsS)
table(complete.cases(valP2vsS))


######-------------make rlogdata and pvals table

rlog=rlogTransformation(dds, blind=TRUE) 
rld=assay(rlog)
head(rld)
length(rld[,1])

rldpvals=cbind(rld,valP1vsP2, valP1vsS, valP2vsS)
head(rldpvals)
dim(rldpvals)
# [1] 20537    51
table(complete.cases(rldpvals))
#FALSE  TRUE 
# 6372 14165 
#write.table(rldpvals, file = "../Data/rldpvals_significance_data2025-03-29.txt", sep = "\t", quote = FALSE, row.names = T)

#################################################################
###################### COMPLEX HEATMAPS #########################
#################################################################

library(ComplexHeatmap)
library(ggplot2)

rldpvals <- read.table(file = "../Data/rldpvals_significance_data2025-03-29.txt", header = TRUE, sep = "\t", check.names = FALSE)

for (file in dataFiles) {
  ### Read data
  data <- read.table(file.path("../Data/gene_lists", file), header = TRUE, sep = "\t")
  
  ### Filter vsd_matrix
  geneCounts <- vsdMatrix[rownames(vsdMatrix) %in% data$target_id, ]
  
  ### Set gene_id as row names
  if (nrow(geneCounts) > 0 && "gene_id" %in% colnames(data)) {
    rownames(geneCounts) <- data$gene_id[match(rownames(geneCounts), data$target_id)]
  }
  
  ### Getting size-ratio data
  growth_data <- read.csv("../Data/Vizer_StrategicGrowth_phenotypicdata_SNPcorrected_20230906.csv", header = T)
  size_ratio <-  growth_data[, c("replicate_ID", "fish_type2", "final_size_ratio")]
  size_ratio
  size_ratio$sampleID <- paste(size_ratio$replicate_ID, size_ratio$fish_type2, sep = "")
  rownames(size_ratio) <- size_ratio$sampleID
  
  
  ### Generate heatmap
  select <- order(rowMeans(geneCounts))
  my_sample_cols <- data.frame(treatment = c(rep("S", 15), rep("P1", 15), rep("P2", 15)))
  row.names(my_sample_cols) <- colnames(geneCounts)
  my_Palette <- c("magenta3", "purple3", "grey50")
  my_colour <- list(treatment = c(P1 = my_Palette[1], P2 = my_Palette[2], S = my_Palette[3]))
  col0 <- colorRampPalette(rev(c("chocolate1", "#FEE090", "grey10", "cyan3", "cyan")))(100)
  
  heatmap_mat <- geneCounts[select, ]
  
  # Get the column names of heatmap_mat
  heatmap_colnames <- colnames(heatmap_mat)
  
  # Reorder rows in size_ratio to match the order of heatmap_mat
  sorted_size_ratio <- size_ratio[heatmap_colnames, "final_size_ratio", drop = FALSE]
  
  
  
  ## Read significance data
  significance_data <- read.table("../Data/rldpvals_significance_data2025-03-29.txt", header = TRUE, row.names = 1)
  colnames(significance_data) <- gsub("^X", "", colnames(significance_data))
  
  ### Filter significance_data
  significant_genes <- significance_data[rownames(significance_data) %in% data$target_id, ]
  
  ### Keep only rows which are significant at pval 0.1 for either P1vsP2 or P1vsS
  filtered_genes <- significant_genes %>%
    filter(padj.P1vsP2 < 0.1 | padj.P1vsS < 0.1)
  
  ### Set gene_id as row names
  if (nrow(filtered_genes) > 0 && "gene_id" %in% colnames(data)) {
    # Match row names of filtered_genes with gene_id column in data
    matched_row_names <- data$gene_id[match(rownames(filtered_genes), data$target_id)]
    
    # Ensure uniqueness of matched_row_names
    unique_row_names <- make.unique(matched_row_names, sep = "_")
    
    # Assign the unique row names to filtered_genes
    rownames(filtered_genes) <- unique_row_names
  }
  
  # Assuming filtered_genes contains the row names that passed the significant threshold
  significant_rows <- rownames(heatmap_mat) %in% rownames(filtered_genes)
  
  # Loop to add '*' for significance between P1vsP2
  for (gene in rownames(heatmap_mat)) {
    if (gene %in% rownames(filtered_genes)) {
      if (!is.na(filtered_genes[gene, "padj.P1vsP2"]) && filtered_genes[gene, "padj.P1vsP2"] < 0.1) {
        # Add double asterisks
        index <- which(rownames(heatmap_mat) == gene)
        rownames(heatmap_mat)[index] <- paste0(gene, "*")
      }
    }
  }
  
  for (gene in rownames(heatmap_mat)) {
    if (gene %in% rownames(filtered_genes)) {
      if (!is.na(filtered_genes[gene, "padj.P1vsS"]) && filtered_genes[gene, "padj.P1vsS"] < 0.1) {
        # Add double asterisks
        index <- which(rownames(heatmap_mat) == gene)
        rownames(heatmap_mat)[index] <- paste0(gene, "**")
      }
    }
  }
  
  
  annotation_col <- data.frame(treatment = c(rep("S", 15), rep("P1", 15), rep("P2", 15)),
                               clutch_ID = c("L6" ,"L2" ,"L3", "L6", "L2", "L3", "L2", "L2", "L2", "L6", "L3", "L3" ,"L3", "L6", 
                                             "L3", "L3", "L6", "L2", "L2", "L3", "L6", "L6", "L6", "L3", "L3", "L6", "L2", "L6", 
                                             "L3", "L6", "L2", "L3", "L6", "L3","L6", "L2", "L3", "L3", "L6", "L2", "L2", "L6", 
                                             "L2", "L2", "L2"))
  
  ann_colors <- list(treatment =  c(P1= "magenta3", P2 ="purple3", S="grey50"),
                     clutch_ID = c(L2 = "brown", L3 = "brown1", L6="burlywood1"))
  
  ha <- HeatmapAnnotation(
    treatment = annotation_col$treatment,
    clutch_ID = annotation_col$clutch_ID,
    col = ann_colors)
  
  # Define custom colors based on the final_size_ratio values
  custom_colors <- ifelse(sorted_size_ratio$final_size_ratio < 0.8, "black", "red")
  
  ha1 <- HeatmapAnnotation(
    ratio = anno_barplot(sorted_size_ratio$final_size_ratio, gp = gpar(fill = custom_colors)))
  
  
  pdf(paste0("../Figures/Heatmap_", gsub("percula\\.|\\.txt", "", file), ".pdf"), width = 8, height = 6)
  print(ComplexHeatmap::pheatmap(geneCounts[select,], 
                                 name="matrix",
                                 cluster_cols= TRUE,
                                 cluster_rows = TRUE,
                                 scale="row", 
                                 color=col0,
                                 top_annotation = ha,
                                 bottom_annotation = ha1,
                                 #annotation_col = annotation_col,
                                 #annotation_colors = ann_colors,
                                 column_split = annotation_col$treatment, 
                                 #column_reorder = list(annotation_col$clutch_ID),
                                 #show_rownames= TRUE, 
                                 row_labels = rownames(heatmap_mat),
                                 show_colnames=T, 
                                 border_color="NA", 
                                 annotation_names_col=F)
  )
  
  dev.off()
}


#################################################################
#~~~~~~~~~~~~~~~COMPLEX HEATMAPS OF SIG GENES ONLY~~~~~~~~~~~~~##
#############only 2 rows of code is changed######################

#running only heatmaps which had the significant rows:

### Read files with gene sets
dataFiles <- c(
  "percula.appetite.regulation.txt",
  "percula.metabolism.glycolysis.txt",
  "percula.metabolism.tca.cycle.txt",
  "percula.ossification.txt",
  "percula.visual.perception.opsins.txt",
  "percula.visual.perception.phototransduction.txt"
)


rldpvals <- read.table(file = "../Data/rldpvals_significance_data2025-03-29.txt", header = TRUE, sep = "\t", check.names = FALSE)

for (file in dataFiles) {
  ### Read data
  data <- read.table(file.path("../Data/gene_lists", file), header = TRUE, sep = "\t")
  
  ### Filter vsd_matrix
  geneCounts <- vsdMatrix[rownames(vsdMatrix) %in% data$target_id, ]
  
  ### Set gene_id as row names
  if (nrow(geneCounts) > 0 && "gene_id" %in% colnames(data)) {
    rownames(geneCounts) <- data$gene_id[match(rownames(geneCounts), data$target_id)]
  }
  
  ### Getting size-ratio data
  growth_data <- read.csv("../Data/Vizer_StrategicGrowth_phenotypicdata_SNPcorrected_20230906.csv", header = T)
  size_ratio <-  growth_data[, c("replicate_ID", "fish_type2", "final_size_ratio")]
  size_ratio
  size_ratio$sampleID <- paste(size_ratio$replicate_ID, size_ratio$fish_type2, sep = "")
  rownames(size_ratio) <- size_ratio$sampleID
  
  
  ### Generate heatmap
  select <- order(rowMeans(geneCounts))
  my_sample_cols <- data.frame(treatment = c(rep("S", 15), rep("P1", 15), rep("P2", 15)))
  row.names(my_sample_cols) <- colnames(geneCounts)
  my_Palette <- c("magenta3", "purple3", "grey50")
  my_colour <- list(treatment = c(P1 = my_Palette[1], P2 = my_Palette[2], S = my_Palette[3]))
  col0 <- colorRampPalette(rev(c("chocolate1", "#FEE090", "grey10", "cyan3", "cyan")))(100)
  
  heatmap_mat <- geneCounts[select, ]
  
  # Get the column names of heatmap_mat
  heatmap_colnames <- colnames(heatmap_mat)
  
  # Reorder rows in size_ratio to match the order of heatmap_mat
  sorted_size_ratio <- size_ratio[heatmap_colnames, "final_size_ratio", drop = FALSE]
  
  
  ## Read significance data
  significance_data <- read.table("../Data/rldpvals_significance_data2025-03-29.txt", header = TRUE, row.names = 1)
  colnames(significance_data) <- gsub("^X", "", colnames(significance_data))
  
  ### Filter significance_data
  significant_genes <- significance_data[rownames(significance_data) %in% data$target_id, ]
  
  ### Keep only rows which are significant at pval 0.1 for either P1vsP2 or P1vsS
  ### Adjust the p-value threshold here as needed - here i`m taking pval < 0.1
  filtered_genes <- significant_genes %>%
    filter(padj.P1vsP2 < 0.1 | padj.P1vsS < 0.1)
  
  ### Set gene_id as row names
  if (nrow(filtered_genes) > 0 && "gene_id" %in% colnames(data)) {
    # Match row names of filtered_genes with gene_id column in data
    matched_row_names <- data$gene_id[match(rownames(filtered_genes), data$target_id)]
    
    # Ensure uniqueness of matched_row_names
    unique_row_names <- make.unique(matched_row_names, sep = "_")
    
    # Assign the unique row names to filtered_genes
    rownames(filtered_genes) <- unique_row_names
  }
  
  # Assuming filtered_genes contains the row names that passed the significant threshold
  significant_rows <- rownames(heatmap_mat) %in% rownames(filtered_genes)
  
  # Loop to add '*' for significance between P1vsP2
  for (gene in rownames(heatmap_mat)) {
    if (gene %in% rownames(filtered_genes)) {
      if (!is.na(filtered_genes[gene, "padj.P1vsP2"]) && filtered_genes[gene, "padj.P1vsP2"] < 0.1) {
        # Add double asterisks
        index <- which(rownames(heatmap_mat) == gene)
        rownames(heatmap_mat)[index] <- paste0(gene, "*")
      }
    }
  }
  
  for (gene in rownames(heatmap_mat)) {
    if (gene %in% rownames(filtered_genes)) {
      if (!is.na(filtered_genes[gene, "padj.P1vsS"]) && filtered_genes[gene, "padj.P1vsS"] < 0.1) {
        # Add double asterisks
        index <- which(rownames(heatmap_mat) == gene)
        rownames(heatmap_mat)[index] <- paste0(gene, "**")
      }
    }
  }
  
  ##**** 2 NEW ROWS below****** ##
  ### Filter heatmap_mat to include only significant genes with asterisks 
  rows_with_asterisks <- rownames(heatmap_mat)[grepl("\\*", rownames(heatmap_mat))] <- rownames(heatmap_mat)[grepl("\\*", rownames(heatmap_mat))]
  heatmap_mat <- heatmap_mat[rownames(heatmap_mat) %in% rows_with_asterisks, ]
  
  # Save the significant heatmap matrix to a CSV file
  #write.csv(heatmap_mat, file = paste0("../Figures/", gsub("\\.|\\.csv$", "", file), "_sig_heatmap_matrix.csv"))
  
  annotation_col <- data.frame(treatment = c(rep("S", 15), rep("P1", 15), rep("P2", 15)),
                               clutch_ID = c("L6" ,"L2" ,"L3", "L6", "L2", "L3", "L2", "L2", "L2", "L6", "L3", "L3" ,"L3", "L6", 
                                             "L3", "L3", "L6", "L2", "L2", "L3", "L6", "L6", "L6", "L3", "L3", "L6", "L2", "L6", 
                                             "L3", "L6", "L2", "L3", "L6", "L3","L6", "L2", "L3", "L3", "L6", "L2", "L2", "L6", 
                                             "L2", "L2", "L2"))
  
  ann_colors <- list(treatment =  c(P1= "magenta3", P2 ="purple3", S="grey50"),
                     clutch_ID = c(L2 = "brown", L3 = "brown1", L6="burlywood1"))
  
  ha <- HeatmapAnnotation(
    treatment = annotation_col$treatment,
    clutch_ID = annotation_col$clutch_ID,
    col = ann_colors)
  
  # Define custom colors based on the final_size_ratio values
  custom_colors <- ifelse(sorted_size_ratio$final_size_ratio < 0.8, "black", "red")
  
  ha1 <- HeatmapAnnotation(
    ratio = anno_barplot(sorted_size_ratio$final_size_ratio, gp = gpar(fill = custom_colors)))
  
  # Define square cell size (you can tweak this)
  cell_size <- 10  # in points, since PDF units are in inches * 72
  
  # Dynamically set PDF width and height
  #pdf_width <- ncol(heatmap_mat) * cell_size / 72  # convert to inches
  #df_height <- nrow(heatmap_mat) * cell_size / 72
  
  pdf(paste0("../Figures/Heatmap_sigDEGs_", gsub("percula\\.|\\.txt", "", file), ".pdf"), width = pdf_width, height = pdf_height)
  print(ComplexHeatmap::pheatmap(heatmap_mat, 
                                 name="matrix",
                                 cluster_cols= TRUE,
                                 cluster_rows = TRUE,
                                 scale="row", 
                                 color=col0,
                                 top_annotation = ha,
                                 bottom_annotation = ha1,
                                 #annotation_col = annotation_col,
                                 #annotation_colors = ann_colors,
                                 column_split = annotation_col$treatment, 
                                 #column_reorder = list(annotation_col$clutch_ID),
                                 #show_rownames= TRUE, 
                                 row_labels = rownames(heatmap_mat),
                                 show_colnames=T, 
                                 border_color="NA", 
                                 annotation_names_col=F,
                                 cellwidth = 8,    # Adjust the width of the entire heatmap
                                 cellheight = 12)  
  )
  
  dev.off()
}



rldpvals <- read.table(file = "../Data/rldpvals_significance_data2025-03-29.txt", header = TRUE, sep = "\t", check.names = FALSE)

for (file in dataFiles) {
  ### Read data
  data <- read.table(file.path("../Data/gene_lists", file), header = TRUE, sep = "\t")
  
  ### Filter vsd_matrix
  geneCounts <- vsdMatrix[rownames(vsdMatrix) %in% data$target_id, ]
  
  ### Set Gene.Name as row names
  if (nrow(geneCounts) > 0 && "Gene.Name" %in% colnames(data)) {
    rownames(geneCounts) <- data$Gene.Name[match(rownames(geneCounts), data$target_id)]
  }
  
  # Ensure unique row names by appending a suffix if duplicates exist
  if (anyDuplicated(rownames(geneCounts)) > 0) {
    rownames(geneCounts) <- make.unique(rownames(geneCounts))
  }
  
  ### Set gene_id as row names
  #if (nrow(geneCounts) > 0 && "gene_id" %in% colnames(data)) {
  #  rownames(geneCounts) <- data$gene_id[match(rownames(geneCounts), data$target_id)]
  #}
  
  ### Getting size-ratio data
  growth_data <- read.csv("../Data/Vizer_StrategicGrowth_phenotypicdata_SNPcorrected_20230906.csv", header = T)
  size_ratio <-  growth_data[, c("replicate_ID", "fish_type2", "final_size_ratio")]
  size_ratio
  size_ratio$sampleID <- paste(size_ratio$replicate_ID, size_ratio$fish_type2, sep = "")
  rownames(size_ratio) <- size_ratio$sampleID
  
  
  ### Generate heatmap
  select <- order(rowMeans(geneCounts))
  my_sample_cols <- data.frame(treatment = c(rep("S", 15), rep("P1", 15), rep("P2", 15)))
  row.names(my_sample_cols) <- colnames(geneCounts)
  my_Palette <- c("magenta3", "purple3", "grey50")
  my_colour <- list(treatment = c(P1 = my_Palette[1], P2 = my_Palette[2], S = my_Palette[3]))
  col0 <- colorRampPalette(rev(c("chocolate1", "#FEE090", "grey10", "cyan3", "cyan")))(100)
  
  heatmap_mat <- geneCounts[select, ]
  
  # Get the column names of heatmap_mat
  heatmap_colnames <- colnames(heatmap_mat)
  
  # Reorder rows in size_ratio to match the order of heatmap_mat
  sorted_size_ratio <- size_ratio[heatmap_colnames, "final_size_ratio", drop = FALSE]
  
  
  ## Read significance data
  significance_data <- read.table("../Data/rldpvals_significance_data2025-03-29.txt", header = TRUE, row.names = 1)
  colnames(significance_data) <- gsub("^X", "", colnames(significance_data))
  
  ### Filter significance_data
  significant_genes <- significance_data[rownames(significance_data) %in% data$target_id, ]
  
  ### Keep only rows which are significant at pval 0.1 for either P1vsP2 or P1vsS
  filtered_genes <- significant_genes %>%
    filter(padj.P1vsP2 < 0.1 | padj.P1vsS < 0.1)
  
  ### Set gene_id as row names
  if (nrow(filtered_genes) > 0 && "gene_id" %in% colnames(data)) {
    # Match row names of filtered_genes with gene_id column in data
    matched_row_names <- data$gene_id[match(rownames(filtered_genes), data$target_id)]
    
    # Ensure uniqueness of matched_row_names
    unique_row_names <- make.unique(matched_row_names, sep = "_")
    
    # Assign the unique row names to filtered_genes
    rownames(filtered_genes) <- unique_row_names
  }
  
  # Assuming filtered_genes contains the row names that passed the significant threshold
  significant_rows <- rownames(heatmap_mat) %in% rownames(filtered_genes)
  
  # Loop to add '*' for significance between P1vsP2
  for (gene in rownames(heatmap_mat)) {
    if (gene %in% rownames(filtered_genes)) {
      if (!is.na(filtered_genes[gene, "padj.P1vsP2"]) && filtered_genes[gene, "padj.P1vsP2"] < 0.1) {
        # Add double asterisks
        index <- which(rownames(heatmap_mat) == gene)
        rownames(heatmap_mat)[index] <- paste0(gene, "*")
      }
    }
  }
  
  for (gene in rownames(heatmap_mat)) {
    if (gene %in% rownames(filtered_genes)) {
      if (!is.na(filtered_genes[gene, "padj.P1vsS"]) && filtered_genes[gene, "padj.P1vsS"] < 0.1) {
        # Add double asterisks
        index <- which(rownames(heatmap_mat) == gene)
        rownames(heatmap_mat)[index] <- paste0(gene, "**")
      }
    }
  }
  
  ##**** 2 NEW ROWS below****** ##
  ### Filter heatmap_mat to include only significant genes with asterisks 
 # rows_with_asterisks <- rownames(heatmap_mat)[grepl("\\*", rownames(heatmap_mat))] <- rownames(heatmap_mat)[grepl("\\*", rownames(heatmap_mat))]
 # heatmap_mat <- heatmap_mat[rownames(heatmap_mat) %in% rows_with_asterisks, ]
  
  rows_with_asterisks <- rownames(heatmap_mat)[grepl("\\*", rownames(heatmap_mat))]
  heatmap_mat <- heatmap_mat[rows_with_asterisks, ]
  
  #ensure row labeles length matches the number of rows
  #row_labels <- rownames(heatmap_mat)
  
  annotation_col <- data.frame(treatment = c(rep("S", 15), rep("P1", 15), rep("P2", 15)),
                               clutch_ID = c("L6" ,"L2" ,"L3", "L6", "L2", "L3", "L2", "L2", "L2", "L6", "L3", "L3" ,"L3", "L6", 
                                             "L3", "L3", "L6", "L2", "L2", "L3", "L6", "L6", "L6", "L3", "L3", "L6", "L2", "L6", 
                                             "L3", "L6", "L2", "L3", "L6", "L3","L6", "L2", "L3", "L3", "L6", "L2", "L2", "L6", 
                                             "L2", "L2", "L2"))
  
  ann_colors <- list(treatment =  c(P1= "magenta3", P2 ="purple3", S="grey50"),
                     clutch_ID = c(L2 = "brown", L3 = "brown1", L6="burlywood1"))
  
  ha <- HeatmapAnnotation(
    treatment = annotation_col$treatment,
    clutch_ID = annotation_col$clutch_ID,
    col = ann_colors)
  
  # Define custom colors based on the final_size_ratio values
  custom_colors <- ifelse(sorted_size_ratio$final_size_ratio < 0.8, "black", "red")
  
  ha1 <- HeatmapAnnotation(
    ratio = anno_barplot(sorted_size_ratio$final_size_ratio, gp = gpar(fill = custom_colors)))
  
  # Check for NA or Inf values in heatmap_mat
  if (any(is.na(heatmap_mat)) || any(is.infinite(heatmap_mat))) {
    cat("Skipping heatmap generation for this file due to NA or Inf values in heatmap_mat.\n")
  } else {
  pdf(paste0("../Figures/Heatmap_sigDEGs_", gsub("percula\\.|\\.txt", "", file), ".pdf"), width = 8, height = 6)
  print(ComplexHeatmap::pheatmap(heatmap_mat, 
                                 name="matrix",
                                 cluster_cols= TRUE,
                                 cluster_rows = TRUE,
                                 scale="row", 
                                 color=col0,
                                 top_annotation = ha,
                                 bottom_annotation = ha1,
                                 #annotation_col = annotation_col,
                                 #annotation_colors = ann_colors,
                                 column_split = annotation_col$treatment, 
                                 #column_reorder = list(annotation_col$clutch_ID),
                                 #show_rownames= TRUE, 
                                 row_labels = rownames(heatmap_mat),
                                 show_colnames=T, 
                                 border_color="NA", 
                                 annotation_names_col=F,
                                 cellwidth = 8,    # Adjust the width of the entire heatmap
                                 cellheight = 12)  
                                 )
  
  dev.off()
}}




#####################################################
################# heatmap of manuscript #############
#####################################################

library(ComplexHeatmap)

#first read in original files
appetite_origin <- read.table("../Data/gene_lists/percula.appetite.regulation.txt", header = TRUE)
krebs_origin <- read.table("../Data/gene_lists/percula.metabolism.tca.cycle.txt", header = TRUE)
glycolysis_origin <- read.table("./Data/gene_lists/percula.metabolism.glycolysis.txt", header = TRUE)


significance_data <- read.table("../Data/rldpvals_significance_data2025-03-29.txt", header = TRUE, row.names = 1)
colnames(significance_data) <- gsub("^X", "", colnames(significance_data))

### Filter based on significance_data
significant_genes_appetite <- significance_data[rownames(significance_data) %in% appetite_origin$target_id, ]
# Assign Gene.Name as row names and ensure uniqueness
rownames(significant_genes_appetite) <- appetite_origin$gene_id[match(rownames(significant_genes_appetite), appetite_origin$target_id)]
#write.csv(significant_genes_appetite, file = "../Data/all_present_genes_appetite.csv", row.names = TRUE, quote = FALSE)

### Filter based on significance_data
significant_genes_krebs <- significance_data[rownames(significance_data) %in% krebs_origin$target_id, ]
# Assign Gene.Name as row names and ensure uniqueness
rownames(significant_genes_krebs) <- krebs_origin$gene_id[match(rownames(significant_genes_krebs), krebs_origin$target_id)]
#write.csv(significant_genes_krebs, file = "../Data/all_present_genes_krebs.csv", row.names = TRUE, quote = FALSE)

### Filter based on significance_data
significant_genes_glycolysis <- significance_data[rownames(significance_data) %in% glycolysis_origin$target_id, ]
# Assign Gene.Name as row names and ensure uniqueness
rownames(significant_genes_glycolysis) <- glycolysis_origin$gene_id[match(rownames(significant_genes_glycolysis), glycolysis_origin$target_id)]
#write.csv(significant_genes_glycolysis, file = "../Data/all_present_genes_glycolysis.csv", row.names = TRUE, quote = FALSE)


### Filter based on significance_data

#Keep only genes which are significant at pval 0.01 for either P1vsP2 or P1vsS
filtered_genes_appetite <- significant_genes_appetite %>%
  filter(pval.P1vsP2 < 0.01 | pval.P1vsS < 0.01)

filtered_genes_krebs <- significant_genes_krebs %>%
  filter(pval.P1vsP2 < 0.01 | pval.P1vsS < 0.01)

filtered_genes_glycolysis <- significant_genes_glycolysis %>%
  filter(pval.P1vsP2 < 0.01 | pval.P1vsS < 0.01)

#Combine the filtered data
merged_data <- rbind(filtered_genes_appetite, filtered_genes_krebs, filtered_genes_glycolysis)
merged_data <- as.matrix(merged_data)
head(merged_data)

# Assume merged_data is a data frame with gene names as rownames
# and columns: pval.P1vsP2, padj.P1vsP2, pval.P1vsS, padj.P1vsS
row_annotations <- lapply(1:nrow(merged_data), function(i) {
  row <- merged_data[i, ]
  pval_p1p2 <- as.numeric(row[["pval.P1vsP2"]])
  pval_p1s  <- as.numeric(row[["pval.P1vsS"]])
  padj_p1p2 <- as.numeric(row[["padj.P1vsP2"]])
  padj_p1s  <- as.numeric(row[["padj.P1vsS"]])
  
  stars <- ""
  if (!is.na(pval_p1p2) && pval_p1p2 < 0.01 && !is.na(pval_p1s) && pval_p1s < 0.01) {
    stars <- "***"
  } else if (!is.na(pval_p1p2) && pval_p1p2 < 0.01) {
    stars <- "*"
  } else if (!is.na(pval_p1s) && pval_p1s < 0.01) {
    stars <- "**"
  }
  
  bold <- (!is.na(padj_p1p2) && padj_p1p2 < 0.1) || (!is.na(padj_p1s) && padj_p1s < 0.1)
  
  list(label = paste0(rownames(merged_data)[i], stars), bold = bold)
})

# Extract new row labels and bold info
row_labels <- sapply(row_annotations, function(x) x$label)
row_fontface <- ifelse(sapply(row_annotations, function(x) x$bold), "bold", "plain")
table(row_labels = row_labels, row_fontface = row_fontface)



annotation_col <- data.frame(social_position = c(rep("S", 15), rep("P1", 15), rep("P2", 15)),
                             clutch_ID = c("L6" ,"L2" ,"L3", "L6", "L2", "L3", "L2", "L2", "L2", "L6", "L3", "L3" ,"L3", "L6", 
                                           "L3", "L3", "L6", "L2", "L2", "L3", "L6", "L6", "L6", "L3", "L3", "L6", "L2", "L6", 
                                           "L3", "L6", "L2", "L3", "L6", "L3","L6", "L2", "L3", "L3", "L6", "L2", "L2", "L6", 
                                           "L2", "L2", "L2"))

ann_colors <- list(social_position =  c(P1= "magenta3", P2 ="purple3", S="grey50"),
                   clutch_ID = c(L2 = "brown", L3 = "brown1", L6="burlywood1"))

ha <- HeatmapAnnotation(
  social_position = annotation_col$social_position,
  clutch_ID = annotation_col$clutch_ID,
  col = ann_colors,
  annotation_name_gp = gpar(fontface = "bold"))

col0 <- colorRampPalette(rev(c("chocolate1", "#FEE090", "grey10", "cyan3", "cyan")))(100)


cell_size <- 14  # in points, since PDF units are in inches * 72

# Dynamically set PDF width and height
pdf_width <- ncol(merged_data) * cell_size / 72  # convert to inches
pdf_height <- nrow(merged_data) * cell_size / 72

pdf(("../Figures/Heatmap_combined_AppetiteKrebsGlycolysis_20250406.pdf"), width = 9, height = 6)

fontsize = 13

ht_map <- ComplexHeatmap::pheatmap(as.matrix(merged_data[ , 1:(ncol(merged_data)-6)]),
                         name="matrix",
                         cluster_cols= TRUE,
                         cluster_rows = TRUE,
                         scale="row", 
                         color=col0,
                         top_annotation = ha,
                         #bottom_annotation = ha1,
                         #annotation_col = annotation_col,
                         #annotation_colors = ann_colors,
                         column_split = annotation_col$social_position, # Split columns by social_position
                         #column_reorder = list(annotation_col$clutch_ID),
                         #show_rownames= TRUE, 
                         row_labels = row_labels,  # Use annotated row names (with *, **, ***)
                        # rownames_gp = row_fontface,  # Bold where needed
                         show_colnames=T, 
                         border_color="NA",
                         cellwidth = cell_size,
                         cellheight = cell_size,
                         annotation_names_col=F,
                         fontsize = fontsize,
                         fontsize_row = fontsize,
                         fontsize_col = fontsize)

draw(ht_map, merge_legend = TRUE)

dev.off()

##########further exploration###################
##### Other metabolism and digestion genes######
################################################

lipase_origin <- read.table("../Data/gene_lists/percula.digestion.digestive.enzyme.lipase.txt", header = TRUE)
protease_origin <- read.table("../Data/gene_lists/percula.digestion.digestive.enzyme.protease.txt", header = TRUE)
b_oxidation_origin <- read.table("../Data/gene_lists/percula.metabolism.beta.oxidation.txt", header = TRUE)
cholesterol_origin <- read.table("../Data/gene_lists/percula.metabolism.cholesterol.biosynthesis.txt", header = TRUE)
fatty_acid_origin <- read.table("../Data/gene_lists/percula.metabolism.fatty.acid.synthesis.txt", header = TRUE)

### Filter based on significance_data
significant_genes_lipase <- significance_data[rownames(significance_data) %in% lipase_origin$target_id, ]
# Assign Gene.Name as row names and ensure uniqueness
rownames(significant_genes_lipase) <- lipase_origin$gene_id[match(rownames(significant_genes_lipase), lipase_origin$target_id)]

### Filter based on significance_data
significant_genes_protease <- significance_data[rownames(significance_data) %in% protease_origin$target_id, ]
# Make the gene IDs unique before assigning as row names
unique_gene_ids <- make.names(protease_origin$gene_id[match(rownames(significant_genes_protease), 
                                                            protease_origin$target_id)], unique = TRUE)
# Assign the unique names as row names
rownames(significant_genes_protease) <- unique_gene_ids

### Filter based on significance_data
significant_genes_b_oxidation <- significance_data[rownames(significance_data) %in% b_oxidation_origin$target_id, ]
# Assign Gene.Name as row names and ensure uniqueness
rownames(significant_genes_b_oxidation) <- b_oxidation_origin$gene_id[match(rownames(significant_genes_b_oxidation), b_oxidation_origin$target_id)]

### Filter based on significance_data
significant_genes_cholesterol <- significance_data[rownames(significance_data) %in% cholesterol_origin$target_id, ]
# Assign Gene.Name as row names and ensure uniqueness
rownames(significant_genes_cholesterol) <- cholesterol_origin$gene_id[match(rownames(significant_genes_cholesterol), cholesterol_origin$target_id)]

### Filter based on significance_data
significant_genes_fatty_acid <- significance_data[rownames(significance_data) %in% fatty_acid_origin$target_id, ]
# Assign Gene.Name as row names and ensure uniqueness
rownames(significant_genes_fatty_acid) <- fatty_acid_origin$gene_id[match(rownames(significant_genes_fatty_acid), fatty_acid_origin$target_id)]

filtered_genes_lipase <- significant_genes_lipase %>%
  filter(pval.P1vsP2 < 0.01 | pval.P1vsS < 0.01)

filtered_genes_protease <- significant_genes_protease %>%
  filter(pval.P1vsP2 < 0.01 | pval.P1vsS < 0.01)

filtered_genes_b_oxidation <- significant_genes_b_oxidation %>%
  filter(pval.P1vsP2 < 0.01 | pval.P1vsS < 0.01)

filtered_genes_cholesterol <- significant_genes_cholesterol %>%
  filter(pval.P1vsP2 < 0.01 | pval.P1vsS < 0.01)

filtered_genes_fatty_acid <- significant_genes_fatty_acid %>%
  filter(pval.P1vsP2 < 0.01 | pval.P1vsS < 0.01)

# Combine the filtered data
merged_metabolism_data <- rbind(filtered_genes_lipase, filtered_genes_b_oxidation, filtered_genes_protease,
                      filtered_genes_cholesterol, filtered_genes_fatty_acid)
merged_metabolism_data <- as.matrix(merged_metabolism_data)
head(merged_metabolism_data)

# Assume merged_data is a data frame with gene names as rownames
# and columns: pval.P1vsP2, padj.P1vsP2, pval.P1vsS, padj.P1vsS
row_annotations_metabol <- lapply(1:nrow(merged_metabolism_data), function(i) {
  row <- merged_metabolism_data[i, ]
  pval_p1p2 <- as.numeric(row[["pval.P1vsP2"]])
  pval_p1s  <- as.numeric(row[["pval.P1vsS"]])
  padj_p1p2 <- as.numeric(row[["padj.P1vsP2"]])
  padj_p1s  <- as.numeric(row[["padj.P1vsS"]])
  
  stars <- ""
  if (!is.na(pval_p1p2) && pval_p1p2 < 0.05 && !is.na(pval_p1s) && pval_p1s < 0.05) {
    stars <- "***"
  } else if (!is.na(pval_p1p2) && pval_p1p2 < 0.05) {
    stars <- "*"
  } else if (!is.na(pval_p1s) && pval_p1s < 0.05) {
    stars <- "**"
  }
  
  bold <- (!is.na(padj_p1p2) && padj_p1p2 < 0.05) || (!is.na(padj_p1s) && padj_p1s < 0.05)
  
  list(label = paste0(rownames(merged_metabolism_data)[i], stars), bold = bold)
})

# Extract new row labels and bold info
row_labels_metabol <- sapply(row_annotations_metabol, function(x) x$label)
row_fontface_metabol <- ifelse(sapply(row_annotations_metabol, function(x) x$bold), "bold", "plain")
table(row_labels_metabol = row_labels_metabol, row_fontface_metabol = row_fontface_metabol)

annotation_col_metabol <- data.frame(social_position = c(rep("S", 15), rep("P1", 15), rep("P2", 15)),
                                      clutch_ID = c("L6" ,"L2" ,"L3", "L6", "L2", "L3", "L2", "L2", "L2", "L6", "L3", "L3" ,"L3", "L6", 
                                                    "L3", "L3", "L6", "L2", "L2", "L3", "L6", "L6", "L6", "L3", "L3", "L6", "L2", "L6", 
                                                    "L3", "L6", "L2", "L3", "L6", "L3","L6", "L2", "L3", "L3", "L6", "L2", "L2", "L6",
                                                    "L2"," L2"," L2"))

ann_colors_metabol <- list(social_position =  c(P1= "magenta3", P2 ="purple3", S="grey50"),
                           clutch_ID = c(L2 = "brown", L3 = "brown1", L6="burlywood1"))

ha_metabol <- HeatmapAnnotation(
  social_position = annotation_col_metabol$social_position,
  clutch_ID = annotation_col_metabol$clutch_ID,
  col = ann_colors_metabol,
  annotation_name_gp = gpar(fontface = "bold"))

col0_metabol <- colorRampPalette(rev(c("chocolate1", "#FEE090", "grey10", "cyan3", "cyan")))(100)
cell_size_metabol <- 14  # in points, since PDF units are in inches * 72

# Dynamically set PDF width and height
pdf_width_metabol <- ncol(merged_metabolism_data) * cell_size_metabol / 72  # convert to inches
pdf_height_metabol <- nrow(merged_metabolism_data) * cell_size_metabol / 72

pdf(("../Figures/Heatmap_combined_metabolism_20250406.pdf"), width = 9, height = 6)
fontsize_metabol = 13
ht_map_metabol <- ComplexHeatmap::pheatmap(as.matrix(merged_metabolism_data[ , 1:(ncol(merged_metabolism_data)-6)]),
                         name="matrix",
                         cluster_cols= TRUE,
                         cluster_rows = TRUE,
                         scale="row", 
                         color=col0_metabol,
                         top_annotation = ha_metabol,
                         #bottom_annotation = ha1,
                         #annotation_col = annotation_col,
                         #annotation_colors = ann_colors,
                         column_split = annotation_col_metabol$social_position, # Split columns by social_position
                         #column_reorder = list(annotation_col$clutch_ID),
                         #show_rownames= TRUE, 
                         row_labels = row_labels_metabol,  # Use annotated row names (with *, **, ***)
                        # rownames_gp = row_fontface,  # Bold where needed
                         show_colnames=T, 
                         border_color="NA",
                         #cellwidth = cell_size_metabol,
                         #cellheight = cell_size_metabol,
                         annotation_names_col=F,
                         fontsize = fontsize_metabol,
                         fontsize_row = fontsize_metabol,
                         fontsize_col = fontsize_metabol)
draw(ht_map_metabol, merge_legend = TRUE)
dev.off()


