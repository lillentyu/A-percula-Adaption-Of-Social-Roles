#WGCNA SCRIPT of **A. percula** gene expression data
#correlating phenotypic traits with gene expression data
#list of phenotypic traits: growth, hue (color), saturation (color), appetite, clutches (L2, L3, L6)


#library used
library(DESeq2)
library(WGCNA)
library(tidyverse)


# Load dds object from previous DESeq2 analyses
load("../Data/transformed_counts_mindepth10_run2.RData") #using transformed counts

# Variance stabilizing transformation, it's important to make sure blind = TRUE 
vst = vst(dds, blind = TRUE)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
# Let WGCNA use multi threads (do only once)
allowWGCNAThreads()  

# transpose data frame
input = t(as.data.frame(assay(vst)))

ncol(input) #number of columns: 20537 number of genes
nrow(input) #number of rows: 45 - matches the number of samples


# filter genes that are useful in WGCNA
#outlier dection and removal
gsg = goodSamplesGenes(input)
gsg$allOK #TRUE all genes have passed, no outlines

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(input)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(input)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  input = input[gsg$goodSamples, gsg$goodGenes]
}
gsg = goodSamplesGenes(input)
gsg$allOK #TRUE

# read in trait spreadsheet
traits = read.csv("../Data/expDesign.csv", sep = ",", row.names = 1)

# cluster samples
sampleTree = hclust(dist(input), method = "average")

# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(traits, signed = TRUE);

# Plot the sample dendrogram and the colors underneath.
#pdf("../Figures/sample_dendro_and_traits_heatmap.pdf")

plotDendroAndColors(sampleTree,
                    traitColors,
                    groupLabels = names(traits),
                    main = "Sample dendrogram and trait heatmap")
#dev.off()

# Choose a set of soft-thresholding powers - using default
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function, having networkType signed!
sft = pickSoftThreshold(input, powerVector = powers, verbose = 5, networkType = "signed")

#Power SFT.R.sq  slope truncated.R.sq  mean.k. median.k. max.k.
#1      1   0.0239  20.60          0.981 10300.00  1.03e+04  10500
#2      2   0.2690 -27.90          0.904  5320.00  5.31e+03   5710
#3      3   0.6160 -20.00          0.921  2850.00  2.82e+03   3310
#4      4   0.7670 -10.40          0.932  1570.00  1.54e+03   2030
#5      5   0.8820  -7.89          0.973   891.00  8.58e+02   1350
#6      6   0.9350  -6.22          0.991   520.00  4.90e+02    958
#7      7   0.9570  -5.09          0.997   312.00  2.85e+02    714
#8      8   0.9670  -4.30          0.998   192.00  1.69e+02    553
#9      9   0.9720  -3.72          0.997   121.00  1.02e+02    443
#10    10   0.9750  -3.29          0.998    78.80  6.27e+01    365
#11    12   0.9720  -2.74          0.995    35.90  2.47e+01    262
#12    14   0.9710  -2.38          0.996    18.10  1.02e+01    199
#13    16   0.9690  -2.13          0.996     9.96  4.43e+00    158
#14    18   0.9660  -1.96          0.995     5.95  2.00e+00    129
#15    20   0.9610  -1.84          0.993     3.80  9.42e-01    108

# Plot the soft threshold scale independence
pdf("../Figures/scale_independence.pdf")
cex1 = 0.9;
plot(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",
     type="n",
     main = paste("Scale independence"));

text(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,
     cex=cex1,
     col="red");
abline(h=0.90,col="red") # this line corresponds to using an R^2 cut-off of h
dev.off()

# Mean connectivity as a function of the soft-thresholding power

pdf("../Figures/Mean_connectivity.pdf")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

#Moving onto automatic module construction and detection

#I choose the power 6, which is the lowest power for which the scale-free topology fit
#index curve flattens out upon reaching a high value (in this case, roughly 0.90).

# automatic network construction and module detection, choosing a soft threshold power of 6;
#I have networkType and TOMType as signed
net = blockwiseModules(input,
                       power = 6,
                       networkType = "signed",
                       TOMType = "signed",
                       minModuleSize = 600,
                       numericLabels = FALSE,
                       mergeCutHeight = 0,
                       saveTOMFileBase = "percula_TOM_net",
                       verbose = 3)

# number of genes associated with each module
table(net$colors)#this gives me all module colors & number of genes after 0.25 mergeHeight
#black      blue     brown     green      grey      pink     red   turquoise    yellow 
#972      1758      1477      1006     10327       781       991      1847      1378 

table(net$unmergedColors) #this gives me all module colors & number of genes without 0.25 mergeHeight (unmerged)
#black      blue     brown     green      grey      pink       red turquoise    yellow 
#1006      1477       781       991     10327       972      1378      1847      1758 

#from this we can see that unmerged and 0.25 merged match in the number of modules and genes in each, note colors do not match up exactly
#I`ll continue to visualise the net$colors 

# Convert labels to colors for plotting
moduleColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]],
                    moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05)

#calculate eigengenes 
MEList = moduleEigengenes(input, colors = net$colors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(12,9)
plot(METree,
     main = "Clustering of module eigengenes",
     xlab = "",
     sub = "") 
abline(h=0, col="red")

table(moduleColors)

# Define numbers of genes and samples
nGenes = ncol(input);
nSamples = nrow(input);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(input, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
# correlations of genes with eigengenes
moduleTraitCor = cor(MEs, traits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(3,3,3,3));
# Display the correlation values within a heatmap plot
#let`s plot only significant correlations
ps600 <- signif(moduleTraitPvalue,1)
cors600  <- signif(moduleTraitCor,2)
textMatrix <-  paste(cors600, "\n(",
                     ps600, ")", sep = "")
textMatrix[ps600>0.05]="-"
dim(textMatrix) <- dim(moduleTraitCor)

#heatmap with only showing the significant correlations
sizeGrWindow(12, 9)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(traits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = TRUE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


#let`s look at the clustering of modules and see which height should I merge
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");

sizeGrWindow(12,9)
plot(METree,
     main = "Clustering of module eigengenes",
     xlab = "",
     sub = "") 
abline(h=0, col="red") #we could carry out a merging at 0.46 height

## Saves unmerged data that you need to create gene scatter plots, heat maps, and GO MWU files from all significant module/trait pairings
save(input, moduleTraitPvalue, moduleTraitCor, moduleColors, MEs, traits, file = "unmergedWGCNAmods600.RData")

#Given clutches - L2, L3, L6 are not an interest of trait I`ll remove these from the heatmap


head(traits)

subset_traits <- traits[, c(-4,-5,-6)]
head(subset_traits)

moduleTraitCor = cor(MEs, subset_traits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(3,3,3,3));
# Display the correlation values within a heatmap plot
#let`s plot only significant correlations
modLabels <- sub("ME","",names(MEs))
ps_sub <- signif(moduleTraitPvalue,1)
cors_sub  <- signif(moduleTraitCor,2)
textMatrix <-  paste(cors_sub, "\n(",
                     ps_sub, ")", sep = "")
textMatrix[ps_sub>0.05]="-"
dim(textMatrix) <- dim(moduleTraitCor)

#heatmap with only showing the significant correlations without clutches displayed
sizeGrWindow(12, 9)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(subset_traits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = TRUE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

#saving all data I need as input for creating heatmaps
save(input, moduleTraitCor,moduleTraitPvalue, moduleColors, MEs, subset_traits, vst, file = "wgcnaMods_fin.RData")

load("./wgcnaMods_fin.RData")

##using some of JK-s code to create WGCNA figure:
mod_df <- data.frame(module = row.names(moduleTraitCor), moduleTraitCor) %>% 
  gather("trait", "corr", -module, factor_key = TRUE) %>% 
  mutate(module = factor(module, levels = rownames(moduleTraitCor))) %>% 
  mutate(module = fct_rev(module))

mod_pval_df <- data.frame(module = row.names(ps_sub), ps_sub) %>% 
  gather("trait", "pval", -module, factor_key = TRUE) %>% 
  mutate(module = factor(module, levels = rownames(ps_sub))) %>% 
  mutate(module = fct_rev(module))

mod_df <- mod_df %>% 
  add_column(lab = mod_pval_df$pval) %>% 
  mutate(lab = ifelse(lab > 0.05, "-", round(corr,2)))

library(ggplot2)
#### Heatmap of modules from WGCNA
theme_bove <- function() {
  # Define the custom theme here
  theme_minimal() + 
    theme(
      panel.grid.major = element_line(color = "grey80"),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold")
    )
}

str(mod_df$corr)
heatmap <- ggplot(data = mod_df, aes(x = trait, y = module, fill = corr)) +
  theme_bove() +
  geom_tile() +
  geom_text(aes(label = lab), color = "#2f2f2f", size = 3) +
  #scale_fill_gradient2(colors = c("dodgerblue2", "skyblue2", "lightcoral", "firebrick1"))
  scale_fill_gradient2(low = "dodgerblue2",  mid = "#FFFFFF",high = "firebrick1") +
  theme(legend.position = "right", legend.title = element_text(size = 12), plot.margin = margin(10, 10, 10, 150)) + 
  theme(axis.ticks = element_blank(), axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, size = 12)) +
  # theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, size = 12)) +
  scale_x_discrete(labels = colnames(moduleTraitCor)) +
  guides(fill = guide_colourbar(barwidth = 0.7, barheight = 10, frame.colour = "black", ticks.colour = "black", title = expression(paste("Pearson", " R"^{2})))) +
  coord_cartesian(expand = FALSE)
heatmap

#### barplot of genes within each module
mct <- table(moduleColors)
modLabels_plot <- as.data.frame(mct[rev(modLabels)])
modLabels_plot2 <- as.character(modLabels_plot$moduleColors)
view(modLabels_plot)
library(cowplot)

barplot <- ggplot(data = modLabels_plot, aes(y = Freq, x = moduleColors, fill = moduleColors)) +
  geom_col(colour = "black", size = 0.1) +
  scale_fill_manual(values = modLabels_plot2) +
  theme_cowplot() +
  coord_flip(expand = FALSE) +
  #scale_y_manual(limits=c(8000, 0), breaks = c(1000, 2000, 3000,4000,5000,6000,7000,8000)) +
  geom_text(aes(label = paste0(rev(modLabels), " (", Freq, ")")), hjust = -0.03, color="black", size = 3) +
  theme(axis.ticks.y = element_blank(), axis.title.y = element_blank(), legend.position = "none", axis.line.y = element_blank(),
        axis.line.x = element_line(size = 0.2), axis.text.y = element_blank(), axis.text.x = element_text(angle = 45,
        hjust = 1, size = 10), axis.title.x = element_text(size = 10)) +
  ylab("# genes\n per module")
  #scale_y_continuous(limits = c(0,4300), breaks = c(200, 400, 600, 800, 1000, 1200, 1400,1800,2500,3000,3500,4300))
barplot

##removing the green module from both plots - green module had over 10K genes; hard to see patterns with such a big module
mod_df <- mod_df %>% filter(module != "MEgreen")
modLabels_plot <- modLabels_plot %>% filter(moduleColors != "MEgreen")


heatmap_nogreen <- ggplot(data = mod_df, aes(x = trait, y = module, fill = corr)) +
  theme_bove() +
  geom_tile() +
  geom_text(aes(label = lab), color = "#2f2f2f", size = 3) +
  #scale_fill_gradient2(colors = c("dodgerblue2", "skyblue2", "lightcoral", "firebrick1"))
  scale_fill_gradient2(low = "dodgerblue2",  mid = "#FFFFFF",high = "firebrick1") +
  theme(legend.position = "right", legend.title = element_text(size = 12), plot.margin = margin(10, 10, 10, 10)) + 
  theme(axis.ticks = element_blank(), axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, size = 12)) +
  # theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, size = 12)) +
  scale_x_discrete(labels = colnames(moduleTraitCor)) +
  guides(fill = guide_colourbar(barwidth = 0.7, barheight = 10, frame.colour = "black", ticks.colour = "black", title = expression(paste("Pearson", " R"^{2})))) +
  coord_cartesian(expand = FALSE)
heatmap_nogreen

#### barplot of genes within each module - no green
moduleColors_nogreen <- moduleColors[moduleColors != "green"]
mct_nogreen <- table(moduleColors_nogreen) 
modLabels_nogreen <- modLabels[modLabels != "green"]
modLabels_plot_nogreen <- as.data.frame(mct_nogreen[rev(modLabels_nogreen)])
modLabels_plot2_nogreen <- as.character(modLabels_plot_nogreen$moduleColors)


barplot_nogreen <- ggplot(data = modLabels_plot_nogreen, aes(y = Freq, x = moduleColors_nogreen, fill = moduleColors_nogreen)) +
  geom_col(colour = "black", size = 0.1) +
  scale_fill_manual(values = modLabels_plot2_nogreen) +
  theme_cowplot() +
  coord_flip(expand = FALSE) +
  #scale_y_manual(limits=c(8000, 0), breaks = c(1000, 2000, 3000,4000,5000,6000,7000,8000)) +
  geom_text(aes(label = paste0(rev(modLabels_nogreen), " (", Freq, ")")), hjust = -0.03, color="black", size = 3) +
  theme(axis.ticks.y = element_blank(), axis.title.y = element_blank(), legend.position = "none", axis.line.y = element_blank(),
        axis.line.x = element_line(size = 0.2), axis.text.y = element_blank(), axis.text.x = element_text(angle = 45,
        hjust = 1, size = 10), axis.title.x = element_text(size = 10)) +
  ylab("# genes\n per module")+ scale_y_continuous(limits = c(0,2200), breaks = c(200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200))

barplot_nogreen

library(cowplot)
# Combine heatmap and barplot side by side
combined_plot <- plot_grid(heatmap_nogreen, barplot_nogreen, ncol = 2, align = "h", rel_widths = c(3, 1.3))
combined_plot
ggsave(file = "../Figures/combined_WGCNA_barplot_no_colour.pdf", height = 5.8, width = 14, units = "in", dpi = 1200)

#wgcna mod/trait plots and files
## A couple things to do before running functions
nGenes <- ncol(input) # extract number of genes; (input = transposed VSD transformed count data)
nSamples <- nrow(input) # extract number of samples; (input = transposed VSD transformed count data)
modNames <- substring(names(MEs), 3) # string of all module names 

geneModuleMembership <- as.data.frame(signedKME(input,MEs));
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) <- paste("MM", modNames, sep="");
names(MMPvalue) <- paste("p.MM", modNames, sep="")
datt=input
traits=subset_traits
vsdWG=vst
write.csv(geneModuleMembership, file = "geneModuleMembership_KME.csv")


# Get unique module names
unique_modules <- unique(moduleColors)
unique_modules
# Create a list of genes by module
genes_by_module <- lapply(unique_modules, function(module) {
  colnames(input)[moduleColors == module]
})

# Name the list elements by their module names
names(genes_by_module) <- unique_modules
for (module in names(genes_by_module)) {
  write.table(
    genes_by_module[[module]],
    file = paste0("Module_", module, "_genes.txt"),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}

#below using steps outlined here
##https://alexslemonade.github.io/refinebio-examples/04-advanced-topics/network-analysis_rnaseq_01_wgcna.html#4_Identifying_co-expression_gene_modules_with_WGCNA_-_RNA-seq
#Which modules have biggest differences across treatment (fish category) groups?

#Run linear model on each module. Limma wants our tests to be per row, so we also need to transpose so the eigengenes are rows
module_eigengenes <- MEs
head(module_eigengenes)
all.equal(rownames(subset_traits), rownames(module_eigengenes))

str(subset_traits)
# Add a new column 'fish_type' based on the row names
subset_traits$fish_type <- ifelse(grepl("S$", rownames(subset_traits)), "S",
                                  ifelse(grepl("P1$", rownames(subset_traits)), "P1",
                                         ifelse(grepl("P2$", rownames(subset_traits)), "P2", NA)))
head(subset_traits)
#make a column that match rownames called Sample_ID
subset_traits$Sample_ID <- rownames(subset_traits)

# Create the design matrix from the `fish-category` variable

des_mat <- model.matrix(~ subset_traits$fish_type)
des_mat
# lmFit() needs a transposed version of the matrix
fit <- limma::lmFit(t(module_eigengenes), design = des_mat)

# Apply empirical Bayes to smooth standard errors
fit <- limma::eBayes(fit)
# Apply multiple testing correction and obtain stats
stats_df <- limma::topTable(fit, number = ncol(module_eigengenes)) %>%
  tibble::rownames_to_column("module")
head(stats_df)

#       module              subset_traits.fish_typeP2 subset_traits.fish_typeS       AveExpr         F      P.Value  adj.P.Val
#1      MEblue                0.15547891               0.19558889 -3.569675e-17 7.4842143 0.0005618845 0.00505696
#2       MEred                0.11234301               0.16050924  9.685539e-19 4.7571899 0.0085897133 0.03865371
#3     MEblack               -0.05671816              -0.12521107 -8.962738e-18 2.7565559 0.0635101254 0.14344433
#4     MEbrown               -0.11425278              -0.10169562  2.397773e-17 2.7527385 0.0637530361 0.14344433
#5    MEyellow                0.05759747               0.10215505 -2.571246e-17 1.8393954 0.1589134789 0.28604426
#6 MEturquoise               -0.03009072              -0.06411815  2.021917e-17 0.7216229 0.4859629572 0.72894444
#significant modules are blue and red for P1 vs P2 vs S 
#almost significant are black and brown - will explore these 2

#let`s explore blue and red mdoules further
#First we need to set up the module eigengene for this module with the sample metadata labels we need.
module_blue_df <- module_eigengenes %>%
  rownames_to_column("accession_code") %>%
  # Here we are performing an inner join with a subset of metadata
  inner_join(subset_traits %>%
               select(Sample_ID, fish_type,),
             by = c("accession_code" = "Sample_ID"))
  
# Plotting the significant modules in treatment (p.adj<0.05)
#significant eigengens
#an eigengene is like an "average expression pattern" for a group of co-expressed genes
library(ggforce)
#MEblue
#shapes = c("P1" = 21, "gen3_D" = 22, "gen4_D"=24)
colvec<- c("P1" = "magenta2","P2" = "purple", "S" = "grey")
# #plot significant cyan,turquoise,pink,black
blue<-ggplot(module_blue_df, aes(x = fish_type, y = MEblue, fill = fish_type, group = fish_type)) +geom_boxplot(outlier.shape = NA) +
  ggforce::geom_sina(maxwidth = 0.3, size = 4, alpha = 0.5) +
  theme_classic() +
  scale_fill_manual(values = colvec) +
  labs(x = "social position", y = "Blue module")+
  guides(fill = guide_legend(override.aes = list(alpha = 1, size = 2))) +  # Adjust legend point size and opacity
  theme(legend.position = "none",                # Position legend (e.g., "bottom", "none")
        legend.title = element_text(size = 12),  # Legend title size
        legend.text = element_text(size = 10))   # Legend text size

blue

#red module
red<-ggplot(module_blue_df, aes(x = fish_type, y = MEred, fill = fish_type, group = fish_type)) +geom_boxplot(outlier.shape = NA) +
  ggforce::geom_sina(maxwidth = 0.3, size = 4, alpha = 0.5) +
  theme_classic() +
  scale_fill_manual(values = colvec) +
  labs(x = "social position", y = "Red module")+
  guides(fill = guide_legend(override.aes = list(alpha = 1, size = 2))) +  # Adjust legend point size and opacity
  theme(legend.position = "none",                # Position legend (e.g., "bottom", "none")
        legend.title = element_text(size = 12),  # Legend title size
        legend.text = element_text(size = 10))   # Legend text size

red

#brown module
brown<-ggplot(module_blue_df, aes(x = fish_type, y = MEbrown, fill = fish_type, group = fish_type)) +geom_boxplot(outlier.shape = NA) +
  ggforce::geom_sina(maxwidth = 0.3, size = 4, alpha = 0.5) +
  theme_classic() +
  scale_fill_manual(values = colvec) +
  labs(x = "social position", y = "Brown module")+
  guides(fill = guide_legend(override.aes = list(alpha = 1, size = 3))) +  # Adjust legend point size and opacity
  theme(legend.position = "none",                # Position legend (e.g., "bottom", "none"),
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 12),  # Legend title size
        legend.text = element_text(size = 10))   # Legend text size

brown

#black module
black<-ggplot(module_blue_df, aes(x = fish_type, y = MEblack, fill = fish_type, group = fish_type)) +geom_boxplot(outlier.shape = NA) +
  ggforce::geom_sina(maxwidth = 0.3, size = 4, alpha = 0.5) +
  theme_classic() +
  scale_fill_manual(values = colvec) +
  labs(x = "social position", y = "Black module")+
  guides(fill = guide_legend(override.aes = list(alpha = 1, size = 2))) +  # Adjust legend point size and opacity
  theme(legend.position = "none",                # Position legend (e.g., "bottom", "none")
        legend.title = element_text(size = 12),  # Legend title size
        legend.text = element_text(size = 10))   # Legend text size

black

#combine the two plots
library(cowplot)
combined_plot <- plot_grid(brown, blue, red, black, ncol = 4, align = "h", rel_widths = c(1, 1))
combined_plot
ggsave(file = "../Figures/WGCNA_modules_boxplots.pdf", height = 5, width = 12, units = "in", dpi = 1200)

#knowing the identified modules of interest, I will proceed to run GO-MWU analysis on them
#to do this I`ll extract module information for GO-MWU
#load("wgcnaMods_fin.RData")

#Looking at the heatmap there are modules of interest:
#I`ll extract modules to run GO-MWU analysis on them
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(input, MEs, use = "p"));
nSamples = nrow(input)
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");

##### generating  GO MWU outputs
modules = sub("MM", "", colnames(geneModuleMembership))
genes = rownames(geneModuleMembership)
inputFiles = c()

for (m in modules){
  moduleGenes = genes[moduleColors==m]
  godf = data.frame(gene=genes,
                    inMod=if_else(genes %in% moduleGenes,
                                  1,
                                  0))
  outname = paste(c("./GO_MWU/", m, "_moduleInput.csv"), collapse="")
  print(paste("module = ", m))
  print(paste("total genes =", length(moduleGenes)))
  write.csv(godf, file=outname, row.names=F, quote=F)
  inputFiles = append(inputFiles, paste(m,"_moduleInput.csv", sep=""))
}
#this writes out all the modules we have; we`ll only use the 4 modules we`re interested for downstream GO-MWU
save(inputFiles, file='./GO_MWU/moduleInputFiles.Rdata')

#[1] "module =  brown"
#[1] "total genes = 1477"
#[1] "module =  green"
#[1] "total genes = 10327"
#[1] "module =  yellow"
#[1] "total genes = 1006"
#[1] "module =  blue"
#[1] "total genes = 1758"
#[1] "module =  red"
#[1] "total genes = 781"
#[1] "module =  magenta"
#[1] "total genes = 1378"
#[1] "module =  turquoise"
#[1] "total genes = 972"
#[1] "module =  black"
#[1] "total genes = 991"
#[1] "module =  pink"
#[1] "total genes = 1847"

# proceed to GO-NWU analysis R script
# I run GO-MWU with Fisher`s exact test (presence/absence) for given gene for given GO term


############################## GO-MWU analysis RESULTS ########################################
#I have already run the GO-MWU analysis and saved the results in the GO_MWU_20250306 folder
#pulling out GO terms and associated genes of interest from WGCNA - GO-MWU modules
#goal is to create heatmaps by combining list of interesting annotated genes with RDL-pvals 

#pulling GO terms and combining with rldpvals and writing csv files to then do heatmaps 

library(tidyverse)
library(dplyr)
#reading in the rldpval file:
rldpvals <- read.csv(file="../Data/percula_RLDandPVALS.csv", row.names=1)
head(rldpvals)
#removing X from sample names
colnames(rldpvals) <- gsub("^X", "", colnames(rldpvals))

#geneID2geneNames
geneID2Name=read.table("./geneID2geneNAME.tab",sep="\t")
head(geneID2Name)
geneID2Name <- geneID2Name[-1,]
colnames(geneID2Name) <- c("GeneID", "GeneName")
#replace empty gg GeneName rows with NA
geneID2Name["GeneName"][geneID2Name["GeneName"] == ''] <- NA


#genes2Go terms
gene2go <- read.table("./DgeneID2GO.tab", header = T)
head(gene2go)

rldpvals_df <- as.data.frame(rldpvals)
colnames(rldpvals_df) 

#changing so the geneIDs (row names) become the first column
rldpvals_df <- cbind(rownames(rldpvals_df),rldpvals_df)
rownames(rldpvals_df) <- NULL
colnames(rldpvals_df) <- c(names(rldpvals_df)) #to not write all the column names
colnames(rldpvals_df)[1] <- "GeneID" 
names(rldpvals_df)

head(rldpvals_df) #rld in the right format

####################GO TERMS##########################
# for all GO terms, only GO terms which had a p-adjusted value < 0.05 were considered significant and were pulled out#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~brown module~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# pulling out GO terms of interest from brown module BP and MF - mostly GO terms associated with growth and development
brownmodBP <- read.table("./GO_MWU/MWU_BP_brown_moduleInput.csv")

#adjust so first row is the header and remove first row
colnames(brownmodBP) <- brownmodBP[1,]
brownmodBP <- brownmodBP[-1,]
head(brownmodBP)
##GO terms of interest: 
# striated muscle cell development - GO:0048741;GO:0014904;GO:0048747;GO:0055002
# muscle cell development - GO:0055001
# skeletal system development - GO:0001501
# dorsal/ventral pattern formation - GO:0009953
# ossification - 	GO:0001503
# cartilage development - GO:0051216
# cartilage morphogenesis - 	GO:0060536
# animal organ morphogenesis - GO:0009887
# tissue morphogenesis - GO:0002009;GO:0048729
# tissue development -GO:0009888
# skeletal muscle tissue development - 	GO:0007519
# muscle tissue development - 	 GO:0014706;GO:0060537
# muscle structure development - GO:0061061
# muscle organ development - 	GO:0007517
# muscle tissue morphogenesis - 	GO:0055008;GO:0060415

#taking all above GOterms and extracting genes associated with them
brown_BP_growthdev <- gene2go %>%
  filter(grepl("GO:0048741|GO:0014904|GO:0048747|GO:0055002|GO:0055001|GO:0001501|GO:0009953|GO:0001503|GO:0051216|GO:0060536|GO:0009887|GO:0002009|GO:0048729|GO:0009888|GO:0007519|GO:0014706|GO:0060537|GO:0061061|GO:0007517|GO:0055008|GO:0060415", GO_Terms)) %>%
  left_join(rldpvals_df)

head(brown_BP_growthdev)

#need to match GO term genes to a table with both ensembl gene ID and gene names
brown_BP_growthdev <- brown_BP_growthdev %>% 
  left_join(geneID2Name, by=c("GeneID"))
colnames(brown_BP_growthdev)
nrow(brown_BP_growthdev)  #full set of genes = 227

#subset by only keeping annotated genes:
brown_BP_growthdev_anno<-subset(brown_BP_growthdev, (!is.na(brown_BP_growthdev$GeneName)))
head(brown_BP_growthdev_anno)
dim(brown_BP_growthdev_anno) #annotated set of genes: 216
#write.csv(brown_BP_growthdev_anno, file.path("./pulledGOterms/brown_BP_growthdevelopment_anno.csv"), row.names = F)

brownmodMF <- read.table("./GO_MWU/MWU_MF_brown_moduleInput.csv")
#adjust so first row is the header and remove first row
colnames(brownmodMF) <- brownmodMF[1,]
brownmodMF <- brownmodMF[-1,]
head(brownmodMF)
##GO terms of interest: 
# growth factor activity GO:0008083
# growth factor activity GO:0008083
# actin binding GO:0003779
# actin filament binding	GO:0051015

#taking all above GOterms and extracting genes associated with them
brown_MF_growth <- gene2go %>%
  filter(grepl("GO:0008083|GO:0003779|GO:0051015", GO_Terms)) %>%
  left_join(rldpvals_df)

head(brown_MF_growth)
#need to match GO term genes to a table with both ensembl gene ID and gene names
brown_MF_growth <- brown_MF_growth %>% 
  left_join(geneID2Name, by=c("GeneID"))
colnames(brown_MF_growth)
nrow(brown_MF_growth)  #full set of genes = 413

#subset by only keeping annotated genes:
brown_MF_growth_anno<-subset(brown_MF_growth, (!is.na(brown_MF_growth$GeneName)))
head(brown_MF_growth_anno)
dim(brown_MF_growth_anno) #annotated set of genes: 275
#write.csv(brown_MF_growth_anno, file.path("./GO_MWU/pulledGOterms/brown_BP_growthdevelopment_anno.csv"), row.names = F)

### combine the two datasets for heatmap

brown_mod_growth_BP <- read.csv("./GO_MWU/pulledGOterms/brown_BP_growthdevelopment_anno.csv", header = T)
brown_mod_growth_MF <- read.csv("./GO_MWU/pulledGOterms/brown_MF_growth_anno.csv", header = T)
head(brown_mod_growth_BP)
dim(brown_mod_growth_BP) #annotated set of genes: 216
dim(brown_mod_growth_MF) #annotated set of genes: 275

# Perform an anti-join to remove matching rows
merged_brown_growth <- anti_join(brown_mod_growth_BP, brown_mod_growth_MF, by = c("GeneID", "GeneName")) %>%
  bind_rows(anti_join(brown_mod_growth_BP, brown_mod_growth_MF, by = c("GeneID", "GeneName")))
dim(merged_brown_growth) #annotated set of genes: 410

# Save the merged dataset
#write.csv(merged_brown_growth, "./GO_MWU/pulledGOterms/mergedbrown_growth_dataset.csv", row.names = FALSE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~blue module~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# pulling out GO terms of interest from blue module BP and MF - mostly GO terms associated with stimulus and behavior
bluemodBP <- read.table("./GO_MWU/MWU_BP_blue_moduleInput.csv")

#adjust so first row is the header and remove first row
colnames(bluemodBP) <- bluemodBP[1,]
bluemodBP <- bluemodBP[-1,]
head(bluemodBP)

#GO terms of interest
# detection of visible light - GO:0009584
# detection of external stimulus - GO:0009583;GO:0009581;GO:0009582
# detection of stimulus - 	GO:0051606
# response to radiation - GO:0009416;GO:0009314
# phototransduction - GO:0007602
# response to abiotic stimulus - GO:0009628
# nervous system process - GO:0007600;GO:0050877
# system process - GO:0003008
# sensory perception of light stimulus  - GO:0007601;GO:0050953
# multicellular organismal process  - GO:0032501
# cone photoresponse recovery - GO:0036368

#taking all above GOterms and extracting genes associated with them
blue_BP_stimulus <- gene2go %>%
  filter(grepl("GO:0009584|GO:0009583|GO:0009581|GO:0009582|GO:0051606|GO:0009416|GO:0009314|GO:0007602|GO:0009628|GO:0007600|GO:0050877|GO:0003008|GO:0007601|GO:0050953|GO:0032501|GO:0036368", GO_Terms)) %>%
  left_join(rldpvals_df)
#need to match GO term genes to a table with both ensembl gene ID and gene names
blue_BP_stimulus <- blue_BP_stimulus %>% 
  left_join(geneID2Name, by=c("GeneID"))
colnames(blue_BP_stimulus)
nrow(blue_BP_stimulus)  #full set of genes = 99

#subset by only keeping annotated genes:
blue_BP_stimulus_anno<-subset(blue_BP_stimulus, (!is.na(blue_BP_stimulus$GeneName)))
head(blue_BP_stimulus_anno)
dim(blue_BP_stimulus_anno) #annotated set of genes: 85
#write.csv(blue_BP_stimulus_anno, file.path("./GO_MWU/pulledGOterms/blue_BP_stimulus_anno.csv"), row.names = F)

#go terms associated with behavior
# locomotory behavior - 	 GO:0007626
# behavior - GO:0007610
# visual behavior - GO:0007634;GO:0007632

#taking all above GOterms and extracting genes associated with them
blue_BP_behavior <- gene2go %>%
  filter(grepl("GO:0007626|GO:0007610|GO:0007634|GO:0007632", GO_Terms)) %>%
  left_join(rldpvals_df)
#need to match GO term genes to a table with both ensembl gene ID and gene names
blue_BP_behavior <- blue_BP_behavior %>% 
  left_join(geneID2Name, by=c("GeneID"))
colnames(blue_BP_behavior)
nrow(blue_BP_behavior)  #full set of genes = 37

#subset by only keeping annotated genes:
blue_BP_behavior_anno<-subset(blue_BP_behavior, (!is.na(blue_BP_behavior$GeneName)))
head(blue_BP_behavior_anno)
dim(blue_BP_behavior_anno) #annotated set of genes: 35
#write.csv(blue_BP_behavior_anno, file.path("./GO_MWU/pulledGOterms/blue_BP_behavior_anno.csv"), row.names = F)


bluemodMF <- read.table("./GO_MWU/MWU_MF_blue_moduleInput.csv")
#adjust so first row is the header and remove first row
colnames(bluemodMF) <- bluemodMF[1,]
bluemodMF <- bluemodMF[-1,]
head(bluemodMF)

##GO terms of interest:
# neurotransmitter receptor - GO:0030594;GO:0022824;GO:0022835
# GABA receptor - GO:0016917;GO:0004890
# hormone - GO:0005179

#taking all above GOterms and extracting genes associated with them
blue_MF_receptor <- gene2go %>%
  filter(grepl("GO:0030594|GO:0022824|GO:0022835|GO:0016917|GO:0004890|GO:0005179", GO_Terms)) %>%
  left_join(rldpvals_df)
#need to match GO term genes to a table with both ensembl gene ID and gene names
blue_MF_receptor <- blue_MF_receptor %>% 
  left_join(geneID2Name, by=c("GeneID"))
colnames(blue_MF_receptor)
nrow(blue_MF_receptor)  #full set of genes = 119

#subset by only keeping annotated genes:
blue_MF_receptor_anno<-subset(blue_MF_receptor, (!is.na(blue_MF_receptor$GeneName)))
head(blue_MF_receptor_anno)
dim(blue_MF_receptor_anno) #annotated set of genes: 88
#write.csv(blue_MF_receptor_anno, file.path("./GO_MWU/pulledGOterms/blue_MF_receptor_anno.csv"), row.names = F)



#########################################################################################################
############## running heatmaps for the above GO terms and genes of interest ############################
library(ComplexHeatmap)
library(pheatmap)
library(ggplot2)


dataFiles <- c(
  "blue_BP_behavior_anno.csv",
  "blue_BP_stimulus_anno.csv",
  "blue_MF_receptor_anno.csv",
  "brown_BP_growthdevelopment_anno.csv",
  "brown_MF_growth_anno.csv"
)


for (file in dataFiles) {
  ### Read data
  data <- read.csv(file.path("./GO_MWU/pulledGOterms", file), header = TRUE)
  
  ### Remove X from column names
  colnames(data) <- gsub("^X", "", colnames(data))
  
  ### Filter data to have vsdMatrix information only columns: 1, 3:47, 54
  geneCounts <- data[, c(1, 3:47, 54)]
  
  ### Remove completely identical rows
  geneCounts <- geneCounts[!duplicated(geneCounts), ]
  
  ### Make gene names unique with _1, _2, etc.
  unique_gene_names <- make.unique(geneCounts$GeneName, sep = "_")
  
  ### Set GeneName as row names, remove GeneID and GeneName from dataset
  rownames(geneCounts) <- unique_gene_names
  geneCounts <- geneCounts[, -c(1, 47)]
  
  ### Remove any rows that are only NA values
  geneCounts <- geneCounts[rowSums(is.na(geneCounts)) != ncol(geneCounts), ]
  
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
  

  ## Get significance data
  #filter - adjust depending on what cut off you want
  p.val=0.05# FDR cutoff - this will be a p-value cutoff of the unadjusted p-values to filter only significant genes
  
  significance_data <- data[data$pval.P1vsP2<=p.val & !is.na(data$pval.P1vsP2) & data$pval.P1vsS<=p.val & !is.na(data$pval.P1vsS),]
 
  for (gene in significance_data$GeneName) {
    if (gene %in% rownames(heatmap_mat)) {
      pval_p1vsp2 <- significance_data$pval.P1vsP2[significance_data$GeneName == gene]
      pval_p1vss <- significance_data$pval.P1vsS[significance_data$GeneName == gene]
      
      label <- gene
      if (!is.na(pval_p1vsp2) && pval_p1vsp2 <= 0.05 && !is.na(pval_p1vss) && pval_p1vss <= 0.05) {
        label <- paste0(gene, "***")
      } else if (!is.na(pval_p1vsp2) && pval_p1vsp2 <= 0.05) {
        label <- paste0(gene, "*")
      } else if (!is.na(pval_p1vss) && pval_p1vss <= 0.05) {
        label <- paste0(gene, "**")
      }
      
      rownames(heatmap_mat)[rownames(heatmap_mat) == gene] <- label
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
  
  ### Only keep significant rows (rows with *)
  heatmap_mat_sig <- heatmap_mat[grep("\\*", rownames(heatmap_mat)), ]
  
  ### Make heatmap_mat_sig a matirix
  heatmap_mat_sig <- as.matrix(heatmap_mat_sig)
  
  # Save the significant heatmap matrix to a CSV file
  write.csv(heatmap_mat_sig, file = paste0("./GO_MWU/pulledGOterms/", gsub("\\.|\\.anno.csv$", "", file), "_sig_heatmap_matrix.csv"))
  
  pdf(paste0("./GO_MWU/pulledGOterms/Heatmap_", gsub("\\.|\\.anno.csv$", "", file), ".pdf"), width = 8, height = 6)
  print(ComplexHeatmap::pheatmap(heatmap_mat_sig, 
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
                                 row_labels = rownames(heatmap_mat_sig),
                                 show_colnames=T, 
                                 border_color="NA", 
                                 annotation_names_col=F)
  )
  
  dev.off()
}

########## manuscript heatmap for brown module growth and development GO terms ##########
########################################################################################

brown_original_BP <- read.csv("./pulledGOterms/brown_BP_growthdevelopment_anno.csv", header = T) 
brown_original_MF <- read.csv("./pulledGOterms/brown_MF_growth_anno.csv", header = T) 
colnames(brown_original_BP) <- gsub("^X", "", colnames(brown_original_BP))
colnames(brown_original_MF) <- gsub("^X", "", colnames(brown_original_MF))

# Filter datasets
#Keep only genes which are significant at pval 0.01 for either P1vsP2 or P1vsS
filtered_brown_original_BP <- brown_original_BP %>%
  filter(pval.P1vsP2 < 0.05 & pval.P1vsS < 0.05)
nrow(filtered_brown_original_BP) # 19


filtered_brown_original_MF <- brown_original_MF %>%
  filter(pval.P1vsP2 < 0.05 & pval.P1vsS < 0.05)
nrow(filtered_brown_original_MF) # 19

#make GeneName rownames
rownames(filtered_brown_original_BP) <- filtered_brown_original_BP$GeneName
#make rownames unique
filtered_brown_original_MF$GeneName <- make.unique(filtered_brown_original_MF$GeneName, sep = "_")
rownames(filtered_brown_original_MF) <- filtered_brown_original_MF$GeneName

#Combine the filtered data
merged_data <- as.data.frame(rbind(filtered_brown_original_BP, filtered_brown_original_MF))

#remove GeneID and GO terms columns (columns 1 and 2)
merged_data <- merged_data[, -c(1, 2)]

head(merged_data)
str(merged_data)
#change data type to numeric



# Assume merged_data is a data frame with gene names as rownames
# and columns: pval.P1vsP2, padj.P1vsP2, pval.P1vsS, padj.P1vsS
row_annotations <- lapply(1:nrow(merged_data), function(i) {
  row <- merged_data[i, ]
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

fontsize = 13

ht_map_growth <- ComplexHeatmap::pheatmap(as.matrix(merged_data[ , 1:(ncol(merged_data)-7)]),
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

draw(ht_map_growth, merge_legend = TRUE)



## reading in the 2 brown module heatmap matrixes and plotting them together
brown_BP_growthdev <- read.csv("./pulledGOterms/brown_BP_growthdevelopment_annocsv_sig_heatmap_matrix.csv", header = TRUE)
brown_MF_growth <- read.csv("./pulledGOterms/brown_MF_growth_annocsv_sig_heatmap_matrix.csv", header = TRUE)
head(brown_BP_growthdev)

#make gene names row names 
rownames(brown_BP_growthdev) <- brown_BP_growthdev$X
#delete that first column then
brown_BP_growthdev <- brown_BP_growthdev %>% 
  select(-X)

rownames(brown_MF_growth) <- brown_MF_growth$X
#delete that first column then
brown_MF_growth <- brown_MF_growth %>% 
  select(-X)

colnames(brown_BP_growthdev) <- gsub("^X", "", colnames(brown_BP_growthdev))
colnames(brown_MF_growth) <- gsub("^X", "", colnames(brown_MF_growth))
head(brown_BP_growthdev)
head(brown_MF_growth)

#merge these 2 datasets
merged_brown_growth <-  rbind(brown_BP_growthdev, brown_MF_growth)

#remove duplicate rows
merged_brown_growth <- merged_brown_growth[!duplicated(merged_brown_growth), ]

#make this a matrix
merged_brown_growth <- as.matrix(merged_brown_growth)

### Getting size-ratio data - no need for size ratio data in this heatmap
#growth_data <- read.csv("../Data/Vizer_StrategicGrowth_phenotypicdata_SNPcorrected_20230906.csv", header = T)
#size_ratio <-  growth_data[, c("replicate_ID", "fish_type2", "final_size_ratio")]
#size_ratio
#size_ratio$sampleID <- paste(size_ratio$replicate_ID, size_ratio$fish_type2, sep = "")
#rownames(size_ratio) <- size_ratio$sampleID
annotation_col <- data.frame(social_position = c(rep("S", 15), rep("P1", 15), rep("P2", 15)),
                             clutch_ID = c("L6" ,"L2" ,"L3", "L6", "L2", "L3", "L2", "L2", "L2", "L6", "L3", "L3" ,"L3", "L6", 
                                           "L3", "L3", "L6", "L2", "L2", "L3", "L6", "L6", "L6", "L3", "L3", "L6", "L2", "L6", 
                                           "L3", "L6", "L2", "L3", "L6", "L3","L6", "L2", "L3", "L3", "L6", "L2", "L2", "L6", 
                                           "L2", "L2", "L2"))

ann_colors <- list(social_position =  c(P1= "magenta3", P2 ="purple3", S="grey50"),
                   clutch_ID = c(L2 = "brown", L3 = "brown1", L6="burlywood1"))
col0 <- colorRampPalette(rev(c("chocolate1", "#FEE090", "grey10", "cyan3", "cyan")))(100)

ha <- HeatmapAnnotation(
  social_position = annotation_col$social_position,
  clutch_ID = annotation_col$clutch_ID,
  col = ann_colors,
  annotation_name_gp = gpar(fontface = "bold"))
  
# Define custom colors based on the final_size_ratio values
#custom_colors <- ifelse(sorted_size_ratio$final_size_ratio < 0.8, "black", "red")

#ha1 <- HeatmapAnnotation(
#  ratio = anno_barplot(sorted_size_ratio$final_size_ratio, gp = gpar(fill = custom_colors)))
fontsize <- 13

cell_size <- 14  # in points, since PDF units are in inches * 72

# Dynamically set PDF width and height
pdf_width <- ncol(merged_data) * cell_size / 72  # convert to inches
pdf_height <- nrow(merged_data) * cell_size / 72

ht <- ComplexHeatmap::pheatmap(merged_brown_growth, 
                          name="matrix",
                          cluster_cols= TRUE,
                          cluster_rows = TRUE,
                          scale="row", 
                          color=col0,
                          top_annotation = ha,
                          #bottom_annotation = ha1,
                          #annotation_col = annotation_col,
                          #annotation_colors = ann_colors,
                          column_split = annotation_col$social_position, 
                          #column_reorder = list(annotation_col$clutch_ID),
                          #show_rownames= TRUE, 
                          row_labels = rownames(merged_brown_growth),
                          show_colnames=T, 
                          border_color="NA", 
                          annotation_names_col=F,
                          cellwidth = cell_size,
                          cellheight = cell_size,
                          fontsize = fontsize,
                          fontsize_row = fontsize,
                          fontsize_col = fontsize)
  
# draw heatmap, merge legends, and increase font size for legends
draw(ht,
     merge_legend = TRUE)

########################################################################################
######################## making bubble plots of pulled GO terms ########################

#Create bubble plots for GO terms which are significant in the brown module

# pulling out GO terms of interest from brown module BP and MF - mostly GO terms associated with growth and development
brownmodBP <- read.table("./GO_MWU/MWU_BP_brown_moduleInput.csv", header = T)

##GO terms of interest: 
# striated muscle cell development - GO:0048741;GO:0014904;GO:0048747;GO:0055002
# muscle cell development - GO:0055001
# skeletal system development - GO:0001501
# dorsal/ventral pattern formation - GO:0009953
# ossification - 	GO:0001503
# cartilage development - GO:0051216
# cartilage morphogenesis - 	GO:0060536
# animal organ morphogenesis - GO:0009887
# tissue morphogenesis - GO:0002009;GO:0048729
# tissue development -GO:0009888
# skeletal muscle tissue development - 	GO:0007519
# muscle tissue development - 	 GO:0014706;GO:0060537
# muscle structure development - GO:0061061
# muscle organ development - 	GO:0007517
# muscle tissue morphogenesis - 	GO:0055008;GO:0060415

# Store your GO terms as a pattern
go_terms_BP <- paste0(
  "GO:0048741|GO:0014904|GO:0048747|GO:0055002|GO:0055001|GO:0001501|GO:0009953|GO:0001503|GO:0051216|GO:0060536|GO:0009887|GO:0002009|GO:0048729|GO:0009888|GO:0007519|GO:0014706|GO:0060537|GO:0061061|GO:0007517|GO:0055008|GO:0060415")

# Filter rows that contain at least one of those GO terms
filtered_brownmodBP <- brownmodBP %>%
  filter(grepl(go_terms_BP, term))
head(filtered_brownmodBP)
nrow(filtered_brownmodBP) # 15

brownmodMF <- read.table("./GO_MWU/MWU_MF_brown_moduleInput.csv", header = T)
head(brownmodMF)

##GO terms of interest: 
# growth factor activity GO:0008083
# growth factor activity GO:0008083
# actin binding GO:0003779
# actin filament binding	GO:0051015

# Store your GO terms as a pattern
go_terms_MF <- paste0(
  "GO:0008083|GO:0003779|GO:0051015")

# Filter rows that contain at least one of those GO terms
filtered_brownmodMF <- brownmodMF %>%
  filter(grepl(go_terms_MF, term))
head(filtered_brownmodMF)
nrow(filtered_brownmodMF) # 3

# Add ontology labels before merging
filtered_brownmodBP$category <- "Biological Process"
filtered_brownmodMF$category <- "Molecular Function"

# Combine the filtered data
merged_brownmod <- rbind(filtered_brownmodBP, filtered_brownmodMF)
head(merged_brownmod)

# Prepare data: reorder GO terms by significance for cleaner y-axis
bubble_data <- merged_brownmod %>%
  arrange(p.adj) %>%
  mutate(name = factor(name, levels = rev(unique(name))))

# Plot
ggplot(bubble_data, aes(x = -log10(p.adj), y = name, size = nseqs, color = category)) +
  geom_point(alpha = 0.8) +
  scale_size(range = c(2, 15),
             breaks = c(10, 20, 50, 100, 200),
             name = "Gene Count") +
  scale_color_manual(values = c(
      "Biological Process" = "brown",
      "Molecular Function" = "lightblue")) +
  labs(
    x = expression(-log[10]~"(adjusted p-value)"),
    y = "GO Term",
    color = "Ontology"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 14),  # Increase legend title size
        legend.text = element_text(size = 14))    # Increase legend label size)

ggsave("../Figures/brown_module_GO_terms_bubble_plot.png", 
       width = 12, height = 8, dpi = 1200, bg="white")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
### Running above code to generate figures one by one - so I can combine them via cowplot

#blue modules
data <- read.csv(file.path("./pulledGOterms/blue_BP_behavior_anno.csv"), header = TRUE)
#go back to line: #819 and re-run all until heatmap
BP_behavior <- ComplexHeatmap::pheatmap(heatmap_mat_sig, 
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
                               row_labels = rownames(heatmap_mat_sig),
                               show_colnames=T, 
                               border_color="NA", 
                               annotation_names_col=F)

BP_behavior
nrow(BP_behavior) #5

data <- read.csv(file.path("./pulledGOterms/blue_BP_stimulus_anno.csv"), header = TRUE)
BP_stimulus <- ComplexHeatmap::pheatmap(heatmap_mat_sig, 
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
                               row_labels = rownames(heatmap_mat_sig),
                               show_colnames=T, 
                               border_color="NA", 
                               annotation_names_col=F)
BP_stimulus
nrow(BP_stimulus) #19

data <- read.csv(file.path("./pulledGOterms/blue_MF_receptor_anno.csv"), header = TRUE)
MF_receptor <- ComplexHeatmap::pheatmap(heatmap_mat_sig, 
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
                               row_labels = rownames(heatmap_mat_sig),
                               show_colnames=T, 
                               border_color="NA", 
                               annotation_names_col=F)
MF_receptor
nrow(MF_receptor) #8

data <- read.csv(file.path("./pulledGOterms/brown_BP_growthdevelopment_anno.csv"), header = TRUE)
BP_growth <- ComplexHeatmap::pheatmap(heatmap_mat_sig, 
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
                               row_labels = rownames(heatmap_mat_sig),
                               show_colnames=T, 
                               border_color="NA", 
                               annotation_names_col=F)
BP_growth
nrow(BP_growth) #19


data <- read.csv(file.path("./pulledGOterms/brown_MF_growth_anno.csv"), header = TRUE)
MF_growth <- ComplexHeatmap::pheatmap(heatmap_mat_sig, 
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
                               row_labels = rownames(heatmap_mat_sig),
                               show_colnames=T, 
                               border_color="NA", 
                               annotation_names_col=F)
MF_growth
nrow(MF_growth) #19

### Combine heatmaps of blue module and heatmaps of brown module
library(grid)

# Capture heatmaps as grobs
BP_behavior_grob <- grid.grabExpr(draw(BP_behavior))
BP_stimulus_grob <- grid.grabExpr(draw(BP_stimulus))
MF_receptor_grob <- grid.grabExpr(draw(MF_receptor))

# Combine the grobs in a grid layout
p1 <- plot_grid(BP_behavior_grob, BP_stimulus_grob, MF_receptor_grob, ncol = 3, nrow = 1)
p1

