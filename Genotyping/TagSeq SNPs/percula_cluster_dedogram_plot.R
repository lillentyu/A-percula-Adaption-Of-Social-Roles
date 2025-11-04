### This script contains code to generate dendogram from the ANGSD matrix data
# Code by Lili Vizer
# Date: 2023 August

#First we start by reading in the bams files.
#Change what files you are working with to plot the cluster dendrogram with and without genotype A

#bams=read.table("bams")[,1] # list of bam files for symbiotic/brown oculina individuals only
bams=read.table("percula_bams")[,1] # list of bam files for all percula individuals
head(bams)
goods=c(1:length(bams))


ma = as.matrix(read.table("./percula_transcriptom_angsd_out.ibsMat")) # identity by state matrix of all individuals

ma=ma[goods,goods]
dimnames(ma)=list(bams[goods],bams[goods])


# trim column and row names to make the cluster dendrogram nicer to read
colnames(ma) = c("1P1_L3", "1P2_L2", "1S_L6", 
                 "10P1_L6",	"10P2_L3",	"10S_L2",	
                 "11P1_L2",	"11P2_L6",	"11S_L3",
                 "13P1_L3", "13P2_L2", "13S_L6",
                 "14P1_L3",	"14P2_L6",	"14S_L2",
                 "15P1_L2", "15P2_L6", "15S_L3",
                 "2P1_L3",	"2P2_L6",	"2S_L2",	
                 "20P1_L6",	"20P2_L3",	"20S_L2",	
                 "21P1_L3", "21P2_L6", "21S_L2",
                 "23P1_L3",	"23P2_L2",	"23S_L6",	
                 "26P1_L6",	"26P2_L2",	"26S_L3",
                 "3P1_L6", "3P2_L2", "3S_L3",
                 "5P1_L6", "5P2_L2", "5S_L3",
                 "6P1_L3", "6P2_L2", "6S_L6",
                 "7P1_L2", "7P2_L6", "7S_L3")

row.names(ma) = c("1P1_L3", "1P2_L2", "1S_L6", 
                  "10P1_L6",	"10P2_L3",	"10S_L2",	
                  "11P1_L2",	"11P2_L6",	"11S_L3",
                  "13P1_L3", "13P2_L2", "13S_L6",
                  "14P1_L3",	"14P2_L6",	"14S_L2",
                  "15P1_L2", "15P2_L6", "15S_L3",
                  "2P1_L3",	"2P2_L6",	"2S_L2",	
                  "20P1_L6",	"20P2_L3",	"20S_L2",	
                  "21P1_L3", "21P2_L6", "21S_L2",
                  "23P1_L3",	"23P2_L2",	"23S_L6",	
                  "26P1_L6",	"26P2_L2",	"26S_L3",
                  "3P1_L6", "3P2_L2", "3S_L3",
                  "5P1_L6", "5P2_L2", "5S_L3",
                  "6P1_L3", "6P2_L2", "6S_L6",
                  "7P1_L2", "7P2_L6", "7S_L3")

# now plot dendrograms
# first with all samples
hc=hclust(as.dist(ma),"ave")
hc
hc.plot = plot(hc,cex=0.7) 

#corrected dendogram
ma_corr <- ma
# trim column and row names to make the cluster dendrogram nicer to read
colnames(ma_corr) = c("1P1_L3", "1P2_L2", "1S_L6", 
                 "10P1_L6",	"10P2_L3",	"10S_L2",	
                 "11P1_L2",	"11P2_L6",	"11S_L3",
                 "13P1_L2", "13P2_L3", "13S_L6",
                 "14P1_L3",	"14P2_L6",	"14S_L2",
                 "15P1_L6", "15P2_L2", "15S_L3",
                 "2P1_L6",	"2P2_L3",	"2S_L2",	
                 "20P1_L6",	"20P2_L3",	"20S_L2",	
                 "21P1_L3", "21P2_L6", "21S_L2",
                 "23P1_L3",	"23P2_L2",	"23S_L6",	
                 "26P1_L6",	"26P2_L2",	"26S_L3",
                 "3P1_L2", "3P2_L6", "3S_L3",
                 "5P1_L6", "5P2_L2", "5S_L3",
                 "6P1_L3", "6P2_L2", "6S_L6",
                 "7P1_L6", "7P2_L2", "7S_L3")

row.names(ma_corr) = c("1P1_L3", "1P2_L2", "1S_L6", 
                  "10P1_L6",	"10P2_L3",	"10S_L2",	
                  "11P1_L2",	"11P2_L6",	"11S_L3",
                  "13P1_L2", "13P2_L3", "13S_L6",
                  "14P1_L3",	"14P2_L6",	"14S_L2",
                  "15P1_L6", "15P2_L2", "15S_L3",
                  "2P1_L6",	"2P2_L3",	"2S_L2",	
                  "20P1_L6",	"20P2_L3",	"20S_L2",	
                  "21P1_L3", "21P2_L6", "21S_L2",
                  "23P1_L3",	"23P2_L2",	"23S_L6",	
                  "26P1_L6",	"26P2_L2",	"26S_L3",
                  "3P1_L2", "3P2_L6", "3S_L3",
                  "5P1_L6", "5P2_L2", "5S_L3",
                  "6P1_L3", "6P2_L2", "6S_L6",
                  "7P1_L6", "7P2_L2", "7S_L3")

# now plot dendrograms
# first with all samples
hc_corr=hclust(as.dist(ma_corr),"ave")
hc_corr
hc.plot = plot(hc_corr,cex=0.7) 


# combine the 2 plots
library(gridGraphics)
library(grid)
library(cowplot)
library(ggplot2)

# Function to capture base plot as a grob
grab_plot <- function(expr) {
  grid.newpage()
  grid.echo(expr)
  grid.grab()
}

# Capture the dendrogram plots
plot1 <- grab_plot(function() plot(hc, main = "", xlab = "Sample", sub = ""))
plot2 <- grab_plot(function() plot(hc_corr, main = "", xlab = "Sample", sub = ""))

# Combine using cowplot
combined_plot <- plot_grid(plot1, plot2, labels = c("A", "B"), ncol = 1, align = "v", rel_heights = c(1, 1))
print(combined_plot)

# Save the combined plot
ggsave("combined_dendrogram_plot.png", plot = combined_plot, width = 8, height = 8, dpi = 1200, bg = "white", units = "in", limitsize = FALSE)

