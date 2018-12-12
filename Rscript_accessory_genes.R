# Set working directory:
setwd('/Users/cynthiarodriguez/Desktop/Final_project_informatics/test_again')

#Read data output we obtained from the previous script"
tablein <- read.table('gene_presence_absence_sorted.fix.txt', 
                      row.names = 1,header = T, sep = '\t', quote = NULL, na.strings=c("","NA"))

# Load gplots and heatmap.plus to make a heatmap.
library(gplots)
library(heatmap.plus)

#Read the data as a dataframe.
tablenum <- as.data.frame(lapply(tablein, function(x) as.integer(x!="0")))
tablenum[is.na(tablenum)] <- 0

#Read the data as a matrix to make heatmap
rownames(tablenum) <- rownames(tablein)
input <-as.matrix(t(tablenum))

#Use the Euclidean distance for clustering of strains in the heatmap.
distfunc <- function(x) dist(x, method="euclidean")

#Make heatmap.
heatmap.2(input, main = "Accessory genes from 129 bifidobacteria strains", distfun=distfunc, 
          trace = "none", margins = c(10,13), cexRow=0.2, cexCol=0.3)
