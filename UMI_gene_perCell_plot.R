#!/usr/bin/Rscript
# plot umi/cell or gene/cell

args <- commandArgs(); # print(args)
dir <- args[6]
Name <- args[7]

library(matrixStats)

print("Plot UMI/gene per cell")
File1 <- paste(Name, ".UMIcounts.csv", sep="")

Df_umi  <- read.csv(paste(dir, File1, sep="/"), header=T, sep = "\t");  # Df_umi[1:3,1:3]
# Df_umi <-read.csv("0807_split_rna.hg19.dedup.UMIcounts.csv", header=T, sep = "\t"); Df_umi[1:3,1:3]

row.names(Df_umi) <- Df_umi$gene; # Df_umi[1:3,1:3]
Df2 <- Df_umi[ ,2:ncol(Df_umi)]; # Df2[1:3,1:3]

umiSum <- colSums(Df2); 
# plot(sort(log10(umiSum)))
geneSum <- as.data.frame(colCounts(Df2>0))
rownames(geneSum) <- colnames(Df2)

# plot(sort(log10(geneSum)))
print("top10 UMI counts: ")
head(sort(umiSum, decreasing = T), 10)

# set cutoff of gene
Cutoff <- 0
Idx <- geneSum > Cutoff; # sum(Idx) # number of cell pass filter
# mean(geneSum[Idx])
colnames(geneSum) <- "Genes"
Df3 <- as.data.frame(sort(geneSum$Genes, decreasing = T))
Df3$Count <- c(1: nrow(Df3)); colnames(Df3) <- c("Gene", "Count")

file3 <- paste(Name,'.detected.genes.pdf', sep="")
pdf(paste(dir,file3, sep="/"))
plot(Df3$Count,log10(Df3$Gene), xlab="Barcode rank", ylab = "log10 (Genes)", main = "Detected Genes per Cell", col="darkblue", pch=16)
garbage <- dev.off()

# set cutoff of UMI
Cutoff2 <- 0
Idx2 <- umiSum > Cutoff2; # sum(Idx) # number of cell pass filter 
# mean(umiSum[Idx2])
Df4 <- as.data.frame(sort(umiSum, decreasing = T))
Df4$Count <- c(1: nrow(Df4)); colnames(Df4) <- c("Umi", "Count")

## UMIs per cell
# head(umiSum); length(umiSum)
umiSum <- as.data.frame(umiSum)
geneSum <- as.data.frame(geneSum)

# plot umi vs genes
Df <- cbind(geneSum$Genes, umiSum[match(rownames(geneSum), rownames(umiSum)), ]); colnames(Df) <- c("Genes","UMIs"); head(Df)
Df <- as.data.frame(Df)

file7 <- paste(Name,'.UMIvsGenes.pdf', sep="")
pdf(paste(dir,file7, sep="/"))
plot(Df$UMIs,Df$Genes, xlab="UMIs", ylab = "Detected Genes", main = "Detected Genes per Cell", col="darkblue", pch=16)
legend("topleft", c(paste("Number of Cells that detected > 0 genes:       ",sum(Df$Genes > 0),sep = ""), 
paste("Number of Cells that detected > 10 genes:     ",sum(Df$Genes > 10), sep = ""),
paste("Number of Cells that detected > 100 genes:   ",sum(Df$Genes > 100),sep = ""),
paste("Number of Cells that detected > 500 genes:   ",sum(Df$Genes > 500),sep = ""),
paste("Number of Cells that detected > 1000 genes: ",sum(Df$Genes > 1000), sep = "")), bty="n") 
garbage <- dev.off()

print(paste("Number of Cells that detected > 0 genes:    ",sum(Df$Genes > 0), sep = ""))
print(paste("Number of Cells that detected > 10 genes:   ",sum(Df$Genes > 10), sep = ""))
print(paste("Number of Cells that detected > 100 genes:  ",sum(Df$Genes > 100), sep = ""))
print(paste("Number of Cells that detected > 500 genes:  ",sum(Df$Genes > 500), sep = ""))
print(paste("Number of Cells that detected > 1000 genes: ",sum(Df$Genes > 1000), sep = ""))