#!/usr/bin/Rscript
# plot umi/cell or gene/cell

args <- commandArgs(); # print(args)
dir <- args[6]
Name <- args[7]

library(matrixStats)

print("Plot UMI/gene per cell")
File1 <- paste(Name, ".UMIcounts.csv", sep="")
File2 <- paste(Name, ".feature.count.percell.txt.summary", sep="")

Df_umi  <- read.csv(paste(dir, File1, sep="/"), header=T, sep = "\t");  # Df_umi[1:3,1:3]
Df_gene <- read.table(paste(dir, File2, sep="/"), header=T, sep = "\t"); # Df_gene[ ,1:3]

# Df_umi <-read.csv("0807_split_rna.hg19.dedup.UMIcounts.csv", header=T, sep = "\t"); Df_umi[1:3,1:3]
# Df_gene <-read.table("0807_split_rna.hg19.feature.count.percell.txt.summary", header=T, sep = "\t"); Df_gene[,1:3]

row.names(Df_umi) <- Df_umi$gene; # Df_umi[1:3,1:3]
Df2 <- Df_umi[ ,2:ncol(Df_umi)]; # Df2[1:3,1:3]

umiSum <- colSums(Df2); 
# plot(sort(log10(umiSum)))
geneSum <-colCounts(Df2>0)
# plot(sort(log10(geneSum)))
print("top10 UMI counts: ")
head(sort(umiSum, decreasing = T), 10)

# set cutoff of gene
Cutoff <- 0
Idx <- geneSum > Cutoff; # sum(Idx) # number of cell pass filter
# mean(geneSum[Idx])
Df3 <- as.data.frame(sort(geneSum, decreasing = T))
Df3$Count <- c(1: nrow(Df3)); colnames(Df3) <- c("Gene", "Count")

file3 <- paste(Name,'.detected.genes.pdf', sep="")
pdf(paste(dir,file3, sep="/"))
plot(Df3$Count,log10(Df3$Gene), xlab="Barcode rank", ylab = "log10 (Genes)", main = "Detected Genes per Cell", col="darkblue", pch=16)
dev.off()

# set cutoff of UMI
Cutoff2 <- 0
Idx2 <- umiSum > Cutoff2; # sum(Idx) # number of cell pass filter 
# mean(umiSum[Idx2])
Df4 <- as.data.frame(sort(umiSum, decreasing = T))
Df4$Count <- c(1: nrow(Df4)); colnames(Df4) <- c("Umi", "Count")

file4 <- paste(Name,'.detected.umi.pdf', sep="")
pdf(paste(dir,file4, sep="/"))
plot(Df4$Count,log10(Df4$Umi), xlab="Barcode rank", ylab = "log10 (UMIs)", main = "Detected UMIs per Cell", col="darkred", pch=16)
dev.off()

# plot reads and UMI
## reads per cell
Df_gene2 <- Df_gene[1, 2:ncol(Df_gene)]
Df_gene3 <- Df_gene2[, Df_gene2 > 0]; dim(Df_gene3)
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}
Newnames <- substrRight(colnames(Df_gene3), 17); # Newnames[2:length(Newnames)]
colnames(Df_gene3) <- Newnames 
Df_gene3 <- t(Df_gene3); colnames(Df_gene3) <- "Reads"; head(Df_gene3)
Df_gene4 <- Df_gene3[order(rownames(Df_gene3)), ]
Df_gene4 <- Df_gene4[Df_gene4 >0]
# head(Df_gene4); length(Df_gene4)

file5 <- paste(Name,'.detected.reads.pdf', sep="")
pdf(paste(dir,file5, sep="/"))
plot(log10(sort(Df_gene4, decreasing = T)), xlab="Barcode rank", ylab = "log10 (Reads)", main = "Detected Reads per Cell", col="orange", pch=16)
dev.off()

## UMIs per cell
# head(umiSum); length(umiSum)
Df_gene4 <- as.data.frame(Df_gene4)
umiSum <- as.data.frame(umiSum)
# combine umi and reads
Df <- cbind(Df_gene4, umiSum[match(rownames(Df_gene4), rownames(umiSum)), ]); colnames(Df) <- c("Reads","UMIs"); head(Df)

file6 <- paste(Name,'.umi.reads.pdf', sep="")
pdf(paste(dir,file6, sep="/"))
plot(Df, main = "Reads vs UMIs per cell", col="forestgreen", pch=16)
dev.off()

