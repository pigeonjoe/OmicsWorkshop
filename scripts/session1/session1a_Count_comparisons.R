##### Read input file #####
df <- read.csv("DESeq2_results_all.csv", header = T, row.names = 1, check.names = F)
dim(df)
View(df)

AluIndVsCon <- as.data.frame(df[,1:3])
AluVsCon <- as.data.frame(df[,4:6])
IndVsCon <- as.data.frame(df[,7:9])


##### Count and subset of Differential Expressed Genes (DEGs) #####
##### Set cutoff to filter #####
pCutoff <- 0.05

### Number of significant DEG IDs by P-value and P-adjusted (Q) cutoffs ###
cat(paste0('significant p-values = ', sum(AluVsCon$pvalue <= pCutoff, na.rm = T), '\n'))
cat(paste0('significant p-adjust = ', sum(AluVsCon$padj <= pCutoff, na.rm = T), '\n'))

cat(paste0('significant p-values = ', sum(IndVsCon$pvalue <= pCutoff, na.rm = T), '\n'))
cat(paste0('significant p-adjust = ', sum(IndVsCon$padj <= pCutoff, na.rm = T), '\n'))

cat(paste0('significant p-values = ', sum(AluIndVsCon$pvalue <= pCutoff, na.rm = T), '\n'))
cat(paste0('significant p-adjust = ', sum(AluIndVsCon$padj <= pCutoff, na.rm = T), '\n'))


##### Subset significant DEG #####
AluVsConSig <- subset(AluVsCon, AluVsCon$padj <= pCutoff)
IndVsConSig <- subset(IndVsCon, IndVsCon$padj <= pCutoff)
AluIndVsConSig <- subset(AluIndVsCon, AluIndVsCon$padj <= pCutoff)


##### Subset up regulated genes #####
AluVsConUp <- subset(AluVsCon, AluVsCon$log2FoldChange > 0 & AluVsCon$padj <= pCutoff)
mean(AluVsConUp$log2FoldChange)

IndVsConUp <- subset(IndVsCon, IndVsCon$log2FoldChange > 0 & IndVsCon$padj <= pCutoff)
mean(IndVsConUp$log2FoldChange)

AluIndVsConUp <- subset(AluIndVsCon, AluIndVsCon$log2FoldChange > 0 & AluIndVsCon$padj <= pCutoff)
mean(AluIndVsConUp$log2FoldChange)


##### Subset down regulated genes #####
AluVsConDown <- subset(AluVsCon, AluVsCon$log2FoldChange < 0 & AluVsCon$padj <= pCutoff)
mean(AluVsConDown$log2FoldChange)

IndVsConDown <- subset(IndVsCon, IndVsCon$log2FoldChange < 0 & IndVsCon$padj <= pCutoff)
mean(IndVsConDown$log2FoldChange)

AluIndVsConDown <- subset(AluIndVsCon, AluIndVsCon$log2FoldChange < 0 & AluIndVsCon$padj <= pCutoff)
mean(AluIndVsConDown$log2FoldChange)

