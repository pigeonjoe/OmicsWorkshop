# load libraries
library(ggplot2)
library(plotly)


# -----------------------------------------------------
# 1. load DESeq2 results
# -----------------------------------------------------
aluminium_indium_deseq_file <- 'aluminium_indium_vs_control_DESeq2_results.csv'
aluminium_indium_dres <- read.table(aluminium_indium_deseq_file, header = T, sep = ',',  row.names = 1, check.names = F)

aluminium_deseq_file <- 'aluminium_vs_control_DESeq2_results.csv'
aluminium_dres <- read.table(aluminium_deseq_file, header = T, sep = ',',  row.names = 1, check.names = F)

indium_deseq_file <- 'indium_vs_control_DESeq2_results.csv'
indium_dres <- read.table(indium_deseq_file, header = T, sep = ',',  row.names = 1, check.names = F)

# remove invalid genes
valid.genes <- !(is.na(aluminium_indium_dres$padj) | is.na(aluminium_dres$padj) | is.na(indium_dres$padj))
aluminium_indium_dres <- aluminium_indium_dres[valid.genes,]
aluminium_dres <- aluminium_dres[valid.genes,]
indium_dres <- indium_dres[valid.genes,]

# check the loaded values
dim(aluminium_indium_dres)
dim(aluminium_dres)
dim(indium_dres)


# -----------------------------------------------------
# 2. calculate observed and estimated log2FC for DEGs 
# -----------------------------------------------------
aluminium_indium_lfc <- aluminium_indium_dres$log2FoldChange
aluminium_lfc <- aluminium_dres$log2FoldChange
indium_lfc <- indium_dres$log2FoldChange

observed_lfc <- aluminium_indium_lfc
expected_lfc <- aluminium_lfc + indium_lfc


# -----------------------------------------------------
# 3. find the overlapping DEGs
# -----------------------------------------------------
padj.cutoff <- 0.05
degs <- (aluminium_indium_dres$padj <= padj.cutoff) & (aluminium_dres$padj <= padj.cutoff) & (indium_dres$padj <= padj.cutoff)

sum(degs)

observed_lfc <- observed_lfc[degs]
expected_lfc <- expected_lfc[degs]


# -----------------------------------------------------
# 4. plot additivity and fit linear model 
# -----------------------------------------------------
plot.data <- data.frame(Observed = observed_lfc, Expected = expected_lfc)

p <- ggplot(plot.data, aes(Observed,Expected)) +
     geom_point() + 
     geom_smooth(method = 'lm', formula = y~x) + 
     geom_abline(intercept = 0, slope = 2, color = 'red', linetype = 'dashed', linewidth = 1.5)
p
#ggplotly(p)
