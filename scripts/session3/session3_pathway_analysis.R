# read DESeq2 results
df <- read.csv("DESeq2_results_all.csv", header = T, row.names = 1, check.names = F)

AluIndVsCon <- as.data.frame(df[,1:3])
AluVsCon <- as.data.frame(df[,4:6])
IndVsCon <- as.data.frame(df[,7:9])

# number of significant DEG IDs by P-value and P-adjusted (Q) cutoff
pcutoff <- 0.05

AluVsConSig <- subset(AluVsCon, AluVsCon$padj <= pcutoff)
cat(paste0('significant p-adjust = ', nrow(AluVsConSig), '\n'))

IndVsConSig <- subset(IndVsCon, IndVsCon$padj <= pcutoff)
cat(paste0('significant p-adjust = ', nrow(IndVsConSig), '\n'))

AluIndVsConSig <- subset(AluIndVsCon, AluIndVsCon$padj <= pcutoff)
cat(paste0('significant p-adjust = ', nrow(AluIndVsConSig), '\n'))

# subset up regulated genes
AluVsConSigUp <- subset(AluVsConSig, AluVsConSig$log2FoldChange > 0)
cat(paste0('significant p-adjust and up-regulated = ', nrow(AluVsConSigUp), '\n'))

IndVsConSigUp <- subset(IndVsConSig, IndVsConSig$log2FoldChange > 0)
cat(paste0('significant p-adjust and up-regulated = ', nrow(IndVsConSigUp), '\n'))

AluIndVsConSigUp <- subset(AluIndVsConSig, AluIndVsConSig$log2FoldChange > 0)
cat(paste0('significant p-adjust and up-regulated = ', nrow(AluIndVsConSigUp), '\n'))

# subset down regulated genes (please complete the line 34-41)
AluVsConSigDown <- subset()
cat(paste0('significant p-adjust and up-regulated = ', nrow(AluVsConSigDown), '\n'))

IndVsConSigDown <- subset()
cat(paste0('significant p-adjust and up-regulated = ', nrow(IndVsConSigDown), '\n'))

AluIndVsConSigDown <- subset()
cat(paste0('significant p-adjust and up-regulated = ', nrow(AluIndVsConSigDown), '\n'))

# mapping genes
# data source: OrthoDB v10.1 (https://www.orthodb.org/)
# dma - Daphnia magna (water flea)
# dme - Drosophila melanogaster (fruit fly)
# hsa - Homo sapiens (human)

# example - map D.magna genes to H.sapiens genes
orthoMap <- read.csv("odb10v1_dma_to_hsa_level_33208_gn.tsv", header = F, row.names = 1, sep = '\t', stringsAsFactors = F)

RefOrtho <- subset(orthoMap, row.names(orthoMap) %in% row.names(df))
RefOrtho <- unique(unlist(strsplit(RefOrtho$V2, split = ';')))
length(RefOrtho)
writeLines(RefOrtho, 'reference_orthologs_hsa.txt')

AluVsConSigOrtho <- subset(orthoMap, row.names(orthoMap) %in% row.names(AluVsConSig))
AluVsConSigOrtho <- unique(unlist(strsplit(AluVsConSigOrtho$V2, split = ';')))
length(AluVsConSigOrtho)
writeLines(AluVsConSigOrtho, 'AluVsConSig_orthologs_hsa.txt')

IndVsConSigOrtho <- subset(orthoMap, row.names(orthoMap) %in% row.names(IndVsConSig))
IndVsConSigOrtho <- unique(unlist(strsplit(IndVsConSigOrtho$V2, split = ';')))
length(IndVsConSigOrtho)
writeLines(IndVsConSigOrtho, 'IndVsConSigOrtho_orthologs_hsa.txt')

AluIndVsConSigOrtho <- subset(orthoMap, row.names(orthoMap) %in% row.names(AluIndVsConSig))
AluIndVsConSigOrtho <- unique(unlist(strsplit(AluIndVsConSigOrtho$V2, split = ';')))
length(AluIndVsConSigOrtho)
writeLines(AluIndVsConSigOrtho, 'AluIndVsConSigOrtho_orthologs_hsa.txt')

# please modify line 50-70 for mapping D.magna genes to D.melanogaster genes



# pathway over-representation analysis
# move on to WEB-based GEne SeT AnaLysis Toolkit (https://www.webgestalt.org/) to complete this step

# If use Homo sapiens genes for mapping
# odb10v1_dma_to_hsa_level_33208_gn.tsv - gene symbol

# If use Drosophila melanogaster for mapping
# odb10v1_dma_to_dme_level_6656.tsv - ensembl gene id



