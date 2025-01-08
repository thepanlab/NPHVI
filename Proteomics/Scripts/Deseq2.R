# module load R/3.6.2-foss-2019b

library('DESeq2')
library('ggplot2')
library('ggrepel')

countsName <- "/ourdisk/hpc/rnafold/dywang/dont_archive/On_process/V45/Sidbers/Proteome/Deseq/combined_proteome.txt"


countData <- read.table(countsName, header = TRUE, sep = "\t",stringsAsFactors=F)  # sort -u file_name >new_file_name
head(countData)

metaDataName <- "/ourdisk/hpc/rnafold/dywang/dont_archive/On_process/V45/Sidbers/Proteome/Deseq/MetaData.txt"


metaData <- read.table(metaDataName, header = TRUE, sep = "\t")
head(metaData)

dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design=~dex, tidy = TRUE)

dds <- DESeq(dds)

# output results

res <- results(dds)
res <- res[order(res$padj),]
head(res)
write.table(res,"/ourdisk/hpc/rnafold/dywang/dont_archive/On_process/V45/Sidbers/Proteome/Deseq/Comp_TPM_Sidbers.txt",quote=F)

# PCA plot
vsdata <- varianceStabilizingTransformation(dds, blind=TRUE)
plotPCA(vsdata, intgroup="dex") 

# volcano plot
res$threhold =  factor(ifelse(res$padj < 0.05 & abs(res$log2FoldChange) >=1, ifelse(res$log2FoldChane>=1, 'Up', 'Down' ), 'NoSignificant'), level = c('Up', 'Down', 'Nosignificant'))
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-30,30)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.05 & log2FoldChange<=0), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.05 & log2FoldChange>=0), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
geom_vline(xintercept=c(-1, 1), lty=3, col = "black", lwd=0.5)
geom_hline(yintercept = -log10(0.05), lty=3, col="black", lwd=0.5)
