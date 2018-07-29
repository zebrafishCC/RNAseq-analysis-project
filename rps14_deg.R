setwd("/home/chengchen/data/RNAseq/countfile/rps14")

library("DESeq2")
#import data with htseq-count method in package DESeq2
directory = "/home/chengchen/data/RNAseq/countfile/rps14"
sampleFiles=grep("count",list.files(directory),value = TRUE)
sampleCondition = c("011Rpm","011Rpm","Ctl","Ctl","Ctl","ctRpm","ctRpm")
sampleName = sub("_count.txt","",sampleFiles)
sampleTable = data.frame(sampleName = sampleName, fileName = sampleFiles,condition=sampleCondition)
sampleTable

ddsHTseq = DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,directory = directory,design = ~condition)
head(counts(ddsHTseq))

ddsHTseq$condition =  relevel(ddsHTseq$condition,ref = "Ctl") #important to define the base level

#vsd normalization for PCA analysis
dds=ddsHTseq[rowSums(counts(ddsHTseq))>=1,]
vsd = vst(dds,blind = FALSE)
pdf("rps14clustering.pdf")
plotPCA(vsd) #PCAplot to find sample clustering
dev.off()

#sample distances
sampleDists = dist(t(assay(vsd)))
library("pheatmap")
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pdf("rps14sample_distances.pdf")
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
dev.off()

#gene clustering
library("genefilter")
library("pheatmap")
topVarGenes = head(order(rowVars(assay(vsd)),decreasing = TRUE),20)
selected = assay(vsd)[topVarGenes,]
selected = selected - rowMeans(selected)
anno = as.data.frame(colData(vsd))
pdf("rps14_top20vsd_heatmap.pdf")
pheatmap(selected,annotation_col = anno)
dev.off()

#data filtering rowSums should be larger than 10 reads, improve analysis speed
nrow(ddsHTseq)
dds=ddsHTseq[rowSums(counts(ddsHTseq))>=10,]
nrow(dds) #check the filtering change

#obtain the DEG results
dds=DESeq(dds)

res1 = results(dds,contrast = c("condition","ctRpm","Ctl")) 
res1
summary(res1)
res05 = res1[res1$padj<0.05&!is.na(res1$padj),] #remove rows with padj less than 0.05
sum(res05$padj<0.05,na.rm = TRUE)
#res05ordered = res05[order(res05$padj),]
keep1 = abs(res05$log2FoldChange) > 1 #keep those genes with reasonable fold change
res05_lfc = res05[keep1,]
nrow(res05_lfc)
head(res05_lfc)
ids = rownames(res05_lfc)
###get zfin_gene_symbol for ensembl genes
library("biomaRt")
ensembl=useMart("ensembl")
ensembl=useDataset("drerio_gene_ensembl",mart=ensembl)
symbols <- getBM(attributes=c('ensembl_gene_id','zfin_id_symbol',"go_id", "name_1006","namespace_1003"),filters='ensembl_gene_id',values=ids, mart=ensembl)
symbols=symbols[!duplicated(symbols$ensembl_gene_id),]
res05_lfc = as.data.frame(res05_lfc)
res05_lfc = cbind(res05_lfc,symbols)
head(res05_lfc)
#res05_lfc=res05_lfc[,-7] #remove ensembl_gene_id column
res05_lfc_ordered = res05_lfc[order(res05_lfc$padj),]
head(res05_lfc_ordered)
#res05ordered_lfc$symbol = mapIds(org.Dr.eg.db,keys = row.names(res05ordered_lfc),column = "SYMBOL",keytype = "ENSEMBL", multiVals="first")
write.csv(res05_lfc_ordered,file="ctRpm_vs_CtlDEG.csv") #save ctRpm control information

##generate report for 536_vs_control DEG
library("ReportingTools")
htmlRep <- HTMLReport(shortName="CtRpm_ctl_report", title="CtRpm_Ctl_report",
                      reportDirectory="./report")
publish(res05_lfc_ordered, htmlRep)
url <- finish(htmlRep)
#browseURL(url)


#011Rpm and control
res2 = results(dds,contrast = c("condition","011Rpm","Ctl"))
summary(res2)
res05 = res2[res2$padj<0.05&!is.na(res2$padj),] #remove rows with padj less than 0.05
sum(res05$padj<0.05,na.rm = TRUE)
#res05ordered = res05[order(res05$padj),]
keep1 = abs(res05$log2FoldChange) > 1 #keep those genes with reasonable fold change
res05_lfc = res05[keep1,]
nrow(res05_lfc)
head(res05_lfc)
ids = rownames(res05_lfc)
###get zfin_gene_symbol for ensembl genes
symbols <- getBM(attributes=c('ensembl_gene_id','zfin_id_symbol',"go_id", "name_1006","namespace_1003"),filters='ensembl_gene_id',values=ids, mart=ensembl)
symbols=symbols[!duplicated(symbols$ensembl_gene_id),]
res05_lfc = as.data.frame(res05_lfc)
res05_lfc = cbind(res05_lfc,symbols)
head(res05_lfc)
#res05_lfc=res05_lfc[,-7] #remove ensembl_gene_id column
res05_lfc_ordered = res05_lfc[order(res05_lfc$padj),]
head(res05_lfc_ordered)
#res05ordered_lfc$symbol = mapIds(org.Dr.eg.db,keys = row.names(res05ordered_lfc),column = "SYMBOL",keytype = "ENSEMBL", multiVals="first")
write.csv(res05_lfc_ordered,file="001Rpm_vs_CtlDEG.csv") #save 011Rpm_vs_Ctl information

##generate report for 536_vs_control DEG
library("ReportingTools")
htmlRep <- HTMLReport(shortName="011Rpm_ctl_report", title="011Rpm_Ctl_report",
                      reportDirectory="./report")
publish(res05_lfc_ordered, htmlRep)
url <- finish(htmlRep)
#browseURL(url)

#011_inject control data quality analysis
##MA plot to see the pattern of DEG genes
library("apeglm")
resultsNames(dds)
res1 = lfcShrink(dds,coef = "condition_011Rpm_vs_Ctl",type = "apeglm")
pdf("011Rpm_controlMAplot.pdf") #save MAplot to pdf
plotMA(res1,ylim=c(-5,5))
dev.off()

##plot pvalue
res1 = results(dds,contrast = c("condition","011Rpm","Ctl"))
pdf("011Rpm_Ctl_pvalue.pdf")
hist(res1$pvalue[res1$baseMean > 1], breaks = 0:20/20,col = "grey50", border = "white")
dev.off()

#0536_inject control data quality analysis
#MAplot
res2 = lfcShrink(dds,coef = "condition_ctRpm_vs_Ctl",type = "apeglm")
pdf("ctRpm_CtlMAplot.pdf") #save MAplot to pdf
plotMA(res2,ylim=c(-5,5))
dev.off()

##plot pvalue
res2 = results(dds,contrast = c("condition","ctRpm","Ctl"))
pdf("ctRpm_Ctl_pvalue.pdf")
hist(res2$pvalue[res2$baseMean > 1], breaks = 0:20/20,col = "grey50", border = "white")
dev.off()