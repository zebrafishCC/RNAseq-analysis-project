#remember to set working directory
#open files containing drug injection count files
setwd("/home/chengchen/data/RNAseq/countfile/drug_injection_only")
dir()
#open summarize count file, count file was obtatined from htseq-count with data integrataion
count = read.csv("drug_final_summarize_count.csv",row.names = 1)
head(count)
#open coldata containing sample information
coldata = read.csv("coldata.csv",row.names = 1)
head(coldata)

#make sure rowname of coldata is consistent with colnames of count data
rownames(coldata)
colnames(count)
colnames(count)=rownames(coldata)
head(count)
all(rownames(coldata)==colnames(count))

#import data for further DESeq2 analyis
library("DESeq2")
dds = DESeqDataSetFromMatrix(countData = count,colData = coldata,design=~condition)
head(coldata)
dds$condition =  relevel(dds$condition,ref = "control") #important to define the base level

#vsd normalization for PCA analysis
#dds=dds[rowSums(counts(dds)>=1),]
#vsd = vst(dds,blind = FALSE)
#plotPCA(vsd) #PCAplot to find sample clustering

#gene clustering
library("genefilter")
library("pheatmap")
topVarGenes = head(order(rowVars(assay(vsd)),decreasing = TRUE),20)
selected = assay(vsd)[topVarGenes,]
selected = selected - rowMeans(selected)
anno = as.data.frame(colData(vsd))
pdf("drug_top20vsd_heatmap.pdf")
pheatmap(selected,annotation_col = anno)
dev.off()

#data filtering rowSums should be larger than 10 reads, improve analysis speed
nrow(dds)
keep = rowSums(counts(dds))>=10 
dds = dds[keep,]
nrow(dds) #check the filtering change

#add gene id
library("biomaRt")
ensembl=useMart("ensembl")
ensembl=useDataset("drerio_gene_ensembl",mart=ensembl)

#obtain the DEG results
dds=DESeq(dds)


#011inject and control
res1 = results(dds,contrast = c("condition","011_inject","control")) 
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
symbols <- getBM(attributes=c('ensembl_gene_id','zfin_id_symbol',"go_id", "name_1006","namespace_1003"),filters='ensembl_gene_id',values=ids, mart=ensembl)
symbols=symbols[!duplicated(symbols$ensembl_gene_id),]
res05_lfc = as.data.frame(res05_lfc)
res05_lfc = cbind(res05_lfc,symbols)
head(res05_lfc)
#res05_lfc=res05_lfc[,-7] #remove ensembl_gene_id column
res05_lfc_ordered = res05_lfc[order(res05_lfc$padj),]
head(res05_lfc_ordered)
#res05ordered_lfc$symbol = mapIds(org.Dr.eg.db,keys = row.names(res05ordered_lfc),column = "SYMBOL",keytype = "ENSEMBL", multiVals="first")
write.csv(res05_lfc_ordered,file="011_vs_controlDEG.csv") #save 011 control information

##html report for 011-control
library("ReportingTools")
htmlRep <- HTMLReport(shortName="011_control_report", title="011_control_report",
                      reportDirectory="./report")
publish(res05_lfc_ordered, htmlRep)
url <- finish(htmlRep)
##GO annotations
#go <- getBM(attributes=c("go_id", "name_1006","namespace_1003"), mart=ensembl)


#536injet and control
res2 = results(dds,contrast = c("condition","536_inject","control"))
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
write.csv(res05_lfc_ordered,file="536_vs_controlDEG.csv") #save 536 control information

##generate report for 536_vs_control DEG
library("ReportingTools")
htmlRep <- HTMLReport(shortName="536_control_report", title="536_control_report",
                      reportDirectory="./report")
publish(res05_lfc_ordered, htmlRep)
url <- finish(htmlRep)
#browseURL(url)

#011_inject control data quality analysis
##MA plot to see the pattern of DEG genes
library("apeglm")
resultsNames(dds)
res1 = lfcShrink(dds,coef = "condition_011_inject_vs_control",type = "apeglm")
pdf("011_controlMAplot.pdf") #save MAplot to pdf
plotMA(res1,ylim=c(-5,5))
dev.off()

##plot pvalue
res1 = results(dds,contrast = c("condition","011_inject","control"))
pdf("011vsControl_pvalue.pdf")
hist(res1$pvalue[res1$baseMean > 1], breaks = 0:20/20,col = "grey50", border = "white")
dev.off()

#0536_inject control data quality analysis
#MAplot
res2 = lfcShrink(dds,coef = "condition_536_inject_vs_control",type = "apeglm")
pdf("536_controlMAplot.pdf") #save MAplot to pdf
plotMA(res2,ylim=c(-5,5))
dev.off()

##plot pvalue
res2 = results(dds,contrast = c("condition","536_inject","control"))
pdf("536vsControl_pvalue.pdf")
hist(res2$pvalue[res2$baseMean > 1], breaks = 0:20/20,col = "grey50", border = "white")
dev.off()