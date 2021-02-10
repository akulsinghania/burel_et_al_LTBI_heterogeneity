################################################################################################
######## Load required libraries
################################################################################################

library(beadarray)
library(lumi)
library(sva)
library(limma)

################################################################################################
######## Read in Raw Illumina files (idat files)
################################################################################################

idatfiles <- dir(pattern=".idat")

x=readIdatFiles(idatFiles = idatfiles)
x<-addFeatureData(x,toAdd = "SYMBOL",annotation = "Humanv4")
save.image("Illumina_idat.RData")

write.table(fData(x),"Annotations.txt",sep="\t")

################################################################################################
######## Process raw data into a normalized matrix
################################################################################################

## create ExpressionSet object
data<-ExpressionSet(assayData=exprs(x),featureData=AnnotatedDataFrame(fData(x)))

## log2 tranform data
t<-lumiT(data,method="log2")

## Quantile normalize the log transformed data
norm<-lumiN(t,method="quantile")

boxplot(norm,las=2,cex.axis=0.4)

write.table(exprs(norm),"Norm_matrix_allprobes.txt",sep="\t")

################################################################################################
######## MA plots to check quality
################################################################################################

library(affyPLM)
pdf("MAplots_log2_q_norm.pdf")
MAplot(norm, cex = 1.0)
dev.off()

################################################################################################
######## Filter data
################################################################################################

## Filter the Norm_matrix_allprobes.txt
## Filtering steps:
## 1) Calculate Median expression value for each gene across all samples, and take bottom 50% probes with lowest median
## 2) Calculate coefficient of variation (COV) expression value for each gene across all samples, and take bottom 50% probes with lowest COV
## 3) Calculate inter quartile range (IQR) expression value for each gene across all samples, and take bottom 50% probes with lowest IQR
## Overlap these 3 lists and take the intersect (those probes that had lowest Median and lowest COV and lowest IQR) - and remove these probes from the main matrix to obtain a filtered matrix with a total of 27960 probes
## Save file as Norm_matrix_filtered.txt


################################################################################################
######## Annotation of Illumina ID probes to Gene symbols
################################################################################################














################################################################################################
######## Principal Component Analysis (PCA)
################################################################################################

## import normalized and filtered expression matrix
data <- read.delim("Norm_matrix_filtered.txt", check.names=FALSE, stringsAsFactors=FALSE)
dim(data)

## separate Gene annotation columns and Expression columns
ann<-data[,1:2]
ann<-ann[order(ann$Probes),]

exprs<-data[,3:155]
rownames(exprs)<-data[,1]

## import patient metadata file
targets<-read.table("MOS_Annotations.txt", header=TRUE, row.names=1, sep="\t")

## make sure samples are in the same order in expression matrix file and metadata file
exprs<-exprs[,sort(colnames(exprs))]
targets<-targets[sort(rownames(targets)),]
colnames(exprs)==rownames(targets)


## calculate principal components (PC) and vairation from each PC
tmatrix<-t(exprs)
pcs<-prcomp(tmatrix)
#summary(pcs)
 
PC1<-pcs$x[,1]
PC2<-pcs$x[,2]
PC3<-pcs$x[,3]
PC4<-pcs$x[,4]

eigs <- pcs$sdev^2
eigs[1] / sum(eigs)

PCA_details <- cbind(PC1, PC2, PC3, PC4)


## Plot PC1 vs PC2
plot(PC1, PC2, col=as.factor(targets$Batch) pch=as.factor(targets$Group), xlab=paste("PC1 (",round((eigs[1]/sum(eigs)*100),2),"%)",sep=""), ylab=paste("PC2 (",round((eigs[2]/sum(eigs)*100),2),"%)",sep=""), cex=1.25, lwd=2,cex.axis=2,cex.lab=1.5)



################################################################################################
######## Batch correction
################################################################################################

combat_correction<-ComBat(as.matrix(exprs),targets$Batch, mod=targets$Group)

write.table(combat_correction, "Combat_corrected_Norm_matrix_filtered.txt",sep="\t")


################################################################################################
######## Differential gene expression analysis (DEG)
################################################################################################

## using limma

Treat <- factor(targets$Group)
Treat

design <- model.matrix(~0+Treat)
colnames(design)

fit <- lmFit(exprs,design)

## figuring in the subject IDs in the design to account for paired samples
corfit <- duplicateCorrelation(exprs,design,block=targets$Subject)
corfit$consensus

fit <- lmFit(exprs,design,block=targets$Subject,correlation=corfit$consensus)
colnames(fit)


## Make contrasts between the comparison groups of interest
## Comparison 1: Active TB (ATB) - Post vs. Pre
## Comparison 2: LTBI-Risk - Post vs. Pre
## Comparison 3: LTBI-Other - Post vs. Pre

cm <- makeContrasts(
  ATB=ATB_Post-ATB_Pre,
  LTBIRisk=LTBIRisk_Post-LTBIRisk_Pre,
  LTBIOther=LTBIOther_Post-LTBIOther_Pre
  levels=design)


fit2 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit2)


## write the DEGs to txt files
ATB<-topTable(fit2, coef="ATB", n=30000)
ATB<-ATB[order(rownames(ATB)),]
ATB<-cbind(ATB,ann)
ATB<-ATB[order(ATB$adj.P.Val),]
write.table(ATB, file="DEG_ATB_Post_vs_Pre.txt", quote=FALSE, row.names=TRUE,sep="\t")

LTBIRisk_Post_Pre<-topTable(fit2, coef="LTBIRisk_Post_Pre", n=30000)
LTBIRisk_Post_Pre<-LTBIRisk_Post_Pre[order(rownames(LTBIRisk_Post_Pre)),]
LTBIRisk_Post_Pre<-cbind(LTBIRisk_Post_Pre,ann)
LTBIRisk_Post_Pre<-LTBIRisk_Post_Pre[order(LTBIRisk_Post_Pre$adj.P.Val),]
write.table(LTBIRisk_Post_Pre, file="DEG_LTBI_Risk_Post_vs_Pre.txt", quote=FALSE, row.names=TRUE,sep="\t")

LTBIOther_Post_Pre<-topTable(fit2, coef="LTBIOther_Post_Pre", n=30000)
LTBIOther_Post_Pre<-LTBIOther_Post_Pre[order(rownames(LTBIOther_Post_Pre)),]
LTBIOther_Post_Pre<-cbind(LTBIOther_Post_Pre,ann)
LTBIOther_Post_Pre<-LTBIOther_Post_Pre[order(LTBIOther_Post_Pre$adj.P.Val),]
write.table(LTBIOther_Post_Pre, file="DEG_LTBI_Other_Post_vs_Pre.txt", quote=FALSE, row.names=TRUE,sep="\t")



################################################################################################
################################################################################################
################################################################################################
