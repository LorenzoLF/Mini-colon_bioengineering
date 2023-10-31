############## Spatiotemporally resolved ex vivo colorectal cancer development in engineered mini-colons##############
############## Mini-colon AKP lines bulk RNA-Seq analysis: differential expression & In vivo AKP signature  ##############



##### ------------------------------------------- 1. Setup -------------------------------------------

##Packages
library(edgeR)
library(limma)
library(org.Mm.eg.db)

##Batch-corrected bulk RNA-Seq data (available in Resources folder)
countdata <- read.table("countdata.txt", sep = "\t", header = T)
#Sample info
sampleinfo <- read.table("sampleinfo.txt", sep ="\t", header = T)



##### ------------------------------------------- 2. Gene filtering -------------------------------------------

##Gene filtering
myCPM <- cpm(countdata)
thresh <- myCPM > 1
keep <- rowSums(thresh) >= 3
counts.keep <- countdata[keep,]

##Conversion to DGEList object
dgeObj <- DGEList(counts.keep, samples = sampleinfo)



##### ------------------------------------------- 3. Normalization -------------------------------------------

##Count normalization with TMM
dgeObj <- calcNormFactors(dgeObj, method = "TMM")
#log2 counts per million
logcounts <- cpm(dgeObj,log=TRUE)
colnames(logcounts) <- dgeObj$samples$Subgroup

##Data exploration
dc = dist (t(logcounts), method = "euclidean")
hc = hclust (dc, method = "average")
plot(hc, labels = sampleinfo$Name)



##### ------------------------------------------- 4. Differential expression -------------------------------------------

##Design matrix creation
group <- dgeObj$samples$Group
group <- factor(group)
design <- model.matrix(~ 0 + group)
#Make the column names of the design matrix nicer
colnames(design) <- gsub("group", "", colnames(design))

##Voom-transformation of the data
v <- voom(dgeObj, design, plot = TRUE)

##Linear model fit
fit <- lmFit(v)

##Contrast matrix creation (for in vivo AKP signature: in vivo AKP vs organoid AKP)
contrast.matrix = makeContrasts(
  vAKPt_vs_org = invivo_AKP - org_AKP,
  levels = design)

##Contrast calculation
fit.cont <- contrasts.fit(fit, contrast.matrix)

##Empirical Bayes
res.limma <- eBayes(fit.cont)

##Gene annotation
id.symbol <- mapIds(org.Mm.eg.db, keys = rownames(fit.cont), keytype = "ENSEMBL", column="SYMBOL")
res.limma$genes <- id.symbol
v$genes <- id.symbol

##Results
results <- topTable(res.limma, p.value=0.05, lfc=1, adjust.method="BH", sort.by="p", number=nrow(res.limma), coef = "vAKPt_vs_org")
head(results)

#*Please see script in Mini-colon.AKP.Gpx2.bulkRNA-Seq.R for examples on heatmap and volcano plot generation
