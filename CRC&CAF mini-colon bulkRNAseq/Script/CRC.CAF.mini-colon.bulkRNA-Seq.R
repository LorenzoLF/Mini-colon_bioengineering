############## Biofabricated patient-derived mini-colons as next-generation colorectal cancer organoids ##############
############## Transcriptomic profiling of CRC mini-colons cultured i) in control medium, ii) in CAF-conditioned medium, and iii) in the presence of CAFs  ##############



##### ------------------------------------------- 1. Setup -------------------------------------------

##Packages
library(edgeR)
library(limma)
library(gplots)
library(heatmap3)
library(org.Hs.eg.db)
library(EnhancedVolcano)

##bulk RNA-Seq data (available in Resources folder)
countdata <- read.table("countdata.txt", sep = "\t", header = T)
#Sample info
sampleinfo <- read.table("sampleinfo.txt", sep ="\t", header = T)
colnames(countdata) <- sampleinfo$Name
head(countdata)



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

##Sample distributions
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)

##Data exploration
dc = dist(t(logcounts), method = "euclidean")
hc = hclust(dc, method="complete")
plot(hc)



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

##Contrast matrix creation
contrast.matrix = makeContrasts(
  cm_vs_ctrl = cm - ctrl, 
  caf_vs_ctrl = caf - ctrl, 
  caf_vs_cm = caf - cm, 
  levels = design)

##Contrast calculation
fit.cont <- contrasts.fit(fit, contrast.matrix)

##Empirical Bayes
res.limma <- eBayes(fit.cont)

##Gene annotation
id.symbol <- mapIds(org.Hs.eg.db, keys = rownames(fit.cont), keytype = "ENSEMBL", column="SYMBOL")
res.limma$genes <- id.symbol
v$genes <- id.symbol

##Results
#For multiple comparisons in the contrast matrix, select comparison (choose between 'coef = "cm_vs_ctrl"', 'coef = "caf_vs_ctrl"',  and 'coef = "caf_vs_cm"':
results <- topTable(res.limma, p.value=0.05, lfc = 1, adjust.method="BH", sort.by="p", number=nrow(res.limma), coef = "cm_vs_ctrl")
head(results)

##Heatmap visualization of differentially expressed genes for each comparison (here "cm_vs_ctrl" as an example)
expM <- logcounts[rownames(results),c(which(dgeObj$samples$Group == "ctrl"), which(dgeObj$samples$Group == "cm"))]
png(filename = "Heatmap.png", width = 900, height = 900)
heatmap3(expM, balanceColor = T, labRow = NA, margins = c(10,25), cexRow = 1.5, cexCol = 1.5)
dev.off()
#Check the file "Heatmap.png" in your working directory to see the heatmap

#Volcano visualization of differentially-expressed genes
results.full <- topTable(res.limma, adjust.method="BH", sort.by="p", number=Inf, coef = "cm_vs_ctrl")
png(filename = "Volcano.png", width = 700, height = 700)
EnhancedVolcano(results.full,
                lab = results.full$ID,
                x = "logFC",
                y = "adj.P.Val",
                pCutoff = 0.05,
                FCcutoff = 1,
                ylab = bquote(~-Log[10] ~ italic(Q)),
                title = 'CAF-conditioned medium vs Control',
                labSize = 4,
                drawConnectors = TRUE, typeConnectors = "open", arrowheads = FALSE, maxoverlapsConnectors = Inf,
                selectLab = c("MMP7","TUBB3","IL1B","CD274"),
                ylim = c(0,15),
                gridlines.major = FALSE,
                gridlines.minor = FALSE)
dev.off()
#Check the file "Volcano.png" in your working directory to see the volcano plot
