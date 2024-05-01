############## Biofabricated patient-derived mini-colons as next-generation colorectal cancer organoids ##############
############## CRC & CAF single-cell RNA-Seq analysis  ##############



##### ------------------------------------------- 1. Setup -------------------------------------------

##Packages
library(Seurat)
library(dplyr)
library(matrixStats)
library(ggplot2)
library(burgertools)
library(plotly)
library(RColorBrewer)
library(Matrix)

##Seed
set.seed(7)

##Theme
theme_set(theme_light())

##Single-cell RNA-Seq data (Cocultures_scRNAseq.zip is available here: https://crc-tme.com/Downloads/)
rawdata <- Read10X(data.dir = "10x input", gene.column = 1)
metadata <- read.table("metadata.tsv", header = T, sep = "\t")
rownames(metadata) <- colnames(rawdata)

##Seurat object
SO <- CreateSeuratObject(counts = rawdata, meta.data = metadata, project = "Coculture", min.cells = 5, min.features = 0, assay = "RNA")
#Subset the patient and conditions of interest
SO <- subset(x = SO, subset = patient == "MS")
Idents(SO) <- "culture"
SO <- subset(x = SO, idents = c("PDO","CAF","PDO+CAF"))
#Add metadata indicating cell type + culture conditions
SO$cellcond <- paste(SO$celltype, SO$culture, sep = ".in.")



##### ------------------------------------------- 2. Pre-processing -------------------------------------------

##Preliminary subsetting, feature selection and normalization for cell QC
SO <- subset(SO, subset = nFeature_RNA > 100)
SO <- ComputeMitoContent(SO)
SO <- FindVariableFeatures(SO,nfeatures = 2000)
AnnotatedVariableFeaturePlot(SO)
SO <- NormalizeData(SO)

##Cell QC and filtering
IdentPlotQC(SO, ident = "celltype")
SO <- subset(SO, subset = nFeature_RNA > 2500 & percent.mito < 20)
table(Idents(SO))
table(SO$celltype)
IdentPlotQC(SO, ident = "celltype")
IdentPlotQC(SO, ident = "culture")



##### ------------------------------------------- 3. Normalization, feature selection and data scaling -------------------------------------------

##Normalization, feature selection and data scaling
SO <- SetIdent(SO, value = "celltype")
DefaultAssay(SO) <- "RNA"
SO <- NormalizeData(object = SO, normalization.method = "LogNormalize", scale.factor = 10000)
SO <- FindVariableFeatures(SO, nfeatures = 2000)
SO <- ScaleData(SO,features = VariableFeatures(SO))



##### ------------------------------------------- 4. Linear reduction and dimensionality determination -------------------------------------------

##Linear reduction and dimensionality determination
SO <- RunPCA(object = SO)
DimPlot(SO, reduction = "pca", cols = c("firebrick2","darkolivegreen3"))
ElbowPlot(SO, 50)



##### ------------------------------------------- 6. Non-linear dimensional reduction -------------------------------------------

##UMAP
SO <- RunUMAP(SO,reduction = "pca",dims = 1:10, n.neighbors = 10, min.dist = 1.0)
DimPlot(SO,pt.size=1, reduction = "umap", cols = c("#FDC086","salmon","#BEAED4","#386CB0"), group.by = "cellcond")

##Feature plots for representative markers
FeaturePlot(SO, "MMP7", pt.size = 2, order=T) & scale_colour_gradientn(colours = (brewer.pal(n = 9, name = "Reds")))
FeaturePlot(SO, "LAMC2", pt.size = 2, order=T) & scale_colour_gradientn(colours = (brewer.pal(n = 9, name = "Reds")))
FeaturePlot(SO, "EMP1", pt.size = 2, order=T) & scale_colour_gradientn(colours = (brewer.pal(n = 9, name = "Reds")))
FeaturePlot(SO, "LGR5", pt.size = 2, order=T) & scale_colour_gradientn(colours = (brewer.pal(n = 9, name = "Reds")))
FeaturePlot(SO, "KRT20", pt.size = 2, order=T) & scale_colour_gradientn(colours = (brewer.pal(n = 9, name = "Reds")))
FeaturePlot(SO, "FABP1", pt.size = 2, order=T) & scale_colour_gradientn(colours = (brewer.pal(n = 9, name = "Reds")))
FeaturePlot(SO, "PHGR1", pt.size = 2, order=T) & scale_colour_gradientn(colours = (brewer.pal(n = 9, name = "Reds")))
FeaturePlot(SO, "MKI67", pt.size = 2, order=T) & scale_colour_gradientn(colours = (brewer.pal(n = 9, name = "Reds")))
FeaturePlot(SO, "NPSR1", pt.size = 2, order=T) & scale_colour_gradientn(colours = (brewer.pal(n = 9, name = "Reds")))
FeaturePlot(SO, "ACTA2", pt.size = 2, order=T) & scale_colour_gradientn(colours = (brewer.pal(n = 9, name = "Reds")))
FeaturePlot(SO, "TAGLN", pt.size = 2, order=T) & scale_colour_gradientn(colours = (brewer.pal(n = 9, name = "Reds")))
FeaturePlot(SO, "IL6", pt.size = 2, order=T) & scale_colour_gradientn(colours = (brewer.pal(n = 9, name = "Reds")))
FeaturePlot(SO, "HGF", pt.size = 2, order=T) & scale_colour_gradientn(colours = (brewer.pal(n = 9, name = "Reds")))
FeaturePlot(SO, "FAP", pt.size = 2, order=T) & scale_colour_gradientn(colours = (brewer.pal(n = 9, name = "Reds")))

##Signature enrichment
signatures.crctme <- list()
signatures.crctme[["Inv.Met"]] <- c("MMP7","LAMB3","LAMC2","EMP1","TSPAN1")
signatures.crctme[["Stemness"]] <- c("LGR5","CD44","SMOC2","LRIG1","EPHB2")
signatures.crctme[["Differentiation"]] <- c("KRT20","FABP1","FABP3","SELENBP1","PHGR1")
signatures.crctme[["Myofibroblast"]] <- c("ACTA2","TAGLN","LAMA2","ACTG1","MYL12A")
signatures.crctme[["Inflammfibroblast"]] <- c("IL6","IL11","IL32","TNFSF18")
signatures.crctme[["Activfibroblast"]] <- c("HGF","FAP","WNT5A","VEGFA")
SO <- Impute(SO, signatures.crctme)
SO <- ScoreSignatures(SO, signatures.crctme)
#Subset epithelial cells and CAFs
SO <- SetIdent(SO, value = "celltype")
SO.ep <- subset(x = SO, idents = "EP")
SO.caf <- subset(x = SO, idents = "CAF")
#Plot signature enrichment in the relevant cell type
FeaturePlot(SO.ep, names(signatures.crctme)[1], pt.size = 2, order=T) & scale_colour_gradientn(colours = (brewer.pal(n = 9, name = "Reds")))
FeaturePlot(SO.ep, names(signatures.crctme)[2], pt.size = 2, order=T) & scale_colour_gradientn(colours = (brewer.pal(n = 9, name = "Reds")))
FeaturePlot(SO.ep, names(signatures.crctme)[3], pt.size = 2, order=T) & scale_colour_gradientn(colours = (brewer.pal(n = 9, name = "Reds")))
FeaturePlot(SO.caf, names(signatures.crctme)[4], pt.size = 2, order=T) & scale_colour_gradientn(colours = (brewer.pal(n = 9, name = "Reds")))
FeaturePlot(SO.caf, names(signatures.crctme)[5], pt.size = 2, order=T) & scale_colour_gradientn(colours = (brewer.pal(n = 9, name = "Reds")))
FeaturePlot(SO.caf, names(signatures.crctme)[6], pt.size = 2, order=T) & scale_colour_gradientn(colours = (brewer.pal(n = 9, name = "Reds")))

##Violin plots
VlnPlot(SO.ep, "Inv.Met", group.by = "cellcond", cols = c("#BEAED4","#386CB0"))
VlnPlot(SO.ep, "MMP7", group.by = "cellcond", cols = c("#BEAED4","#386CB0"))
