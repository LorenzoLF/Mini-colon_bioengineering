############## Spatiotemporally resolved ex vivo colorectal cancer development in engineered mini-colons##############
############## Mini-colon single-cell RNA-Seq analysis  ##############



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

##Single-cell RNA-Seq data (available in Resources folder)
rawdata <- Read10X(data.dir = "raw_feature_bc_matrix", gene.column = 2)

##Seurat object
SO <- CreateSeuratObject(counts = rawdata, project = "Mini-colon", min.cells = 5, min.features = 0, assay = "RNA")

##Gene signatures (will be used for cell classification)
#Signatures for mini-gut cell types from PMID: 32939089
stem_prog_ref <- c("Olfm4","Lgr5","Cdca7","Smoc2","Slc21a2","Clca3b","Ifitm3","Soat1","Axin2","Igfbp4","Shmt1","Rnf43","Rgmb","Gstt2","Gas6","Aqp1","Tnfrsf19","Cftr","Ccdc34","Apex1")
dividing_ref <- c("Ube2c","Mki67","Birc5","Cenpa","Cdca3","Cdca8","Cenpf","Smc4","Top2a","Ccnb2","Cks2","Pclaf","Cdc20","Spc24","Cenpe","Tuba1b","Cenpw","Ccnb1","Smc2")
enterocyte_ref <- c("Aldob","Fabp1","Rbp2","Apoa1","Apoa4","Ces2a","Reg1","Reg3a","Prap1","Gstm3","Aldh1a1","Spink1","Gsta1","Tkfc","Cyp2b10","Maoa","Cyp4f14","Cbr3","Cyp3a11","Cyp2c29","Mttp","Ada","Cyp3a13","Clca4a","Ces2c","Clca4b","Prl2c3","Hsd17b11","Krt20","Guca2b","Muc3","Cyp2c66","Dgat2","Me1")
paneth_ref <- c("Defa30","Defa29","Ang4","Defa22","Defa17","Defa34","Defa24","Defa26","Defa21","Lyz1","Defa23","Defa36","Mmp7","Defa35","Defa3","Defa27")
goblet_ref <- c("Fcgbp","Txndc5","Agr2","Creb3l1","Rep15","Tff3","Spdef","Galnt7","Galnt12","Sdf2l1","Asph","Klf4","Kctd14","Yipf1","Ern2","Pdia5","Atoh1","Stx17","Ang")
tuft_ref <- c("Rgs13","Dclk1","Alox5ap","Ltc4s","Lrmp","Avil","Fyb","Kctd12","Gng13","Sh2d6","Hck","Nrep","Trpm5","Pygl","Atp2a3","Bpgm","Pik3r5","Hpgds","Pou2f3")
enteroendocrine_ref <- c("Chgb","Gcg","Chga","Tac1","Ghrl","Sct","Cck","Resp18", "Hmgn3","Tph1","Scg2","Neb","Etv1","Cpe","Bex2","Pcsk1n","Scn3a","Cadps","Slc18a1","Rasd1")
signatures.celltype <- list()
signatures.celltype[["Stem_Prog_ref"]] <- stem_prog_ref
signatures.celltype[["Dividing_ref"]] <- dividing_ref
signatures.celltype[["Enterocyte_ref"]] <- enterocyte_ref
signatures.celltype[["Paneth_ref"]] <- paneth_ref
signatures.celltype[["Goblet_ref"]] <- goblet_ref
signatures.celltype[["Tuft_ref"]] <- tuft_ref
signatures.celltype[["Enteroendocrine_ref"]] <- enteroendocrine_ref
#Signatures for cancer cell stemness evaluation
signatures.csc <- list()
signatures.csc[["CancerStemCell"]] <- c("Cd44", "Lgr5", "Sox9")
signatures.csc[["CancerDiffCell"]] <- c("Krt20", "Apoc2", "Fabp2")
#Signatures for enterocyte zonation from PMID: 30270040 (available in Resources folder)
signatures <- read.table("PMID30270040_cluster_markers.txt", sep = "\t", fill = T)
signatures.zonation <- list()
for(i in 1:nrow(signatures)) {
  name <- signatures[i,1]
  genes <- as.character(signatures[i,c(3:ncol(signatures))])
  signatures.zonation[[name]] <- genes
}
#Signatures for CMS scoring from PMID: 35773407 (available in Resources folder)
signatures <- read.table("PMID35773407_iCMS_signatures.txt", sep = "\t", fill = T)
signatures.cms <- list()
for(i in 1:nrow(signatures)) {
  name <- signatures[i,1]
  genes <- as.character(signatures[i,c(3:ncol(signatures))])
  signatures.cms[[name]] <- genes
}



##### ------------------------------------------- 2. Pre-processing -------------------------------------------

##Preliminary subsetting, feature selection and normalization for cell QC
SO <- subset(SO, subset = nFeature_RNA > 100)
SO <- ComputeMitoContent(SO)
SO <- FindVariableFeatures(SO, nfeatures = 2000)
AnnotatedVariableFeaturePlot(SO)
SO <- NormalizeData(SO)

##Preliminary classification for cell QC
#Cell type
SO <- Impute(SO, signatures.celltype)
SO <- ScoreSignatures(SO, signatures.celltype)
SO <- Classify(SO,signatures.celltype, metadata.name = "celltype")
#Cell stemness
SO <- Impute(SO,signatures.csc)
SO <- ScoreSignatures(SO,signatures.csc)
SO <- Classify(SO,signatures.celltype, metadata.name = "stemness")
#Enterocyte zonation
SO <- Impute(SO,signatures.zonation)
SO <- ScoreSignatures(SO,signatures.zonation)
SO <- Classify(SO,signatures.zonation, metadata.name = "enterocytezonation")

##Cell QC and filtering
IdentPlotQC(SO, ident = "celltype")
SO <- subset(SO, subset = nFeature_RNA > 3000 & percent.mito < 20)
IdentPlotQC(SO, ident = "celltype")



##### ------------------------------------------- 3. Normalization, feature selection and data scaling -------------------------------------------

##Normalization, feature selection and data scaling
SO <- NormalizeData(object = SO, normalization.method = "LogNormalize", scale.factor = 10000)
SO <- FindVariableFeatures(SO, nfeatures = 2000)
SO <- ScaleData(SO,features = VariableFeatures(SO))



##### ------------------------------------------- 4. Linear reduction and dimensionality determination -------------------------------------------

##Linear reduction and dimensionality determination
SO <- RunPCA(object = SO)
DimPlot(SO, reduction = "pca")
ElbowPlot(SO, 50)



##### ------------------------------------------- 5. Cell clustering -------------------------------------------

##Cell clustering
#Cluster identification
SO <- FindNeighbors(SO, dims = 1:10)
SO <- FindClusters(SO, resolution = 1.7)
#Cluster labeling based on gene expression
new.cluster.ids <- c("ISC", "Middle_Enterocyte", "Bottom_Enterocyte", "Bottom_Enterocyte", "Top_Enterocyte", "Middle_Enterocyte", "Bottom_Enterocyte", "Top_Enterocyte", "ISC", "Proliferative", "Middle_Enterocyte", "Prog/Immature", "Prog/Immature", "Proliferative", "Goblet", "Prog/Immature", "Enteroendocrine")
names(new.cluster.ids) <- levels(SO)
SO <- SetIdent(SO, value = "seurat_clusters")
SO <- RenameIdents(SO, new.cluster.ids)
SO$celltype_curated <- Idents(SO)
#DotPlot showing the expression of representative markers
SO@active.ident <- factor(SO@active.ident, 
                          levels=c("ISC", 
                                   "Proliferative",
                                   "Prog/Immature", 
                                   "Bottom_Enterocyte", 
                                   "Middle_Enterocyte", 
                                   "Top_Enterocyte", 
                                   "Goblet", 
                                   "Enteroendocrine"))
DotPlot(SO, features = c("Lgr5","Mki67","Cd44","Krt20","Aldob","Clca4a","Muc2","Tff3","Neurod1","Chga"), cols = c("lightgrey", "red"), dot.scale = 8) + RotatedAxis()



##### ------------------------------------------- 6. Non-linear dimensional reduction -------------------------------------------

##UMAP
SO <- RunUMAP(SO,reduction = "pca",dims = 1:10, n.neighbors=10, min.dist = 1.5)
col_custom <- c("#7FC97F","#BEAED4","#386CB0","#6CA6CD", "#666666","#458B00","#FA8072","#FDC086")
DimPlot(SO,pt.size=3, reduction = "umap", cols = col_custom, group.by = "celltype_curated")

##Feature plots for representative markers
FeaturePlot(SO,c("Lgr5"), cols = c("lightgrey","red"), reduction = "umap", pt.size = 2)
FeaturePlot(SO,c("Mki67"), cols = c("lightgrey","red"), reduction = "umap", pt.size = 2)
FeaturePlot(SO,c("Cd44"), cols = c("lightgrey","red"), reduction = "umap", pt.size = 2)
FeaturePlot(SO,c("Muc2"), cols = c("lightgrey","red"), reduction = "umap", pt.size = 2)
FeaturePlot(SO,c("Neurod1"), cols = c("lightgrey","red"), reduction = "umap", pt.size = 2)
FeaturePlot(SO,c("Krt20"), cols = c("lightgrey","red"), reduction = "umap", pt.size = 2)
FeaturePlot(SO,c("Aldob"), cols = c("lightgrey","red"), reduction = "umap", pt.size = 2)
FeaturePlot(SO,c("Clca4a"), cols = c("lightgrey","red"), reduction = "umap", pt.size = 2)



##### ------------------------------------------- 7. Celltag-based lineage tracing -------------------------------------------


##Celltag information loading (available in Resources folder, obtained following the pipeline from https://github.com/morris-lab/CellTagR)
celltags <- as.data.frame(read.table("celltags.txt",header = T, sep = "\t"))
rownames(celltags) <- celltags$cell.barcode
celltags$cell.barcode <- NULL
SO$clones <- celltags
SO <- SetIdent(SO, value = "celltype_curated")
DimPlot(SO, pt.size=1, reduction = "umap", split.by = "clones",cols = col_custom, ncol = 12)

##Tumor clones analyses based on the parameters described in Materials and Methods
SO <- SetIdent(SO, value = "clones")
SO.clone.tumor <- subset(x = SO, idents = c("1","2","6","8","10","11","13","14","17","20","24","25","32","33","38","48","50","51"))
SO.clone.tumor <- SetIdent(SO.clone.tumor, value = "celltype_curated")
DimPlot(SO.clone.tumor,pt.size=2, reduction = "umap", split.by = "clones",cols = col_custom, ncol = 6)
DimPlot(SO.clone.tumor,pt.size=3, reduction = "umap",cols = col_custom)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
DimPlot(SO.clone.tumor,pt.size=3, reduction = "umap", group.by = "clones", cols = col_vector)
#Heatmap for inter-tumor markers as described in Materials and Methods (requires to scale all features of SO). Gene list available in Resources folder
SO.scale.all <- SO
SO.scale.all <- ScaleData(SO.scale.all, features = rownames(SO.scale.all))
SO.scale.all <- SetIdent(SO.scale.all, value = "clones")
SO.clone.tumor <- subset(x = SO.scale.all, idents = c("1","2","6","8","10","11","13","14","17","20","24","25","32","33","38","48","50","51"))
genes.heatmap <- read.table("genes_heatmap.txt")
genes.heatmap <- as.matrix(genes.heatmap)[,1]
DoHeatmap(SO.clone.tumor, features = genes.heatmap, draw.lines=F) + NoLegend() + scale_fill_gradientn(colors = colorRampPalette(c("navyblue", "white", "white", "firebrick"))(256))
#Violin plots
VlnPlot(SO.clone.tumor, c("Cdkn2a"), pt.size = 2, group.by = "clones")
VlnPlot(SO.clone.tumor, c("Il1a"), pt.size = 2, group.by = "clones")
VlnPlot(SO.clone.tumor, c("Slpi"), pt.size = 2, group.by = "clones")
VlnPlot(SO.clone.tumor, c("Aqp5"), pt.size = 2, group.by = "clones")
VlnPlot(SO.clone.tumor, c("Prdm16"), pt.size = 2, group.by = "clones")
VlnPlot(SO.clone.tumor, c("Lgr5"), pt.size = 2, group.by = "clones")
VlnPlot(SO.clone.tumor, c("Cd44"), pt.size = 2, group.by = "clones")
#CMS scoring (global and in several clones)
SO.clone.tumor <- Impute(SO.clone.tumor,signatures.cms)
SO.clone.tumor <- ScoreSignatures(SO.clone.tumor,signatures.cms)
SO.clone.tumor <- Classify(SO.clone.tumor,signatures.cms, metadata.name = "CMS")
table(Idents(SO.clone.tumor))
FeaturePlot(SO.clone.tumor, features = c("iCMS2"), cols = c("lightgrey","red"), pt.size = 3)
FeaturePlot(SO.clone.tumor, features = c("iCMS3"), cols = c("lightgrey","red"), pt.size = 3)
#Highlight some clusters in DimPlot
SO <- SetIdent(SO, value = "clones")
cluster1.cells <- WhichCells(SO, idents = c("1"))
cluster8.cells <- WhichCells(SO, idents = c("8"))
cluster14.cells <- WhichCells(SO, idents = c("14"))
DimPlot(SO.clone.tumor, cells.highlight= list(cluster1.cells, cluster14.cells, cluster8.cells), cols.highlight = c("#BEAED4", "#FDC086", "#7FC97F"), cols= "gray90",pt.size=3, sizes.highlight = 3)
#Explore individual clones
SO.clone10 <- subset(x = SO, subset = clones == "10")
FeaturePlot(SO.clone10,c("Mki67"), cols = c("lightgrey","red"), reduction = "umap", pt.size = 4)
FeaturePlot(SO.clone10,c("Cd44"), cols = c("lightgrey","red"), reduction = "umap", pt.size = 4)
FeaturePlot(SO.clone10,c("Krt20"), cols = c("lightgrey","red"), reduction = "umap", pt.size = 4)
#Scatter plot to compare Gpx2 expression and cancer stem cell signature
SO.clone.tumor <- SetIdent(SO.clone.tumor, value = "clones")
FeatureScatter(object = SO.clone.tumor, feature1 = 'CancerStemCell', feature2 = 'Gpx2', pt.size = 3, cols = col_vector)

##Healthy clones analyses based on the parameters described in Materials and Methods
SO <- SetIdent(SO, value = "clones")
SO.clone.healthy <- subset(x = SO, idents = c("4","9","15","16","19","21","29","31","34","35","42","46","52","56","61","62"))
SO.clone.healthy <- SetIdent(SO.clone.healthy, value = "celltype_curated")
DimPlot(SO.clone.healthy,pt.size=3, reduction = "umap",cols = col_custom)
DimPlot(SO.clone.healthy,pt.size=3, reduction = "umap", group.by = "clones", cols = col_vector)

##Combined plot with healthy and tumor clones
SO <- SetIdent(SO, value = "clones")
SO.clone <- subset(x = SO, idents = c("4","9","15","16","19","21","29","31","34","35","42","46","52","56","61","62","1","2","6","8","10","11","13","14","17","20","24","25","32","33","38","48","50","51"))
SO.clone <- SetIdent(SO.clone, value = "celltype_curated")
DimPlot(SO.clone,pt.size=2, reduction = "umap", split.by = "clones",cols = col_custom, ncol = 9)

