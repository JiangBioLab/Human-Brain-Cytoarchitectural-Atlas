library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)


setwd('/home/wangpp/DATA/Brain/spatial/SpatialAnalysis')

######################################################
V1 <- Load10X_Spatial('V1-17') 
V1 <- SCTransform(V1, assay = "Spatial", verbose = FALSE)
V1 <- RunPCA(V1, assay = "SCT", verbose = FALSE)
V1 <- FindNeighbors(V1, reduction = "pca", dims = 1:30)
V1 <- FindClusters(V1, verbose = FALSE)
V1 <- RunUMAP(V1, reduction = "pca", dims = 1:30)
saveRDS(V1,'V1.rds')

MTG <- Load10X_Spatial('MTG-32') 
MTG <- SCTransform(MTG, assay = "Spatial", verbose = FALSE)
MTG <- RunPCA(MTG, assay = "SCT", verbose = FALSE)
MTG <- FindNeighbors(MTG, reduction = "pca", dims = 1:30)
MTG <- FindClusters(MTG, verbose = FALSE)
MTG <- RunUMAP(MTG, reduction = "pca", dims = 1:30)
saveRDS(MTG,'MTG.rds')

M1 <- Load10X_Spatial('M1-5') 
M1 <- SCTransform(M1, assay = "Spatial", verbose = FALSE)
M1 <- RunPCA(M1, assay = "SCT", verbose = FALSE)
M1 <- FindNeighbors(M1, reduction = "pca", dims = 1:30)
M1 <- FindClusters(M1, verbose = FALSE)
M1 <- RunUMAP(M1, reduction = "pca", dims = 1:30)
saveRDS(M1,'M1.rds')

OFC <- Load10X_Spatial('OFC-12') 
OFC <- SCTransform(OFC, assay = "Spatial", verbose = FALSE)
OFC <- RunPCA(OFC, assay = "SCT", verbose = FALSE)
OFC <- FindNeighbors(OFC, reduction = "pca", dims = 1:30)
OFC <- FindClusters(OFC, verbose = FALSE)
OFC <- RunUMAP(OFC, reduction = "pca", dims = 1:30)
saveRDS(OFC,'OFC.rds')

OFC4 <- Load10X_Spatial('OFC4') 
OFC4 <- SCTransform(OFC4, assay = "Spatial", verbose = FALSE)
OFC4 <- RunPCA(OFC4, assay = "SCT", verbose = FALSE)
OFC4 <- FindNeighbors(OFC4, reduction = "pca", dims = 1:30)
OFC4 <- FindClusters(OFC4, verbose = FALSE)
OFC4 <- RunUMAP(OFC4, reduction = "pca", dims = 1:30)
saveRDS(OFC4,'OFC4.rds')

V18 <- Load10X_Spatial('V18') 
V18 <- SCTransform(V18, assay = "Spatial", verbose = FALSE)
V18 <- RunPCA(V18, assay = "SCT", verbose = FALSE)
V18 <- FindNeighbors(V18, reduction = "pca", dims = 1:30)
V18 <- FindClusters(V18, verbose = FALSE)
V18 <- RunUMAP(V18, reduction = "pca", dims = 1:30)
saveRDS(V18,'V18.rds')

V19 <- Load10X_Spatial('V19') 
V19 <- SCTransform(V19, assay = "Spatial", verbose = FALSE)
V19 <- RunPCA(V19, assay = "SCT", verbose = FALSE)
V19 <- FindNeighbors(V19, reduction = "pca", dims = 1:30)
V19 <- FindClusters(V19, verbose = FALSE)
V19 <- RunUMAP(V19, reduction = "pca", dims = 1:30)
saveRDS(V19,'V19.rds')


DLPFC <- Load10X_Spatial('DLPFC-8') 
DLPFC <- SCTransform(DLPFC, assay = "Spatial", verbose = FALSE)
DLPFC <- RunPCA(DLPFC, assay = "SCT", verbose = FALSE)
DLPFC <- FindNeighbors(DLPFC, reduction = "pca", dims = 1:30)
DLPFC <- FindClusters(DLPFC, verbose = FALSE)
DLPFC <- RunUMAP(DLPFC, reduction = "pca", dims = 1:30)
saveRDS(DLPFC,'DLPFC.rds')


DLPFC <- Load10X_Spatial('DLPFC-9') 
DLPFC <- SCTransform(DLPFC, assay = "Spatial", verbose = FALSE)
DLPFC <- RunPCA(DLPFC, assay = "SCT", verbose = FALSE)
DLPFC <- FindNeighbors(DLPFC, reduction = "pca", dims = 1:30)
DLPFC <- FindClusters(DLPFC, verbose = FALSE)
DLPFC <- RunUMAP(DLPFC, reduction = "pca", dims = 1:30)
saveRDS(DLPFC,'DLPFC2.rds')


S128 <- Load10X_Spatial('S128-1') 
S128 <- SCTransform(S128, assay = "Spatial", verbose = FALSE)
S128 <- RunPCA(S128, assay = "SCT", verbose = FALSE)
S128 <- FindNeighbors(S128, reduction = "pca", dims = 1:30)
S128 <- FindClusters(S128, verbose = FALSE)
S128 <- RunUMAP(S128, reduction = "pca", dims = 1:30)
saveRDS(S128,'S1.rds')


S128 <- Load10X_Spatial('S128-2') 
S128 <- SCTransform(S128, assay = "Spatial", verbose = FALSE)
S128 <- RunPCA(S128, assay = "SCT", verbose = FALSE)
S128 <- FindNeighbors(S128, reduction = "pca", dims = 1:30)
S128 <- FindClusters(S128, verbose = FALSE)
S128 <- RunUMAP(S128, reduction = "pca", dims = 1:30)
saveRDS(S128,'S1-2.rds')


WZHT7 <- Load10X_Spatial('WZHT-7') 
WZHT7 <- SCTransform(WZHT7, assay = "Spatial", verbose = FALSE)
WZHT7 <- RunPCA(WZHT7, assay = "SCT", verbose = FALSE)
WZHT7 <- FindNeighbors(WZHT7, reduction = "pca", dims = 1:30)
WZHT7 <- FindClusters(WZHT7, verbose = FALSE)
WZHT7 <- RunUMAP(WZHT7, reduction = "pca", dims = 1:30)
saveRDS(WZHT7,'WZHT7.rds')


XRH1 <- Load10X_Spatial('XRH-1') 
XRH1 <- SCTransform(XRH1, assay = "Spatial", verbose = FALSE)
XRH1 <- RunPCA(XRH1, assay = "SCT", verbose = FALSE)
XRH1 <- FindNeighbors(XRH1, reduction = "pca", dims = 1:30)
XRH1 <- FindClusters(XRH1, verbose = FALSE)
XRH1 <- RunUMAP(XRH1, reduction = "pca", dims = 1:30)
saveRDS(XRH1,'XRH1.rds')



######################################################



V1 <- Load10X_Spatial('V1') 

plot1 <- VlnPlot(V1, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(V1, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

#saveRDS(plot2,'V1.nCount.spacial.rds')

V1 <- SCTransform(V1, assay = "Spatial", verbose = FALSE)

SpatialFeaturePlot(V1, features = c("HPCA", "TTR"))

p1 <- SpatialFeaturePlot(V1, features = "HPCA", pt.size.factor = 1)
p2 <- SpatialFeaturePlot(V1, features = "HPCA", alpha = c(0.1, 1))
p1 + p2


V1 <- RunPCA(V1, assay = "SCT", verbose = FALSE)
V1 <- FindNeighbors(V1, reduction = "pca", dims = 1:30)
V1 <- FindClusters(V1, verbose = FALSE)
V1 <- RunUMAP(V1, reduction = "pca", dims = 1:30)

#saveRDS(V1,'V1.rds')

p1 <- DimPlot(V1, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(V1, label = TRUE, label.size = 3)
p1 + p2

SpatialDimPlot(V1, cells.highlight = CellsByIdentities(object = V1, idents = c(2, 1, 4, 3,
                                                                                     5, 8)), facet.highlight = TRUE, ncol = 3)
SpatialDimPlot(V1, interactive = TRUE)

de_markers <- FindMarkers(V1, ident.1 = 5, ident.2 = 6)
SpatialFeaturePlot(object = V1, features = rownames(de_markers)[1:3], alpha = c(0.1, 1), ncol = 3)

V1 <- FindSpatiallyVariableFeatures(V1, assay = "SCT", features = VariableFeatures(V1)[1:1000],
                                       selection.method = "markvariogram")

top.features <- head(SpatiallyVariableFeatures(V1, selection.method = "markvariogram"), 6)
SpatialFeaturePlot(V1, features = top.features, ncol = 3, alpha = c(0.1, 1))

cortex <- subset(V1, idents = c(1, 2, 3, 4, 6, 7))
# now remove additional cells, use SpatialDimPlots to visualize what to remove
# SpatialDimPlot(cortex,cells.highlight = WhichCells(cortex, expression = image_imagerow > 400
# | image_imagecol < 150))
cortex <- subset(cortex, anterior1_imagerow > 400 | anterior1_imagecol < 150, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 275 & anterior1_imagecol > 370, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 250 & anterior1_imagecol > 440, invert = TRUE)

p1 <- SpatialDimPlot(cortex, crop = TRUE, label = TRUE)
p2 <- SpatialDimPlot(cortex, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 3)
p1 + p2






SpatialFeaturePlot(brain, features = "Ttr", interactive = TRUE)

LinkedDimPlot(brain)










library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)

MTG <- Load10X_Spatial('/home/wangpp/DATA/Brain/spatial/R21009538_out/outs') 

plot1 <- VlnPlot(MTG, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(MTG, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

MTG <- SCTransform(MTG, assay = "Spatial", verbose = FALSE)

SpatialFeaturePlot(MTG, features = c("HPCA", "TTR"))

p1 <- SpatialFeaturePlot(MTG, features = "HPCA", pt.size.factor = 1)
p2 <- SpatialFeaturePlot(MTG, features = "HPCA", alpha = c(0.1, 1))
p1 + p2


MTG <- RunPCA(MTG, assay = "SCT", verbose = FALSE)
MTG <- FindNeighbors(MTG, reduction = "pca", dims = 1:30)
MTG <- FindClusters(MTG, verbose = FALSE)
MTG <- RunUMAP(MTG, reduction = "pca", dims = 1:30)

p1 <- DimPlot(MTG, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(MTG, label = TRUE, label.size = 3)
p1 + p2

SpatialDimPlot(MTG, cells.highlight = CellsByIdentities(object = MTG, idents = c(2, 1, 4, 3,
                                                                                 5, 8)), facet.highlight = TRUE, ncol = 3)
SpatialDimPlot(MTG, interactive = TRUE)

de_markers <- FindMarkers(MTG, ident.1 = 5, ident.2 = 6)
SpatialFeaturePlot(object = MTG, features = rownames(de_markers)[1:3], alpha = c(0.1, 1), ncol = 3)

MTG <- FindSpatiallyVariableFeatures(MTG, assay = "SCT", features = VariableFeatures(MTG)[1:1000],
                                     selection.method = "markvariogram")

top.features <- head(SpatiallyVariableFeatures(MTG, selection.method = "markvariogram"), 6)
SpatialFeaturePlot(MTG, features = top.features, ncol = 3, alpha = c(0.1, 1))

cortex <- subset(MTG, idents = c(1, 2, 3, 4, 6, 7))
# now remove additional cells, use SpatialDimPlots to visualize what to remove
# SpatialDimPlot(cortex,cells.highlight = WhichCells(cortex, expression = image_imagerow > 400
# | image_imagecol < 150))
cortex <- subset(cortex, anterior1_imagerow > 400 | anterior1_imagecol < 150, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 275 & anterior1_imagecol > 370, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 250 & anterior1_imagecol > 440, invert = TRUE)

p1 <- SpatialDimPlot(cortex, crop = TRUE, label = TRUE)
p2 <- SpatialDimPlot(cortex, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 3)
p1 + p2






SpatialFeaturePlot(brain, features = "Ttr", interactive = TRUE)

LinkedDimPlot(brain)










OFC_tile = "-A1_"
M1_tile = "-A4_"
V1_title = "-A11_"
MTG_title = "-A12_"

human_brain_exp = readRDS("counts.v2.rds")
human_meta = read.table("/data/share/Brain/Brain_subClusteringID_seurat_merge_v6.txt",header = T,sep='\t')

###
ofc <- human_meta[grep(OFC_tile,human_meta$cells),]
write.table(ofc,file = "OFC_cell_meta.txt",sep = "\t",row.names = F,quote = F)

ofc_count <- human_brain_exp[,ofc$cells]
write.csv(ofc_count, file='OFC_count.csv')

###
m1 <- human_meta[grep(M1_tile,human_meta$cells),]
write.table(m1,file = "M1_cell_meta.txt",sep = "\t",row.names = F,quote = F)

m1_count <- human_brain_exp[,m1$cells]
write.csv(m1_count, file='M1_count.csv')


###
v1 <- human_meta[grep(V1_title,human_meta$cells),]
write.table(v1,file = "V1_cell_meta.txt",sep = "\t",row.names = F,quote = F)

v1_count <- human_brain_exp[,v1$cells]
write.csv(v1_count, file='V1_count.csv')


###
mtg <- human_meta[grep(MTG_tile,human_meta$cells),]
write.table(mtg,file = "MTG_cell_meta.txt",sep = "\t",row.names = F,quote = F)

mtg_count <- human_brain_exp[,mtg$cells]
write.csv(mtg_count, file='MTG_count.csv')


















