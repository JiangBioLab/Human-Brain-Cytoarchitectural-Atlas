library(scrattch.hicat)
library(dendextend)
library(Seurat)

#seu0 <- readRDS('seu.harmony.anno.rds')
#result <- readRDS('seu.hicat.de80.niter100.rds')
#dend.result <- readRDS('seu.hicat.de80.niter100.dend.rds')
#
#c <- c("25","30","29","31","58","425","68","70","40","36","39","51","42","48","49","72","93","125","111","104","109",
#"98","115","123","124","88","65","77","86","89")
#
####################### XXX markers ##############################################
#clusters1 <- c("40","36","39","51","42","48","49")
#clusters2 <- c("25","30","29","31","58","425","68","70",,"72","93","125","111","104","109",
#               "98","115","123","124","88","65","77","86","89")
#cells1 <- colnames(seu0)[seu0$hicat_cluster_merge %in% clusters1]
#cells2 <- colnames(seu0)[seu0$hicat_cluster_merge %in% clusters2]
#
#XXX <- FindMarkers(seu0,
#                   ident.1 = clusters1,
#                   ident.2 = clusters2)
#saveRDS(XXX,'Exc_FEZF2.1_markers.rds')
#
XXX<-readRDS('Exc_FEZF2.1_markers.rds')

top5 <- XXX[XXX$avg_log2FC>0,][1:50,]
genes <- rownames(top5)

#FEZF2.1_up-regulated.pdf 30x20
DotPlot(seu0, features=unique(genes)) + RotatedAxis()
#
#################################################################################

#
#clusters <- c("40","36","39","51","42","48","49")
#
#counts <- seu0@assays$RNA@counts[,seu0$hicat_cluster_merge %in% clusters]
#
#seu <- CreateSeuratObject(
#  counts,
#  min.cells = 0,
#  min.features = 0
#)
## 33014 36013
#
#meta <- data.frame(cells = colnames(seu0),
#                   sample = seu0$sample,
#                   individual = seu0$individual,
#                   region = seu0$region,
#                   hicat_cluster = seu0$hicat_cluster,
#                   hicat_cluster_merge = seu0$hicat_cluster_merge)
#rownames(meta) <- meta$cells
#
#seu$sample <- meta[colnames(seu),'sample']
#seu$individual <- meta[colnames(seu),'individual']
#seu$region <- meta[colnames(seu),'region']
#seu$hicat_cluster <- meta[colnames(seu),'hicat_cluster']
#seu$hicat_cluster_merge <- meta[colnames(seu),'hicat_cluster_merge']
#
#
#seu$hicat_cluster_merge <- factor(seu$hicat_cluster_merge, levels = clusters)
#
#saveRDS(seu,'seu.harmony.FEZF2.1.rds')

seu <- readRDS('seu.harmony.FEZF2.1.rds')

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu <- NormalizeData(object = seu, normalization.method = "LogNormalize", scale.factor = 1e4)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu, features = rownames(seu))
seu <- RunPCA(seu, features = VariableFeatures(object = seu))
seu <- RunHarmony(seu, "sample")
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:20)
seu <- FindClusters(seu, resolution = 1)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:40)

seu@active.ident <- seu$hicat_cluster_merge
markers <- FindAllMarkers(seu)

saveRDS(markers, 'seu.harmony.FEZF2.1.markers.rds')
saveRDS(seu,'seu.harmony.FEZF2.1.rds')



################################################
markers <- readRDS('seu.harmony.FEZF2.1.markers.rds')
seu <- readRDS('seu.harmony.FEZF2.1.rds')

DimPlot(seu,label=T)

library(dplyr)
top5 <- top_n(group_by(markers,cluster), 5, avg_log2FC)
genes <- top5$gene

DotPlot(seu, features=unique(genes)) + RotatedAxis()

markers2 <- markers[,c("cluster", "gene","avg_log2FC", "pct.1", "pct.2","p_val",  "p_val_adj")]
write.csv(markers2, 'FEZF2.1.markers.csv',row.names = FALSE, quote = F)




############### K means #######################
seu0 <- readRDS('seu.harmony.FEZF2.1.rds')
clusters <- c("51","42","48","49")

counts <- seu0@assays$RNA@counts[,seu0$hicat_cluster_merge %in% clusters]

seu <- CreateSeuratObject(
  counts,
  min.cells = 0,
  min.features = 0
)

meta <- data.frame(cells = colnames(seu0),
                   sample = seu0$sample,
                   individual = seu0$individual,
                   region = seu0$region,
                   hicat_cluster = seu0$hicat_cluster,
                   hicat_cluster_merge = seu0$hicat_cluster_merge)
rownames(meta) <- meta$cells

seu$sample <- meta[colnames(seu),'sample']
seu$individual <- meta[colnames(seu),'individual']
seu$region <- meta[colnames(seu),'region']
seu$hicat_cluster <- meta[colnames(seu),'hicat_cluster']
seu$hicat_cluster_merge <- meta[colnames(seu),'hicat_cluster_merge']

seu$hicat_cluster_merge <- factor(seu$hicat_cluster_merge, levels = clusters)

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu <- NormalizeData(object = seu, normalization.method = "LogNormalize", scale.factor = 1e4)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu, features = rownames(seu))
seu <- RunPCA(seu, features = VariableFeatures(object = seu))
seu <- RunHarmony(seu, "sample")
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:20)
seu <- FindClusters(seu, resolution = 1)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:20)

seu@active.ident <- seu$hicat_cluster_merge
markers <- FindAllMarkers(seu)

saveRDS(markers, 'seu.harmony.FEZF2.4.markers.rds')
saveRDS(seu,'seu.harmony.FEZF2.4.rds')



seu<-readRDS('seu.harmony.FEZF2.4.rds')

DimPlot(seu, group.by = 'hicat_cluster')

DimPlot(seu, group.by = 'region')
DimPlot(seu, group.by = 'individual')


genes <- c('FEZF2','PDYN','FFAR4','PROKR2','CFTR','KLK7','POGK')

FeaturePlot(seu,
            ncol = 4,
            features=genes)

data <- as.matrix(seu@assays$RNA@data)
km_result <- kmeans(t(data), 6)
table(km_result$cluster)

seu$kmeans <- km_result$cluster[colnames(seu)]


DimPlot(seu, group.by = 'kmeans')

DimPlot(seu, group.by = 'seurat_clusters')

DotPlot(seu,features = genes,group.by = 'kmeans') + RotatedAxis()
DotPlot(seu,features = genes,group.by = 'seurat_clusters') + RotatedAxis()


seu$kmeans <- factor(seu$kmeans)
seu@active.ident <- seu$kmeans
markers <- FindAllMarkers(seu)
top5 <- top_n(group_by(markers,cluster), 5, avg_log2FC)
genes <- top5$gene
DotPlot(seu, features=unique(genes)) + RotatedAxis()


seu@active.ident <- seu$seurat_clusters
markers <- FindAllMarkers(seu)
top5 <- top_n(group_by(markers,cluster), 5, avg_log2FC)
genes <- top5$gene
DotPlot(seu, features=unique(genes)) + RotatedAxis()



############### K means 2 #######################
seu0 <- readRDS('seu.harmony.FEZF2.1.rds')
clusters <- c("40","36","39")

counts <- seu0@assays$RNA@counts[,seu0$hicat_cluster_merge %in% clusters]

seu <- CreateSeuratObject(
  counts,
  min.cells = 0,
  min.features = 0
)

meta <- data.frame(cells = colnames(seu0),
                   sample = seu0$sample,
                   individual = seu0$individual,
                   region = seu0$region,
                   hicat_cluster = seu0$hicat_cluster,
                   hicat_cluster_merge = seu0$hicat_cluster_merge)
rownames(meta) <- meta$cells

seu$sample <- meta[colnames(seu),'sample']
seu$individual <- meta[colnames(seu),'individual']
seu$region <- meta[colnames(seu),'region']
seu$hicat_cluster <- meta[colnames(seu),'hicat_cluster']
seu$hicat_cluster_merge <- meta[colnames(seu),'hicat_cluster_merge']

seu$hicat_cluster_merge <- factor(seu$hicat_cluster_merge, levels = clusters)

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu <- NormalizeData(object = seu, normalization.method = "LogNormalize", scale.factor = 1e4)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu, features = rownames(seu))
seu <- RunPCA(seu, features = VariableFeatures(object = seu))
seu <- RunHarmony(seu, "sample")
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:20)
seu <- FindClusters(seu, resolution = 1)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:20)

seu@active.ident <- seu$hicat_cluster_merge
markers <- FindAllMarkers(seu)

saveRDS(markers, 'seu.harmony.FEZF2.5.markers.rds')
saveRDS(seu,'seu.harmony.FEZF2.5.rds')



seu<-readRDS('seu.harmony.FEZF2.5.rds')

DimPlot(seu, group.by = 'hicat_cluster')

DimPlot(seu, group.by = 'region')
DimPlot(seu, group.by = 'individual')


genes <- c('FEZF2','THEMIS','OR1L8','LINC01116','RGPD6','CFTR','C9orf135-AS1','FILIP1L','SH2D1B','SMYD1')

FeaturePlot(seu,
            ncol = 3,
            features=genes)

data <- as.matrix(seu@assays$RNA@data)
km_result <- kmeans(t(data), 8)
table(km_result$cluster)

seu$kmeans <- km_result$cluster[colnames(seu)]


DimPlot(seu, group.by = 'kmeans')

DimPlot(seu, group.by = 'seurat_clusters')

DotPlot(seu,features = genes,group.by = 'kmeans') + RotatedAxis()
DotPlot(seu,features = genes,group.by = 'seurat_clusters') + RotatedAxis()


seu$kmeans <- factor(seu$kmeans)
seu@active.ident <- seu$kmeans
markers <- FindAllMarkers(seu)
top5 <- top_n(group_by(markers,cluster), 5, avg_log2FC)
genes <- top5$gene
DotPlot(seu, features=unique(genes)) + RotatedAxis()


seu@active.ident <- seu$seurat_clusters
markers <- FindAllMarkers(seu)
top5 <- top_n(group_by(markers,cluster), 5, avg_log2FC)
genes <- top5$gene
DotPlot(seu, features=unique(genes)) + RotatedAxis()


seu <- readRDS('seu.harmony.FEZF2.4.rds')
DimPlot(seu,group.by = 'seurat_clusters')

genes0 <- c('MCTP2', 'SEMA3D', 'AC022126.1', 'DLC1', 'NFIA', 'AC068587.4', 'ZFHX3', 'ASTN2', 'SEZ6L', 'MDFIC', 'KIAA1217', 'LRMDA', 'LINC01479', 'OR3A2', 'AC092114.1')
genes1 <- c('AL136456.1', 'PCOLCE2', 'PIEZO2', 'KIAA1217', 'ROBO2', 'KCNK13', 'ADD3-AS1', 'ZFHX3', 'PCSK5', 'ZBTB7C', 'CCN2', 'PCDH11X', 'KCNK2', 'MARCH1')
genes2 <- c('NREP', 'MDFIC', 'AL138927.1', 'NPFFR2', 'MARCH1', 'LINC02718', 'AC022126.1', 'KIAA1217', 'ZFHX3', 'XIST', 'MIAT', 'ASTN2', 'LINC01060', 'LINC02008', 'CELSR1')
genes3 <- c('ZFHX3', 'DLC1', 'KIAA1217', 'ASTN2', 'HS3ST4')
genes4 <- c('ITPR2', 'AC023503.1', 'MARCH1', 'KIAA1217', 'MT-ND3', 'PCSK5', 'NFIA', 'ZFHX3', 'NPFFR2', 'RASEF', 'SEMA3D', 'MT-ND1', 'MCTP2', 'MSC-AS1')
genes5 <- c('NREP-AS1', 'NREP', 'AC008696.2', 'MARCH1', 'ZFHX3', 'MDFIC', 'AC091885.2', 'AC016598.2', 'CDH9', 'PCDH9-AS2', 'NPFFR2', 'DPP10-AS3', 'LINC02718')
genes6 <- c('HS3ST4', 'UTRN', 'HCRTR2', 'ATP8B4', 'ADAMTSL3', 'RNF220', 'LHFPL3', 'SCUBE1', 'PCDH11X', 'KIAA1217', 'ERG', 'SORCS3', 'PCSK5', 'INPP4B', 'MGAT4C')
genes7 <- c('ROBO2', 'SCUBE1', 'AC092131.1', 'SGCZ', 'AL391117.1', 'SDK2', 'MGAT4C', 'SLC15A5', 'VCAN', 'PIK3C2G', 'THSD7B', 'HS3ST5', 'ADD3-AS1', 'CCN2')
genes8 <- c('AC068587.4', 'ZFHX3', 'WIPF1', 'LINC01568', 'DLC1', 'SEMA3D', 'NFIA', 'MRC2', 'AC024145.1', 'CDC14A', 'CTXND1', 'MDFIC', 'SCML4')
genes9 <- c('MDFIC', 'MCTP2', 'KIAA1217', 'AC022126.1', 'MARCH1', 'FAM198B-AS1', 'ASTN2', 'NPFFR2', 'BCAS1', 'AC008571.2')

genes <- c(genes0,genes1,genes2,genes3,genes4,genes5,genes6,genes7,genes8,genes9)
#FEZF2.4_marker.pdf
DotPlot(seu, features=unique(genes), group.by = 'seurat_clusters') + RotatedAxis()


pdf('FEZF2.4_cluster9_marker.pdf',width = 6,height = 4)
DotPlot(seu, features=unique(genes9), group.by = 'seurat_clusters') + RotatedAxis()
dev.off()


seu <- readRDS('seu.harmony.FEZF2.5.rds')
DimPlot(seu,group.by = 'seurat_clusters')

genes0 <- c('TRPM3', 'SULF1', 'RYR3', 'ASIC2', 'TSPAN18', 'LINC02082', 'CNTNAP3', 'SRGAP1', 'TMEFF2', 'SEMA3E', 'ADAMTSL1', 'RALYL', 'DLC1', 'LRP1B')
genes1 <- c('ANKRD18A', 'ADAMTSL1', 'SEMA3E', 'SLC24A2', 'KIAA1217', 'HS3ST4', 'TRPM3', 'ASIC2', 'RYR3', 'DLC1', 'CCDC88A', 'AC092683.1', 'LINC00342')
genes2 <- c('FRMPD4', 'XIST', 'LINC02718', 'SEMA3E', 'SGCG', 'DLC1', 'KIAA1217', 'TRPM3', 'HS3ST4', 'LINC02082', 'ADAMTSL1', 'DPP10-AS3', 'FOXP2', 'RALYL', 'UNC80')
genes3 <- c('AC046195.2', 'AC046195.1', 'SEMA5A', 'DDR2', 'BAALC-AS1', 'TMC3', 'IL4R', 'AC008696.2', 'LUZP2', 'IQCJ-SCHIP1', 'SLC8A1-AS1', 'GNG12-AS1', 'SULF1')
genes4 <- c('LGR6', 'SEMA5A', 'SULF1', 'KIAA1217', 'CSMD1', 'ROBO2', 'RYR3', 'TSPAN18', 'LINC02082', 'HS3ST2', 'FRMPD4', 'AC090578.1')
genes5 <- c('ADAMTSL1', 'ASIC2', 'ROBO2', 'CSMD1', 'DLC1', 'AC021613.1', 'SEMA5A', 'TRPM3', 'PDZRN4', 'RBM20', 'MCTP1', 'MEGF11', 'SULF1', 'AC019211.1')
genes6 <- c('S100Z', 'AC090578.1', 'SGCZ', 'KIAA1217', 'GRM8', 'LINC01435', 'RALYL', 'PALMD', 'CAV1', 'TRPM3', 'ARHGAP6', 'CLSTN2', 'ABO', 'SRGAP1')
genes7 <- c('LINC01194', 'AC046195.1', 'S100Z', 'SULF1', 'TRPM3', 'SRGAP1', 'RALYL', 'ASIC2', 'AC046195.2', 'PCSK5', 'LINC02082', 'SEMA3E', 'CNTNAP3', 'SGK1')
genes8 <- c('ADAMTSL1', 'SLC8A1', 'TEAD1', 'IL4R', 'MCC', 'SEMA3E', 'SLIT1', 'DLC1', 'ARHGAP10', 'KCNH7', 'PTPRG', 'AC068587.4', 'RHOJ', 'LUZP2')
genes9 <- c('MCC', 'ANTXR1', 'ADAMTSL1', 'TMEFF2', 'VWA2', 'AC092969.1', 'PCDH11Y', 'RYR3', 'NRXN3', 'DLC1', 'TRPM3')
genes10 <- c('SYT6', 'VWA2', 'ITGA11', 'SEMA5A', 'DDR2', 'DLC1', 'SERPINE2', 'HS3ST4', 'IQCJ-SCHIP1', 'KIAA1217')
genes11 <- c('ADAMTSL1', 'TRPM3')
genes12 <- c('FRMPD4', 'MOG', 'CNDP1', 'LINC00639', 'LINC01170', 'FA2H', 'PDE8A', 'MYRF', 'ST18', 'AC079352.1')

genes <- c(genes0,genes1,genes2,genes3,genes4,genes5,genes6,genes7,genes8,genes9,genes10,genes11,genes12)

#FEZF2.5_marker.pdf
DotPlot(seu, features=unique(genes), group.by = 'seurat_clusters') + RotatedAxis()


pdf('FEZF2.5_cluster12_marker.pdf',width = 8,height = 6)
DotPlot(seu, features=unique(genes12), group.by = 'seurat_clusters') + RotatedAxis()
dev.off()








############### cluster merge ####################################################
seu <- readRDS('seu.harmony.FEZF2.4.rds')

seu$seurat_clusters_merge <- as.vector(seu$seurat_clusters)
seu$seurat_clusters_merge[seu$seurat_clusters %in% c('0','9')] <- '0,9'
seu$seurat_clusters_merge[seu$seurat_clusters %in% c('1','2','5')] <- '1,2,5'

seu$seurat_clusters_merge <- factor(seu$seurat_clusters_merge)

DimPlot(seu, group.by = 'seurat_clusters', label = T)
DimPlot(seu, group.by = 'seurat_clusters_merge', label = T)

seu@active.ident <- seu$seurat_clusters
markers <- FindAllMarkers(seu)
saveRDS(markers, 'seu.harmony.FEZF2.4.markers.rds')

seu@active.ident <- seu$seurat_clusters_merge
markers <- FindAllMarkers(seu)
saveRDS(markers, 'seu.harmony.FEZF2.4.merge.markers.rds')


markers1 <- readRDS('seu.harmony.FEZF2.4.markers.rds')
markers2 <- readRDS('seu.harmony.FEZF2.4.merge.markers.rds')

library(dplyr)
top1 <- top_n(group_by(markers1,cluster), 5, avg_log2FC)
genes1 <- top1$gene

top2 <- top_n(group_by(markers2,cluster), 5, avg_log2FC)
genes2 <- top2$gene

DotPlot(seu, features=unique(genes1), group.by = 'seurat_clusters') + RotatedAxis()
DotPlot(seu, features=unique(genes2), group.by = 'seurat_clusters_merge') + RotatedAxis()





############### cluster merge ####################################################
seu <- readRDS('seu.harmony.FEZF2.5.rds')

seu$seurat_clusters_merge <- as.vector(seu$seurat_clusters)
seu$seurat_clusters_merge[seu$seurat_clusters %in% c('0','4','5','11')] <- '0,4,5,11'
seu$seurat_clusters_merge[seu$seurat_clusters %in% c('2','9','10')] <- '2,9,10'
seu$seurat_clusters_merge[seu$seurat_clusters %in% c('6','7')] <- '6,7'

seu$seurat_clusters_merge <- factor(seu$seurat_clusters_merge)

DimPlot(seu, group.by = 'seurat_clusters', label = T)
DimPlot(seu, group.by = 'seurat_clusters_merge', label = T)

seu@active.ident <- seu$seurat_clusters
markers <- FindAllMarkers(seu)
saveRDS(markers, 'seu.harmony.FEZF2.5.markers.rds')

seu@active.ident <- seu$seurat_clusters_merge
markers <- FindAllMarkers(seu)
saveRDS(markers, 'seu.harmony.FEZF2.5.merge.markers.rds')


markers1 <- readRDS('seu.harmony.FEZF2.5.markers.rds')
markers2 <- readRDS('seu.harmony.FEZF2.5.merge.markers.rds')

library(dplyr)
top1 <- top_n(group_by(markers1,cluster), 5, avg_log2FC)
genes1 <- top1$gene

top2 <- top_n(group_by(markers2,cluster), 5, avg_log2FC)
genes2 <- top2$gene

DotPlot(seu, features=unique(genes1), group.by = 'seurat_clusters') + RotatedAxis()
DotPlot(seu, features=unique(genes2), group.by = 'seurat_clusters_merge') + RotatedAxis()

markers2 <- markers[,c("cluster", "gene","avg_log2FC", "pct.1", "pct.2","p_val",  "p_val_adj")]
write.csv(markers2, 'FEZF2.5.markers.csv',row.names = FALSE, quote = F)




