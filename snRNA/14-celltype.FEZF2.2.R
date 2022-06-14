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
#clusters1 <- c("68","70")
#clusters2 <- c("25","30","29","31","58","425","40","36","39","51","42","48","49","72","93","125","111","104","109",
#               "98","115","123","124","88","65","77","86","89")
#cells1 <- colnames(seu0)[seu0$hicat_cluster_merge %in% clusters1]
#cells2 <- colnames(seu0)[seu0$hicat_cluster_merge %in% clusters2]
#
#XXX <- FindMarkers(seu0,
#                   ident.1 = clusters1,
#                   ident.2 = clusters2)
#saveRDS(XXX,'Exc_FEZF2.2_markers.rds')
#
#top5 <- XXX[XXX$avg_log2FC>0,][1:50,]
#genes <- rownames(top5)
#
##FEZF2.2_up-regulated.pdf 30x20
#DotPlot(seu0, features=unique(genes)) + RotatedAxis()
#
#################################################################################

#
#clusters <- c("68","70")
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
#saveRDS(seu,'seu.harmony.FEZF2.2.rds')

seu <- readRDS('seu.harmony.FEZF2.2.rds')

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

saveRDS(markers, 'seu.harmony.FEZF2.2.markers.rds')
saveRDS(seu,'seu.harmony.FEZF2.2.rds')



################################################
markers <- readRDS('seu.harmony.FEZF2.2.markers.rds')
seu <- readRDS('seu.harmony.FEZF2.2.rds')

DimPlot(seu,label=T)

library(dplyr)
top5 <- top_n(group_by(markers,cluster), 5, avg_log2FC)
genes <- top5$gene

DotPlot(seu, features=unique(genes)) + RotatedAxis()

markers2 <- markers[,c("cluster", "gene","avg_log2FC", "pct.1", "pct.2","p_val",  "p_val_adj")]
write.csv(markers2, 'FEZF2.2.markers.csv',row.names = FALSE, quote = F)


############### K means #######################
seu <- readRDS('seu.harmony.FEZF2.2.rds')

seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:20)
seu <- FindClusters(seu, resolution = 1)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:20)

DimPlot(seu, group.by = 'hicat_cluster')

DimPlot(seu, group.by = 'region')
DimPlot(seu, group.by = 'individual')


genes <- c('FEZF2','CSN1S1','ASGR2','LINC01107')

FeaturePlot(seu,
            ncol = 2,
            features=genes)

data <- as.matrix(seu@assays$RNA@data)
km_result <- kmeans(t(data), 4)
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

DimPlot(seu,group.by = 'seurat_clusters')


seu <- readRDS('seu.harmony.FEZF2.2.rds')
DimPlot(seu,group.by = 'seurat_clusters')

genes0 <- c('VAT1L', 'AC016074.2', 'GRIK2', 'ADRA1A', 'ESRRG', 'CHST8', 'GRM8', 'RNF152', 'C14orf39', 'MIR34AHG', 'SULF2', 'ERG', 'PDE1C')
genes1 <- c('GRIK2', 'ADRA1A', 'ESRRG', 'VAT1L', 'RAB37', 'MAST4', 'CDH20', 'CABP1')
genes2 <- c('AC005100.1', 'AC016074.2', 'LINC00922', 'LINC02424', 'OR51E2', 'PCDH15', 'LINC00943', 'ANGPT1', 'ADRA1A', 'ADAMTSL3', 'ROBO2', 'AL589693.1', 'LINC00861')
genes3 <- c('GRIK2', 'LINC00861', 'ADRA1A', 'AC016074.2', 'VAT1L', 'BCL11B', 'KCNN2', 'SLC66A1L', 'GRM8', 'AC015522.1', 'MYO16', 'TENM3-AS1', 'ESRRG')

genes <- c(genes0,genes1,genes2,genes3)

#FEZF2.2_marker.pdf
DotPlot(seu, features=unique(genes), group.by = 'seurat_clusters') + RotatedAxis()


pdf('FEZF2.2_cluster3_marker.pdf',width = 6,height = 4)
DotPlot(seu, features=unique(genes3), group.by = 'seurat_clusters') + RotatedAxis()
dev.off()

#### hicat ##################
seu <- readRDS('seu.harmony.FEZF2.2.rds')
result <- readRDS('seu.hicat.FEZF2.2.rds')
seu$re_hicat <- result$cl.result$cl[colnames(seu)]
seu$re_hicat <- factor(seu$re_hicat)
DimPlot(seu,group.by = 're_hicat')

genes <- c('FEZF2','CSN1S1','ASGR2','LINC01107')
DotPlot(seu,features = genes,group.by = 're_hicat') + RotatedAxis()

seu@active.ident <- seu$re_hicat
markers <- FindAllMarkers(seu)
top5 <- top_n(group_by(markers,cluster), 5, avg_log2FC)
genes <- top5$gene
DotPlot(seu, features=unique(genes)) + RotatedAxis()





############### cluster merge ####################################################
seu <- readRDS('seu.harmony.FEZF2.2.rds')

seu$seurat_clusters_merge <- as.vector(seu$seurat_clusters)
seu$seurat_clusters_merge[seu$seurat_clusters %in% c('3','0')] <- '0,3'

seu$seurat_clusters_merge <- factor(seu$seurat_clusters_merge)

DimPlot(seu, group.by = 'seurat_clusters', label = T)
DimPlot(seu, group.by = 'seurat_clusters_merge', label = T)

seu@active.ident <- seu$seurat_clusters
markers <- FindAllMarkers(seu)
saveRDS(markers, 'seu.harmony.FEZF2.2.markers.rds')

seu@active.ident <- seu$seurat_clusters_merge
markers <- FindAllMarkers(seu)
saveRDS(markers, 'seu.harmony.FEZF2.2.merge.markers.rds')


markers1 <- readRDS('seu.harmony.FEZF2.2.markers.rds')
markers2 <- readRDS('seu.harmony.FEZF2.2.merge.markers.rds')

library(dplyr)
top1 <- top_n(group_by(markers1,cluster), 5, avg_log2FC)
genes1 <- top1$gene

top2 <- top_n(group_by(markers2,cluster), 5, avg_log2FC)
genes2 <- top2$gene

DotPlot(seu, features=unique(genes1), group.by = 'seurat_clusters') + RotatedAxis()
DotPlot(seu, features=unique(genes2), group.by = 'seurat_clusters_merge') + RotatedAxis()



