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
clusters1 <- c("58","425")
clusters2 <- c("25","30","29","31","68","70","40","36","39","51","42","48","49","72","93","125","111","104","109",
               "98","115","123","124","88","65","77","86","89")
cells1 <- colnames(seu0)[seu0$hicat_cluster_merge %in% clusters1]
cells2 <- colnames(seu0)[seu0$hicat_cluster_merge %in% clusters2]

XXX <- FindMarkers(seu0,
                   ident.1 = clusters1,
                   ident.2 = clusters2)
saveRDS(XXX,'Exc_THEMIS_markers.rds')


XXX<-readRDS('Exc_THEMIS_markers.rds')

top5 <- XXX[XXX$avg_log2FC>0,][1:50,]
genes <- rownames(top5)

#THEMIS_up-regulated.pdf 30x20
DotPlot(seu0, features=unique(genes)) + RotatedAxis()

#################################################################################


#clusters <- c("58","425")
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
#saveRDS(seu,'seu.harmony.THEMIS.rds')

seu <- readRDS('seu.harmony.THEMIS.rds')

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu <- NormalizeData(object = seu, normalization.method = "LogNormalize", scale.factor = 1e4)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu, features = rownames(seu))
seu <- RunPCA(seu, features = VariableFeatures(object = seu))
seu <- RunHarmony(seu, "sample")
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:10)
seu <- FindClusters(seu, resolution = 1)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:20)

seu@active.ident <- seu$hicat_cluster_merge
markers <- FindAllMarkers(seu)

saveRDS(markers, 'seu.harmony.THEMIS.markers.rds')
saveRDS(seu,'seu.harmony.THEMIS.rds')



################################################
markers <- readRDS('seu.harmony.THEMIS.markers.rds')
seu <- readRDS('seu.harmony.THEMIS.rds')

DimPlot(seu,label=T)

library(dplyr)
top5 <- top_n(group_by(markers,cluster), 5, avg_log2FC)
genes <- top5$gene

DotPlot(seu, features=unique(genes)) + RotatedAxis()

markers2 <- markers[,c("cluster", "gene","avg_log2FC", "pct.1", "pct.2","p_val",  "p_val_adj")]
write.csv(markers2, 'THEMIS.markers.csv',row.names = FALSE, quote = F)




############### K means #######################
seu <- readRDS('seu.harmony.THEMIS.rds')

seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:20)
seu <- FindClusters(seu, resolution = 1)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:20)

DimPlot(seu, group.by = 'hicat_cluster')

DimPlot(seu, group.by = 'region')
DimPlot(seu, group.by = 'individual')

DimPlot(seu, group.by = 'seurat_clusters')

data <- as.matrix(seu@assays$RNA@data)
km_result <- kmeans(t(data), 4)
table(km_result$cluster)

seu$kmeans <- km_result$cluster[colnames(seu)]


DimPlot(seu, group.by = 'kmeans')

genes <- c('THEMIS','LINC00343','SLN','SNTG2','SMYD1')

FeaturePlot(seu,
            ncol = 3,
            features=genes)

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




seu <- readRDS('seu.harmony.THEMIS.rds')
DimPlot(seu,group.by = 'seurat_clusters')

genes0 <- c()
genes1 <- c()
genes2 <- c()
genes3 <- c()
genes4 <- c()
genes5 <- c()
genes6 <- c()
genes7 <- c()
genes8 <- c()

genes <- c(genes0,genes1,genes2,genes3,genes4,genes5,genes6,genes7,genes8)

#THEMIS_marker.pdf
DotPlot(seu, features=unique(genes), group.by = 'seurat_clusters') + RotatedAxis()


pdf('THEMIS_cluster8_marker.pdf',width = 8,height = 6)
DotPlot(seu, features=unique(genes8), group.by = 'seurat_clusters') + RotatedAxis()
dev.off()






############### cluster merge ####################################################
seu <- readRDS('seu.harmony.THEMIS.rds')

seu$seurat_clusters_merge <- as.vector(seu$seurat_clusters)
seu$seurat_clusters_merge[seu$seurat_clusters %in% c('1','5')] <- '1,5'

seu$seurat_clusters_merge <- factor(seu$seurat_clusters_merge)

DimPlot(seu, group.by = 'seurat_clusters', label = T)
DimPlot(seu, group.by = 'seurat_clusters_merge', label = T)

seu@active.ident <- seu$seurat_clusters
markers <- FindAllMarkers(seu)
saveRDS(markers, 'seu.harmony.THEMIS.markers.rds')

seu@active.ident <- seu$seurat_clusters_merge
markers <- FindAllMarkers(seu)
saveRDS(markers, 'seu.harmony.THEMIS.merge.markers.rds')


markers1 <- readRDS('seu.harmony.THEMIS.markers.rds')
markers2 <- readRDS('seu.harmony.THEMIS.merge.markers.rds')

library(dplyr)
top1 <- top_n(group_by(markers1,cluster), 5, avg_log2FC)
genes1 <- top1$gene

top2 <- top_n(group_by(markers2,cluster), 5, avg_log2FC)
genes2 <- top2$gene

DotPlot(seu, features=unique(genes1), group.by = 'seurat_clusters') + RotatedAxis()
DotPlot(seu, features=unique(c('CAR3',genes2)), group.by = 'seurat_clusters_merge') + RotatedAxis()


markers <- readRDS('seu.harmony.THEMIS.merge.markers.rds')
markers2 <- markers[,c("cluster", "gene","avg_log2FC", "pct.1", "pct.2","p_val",  "p_val_adj")]
write.csv(markers2, 'THEMIS.markers.csv',row.names = FALSE, quote = F)




