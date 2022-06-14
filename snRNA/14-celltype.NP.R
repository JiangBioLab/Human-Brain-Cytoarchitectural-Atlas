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
#clusters1 <- c("25","30","29","31")
#clusters2 <- c("58","425","68","70","40","36","39","51","42","48","49","72","93","125","111","104","109",
#               "98","115","123","124","88","65","77","86","89")
#cells1 <- colnames(seu0)[seu0$hicat_cluster_merge %in% clusters1]
#cells2 <- colnames(seu0)[seu0$hicat_cluster_merge %in% clusters2]
#
#XXX <- FindMarkers(seu0,
#                   ident.1 = clusters1,
#                   ident.2 = clusters2)
#saveRDS(XXX,'Exc_FEZF2.3_markers.rds')
#
#top5 <- XXX[XXX$avg_log2FC>0,][1:50,]
#genes <- rownames(top5)
#
##FEZF2.3_up-regulated.pdf 30x20
#DotPlot(seu0, features=unique(genes)) + RotatedAxis()
#
#################################################################################
#
#
#clusters <- c("25","30","29","31")
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
#saveRDS(seu,'seu.harmony.FEZF2.3.rds')

seu <- readRDS('seu.harmony.FEZF2.3.rds')

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

saveRDS(markers, 'seu.harmony.FEZF2.3.markers.rds')
saveRDS(seu,'seu.harmony.FEZF2.3.rds')



################################################
markers <- readRDS('seu.harmony.FEZF2.3.markers.rds')
seu <- readRDS('seu.harmony.FEZF2.3.rds')

DimPlot(seu,label=T)

library(dplyr)
top5 <- top_n(group_by(markers,cluster), 5, avg_log2FC)
genes <- top5$gene

DotPlot(seu, features=unique(genes)) + RotatedAxis()

markers2 <- markers[,c("cluster", "gene","avg_log2FC", "pct.1", "pct.2","p_val",  "p_val_adj")]
write.csv(markers2, 'FEZF2.3.markers.csv',row.names = FALSE, quote = F)




############### K means #######################
seu <- readRDS('seu.harmony.FEZF2.3.rds')

seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:20)
seu <- FindClusters(seu, resolution = 1)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:20)

DimPlot(seu, group.by = 'hicat_cluster')

DimPlot(seu, group.by = 'region')
DimPlot(seu, group.by = 'individual')

DimPlot(seu, group.by = 'seurat_clusters')

data <- as.matrix(seu@assays$RNA@data)
km_result <- kmeans(t(data), 6)
table(km_result$cluster)

seu$kmeans <- km_result$cluster[colnames(seu)]


DimPlot(seu, group.by = 'kmeans')

genes <- c('FEZF2','IFNG-AS1','RNF144A-AS1','NREP-AS1','PKD2L1')

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


seu <- readRDS('seu.harmony.FEZF2.3.rds')
DimPlot(seu,group.by = 'seurat_clusters')

genes0 <- c('ZNF385D', 'IFNG-AS1', 'NPSR1-AS1', 'KCNIP1', 'AC068587.4', 'ADAMTSL1', 'PDE9A', 'LHFPL3', 'AL355076.2', 'HTR2C', 'PDLIM1', 'LINC01811', 'SPHKAP')
genes1 <- c('IFNG-AS1', 'SPON1', 'AC106894.1', 'ADAMTS12', 'ZNF385D', 'CCDC80', 'NPSR1-AS1', 'PWRN1', 'MIR3681HG', 'KCNIP1', 'LINC02218', 'SAMD3', 'GDNF-AS1', 'NXPH2')
genes2 <- c('ITGA8', 'AFDN', 'HTR2C', 'SORCS2', 'KCNIP4')
genes3 <- c('HTR2C', 'AC009081.2', 'AC107208.1', 'FSHR', 'DAPL1', 'GHR', 'KCNIP4', 'KCNIP1', 'SPATA22', 'LHFPL3', 'GDNF-AS1', 'MYOCD')
genes4 <- c('LRRC4C', 'AC119674.1', 'GRIK1', 'SGO1-AS1', 'IFNG-AS1', 'NPSR1', 'GDNF-AS1', 'RUNX2', 'GPC6', 'ITGA8', 'TSHZ2', 'NPSR1-AS1', 'LRRC1')
genes5 <- c('IFNG-AS1', 'XIST', 'AC068587.4', 'HTR2C', 'KCNIP1', 'NPSR1-AS1', 'NXPH2', 'MYRIP', 'ZNF385D', 'LHFPL3', 'KCNIP4')
genes6 <- c('ZNF385D-AS2', 'NPSR1', 'AC068587.4', 'NPSR1-AS1', 'AC009264.1', 'CD200R1L', 'ADAMTS12', 'LINC02042', 'TOX3', 'IFNG-AS1', 'SORCS2', 'CARD11', 'KCNIP1', 'ZNF385D')
genes7 <- c('ZNF385D-AS2', 'CD200R1L', 'ZNF385D', 'IFNG-AS1', 'ROBO3', 'KCNIP4', 'HTR2C', 'AC093610.1', 'NXPH2', 'AC068587.4', 'TMEM155')
genes8 <- c('GPC5', 'GRIK1', 'CDH6', 'KANK4', 'SLC24A3', 'ASIC2', 'CARD11', 'AC009081.2', 'NPSR1-AS1', 'GALNT10', 'CD36', 'DAPL1')

genes <- c(genes0,genes1,genes2,genes3,genes4,genes5,genes6,genes7,genes8)

#FEZF2.3_marker.pdf
DotPlot(seu, features=unique(genes), group.by = 'seurat_clusters') + RotatedAxis()


pdf('FEZF2.3_cluster8_marker.pdf',width = 8,height = 6)
DotPlot(seu, features=unique(genes8), group.by = 'seurat_clusters') + RotatedAxis()
dev.off()



#### hicat ##################
seu <- readRDS('seu.harmony.FEZF2.3.rds')
result <- readRDS('seu.hicat.FEZF2.3.rds')
seu$re_hicat <- result$cl.result$cl[colnames(seu)]
seu$re_hicat <- factor(seu$re_hicat)
DimPlot(seu,group.by = 're_hicat')

genes <- c('FEZF2','IFNG-AS1','RNF144A-AS1','NREP-AS1','PKD2L1')
DotPlot(seu,features = genes,group.by = 're_hicat') + RotatedAxis()

seu@active.ident <- seu$re_hicat
markers <- FindAllMarkers(seu)
top5 <- top_n(group_by(markers,cluster), 5, avg_log2FC)
genes <- top5$gene
DotPlot(seu, features=unique(genes)) + RotatedAxis()





############### cluster merge ####################################################
seu <- readRDS('seu.harmony.FEZF2.3.rds')

seu$seurat_clusters_merge <- as.vector(seu$seurat_clusters)
seu$seurat_clusters_merge[seu$seurat_clusters %in% c('0','5','6','7')] <- '0,5,6,7'

seu$seurat_clusters_merge <- factor(seu$seurat_clusters_merge)

DimPlot(seu, group.by = 'seurat_clusters', label = T)
DimPlot(seu, group.by = 'seurat_clusters_merge', label = T)

seu@active.ident <- seu$seurat_clusters
markers <- FindAllMarkers(seu)
saveRDS(markers, 'seu.harmony.FEZF2.3.markers.rds')

seu@active.ident <- seu$seurat_clusters_merge
markers <- FindAllMarkers(seu)
saveRDS(markers, 'seu.harmony.FEZF2.3.merge.markers.rds')


markers1 <- readRDS('seu.harmony.FEZF2.3.markers.rds')
markers2 <- readRDS('seu.harmony.FEZF2.3.merge.markers.rds')

library(dplyr)
top1 <- top_n(group_by(markers1,cluster), 5, avg_log2FC)
genes1 <- top1$gene

top2 <- top_n(group_by(markers2,cluster), 5, avg_log2FC)
genes2 <- top2$gene

DotPlot(seu, features=unique(genes1), group.by = 'seurat_clusters') + RotatedAxis()
DotPlot(seu, features=unique(genes2), group.by = 'seurat_clusters_merge') + RotatedAxis()





