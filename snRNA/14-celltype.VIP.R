library(scrattch.hicat)
library(dendextend)
library(Seurat)

seu0 <- readRDS('seu.harmony.anno.rds')
result <- readRDS('seu.hicat.de80.niter100.rds')
dend.result <- readRDS('seu.hicat.de80.niter100.dend.rds')




###################### XXX markers ##############################################
clusters1 <- c("185","180","175","190","173","172")
clusters2 <- c("226","192","199","198","202","201","204_211","206",
               "146_157","144_156","148","160_162","127",
               "233","232","247","229","279_280","222","282","264","235_243","245","256","270","224","220","217", "214",
               "384","267","339","373","334","368","342","380","284","297",
               "336","268","378","375","377","383","311", "303","309","301","299","329_333","327")
cells1 <- colnames(seu0)[seu0$hicat_cluster_merge %in% clusters1]
cells2 <- colnames(seu0)[seu0$hicat_cluster_merge %in% clusters2]

XXX <- FindMarkers(seu0,
                   ident.1 = clusters1,
                   ident.2 = clusters2)
saveRDS(XXX,'Inh_XXX_markers.rds')

top5 <- XXX[XXX$avg_log2FC>0,][1:50,]
genes <- rownames(top5)

#XXX_up-regulated.pdf 30x20
DotPlot(seu0, features=unique(genes)) + RotatedAxis()

#################################################################################



#
#
#clusters <- c("233","232","247","229","279_280","222","282",
#              "264","235_243","245","256","270","224","220",
#              "217", "214")
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
#saveRDS(seu,'seu.harmony.VIP.rds')

seu <- readRDS('seu.harmony.VIP.rds')

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

saveRDS(markers, 'seu.harmony.VIP.markers.rds')
saveRDS(seu,'seu.harmony.VIP.rds')



################################################
markers <- readRDS('seu.harmony.VIP.markers.rds')
seu <- readRDS('seu.harmony.VIP.rds')

DimPlot(seu,label=T)

library(dplyr)
top5 <- top_n(group_by(markers,cluster), 5, avg_log2FC)
genes <- top5$gene

DotPlot(seu, features=unique(genes)) + RotatedAxis()

markers2 <- markers[,c("cluster", "gene","avg_log2FC", "pct.1", "pct.2","p_val",  "p_val_adj")]
write.csv(markers2, 'VIP.markers.csv',row.names = FALSE, quote = F)



