library(Seurat)
library(dendextend)
library(harmony)

seu0 <- readRDS('seu.harmony.anno.rds')


clusters <- c("146_157","144_156","148","160_162","127")

counts <- seu0@assays$RNA@counts[,seu0$hicat_cluster_merge %in% clusters]

seu <- CreateSeuratObject(
  counts,
  min.cells = 0,
  min.features = 0
)
# 33014 36013

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

saveRDS(seu,'seu.harmony.SV2C.rds')

seu <- readRDS('seu.harmony.SV2C.rds')

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

saveRDS(markers, 'seu.harmony.SV2C.markers.rds')
saveRDS(seu,'seu.harmony.SV2C.rds')





################################################
markers <- readRDS('seu.harmony.SV2C.markers.rds')
seu <- readRDS('seu.harmony.SV2C.rds')

DimPlot(seu,label=T)

library(dplyr)
top5 <- top_n(group_by(markers,cluster), 5, avg_log2FC)
genes <- top5$gene

DotPlot(seu, features=unique(genes)) + RotatedAxis()

markers2 <- markers[,c("cluster", "gene","avg_log2FC", "pct.1", "pct.2","p_val",  "p_val_adj")]
write.csv(markers2, 'SV2C.markers.csv',row.names = FALSE, quote = F)



