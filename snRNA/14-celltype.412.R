library(scrattch.hicat)
library(dendextend)
library(Seurat)
library(harmony)

seu0 <- readRDS('seu.harmony.anno.rds')
result <- readRDS('seu.hicat.de80.niter100.rds')
dend.result <- readRDS('seu.hicat.de80.niter100.dend.rds')

c <- c("386","389","392","395","478","476","477","21","1","17","19","22","401","403","404","453","446","454","456","461","463","471","466","469","448","412","413")

###################### XXX markers ##############################################
clusters1 <- c("412")
clusters2 <- c("386","389","392","395","478","476","477","21","1","17","19","22","401","403","404","453","446","454","456","461","463","471","466","469","448","413")
cells1 <- colnames(seu0)[seu0$hicat_cluster_merge %in% clusters1]
cells2 <- colnames(seu0)[seu0$hicat_cluster_merge %in% clusters2]

XXX <- FindMarkers(seu0,
                   ident.1 = clusters1,
                   ident.2 = clusters2)
saveRDS(XXX,'NonN_X1.1_markers.rds')

XXX<-readRDS('NonN_X1.1_markers.rds')

top5 <- XXX[XXX$avg_log2FC>0,][1:50,]
genes <- rownames(top5)

#412_up-regulated.pdf 30x20
DotPlot(seu0, features=unique(genes)) + RotatedAxis()


######################################################################################

seu0 <- readRDS('seu.harmony.X1.rds')

clusters <- c("412")

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

saveRDS(seu,'seu.harmony.X1.1.rds')

seu <- readRDS('seu.harmony.X1.1.rds')

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu <- NormalizeData(object = seu, normalization.method = "LogNormalize", scale.factor = 1e4)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu, features = rownames(seu))
seu <- RunPCA(seu, features = VariableFeatures(object = seu))
seu <- RunHarmony(seu, "sample")
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:10)
seu <- FindClusters(seu, resolution = 0.5)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:10)

#seu@active.ident <- seu$hicat_cluster_merge
markers <- FindAllMarkers(seu)

saveRDS(markers, 'seu.harmony.X1.1.markers.rds')
saveRDS(seu,'seu.harmony.X1.1.rds')



################################################
markers <- readRDS('seu.harmony.X1.1.markers.rds')
seu <- readRDS('seu.harmony.X1.1.rds')

DimPlot(seu,label=T)
DimPlot(seu,group.by = 'individual')
DimPlot(seu,group.by = 'region')

library(dplyr)
top5 <- top_n(group_by(markers,cluster), 20, avg_log2FC)
genes <- top5$gene

genes <- c(genes,
           'SOX2', 'PAX6', 'HES5',
           'EOMES', 'NEUROG2', 'BTG2',
           'NEUROD2', 'TUBB3', 'NEUROD6',
           'DLX2', 'GAD1', 'GAD2',
           'SST', 'NPY', 'LHX6', 'NXPH1', 'NXPH2',
           'PAX6', 'SP8', 'CXCL14', 'HTR3A',
           'MEIS2', 'ETV1', 'SP8',
           'CCK', 'PCP4',
           'BTG2', 'NEUROG2', 'HES6', 'FABP7', 'DBI', 'SLC1A3',
           'FEZF2', 'TLE4', 'BCL11B', 'CUX1', 'POU3F3', 'SATB2',
           'GADD45G', 'NEUROG2', 'SSTR2', 'NEUROD1'
)

DotPlot(seu, features=unique(genes)) + RotatedAxis()

markers2 <- markers[,c("cluster", "gene","avg_log2FC", "pct.1", "pct.2","p_val",  "p_val_adj")]
write.csv(markers2, 'X1.1.markers.csv',row.names = FALSE, quote = F)

FeaturePlot(seu, 
            features = unique(genes),
            ncol = 5)



