library(scrattch.hicat)
library(dendextend)
library(Seurat)

seu <- readRDS('seu.v2.rds')
meta <- readRDS('seu.harmony.anno.meta.v2.rds')

seu$hicat_cluster_newID <- meta[colnames(seu),'hicat_cluster_classes']
seu$hicat_cluster_newID <- factor(seu$hicat_cluster_newID)
seu@active.ident <- seu$hicat_cluster_newID

clusters <- levels(seu$hicat_cluster_newID)
length(clusters)
#4

markers <- c()
for(i in 1:4){
  tmp <- FindMarkers(seu, ident.1 = clusters[i])
  tmp$gene <- rownames(tmp)
  tmp$cluster <- rep(clusters[i],dim(tmp)[1])
  markers <- rbind(markers,tmp)
  saveRDS(markers,'markers_classes.rds')
}

saveRDS(markers,'markers_classes.rds')


