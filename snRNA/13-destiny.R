library(destiny)
library(ggthemes)
library(Seurat)

seu <- readRDS('seu.harmony.anno.Glu.rds')
norm.dat<-seu@assays$RNA@data

Glu <- seu@assays$RNA@data
cellLabels <- seu$hicat_cluster_newID
colnames(Glu) <- cellLabels

Glu<-as.matrix(Glu)

# Make a diffusion map.
Glu_dm <- DiffusionMap(t(Glu)) #,n_pcs = 50
# Start:14:54
# End: 20:00-23:00?
saveRDS(Glu_dm,'Glu_destiny.rds')