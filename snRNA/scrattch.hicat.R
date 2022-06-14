#library(tasic2016data)
library(Seurat)
library(scrattch.hicat)
library(dendextend)
library(dplyr)
library(matrixStats)
library(Matrix)

seu <- readRDS('seu.rds')
norm.dat<-seu@assays$RNA@data

de.param <- de_param(padj.th     = 0.05, 
                     lfc.th      = 1, 
                     low.th      = 1, 
                     q1.th       = 0.3,
                     q2.th       = NULL,
                     q.diff.th   = 0.7, 
                     de.score.th = 100,
                     min.cells = 20)

### 跑一次聚类 ###########
select.cells <- colnames(norm.dat)
result <- scrattch.hicat::iter_clust(norm.dat=norm.dat, 
                                     dim.method = "pca", 
                                     select.cells=select.cells, 
                                     de.param = de.param)




set.seed(123)
### 跑100次聚类计算consensus #############
result <- run_consensus_clust(norm.dat, 
                              niter = 100, 
                              #select.cells = select.cells,
                              de.param = de.param, 
                              #rm.eigen = rm.eigen, 
                              dim.method = "pca", 
                              output_dir = "hicat_out", 
                              mc.cores = 10)


###### Tree ####################
cl.clean <- result$cl.result$cl
de.genes <- result$cl.result$de.genes

select.markers <- select_markers(norm.dat, 
                                 cl.clean, 
                                 de.genes = de.genes,
                                 n.markers = 50)

marker.genes <- select.markers$markers
de.genes <- select.markers$de.genes

cl.med <- get_cl_medians(norm.dat[marker.genes,], 
                         cl.clean)

# Build the dendrogram
dend.result <- build_dend(cl.med,
                          #l.rank, 
                          #l.color,
                          nboot = 100)
plot(dend.result$dend)


