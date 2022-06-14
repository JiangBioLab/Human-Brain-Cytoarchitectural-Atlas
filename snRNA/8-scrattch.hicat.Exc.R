#library(tasic2016data)
library(Seurat)
library(scrattch.hicat)
library(dendextend)
library(dplyr)
library(matrixStats)
library(Matrix)

seu <- readRDS('seu.harmony.anno.Glu.rds')
norm.dat<-seu@assays$RNA@data

de.param <- de_param(padj.th     = 0.05, 
                     lfc.th      = 1, 
                     low.th      = 1, 
                     q1.th       = 0.3,
                     q2.th       = NULL,
                     q.diff.th   = 0.7, 
                     de.score.th = 40,
                     min.cells = 20)


load('subsample_PCA_de40_niter100_Exc/cells.100.rda')
result <- scrattch.hicat::iter_clust(norm.dat=norm.dat, dim.method = "pca", select.cells=select.cells, de.param = de.param)
save(result, file='subsample_PCA_de40_niter100_Exc/result.100.rda')

load('subsample_PCA_de40_niter100_Exc/cells.10.rda')
result <- scrattch.hicat::iter_clust(norm.dat=norm.dat, dim.method = "pca", select.cells=select.cells, de.param = de.param)
save(result, file='subsample_PCA_de40_niter100_Exc/result.10.rda')

load('subsample_PCA_de40_niter100_Exc/cells.1.rda')
result <- scrattch.hicat::iter_clust(norm.dat=norm.dat, dim.method = "pca", select.cells=select.cells, de.param = de.param)
save(result, file='subsample_PCA_de40_niter100_Exc/result.1.rda')

load('subsample_PCA_de40_niter100_Exc/cells.2.rda')
result <- scrattch.hicat::iter_clust(norm.dat=norm.dat, dim.method = "pca", select.cells=select.cells, de.param = de.param)
save(result, file='subsample_PCA_de40_niter100_Exc/result.2.rda')

load('subsample_PCA_de40_niter100_Exc/cells.3.rda')
result <- scrattch.hicat::iter_clust(norm.dat=norm.dat, dim.method = "pca", select.cells=select.cells, de.param = de.param)
save(result, file='subsample_PCA_de40_niter100_Exc/result.3.rda')

load('subsample_PCA_de40_niter100_Exc/cells.4.rda')
result <- scrattch.hicat::iter_clust(norm.dat=norm.dat, dim.method = "pca", select.cells=select.cells, de.param = de.param)
save(result, file='subsample_PCA_de40_niter100_Exc/result.4.rda')

load('subsample_PCA_de40_niter100_Exc/cells.5.rda')
result <- scrattch.hicat::iter_clust(norm.dat=norm.dat, dim.method = "pca", select.cells=select.cells, de.param = de.param)
save(result, file='subsample_PCA_de40_niter100_Exc/result.5.rda')

load('subsample_PCA_de40_niter100_Exc/cells.6.rda')
result <- scrattch.hicat::iter_clust(norm.dat=norm.dat, dim.method = "pca", select.cells=select.cells, de.param = de.param)
save(result, file='subsample_PCA_de40_niter100_Exc/result.6.rda')

load('subsample_PCA_de40_niter100_Exc/cells.7.rda')
result <- scrattch.hicat::iter_clust(norm.dat=norm.dat, dim.method = "pca", select.cells=select.cells, de.param = de.param)
save(result, file='subsample_PCA_de40_niter100_Exc/result.7.rda')

load('subsample_PCA_de40_niter100_Exc/cells.8.rda')
result <- scrattch.hicat::iter_clust(norm.dat=norm.dat, dim.method = "pca", select.cells=select.cells, de.param = de.param)
save(result, file='subsample_PCA_de40_niter100_Exc/result.8.rda')

load('subsample_PCA_de40_niter100_Exc/cells.9.rda')
result <- scrattch.hicat::iter_clust(norm.dat=norm.dat, dim.method = "pca", select.cells=select.cells, de.param = de.param)
save(result, file='subsample_PCA_de40_niter100_Exc/result.9.rda')





#set.seed(12345)
#
#result <- run_consensus_clust_wpp(norm.dat, 
#                              niter = 100, 
#                              #select.cells = select.cells,
#                              de.param = de.param, 
#                              #rm.eigen = rm.eigen, 
#                              dim.method = "pca", 
#                              output_dir = "subsample_PCA_de80_niter100_Exc", 
#                              mc.cores = 10)
#saveRDS(result,'seu.hicat.run_consensus_clust.de80.niter100.Exc.rds')



#set.seed(12345)
#result <- run_consensus_clust_wpp(norm.dat, 
#                                  niter = 100, 
#                                  #select.cells = select.cells,
#                                  de.param = de.param, 
#                                  #rm.eigen = rm.eigen, 
#                                  dim.method = "pca", 
#                                  output_dir = "subsample_PCA_de80_niter100", 
#                                  mc.cores = 10)
#saveRDS(result,'seu.hicat.run_consensus_clust.de80.niter100.rds')




#result <- readRDS('seu.hicat.run_consensus_clust.de80.niter100.rds')
#result.dend <- readRDS('seu.hicat.run_consensus_clust.de80.niter100.dend.rds')
#
#
#dend <- result.dend$dend
#
#seu$tmp <- result$cl.result$cl[colnames(seu)]
#clusters <- unique(seu$tmp)
#anno <- cbind(clusters,clusters)
#for(i in 1:length(clusters)){
#  a <- unique(seu$hicat_cluster_subclasses[seu$tmp==clusters[i]])
#  b <- ''
#  for(j in 1:length(a)){
#    b <- paste0(b,a[j],';')
#  }
#  anno[i,2] <- b
#}
#
#rownames(anno) <- anno[,1]
#anno[,2] <- paste0(anno[,1],':',anno[,2])
#
#dend.label <- dend
#labels(dend.label) <- anno[labels(dend),2]
##dend.de80.pdf 20x4
#plot(dend,horiz = T)
#plot(dend.label,horiz = T)
#
#
#
#
#
#
#