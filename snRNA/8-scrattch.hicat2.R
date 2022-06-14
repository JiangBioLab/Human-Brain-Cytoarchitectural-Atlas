#library(tasic2016data)
library(Seurat)
library(scrattch.hicat)
library(dendextend)
library(dplyr)
library(matrixStats)
library(Matrix)
source('consensusCluster.R')

norm.dat<-readRDS('norm.dat.rds')

de.param <- de_param(padj.th     = 0.05, 
                     lfc.th      = 1, 
                     low.th      = 1, 
                     q1.th       = 0.3,
                     q2.th       = NULL,
                     q.diff.th   = 0.7, 
                     de.score.th = 120,
                     min.cells = 20)

set.seed(12345)

result <- run_consensus_clust_wpp(norm.dat, 
                              niter = 100, 
                              #select.cells = select.cells,
                              de.param = de.param, 
                              #rm.eigen = rm.eigen, 
                              dim.method = "pca", 
                              output_dir = "subsample_PCA_de120_niter100", 
                              mc.cores = 1)
saveRDS(result,'seu.hicat.de120.niter100.rds')



############ Tree #########################################
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
saveRDS(dend.result, 'seu.hicat.de120.niter100.dend.rds')
#dend <- dend.result$dend
####attach cluster labels to the leaves of the tree 
#dend.labeled <- dend
##labels(dend.labeled) <- consensus.cl.df[labels(dend), "cluster_label"]
#plot(dend.labeled)










