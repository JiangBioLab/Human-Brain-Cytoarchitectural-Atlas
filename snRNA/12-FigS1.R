library(scrattch.hicat)
library(dendextend)
library(Seurat)

seu <- readRDS('seu.harmony.anno.rds')

result<-readRDS('seu.hicat.run_consensus_clust.de100.niter100.rds')

result$co.result$cl.list


