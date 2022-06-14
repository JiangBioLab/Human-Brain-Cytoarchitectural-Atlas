library(Seurat)
library(scrattch.hicat)
library(dendextend)
library(dplyr)
library(matrixStats)
library(Matrix)

norm.dat <- readRDS('norm.dat.rds')
anno <- readRDS('anno.rds')


#celltype <- unique(anno$celltype)
#DE_genes <- c()
#for(i in 1:length(celltype)){
#  cl <- as.vector(anno$celltype)
#  names(cl) <- anno$cells
#  cl[!cl == celltype[i]] <- 'Other'
#  de <- de_all_pairs(norm.dat, cl)
#  de2 <- cbind(rep(names(de)[1],dim(de[[1]])[1]),rownames(de[[1]]),de[[1]])
#  colnames(de2)[1:2] <- c('celltype','gene')
#  DE_genes <- rbind(DE_genes, de2)
#  print(paste0("Done:",celltype[i]))
#}
#saveRDS(DE_genes, 'celltype_DE_genes.rds')



celltype <- unique(anno$celltype)
DE_genes <- c()
for(i in 1:length(celltype)){
  anno.tmp <- anno[anno$celltype==celltype[i],]
  norm.dat.tmp <- norm.dat[,rownames(anno.tmp)]
  
  cluster <- unique(anno.tmp$cluster)
  for(j in 1:length(cluster)){
    cl <- as.vector(anno.tmp$cluster)
    names(cl) <- anno.tmp$cells
    cl[!cl == cluster[j]] <- 'Other'
    de <- de_all_pairs(norm.dat.tmp, cl)
    de2 <- cbind(rep(celltype[i],dim(de[[1]])[1]),rep(names(de)[1],dim(de[[1]])[1]),rownames(de[[1]]),de[[1]])
    colnames(de2)[1:3] <- c('celltype','cluster','gene')
    DE_genes <- rbind(DE_genes, de2)
    print(paste0("Done:",celltype[i]," ",cluster[j]))
  }
  
}
saveRDS(DE_genes, 'cluster_DE_genes.rds')




