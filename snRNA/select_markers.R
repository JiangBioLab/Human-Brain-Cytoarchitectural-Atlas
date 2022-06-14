library(Seurat)
library(scrattch.hicat)
library(dendextend)
library(dplyr)
library(matrixStats)
library(Matrix)

barcode_df2 <- readRDS('Brain_subClusteringID_seurat_merge.rds')
norm.dat <- readRDS('norm.dat.rds')

cl.clean <- barcode_df2$subClusteringID_merge
names(cl.clean) <- barcode_df2$barcode

cl <- gsub(',','p',cl.clean)
cl <- gsub('_','n',cl)
cl <- gsub('\\.','d',cl)
cn <- as.character(sort(unique(cl)))
pairs= create_pairs(cn)
row.names(pairs) <- paste(pairs[, 1], pairs[, 2], sep = "_")
DElist <- c()
for(i in 1:dim(pairs)[1]){ #
  print(i)
  cl0 <- cl[cl %in% c(pairs[i,1],pairs[i,2])]
  norm.dat0 <- norm.dat[,names(cl0)]
  de <- select_markers(norm.dat0, 
                       cl0, 
                       n.markers = 50)
  DElist$markers <- unique(c(DElist$markers,de$markers))
  DElist$de.genes[[rownames(pairs)[i]]] <- de$de.genes[[1]]
}

saveRDS(DElist,'DElist.rds')


