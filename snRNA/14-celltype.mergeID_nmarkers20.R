library(Seurat)
library(scrattch.hicat)
library(dendextend)
library(dplyr)
library(matrixStats)
library(Matrix)

barcode_df2 <- readRDS('Brain_subClusteringID_seurat_merge_v2.rds')
norm.dat <- readRDS('norm.dat.v2.rds')

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
  if(i == 7175){
    next
  }
  print(i)
  cl0 <- cl[cl %in% c(pairs[i,1],pairs[i,2])]
  norm.dat0 <- norm.dat[,names(cl0)]
  de <- select_markers(norm.dat0, 
                      cl0, 
                      n.markers = 20)
  DElist$markers <- unique(c(DElist$markers,de$markers))
  DElist$de.genes[[rownames(pairs)[i]]] <- de$de.genes[[1]]
  
  if(!(i %% 10)){
    saveRDS(DElist,'DElist.rds')
  }
}

saveRDS(DElist,'Brain_subClusteringID_seurat_merge_DE20_2.rds')


DE <- readRDS('Brain_subClusteringID_seurat_merge_DE20_2.rds')

marker.genes <- DE$markers
de.genes <- DE$de.genes

cl.med <- get_cl_medians(norm.dat[marker.genes,], 
                         cl)

dend.result <- build_dend(cl.med,
                          nboot = 100)
saveRDS(dend.result, 'Brain_subClusteringID_seurat_merge_dend2.rds')


aa <- labels(dend.result$dend)

aa <- gsub('p',',',aa)
aa <- gsub('n','_',aa)
aa <- gsub('d','\\.',aa)

labels(dend.result$dend) <- aa

saveRDS(dend.result, 'Brain_subClusteringID_seurat_merge_dend_renames2.rds')

### dend.pdf 30x4
plot(dend.result$dend,
     horiz = T)

