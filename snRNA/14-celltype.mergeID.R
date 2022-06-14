library(Seurat)
library(scrattch.hicat)
library(dendextend)
library(dplyr)
library(matrixStats)
library(Matrix)
#source('wpp_de.genes.R')
#source('wpp_markers.R')

#seu <- readRDS('seu.harmony.anno.rds')
#
#
#barcode_df <- FetchData(seu, vars = c("disease","hicat_cluster","hicat_cluster_merge","hicat_cluster_anno"))
#barcode_df$barcode <- rownames(barcode_df)
#
#barcode_df$subClusteringID <- as.vector(barcode_df$hicat_cluster_merge)
#
#DimPlot(seu,group.by='kmeans')
#DimPlot(seu,group.by='seurat_clusters')
#
#
#data <- as.matrix(seu@assays$RNA@data)
#km_result <- kmeans(t(data), 8)
#table(km_result$cluster)
#seu$kmeans <- km_result$cluster[colnames(seu)]
#
#
#
#seu <- readRDS('seu.harmony.LINC00507.rds')
#barcode_df[colnames(seu),'subClusteringID'] <- paste0('LINC00507_',seu$seurat_clusters)#
#
#seu <- readRDS('seu.harmony.THEMIS.rds')
#barcode_df[colnames(seu),'subClusteringID'] <- paste0('THEMIS_',seu$seurat_clusters)#
#
#seu <- readRDS('seu.harmony.RORB.rds')
#barcode_df[colnames(seu),'subClusteringID'] <- paste0('RORB_',seu$seurat_clusters)#
#
#seu <- readRDS('seu.harmony.FEZF2.2.rds')
#barcode_df[colnames(seu),'subClusteringID'] <- paste0('FEZF2.2_',seu$seurat_clusters)#
#
#seu <- readRDS('seu.harmony.FEZF2.3.rds')
#barcode_df[colnames(seu),'subClusteringID'] <- paste0('FEZF2.3_',seu$seurat_clusters)#
#
#seu <- readRDS('seu.harmony.FEZF2.4.rds')
#barcode_df[colnames(seu),'subClusteringID'] <- paste0('FEZF2.4_',seu$seurat_clusters)#
#
#seu <- readRDS('seu.harmony.FEZF2.5.rds')
#barcode_df[colnames(seu),'subClusteringID'] <- paste0('FEZF2.5_',seu$seurat_clusters)#
#
#
#barcode_df2 <- barcode_df[,c("barcode","hicat_cluster","hicat_cluster_merge","hicat_cluster_anno","subClusteringID" )]
#
#write.table(barcode_df2,'Brain_subClusteringID_seurat.txt', col.names = T,row.names = F,sep='\t',quote = F)
#
#
#
#barcode_df2$subClusteringID_merge <- barcode_df2$subClusteringID
#
#barcode_df2$subClusteringID_merge[barcode_df2$subClusteringID %in% c('FEZF2.2_3','FEZF2.2_0')] <- 'FEZF2.2_0,3'
#barcode_df2$subClusteringID_merge[barcode_df2$subClusteringID %in% c('FEZF2.3_0','FEZF2.3_5','FEZF2.3_6','FEZF2.3_7')] <- 'FEZF2.3_0,5,6,7'
#barcode_df2$subClusteringID_merge[barcode_df2$subClusteringID %in% c('FEZF2.4_0','FEZF2.4_9')] <- 'FEZF2.4_0,9'
#barcode_df2$subClusteringID_merge[barcode_df2$subClusteringID %in% c('FEZF2.4_1','FEZF2.4_2','FEZF2.4_5')] <- 'FEZF2.4_1,2,5'
#barcode_df2$subClusteringID_merge[barcode_df2$subClusteringID %in% c('FEZF2.5_0','FEZF2.5_4','FEZF2.5_5','FEZF2.5_11')] <- 'FEZF2.5_0,4,5,11'
#barcode_df2$subClusteringID_merge[barcode_df2$subClusteringID %in% c('FEZF2.5_2','FEZF2.5_9','FEZF2.5_10')] <- 'FEZF2.5_2,9,10'
#barcode_df2$subClusteringID_merge[barcode_df2$subClusteringID %in% c('FEZF2.5_6','FEZF2.5_7')] <- 'FEZF2.5_6,7'
#barcode_df2$subClusteringID_merge[barcode_df2$subClusteringID %in% c('LINC00507_0','LINC00507_6','LINC00507_16')] <- 'LINC00507_0,6,16'
#barcode_df2$subClusteringID_merge[barcode_df2$subClusteringID %in% c('LINC00507_3','LINC00507_1','LINC00507_11')] <- 'LINC00507_1,3,11'
#barcode_df2$subClusteringID_merge[barcode_df2$subClusteringID %in% c('LINC00507_12','LINC00507_4','LINC00507_7')] <- 'LINC00507_4,7,12'
#barcode_df2$subClusteringID_merge[barcode_df2$subClusteringID %in% c('LINC00507_8','LINC00507_13','LINC00507_19')] <- 'LINC00507_8,13,19'
#barcode_df2$subClusteringID_merge[barcode_df2$subClusteringID %in% c('LINC00507_5','LINC00507_15')] <- 'LINC00507_5,15'
#barcode_df2$subClusteringID_merge[barcode_df2$subClusteringID %in% c('RORB_1','RORB_17')] <- 'RORB_1,17'
#barcode_df2$subClusteringID_merge[barcode_df2$subClusteringID %in% c('RORB_10','RORB_19')] <- 'RORB_10_19'
#barcode_df2$subClusteringID_merge[barcode_df2$subClusteringID %in% c('RORB_14','RORB_8','RORB_18','RORB_4')] <- 'RORB_4,8,14,18'
#barcode_df2$subClusteringID_merge[barcode_df2$subClusteringID %in% c('RORB_2','RORB_9','RORB_20')] <- 'RORB_2,9,20'
#barcode_df2$subClusteringID_merge[barcode_df2$subClusteringID %in% c('RORB_3','RORB_6','RORB_7')] <- 'RORB_3,6,7'
#barcode_df2$subClusteringID_merge[barcode_df2$subClusteringID %in% c('THEMIS_1','THEMIS_5')] <- 'THEMIS_1,5'
#
##x1 <- readRDS('seu.harmony.X1.1.rds')
##barcode_df2[colnames(x1),'subClusteringID_merge'] <- paste0('412_',x1$seurat_clusters)
#
#barcode_df2$subclass <- barcode_df2$hicat_cluster_anno
#barcode_df2$subclass[barcode_df2$hicat_cluster_anno %in% c('Exc FEZF2.2')] <- 'ET'
#barcode_df2$subclass[barcode_df2$hicat_cluster_anno %in% c('Exc FEZF2.3')] <- 'NP'
#barcode_df2$subclass[barcode_df2$hicat_cluster_anno %in% c('Exc LINC00507')] <- 'L2/3_IT'
#barcode_df2$subclass[barcode_df2$hicat_cluster_anno %in% c('Exc RORB')] <- 'L5_IT'
#barcode_df2$subclass[barcode_df2$hicat_cluster_anno %in% c('Exc THEMIS')] <- 'ET'
#barcode_df2$subclass[barcode_df2$subClusteringID %in% c('Exc THEMIS')] <- 'ET'
#
#write.table(barcode_df2,'Brain_subClusteringID_seurat_merge_v1.txt', col.names = T,row.names = F,sep='\t',quote = F)
#saveRDS(barcode_df2, 'Brain_subClusteringID_seurat_merge_v1.rds')

#barcode_df2 <- barcode_df2[!barcode_df2$subClusteringID_merge %in% c('412','448'),]
#write.table(barcode_df2,'Brain_subClusteringID_seurat_merge_v2.txt', col.names = T,row.names = F,sep='\t',quote = F)
#saveRDS(barcode_df2, 'Brain_subClusteringID_seurat_merge_v2.rds')



######## Tree ########################
#seu <- readRDS('seu.harmony.anno.rds')
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

#DElist <- c()
#for(i in 1:dim(pairs)[1]){ #
#  print(i)
#  cl0 <- cl[cl %in% c(pairs[i,1],pairs[i,2])]
#  norm.dat0 <- norm.dat[,names(cl0)]
#  de <- select_markers(norm.dat0, 
#                      cl0, 
#                      n.markers = 20)
#  DElist$markers <- unique(c(DElist$markers,de$markers))
#  DElist$de.genes[[rownames(pairs)[i]]] <- de$de.genes[[1]]
#}
#
#saveRDS(DElist,'DElist.rds')
#


DElist0<-readRDS('DElist.rds')

#pairs <- pairs[!row.names(pairs) %in% names(DElist0$de.genes),]
#saveRDS(pairs,'pairs.rds')
table(row.names(pairs) %in% names(DElist0$de.genes))
#FALSE  TRUE
#2315  7696


DE1 <- readRDS('DElist1.rds')
DE2 <- readRDS('DElist2.rds')
DE3 <- readRDS('DElist3.rds')
DE4 <- readRDS('DElist4.rds')
DE5 <- readRDS('DElist5.rds')
DE6 <- readRDS('DElist6.rds')
DE7 <- readRDS('DElist7.rds')
DE8 <- readRDS('DElist8.rds')
DE9 <- readRDS('DElist9.rds')
DE10 <- readRDS('DElist10.rds')
DE11 <- readRDS('DElist1-1.rds')

DE <- c()
DE$de.genes <- c(DE1$de.genes,
                 DE2$de.genes,
                 DE3$de.genes,
                 DE4$de.genes,
                 DE5$de.genes,
                 DE6$de.genes,
                 DE7$de.genes,
                 DE8$de.genes,
                 DE9$de.genes,
                 DE10$de.genes,
                 DE11$de.genes)

DE$markers <- unique(c(DE1$markers,
                       DE2$markers,
                       DE3$markers,
                       DE4$markers,
                       DE5$markers,
                       DE6$markers,
                       DE7$markers,
                       DE8$markers,
                       DE9$markers,
                       DE10$markers,
                       DE11$markers))

table(rownames(pairs) %in% names(DE$de.genes))
#FALSE  TRUE 
#7379  2632 

pairs2 <- pairs[!rownames(pairs) %in% names(DE$de.genes),]
for(j in 1:dim(pairs2)[1]){
  DE$de.genes[[rownames(pairs2)[j]]] <- DElist0$de.genes[[rownames(pairs2)[j]]]
}
DE$markers <- unique(c(DE$markers,DElist0$markers))

table(rownames(pairs) %in% names(DE$de.genes))
#FALSE  TRUE 
#    1 10010 

rownames(pairs)[!(rownames(pairs) %in% names(DE$de.genes))]
#"395_476"

saveRDS(DE,'Brain_subClusteringID_seurat_merge_DE20.rds')


marker.genes <- DE$markers
de.genes <- DE$de.genes

cl.med <- get_cl_medians(norm.dat[marker.genes,], 
                         cl)

dend.result <- build_dend(cl.med,
                          nboot = 100)
saveRDS(dend.result, 'Brain_subClusteringID_seurat_merge_dend.rds')


aa <- labels(dend.result$dend)

aa <- gsub('p',',',aa)
aa <- gsub('n','_',aa)
aa <- gsub('d','\\.',aa)

labels(dend.result$dend) <- aa

saveRDS(dend.result, 'Brain_subClusteringID_seurat_merge_dend_renames.rds')

### dend.pdf 30x4
plot(dend.result$dend,
     horiz = T)
