library(scrattch.hicat)
library(dendextend)
library(Seurat)

seu0 <- readRDS('seu.harmony.anno.rds')

#
#clusters1 <- c("448","412","413")
#clusters2 <- c("386","389","392","395","478","476","477","21","1","17","19","22","401","403","404","453","446","454","456","461","463","471","466","469")
#cells1 <- colnames(seu0)[seu0$hicat_cluster_merge %in% clusters1]
#cells2 <- colnames(seu0)[seu0$hicat_cluster_merge %in% clusters2]
#
#XXX <- FindMarkers(seu0,
#                   ident.1 = clusters1,
#                   ident.2 = clusters2)
#saveRDS(XXX,'NonN_X1_markers.rds')
#
#
#clusters1 <- c("461","463","471","466","469")
#clusters2 <- c("386","389","392","395","478","476","477","21","1","17","19","22","401","403","404","453","446","454","456","448","412","413")
#cells1 <- colnames(seu0)[seu0$hicat_cluster_merge %in% clusters1]
#cells2 <- colnames(seu0)[seu0$hicat_cluster_merge %in% clusters2]
#
#XXX <- FindMarkers(seu0,
#                   ident.1 = clusters1,
#                   ident.2 = clusters2)
#saveRDS(XXX,'NonN_X2_markers.rds')
#
#
#clusters1 <- c("453","446","454","456")
#clusters2 <- c("386","389","392","395","478","476","477","21","1","17","19","22","401","403","404","461","463","471","466","469","448","412","413")
#cells1 <- colnames(seu0)[seu0$hicat_cluster_merge %in% clusters1]
#cells2 <- colnames(seu0)[seu0$hicat_cluster_merge %in% clusters2]
#
#XXX <- FindMarkers(seu0,
#                   ident.1 = clusters1,
#                   ident.2 = clusters2)
#saveRDS(XXX,'NonN_X3_markers.rds')
#
#
#
#clusters1 <- c("401","403","404")
#clusters2 <- c("386","389","392","395","478","476","477","21","1","17","19","22","453","446","454","456","461","463","471","466","469","448","412","413")
#cells1 <- colnames(seu0)[seu0$hicat_cluster_merge %in% clusters1]
#cells2 <- colnames(seu0)[seu0$hicat_cluster_merge %in% clusters2]
#
#XXX <- FindMarkers(seu0,
#                   ident.1 = clusters1,
#                   ident.2 = clusters2)
#saveRDS(XXX,'NonN_X4_markers.rds')
#
#
#clusters1 <- c("21","1","17","19","22")
#clusters2 <- c("386","389","392","395","478","476","477","401","403","404","453","446","454","456","461","463","471","466","469","448","412","413")
#cells1 <- colnames(seu0)[seu0$hicat_cluster_merge %in% clusters1]
#cells2 <- colnames(seu0)[seu0$hicat_cluster_merge %in% clusters2]
#
#XXX <- FindMarkers(seu0,
#                   ident.1 = clusters1,
#                   ident.2 = clusters2)
#saveRDS(XXX,'NonN_X5_markers.rds')


clusters1 <- c("478","476","477")
clusters2 <- c("386","389","392","395","21","1","17","19","22","401","403","404","453","446","454","456","461","463","471","466","469","448","412","413")
cells1 <- colnames(seu0)[seu0$hicat_cluster_merge %in% clusters1]
cells2 <- colnames(seu0)[seu0$hicat_cluster_merge %in% clusters2]

XXX <- FindMarkers(seu0,
                   ident.1 = clusters1,
                   ident.2 = clusters2)
saveRDS(XXX,'NonN_X6_markers.rds')



clusters1 <- c("386","389","392","395")
clusters2 <- c("478","476","477","21","1","17","19","22","401","403","404","453","446","454","456","461","463","471","466","469","448","412","413")
cells1 <- colnames(seu0)[seu0$hicat_cluster_merge %in% clusters1]
cells2 <- colnames(seu0)[seu0$hicat_cluster_merge %in% clusters2]

XXX <- FindMarkers(seu0,
                   ident.1 = clusters1,
                   ident.2 = clusters2)
saveRDS(XXX,'NonN_X7_markers.rds')



