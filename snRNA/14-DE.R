library(scrattch.hicat)
library(dendextend)
library(Seurat)

seu0 <- readRDS('seu.harmony.anno.rds')


clusters1 <- c("88","65","77","86","89")
clusters2 <- c("25","30","29","31","58","425","68","70","40","36","39","51","42","48","49","72","93","125","111","104","109",
               "98","115","123","124")
cells1 <- colnames(seu0)[seu0$hicat_cluster_merge %in% clusters1]
cells2 <- colnames(seu0)[seu0$hicat_cluster_merge %in% clusters2]

XXX <- FindMarkers(seu0,
                   ident.1 = clusters1,
                   ident.2 = clusters2)
saveRDS(XXX,'Exc_LINC00507_markers.rds')


clusters1 <- c("40","36","39","51","42","48","49")
clusters2 <- c("25","30","29","31","58","425","68","70","72","93","125","111","104","109",
               "98","115","123","124","88","65","77","86","89")
cells1 <- colnames(seu0)[seu0$hicat_cluster_merge %in% clusters1]
cells2 <- colnames(seu0)[seu0$hicat_cluster_merge %in% clusters2]

XXX <- FindMarkers(seu0,
                   ident.1 = clusters1,
                   ident.2 = clusters2)
saveRDS(XXX,'Exc_FEZF2.1_markers.rds')



clusters1 <- c("68","70")
clusters2 <- c("25","30","29","31","58","425","40","36","39","51","42","48","49","72","93","125","111","104","109",
               "98","115","123","124","88","65","77","86","89")
cells1 <- colnames(seu0)[seu0$hicat_cluster_merge %in% clusters1]
cells2 <- colnames(seu0)[seu0$hicat_cluster_merge %in% clusters2]

XXX <- FindMarkers(seu0,
                   ident.1 = clusters1,
                   ident.2 = clusters2)
saveRDS(XXX,'Exc_FEZF2.2_markers.rds')



clusters1 <- c("25","30","29","31")
clusters2 <- c("58","425","68","70","40","36","39","51","42","48","49","72","93","125","111","104","109",
               "98","115","123","124","88","65","77","86","89")
cells1 <- colnames(seu0)[seu0$hicat_cluster_merge %in% clusters1]
cells2 <- colnames(seu0)[seu0$hicat_cluster_merge %in% clusters2]

XXX <- FindMarkers(seu0,
                   ident.1 = clusters1,
                   ident.2 = clusters2)
saveRDS(XXX,'Exc_FEZF2.3_markers.rds')


clusters1 <- c("72","93","125","111","104","109",
               "98","115","123","124")
clusters2 <- c("25","30","29","31","58","425","68","70","40","36","39","51","42","48","49","88","65","77","86","89")
cells1 <- colnames(seu0)[seu0$hicat_cluster_merge %in% clusters1]
cells2 <- colnames(seu0)[seu0$hicat_cluster_merge %in% clusters2]

XXX <- FindMarkers(seu0,
                   ident.1 = clusters1,
                   ident.2 = clusters2)
saveRDS(XXX,'Exc_RORB_markers.rds')


