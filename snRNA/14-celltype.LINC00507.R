library(scrattch.hicat)
library(dendextend)
library(Seurat)
library(harmony)

#seu0 <- readRDS('seu.harmony.v2.rds')
#result <- readRDS('seu.hicat.de80.niter100.rds')
#dend.result <- readRDS('Brain_subClusteringID_seurat_merge_dend_renames.rds')
#barcode_df2 <- readRDS('Brain_subClusteringID_seurat_merge_v2.rds')
#
#seu0$hicat_cluster_merge <- barcode_df2[colnames(seu0),'subClusteringID_merge']
#seu0$hicat_cluster_merge <- factor(seu0$hicat_cluster_merge, 
#                                  levels = rev(labels(dend.result$dend)))
#
#
#c <- c("25","30","29","31","58","425","68","70","40","36","39","51","42","48","49","72","93","125","111","104","109",
#"98","115","123","124","88","65","77","86","89")
#
####################### XXX markers ##############################################
#clusters1 <- c("88","65","77","86","89")
#clusters2 <- c("25","30","29","31","58","425","68","70","40","36","39","51","42","48","49","72","93","125","111","104","109",
#               "98","115","123","124")
#cells1 <- colnames(seu0)[seu0$hicat_cluster_merge %in% clusters1]
#cells2 <- colnames(seu0)[seu0$hicat_cluster_merge %in% clusters2]
#
#XXX <- FindMarkers(seu0,
#                   ident.1 = clusters1,
#                   ident.2 = clusters2)
#saveRDS(XXX,'Exc_LINC00507_markers.rds')
#
#XXX<-readRDS('Exc_LINC00507_markers.rds')
#top5 <- XXX[XXX$avg_log2FC>0,][1:50,]
#genes <- rownames(top5)
#
##LINC00507_up-regulated.pdf 30x20
#DotPlot(seu0, features=unique(genes)) + RotatedAxis()

#################################################################################

#
#clusters <- c("LINC00507_18", "LINC00507_9", "LINC00507_8,13,19", "LINC00507_0,6,16", "LINC00507_2", "LINC00507_17", "LINC00507_1,3,11", 
#              "LINC00507_5,15", "LINC00507_10", "LINC00507_4,7,12")
#
#counts <- seu0@assays$RNA@counts[,seu0$hicat_cluster_merge %in% clusters]
#
#seu <- CreateSeuratObject(
#  counts,
#  min.cells = 0,
#  min.features = 0
#)
## 33014 36013
#
#meta <- data.frame(cells = colnames(seu0),
#                   sample = seu0$sample,
#                   individual = seu0$individual,
#                   region = seu0$region,
#                   hicat_cluster = seu0$hicat_cluster,
#                   hicat_cluster_merge = seu0$hicat_cluster_merge)
#rownames(meta) <- meta$cells
#
#seu$sample <- meta[colnames(seu),'sample']
#seu$individual <- meta[colnames(seu),'individual']
#seu$region <- meta[colnames(seu),'region']
#seu$hicat_cluster <- meta[colnames(seu),'hicat_cluster']
#seu$hicat_cluster_merge <- meta[colnames(seu),'hicat_cluster_merge']
#
#
#seu$hicat_cluster_merge <- factor(seu$hicat_cluster_merge, levels = clusters)
#
#saveRDS(seu,'seu.harmony.LINC00507.rds')

seu <- readRDS('seu.harmony.LINC00507.rds')

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu <- NormalizeData(object = seu, normalization.method = "LogNormalize", scale.factor = 1e4)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu, features = rownames(seu))
seu <- RunPCA(seu, features = VariableFeatures(object = seu))
seu <- RunHarmony(seu, "sample")
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:20)
seu <- FindClusters(seu, resolution = 1)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:20)

seu@active.ident <- seu$hicat_cluster_merge
markers <- FindAllMarkers(seu)

saveRDS(markers, 'seu.harmony.LINC00507.markers.rds')
saveRDS(seu,'seu.harmony.LINC00507.rds')



################################################
#markers <- readRDS('seu.harmony.LINC00507.markers.rds')
#seu <- readRDS('seu.harmony.LINC00507.rds')
#
#DimPlot(seu,label=T)
#
#library(dplyr)
#top5 <- top_n(group_by(markers,cluster), 5, avg_log2FC)
#genes <- top5$gene
#
#DotPlot(seu, features=unique(genes)) + RotatedAxis()
#
#markers2 <- markers[,c("cluster", "gene","avg_log2FC", "pct.1", "pct.2","p_val",  "p_val_adj")]
#write.csv(markers2, 'LINC00507.markers.csv',row.names = FALSE, quote = F)
#
#
#
#
#
################ K means #######################
#seu <- readRDS('seu.harmony.LINC00507.rds')
#
#seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:10)
#seu <- FindClusters(seu, resolution = 1)
#seu <- RunUMAP(seu, reduction = "harmony", dims = 1:10)
#
#DimPlot(seu, group.by = 'hicat_cluster')
#
#DimPlot(seu, group.by = 'region')
#DimPlot(seu, group.by = 'individual')
#
#
#genes <- c('LINC00507','LAMP5','KCNG3','ATP7B','GLRA3','RORB','RTKN2','DSG3','CARM1P1','OPALIN','MAP6D1')
#
#FeaturePlot(seu,
#            ncol = 4,
#            features=genes)
#
#data <- as.matrix(seu@assays$RNA@data)
#km_result <- kmeans(t(data), 8)
#table(km_result$cluster)
#
#seu$kmeans <- km_result$cluster[colnames(seu)]
#
#
#DimPlot(seu, group.by = 'kmeans')
#
#DimPlot(seu, group.by = 'seurat_clusters', label=T)
#
#DotPlot(seu,features = genes,group.by = 'kmeans') + RotatedAxis()
#DotPlot(seu,features = genes,group.by = 'seurat_clusters') + RotatedAxis()
#
#
#seu$kmeans <- factor(seu$kmeans)
#seu@active.ident <- seu$kmeans
#saveRDS(seu, 'seu.harmony.LINC00507.rds')
#
#markers <- FindAllMarkers(seu)
#saveRDS(markers, 'seu.harmony.LINC00507.kmeans.markers.rds')
#top5 <- top_n(group_by(markers,cluster), 5, avg_log2FC)
#genes <- top5$gene
#DotPlot(seu, features=unique(genes)) + RotatedAxis()
#
#
#seu@active.ident <- seu$seurat_clusters
#markers <- FindAllMarkers(seu)
#saveRDS(markers, 'seu.harmony.LINC00507.seuCluster.markers.rds')
#markers <-readRDS('seu.harmony.LINC00507.seuCluster.markers.rds')
#top5 <- top_n(group_by(markers,cluster), 5, avg_log2FC)
#genes <- top5$gene
#DotPlot(seu, features=unique(genes)) + RotatedAxis()
#
#
#
#
#
#seu <- readRDS('seu.harmony.LINC00507.rds')
#DimPlot(seu,group.by = 'seurat_clusters')
#
#genes0 <- c('PCDH11Y', 'AC021613.1', 'PTPRD', 'GREB1L', 'CDH13', 'PDGFD', 'LINC01500', 'NLGN4Y', 'ACVR1C', 'MGAT4C', 'USP9Y', 'ZNF804A')
#genes1 <- c('CNTN5', 'CCBE1', 'GRM1', 'EPHA6', 'SAMD3', 'SGCD', 'DOK6', 'GRID2', 'EFNA5', 'CADPS2')
#genes2 <- c('XIST', 'CNTNAP5', 'NKAIN2', 'MAML3', 'KCNH4', 'SNTG1', 'ABLIM1', 'HHAT', 'UNC13C')
#genes3 <- c('KAZN', 'LINC01331', 'GPC6', 'AC110023.1', 'ADAMTS3', 'MARCH1', 'PDE1C', 'PDZD2', 'AC091078.1', 'PALMD', 'CTNNA2', 'L3MBTL4', 'AC007368.1')
#genes4 <- c('GLIS3', 'GRID2', 'RORB', 'NTNG1', 'PDZD2', 'LINC02296', 'GRIK1', 'CDH13', 'GPC6', 'AC007368.1', 'ROBO2', 'SGCZ', 'COL5A2')
#genes5 <- c('SV2C', 'SHISA9', 'GLIS3', 'ROBO2', 'VAV3', 'NRG1', 'GALNTL6', 'SYT2', 'TMEM132C', 'KAZN', 'MEG8', 'SCN1A-AS1', 'MIR34AHG', 'CDH18', 'TOX2')
#genes6 <- c('LINC02263', 'AC110614.1', 'LRP1B', 'PDGFD', 'AC022325.2', 'CNTN5', 'LSAMP', 'LINC01378', 'SDK1', 'AC097512.1', 'AC093607.1', 'SEMA6D')
#genes7 <- c('SYT10', 'AC110014.1', 'AC117453.1', 'NTNG1', 'LINC02822', 'AJ009632.2', 'HDAC9', 'LINC01478', 'RPS6KA2', 'ADGRF5', 'ADCY8', 'LRRC4C', 'LINC02301', 'COL5A2')
#genes8 <- c('TRHDE', 'PLCH1', 'NTNG1', 'COL5A2', 'LRRC4C', 'COBLL1', 'CLMN', 'LPP', 'AL157769.1', 'CNTNAP5', 'PCDH7', 'SMARCA2', 'PTPRT')
#genes9 <- c('CTNNA2', 'MCTP1', 'NFIA', 'CNTNAP5', 'LRP1B', 'THEMIS', 'NCAM2', 'GRIK3', 'LINC00343', 'TAFA2', 'SLC35F3', 'HS3ST5', 'SOX5', 'PXDN', 'TMEM233')
#genes10 <- c('ERBB4', 'GRIK1', 'RBMS3')
#genes11 <- c('CCBE1', 'TRPC6', 'LINC02822', 'PKN2-AS1', 'ITPR2', 'CRYBG3', 'SERPINE2', 'CNTN6', 'PRKD1', 'SFMBT2', 'LINC02055', 'LIPC', 'ADAMTS9-AS2')
#genes12 <- c('COL5A2', 'NTNG1', 'LINC00326', 'GRB14', 'NLGN4Y', 'MAGI2', 'ADCY8', 'PRKY', 'USP9Y', 'TTTY14', 'SYNPR', 'PCDH7', 'ADGRF5', 'KCNB2')
#genes13 <- c('COL5A2', 'NTNG1', 'DCC', 'RMST', 'RIT2', 'COBLL1', 'RORB', 'KCNH5', 'CNTNAP2', 'LRRC4C', 'GLIS3', 'GRIA1', 'SNTG1', 'LINC02296', 'UNC5D')
#genes14 <- c('XKR6', 'AGAP1', 'ATXN1', 'RORA', 'KHDRBS3', 'FAF1', 'SLC8A1', 'ATP8A2')
#genes15 <- c('PRKD1', 'LAMA2', 'PKN2-AS1', 'AL138701.2', 'SV2C', 'AC060765.2', 'SPOCK3', 'CDH18', 'NCALD', 'LINC02552')
#genes16 <- c('AC060765.1', 'MFGE8', 'AC067956.1', 'ESRRG', 'CNTN3', 'ENOX1', 'AC019211.1', 'LAMA2', 'AC110614.1', 'AC093843.1', 'KCNMB2', 'ALCAM')
#genes17 <- c('ST18', 'LINC01608', 'DOCK5', 'FRMD4B', 'LPAR1', 'RNF220', 'BCAS1', 'FA2H', 'LINC00639', 'CNDP1', 'CERCAM', 'SHROOM4', 'ZNF536', 'UGT8', 'ENPP2')
#genes18 <- c('SEMA5A', 'SORCS3', 'ZNF804A', 'CSMD1', 'CNTN5', 'AC060765.1', 'ACVR1C', 'CHSY3', 'AL359636.2', 'CDH13', 'EPHB1', 'ADAMTS19', 'FRMD4A')
#genes19 <- c('AC092448.1', 'LINC02388', 'KCNH8', 'VAV3', 'CUX1', 'AC112770.1', 'ADGRL4', 'AC020637.1', 'POU6F2', 'TAFA2', 'PLCXD2', 'ADAMTS17', 'PLEKHH2', 'LRRC4C')
#
#
#genes <- c(genes0,genes1,genes2,genes3,genes4,genes5,genes6,genes7,genes8,genes9,
#           genes10,genes11,genes12,genes13,genes14,genes15,genes16,genes17,genes18,genes19)
##LINC00507_marker.pdf 8x60
#DotPlot(seu, features=unique(genes), group.by = 'seurat_clusters') + RotatedAxis()
#
#pdf('LINC00507_cluster11_marker.pdf',width = 8,height = 6)
#DotPlot(seu, features=unique(genes11), group.by = 'seurat_clusters') + RotatedAxis()
#dev.off()
#
#
#
#
################# COMET ########################
#seu <- readRDS('seu.harmony.LINC00507.rds')
#marker <- as.matrix(seu@assays$RNA@data)
#write.table(marker, 'COMET/LINC00507_marker.txt', col.names=T, row.names=T, sep='\t',quote=F)
#
#vis <- seu@reductions$umap@cell.embeddings
#write.table(vis, 'COMET/LINC00507_vis.txt', col.names=F, row.names=T, sep='\t',quote=F)
#
#cluster <- cbind(colnames(seu),paste0('LINC00507_',seu$seurat_clusters))
#write.table(cluster, 'COMET/LINC00507_cluster.txt', col.names=F, row.names=F, sep='\t',quote=F)
#
### run COMET
### 命令行下输入
### Comet marker.txt vis.txt cluster.txt output/
#
#
#
#
#
################ cluster merge ####################################################
#seu <- readRDS('seu.harmony.LINC00507.rds')
#
#seu$seurat_clusters_merge <- as.vector(seu$seurat_clusters)
#seu$seurat_clusters_merge[seu$seurat_clusters %in% c('0','6','16')] <- '0,6,16'
#seu$seurat_clusters_merge[seu$seurat_clusters %in% c('1','3','11')] <- '1,3,11'
#seu$seurat_clusters_merge[seu$seurat_clusters %in% c('4','7','12')] <- '4,7,12'
#seu$seurat_clusters_merge[seu$seurat_clusters %in% c('8','13','19')] <- '8,13,19'
#seu$seurat_clusters_merge[seu$seurat_clusters %in% c('5','15')] <- '5,15'
#
#seu$seurat_clusters_merge <- factor(seu$seurat_clusters_merge)
#
#DimPlot(seu, group.by = 'seurat_clusters', label = T)
#DimPlot(seu, group.by = 'seurat_clusters_merge', label = T)
#
#seu@active.ident <- seu$seurat_clusters
#markers <- FindAllMarkers(seu)
#saveRDS(markers, 'seu.harmony.LINC00507.markers.rds')
#
#seu@active.ident <- seu$seurat_clusters_merge
#markers <- FindAllMarkers(seu)
#saveRDS(markers, 'seu.harmony.LINC00507.merge.markers.rds')
#
#
#markers1 <- readRDS('seu.harmony.LINC00507.markers.rds')
#markers2 <- readRDS('seu.harmony.LINC00507.merge.markers.rds')
#
#library(dplyr)
#top1 <- top_n(group_by(markers1,cluster), 5, avg_log2FC)
#genes1 <- top1$gene
#
#top2 <- top_n(group_by(markers2,cluster), 5, avg_log2FC)
#genes2 <- top2$gene
#
#DotPlot(seu, features=unique(genes1), group.by = 'seurat_clusters') + RotatedAxis()
#DotPlot(seu, features=unique(genes2), group.by = 'seurat_clusters_merge') + RotatedAxis()
#
#
#
#
markers <- readRDS('seu.harmony.LINC00507.markers.rds')
seu <- readRDS('seu.harmony.LINC00507.rds')

DimPlot(seu, group.by = 'hicat_cluster_merge', label = T)

top <- top_n(group_by(markers,cluster), 5, avg_log2FC)
genes <- top$gene

DotPlot(seu, features=unique(c(genes)), group.by = 'hicat_cluster_merge') + RotatedAxis()


markers2 <- markers[,c("cluster", "gene","avg_log2FC", "pct.1", "pct.2","p_val",  "p_val_adj")]
write.csv(markers2, 'LINC00507.markers.csv',row.names = FALSE, quote = F)


