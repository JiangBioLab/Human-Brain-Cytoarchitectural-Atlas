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


#c <- c("25","30","29","31","58","425","68","70","40","36","39","51","42","48","49","72","93","125","111","104","109",
#"98","115","123","124","88","65","77","86","89")
#
####################### XXX markers ##############################################
#clusters1 <- c("72","93","125","111","104","109",
#               "98","115","123","124")
#clusters2 <- c("25","30","29","31","58","425","68","70","40","36","39","51","42","48","49","88","65","77","86","89")
#cells1 <- colnames(seu0)[seu0$hicat_cluster_merge %in% clusters1]
#cells2 <- colnames(seu0)[seu0$hicat_cluster_merge %in% clusters2]
#
#XXX <- FindMarkers(seu0,
#                   ident.1 = clusters1,
#                   ident.2 = clusters2)
#saveRDS(XXX,'Exc_RORB_markers.rds')
#
#
#XXX<-readRDS('Exc_RORB_markers.rds')
#top5 <- XXX[XXX$avg_log2FC>0,][1:50,]
#genes <- rownames(top5)
#
##RORB_up-regulated.pdf 30x20
#DotPlot(seu0, features=unique(genes)) + RotatedAxis()

#################################################################################
#
#clusters0 <- c("446", "453", "454", "456", "463", "471", "461", 
#               "466", "469", "401", "403", "404", "413", "297", 
#               "284", "380", "342", "368", "334", "373", "339", 
#               "267", "384", "299", "327", "301", "329_333", "309", 
#               "303", "311", "377", "375", "378", "383", "268", 
#               "336", "214", "217", "220", "224", "270", "256", 
#               "235_243", "245", "264", "282", "222", "279_280", "229", 
#               "247", "232", "233", "190", "172", "173", "175", 
#               "180", "185", "201", "202", "198", "199", "192", 
#               "226", "204_211", "206", "127", "160_162", "148", "144_156", 
#               "146_157", "THEMIS_1,5", "THEMIS_7", "THEMIS_9", "THEMIS_2", "THEMIS_3", "THEMIS_10", 
#               "THEMIS_0", "THEMIS_4", "THEMIS_6", "THEMIS_8", "FEZF2.2_1", "FEZF2.2_0,3", "FEZF2.2_2", 
#               "LINC00507_18", "LINC00507_9", "LINC00507_8,13,19", "LINC00507_0,6,16", "LINC00507_2", "LINC00507_17", "LINC00507_1,3,11", 
#               "LINC00507_5,15", "LINC00507_10", "LINC00507_4,7,12", "RORB_0", "RORB_11", "RORB_16", "RORB_2,9,20", 
#               "RORB_3,6,7", "RORB_21", "RORB_22", "RORB_5", "RORB_10_19", "RORB_1,17", "RORB_23", 
#               "RORB_12", "RORB_4,8,14,18", "LINC00507_14", "RORB_13", "RORB_15", "FEZF2.4_0,9", "FEZF2.4_4", 
#               "FEZF2.4_1,2,5", "FEZF2.4_8", "FEZF2.4_6", "FEZF2.4_7", "FEZF2.5_12", "FEZF2.5_6,7", "FEZF2.5_8", 
#               "FEZF2.5_2,9,10", "FEZF2.5_0,4,5,11", "FEZF2.5_3", "FEZF2.4_3", "FEZF2.5_1", "FEZF2.3_2", "FEZF2.3_4", 
#               "FEZF2.3_8", "FEZF2.3_3", "FEZF2.3_0,5,6,7", "FEZF2.3_1", "478", "476", "477", 
#               "386", "389", "392", "395", "21", "1", "17", 
#               "19", "22")
#
#clusters <- c("RORB_0", "RORB_11", "RORB_16", "RORB_2,9,20", 
#              "RORB_3,6,7", "RORB_21", "RORB_22", "RORB_5", "RORB_10_19", "RORB_1,17", "RORB_23", 
#              "RORB_12", "RORB_4,8,14,18", "LINC00507_14", "RORB_13", "RORB_15")
#
#counts <- seu0@assays$RNA@counts[,seu0$hicat_cluster_merge %in% clusters]
#
#seu <- CreateSeuratObject(
#  counts,
#  min.cells = 0,
#  min.features = 0
#)
#dim(seu)
## 33261 24244
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
#saveRDS(seu,'seu.harmony.RORB.rds')

seu <- readRDS('seu.harmony.RORB.rds')

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

saveRDS(markers, 'seu.harmony.RORB.markers.rds')
saveRDS(seu,'seu.harmony.RORB.rds')



#################################################
markers <- readRDS('seu.harmony.RORB.markers.rds')
seu <- readRDS('seu.harmony.RORB.rds')
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
#write.csv(markers2, 'RORB.markers.csv',row.names = FALSE, quote = F)
#
#
################ K means #######################
#seu <- readRDS('seu.harmony.RORB.rds')
#
#seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:20)
#seu <- FindClusters(seu, resolution = 1)
#seu <- RunUMAP(seu, reduction = "harmony", dims = 1:20)
#
#DimPlot(seu, group.by = 'hicat_cluster')
#
#DimPlot(seu, group.by = 'region')
#DimPlot(seu, group.by = 'individual')
#
#
#genes <- c('RORB','CCDC68','PTPN3','THEMIS','ENPEP','TNNT2','LAMA4','VILL','OTOGL','LINC01202','SLC22A18','THFAIP6','LNX2','RPRM','MED8','FGF10')
#
#FeaturePlot(seu,
#            ncol = 4,
#            features=genes)
#
#data <- as.matrix(seu@assays$RNA@data)
#km_result <- kmeans(t(data), 4)
#table(km_result$cluster)
#
#seu$kmeans <- km_result$cluster[colnames(seu)]
#
#
#DimPlot(seu, group.by = 'kmeans')
#
#DimPlot(seu, group.by = 'seurat_clusters',label=T)
#
#DotPlot(seu,features = genes,group.by = 'kmeans') + RotatedAxis()
#DotPlot(seu,features = genes,group.by = 'seurat_clusters') + RotatedAxis()
#
#
#seu$kmeans <- factor(seu$kmeans)
#seu@active.ident <- seu$kmeans
#markers <- FindAllMarkers(seu)
#saveRDS(markers,'seu.harmony.RORB.kmeans.markers.rds')
#markers <- readRDS('seu.harmony.RORB.kmeans.markers.rds')
#top5 <- top_n(group_by(markers,cluster), 5, avg_log2FC)
#genes <- top5$gene
#DotPlot(seu, features=unique(genes)) + RotatedAxis()
#
#
#seu@active.ident <- seu$seurat_clusters
#markers <- FindAllMarkers(seu)
#saveRDS(markers,'seu.harmony.RORB.seuCluster.markers.rds')
#markers <- readRDS('seu.harmony.RORB.seuCluster.markers.rds')
#top5 <- top_n(group_by(markers,cluster), 5, avg_log2FC)
#genes <- top5$gene
#DotPlot(seu, features=unique(genes)) + RotatedAxis()
#
#
#seu <- readRDS('seu.harmony.RORB.rds')
#DimPlot(seu,group.by = 'seurat_clusters')
#
#genes0 <- c('RMST', 'PLCH1', 'COL22A1', 'MME', 'LINC01915', 'PCBP3', 'OTOGL', 'LCORL', 'F11-AS1', 'AC007656.2', 'AL392023.2', 'FOXO1', 'ITGA11')
#genes1 <- c('AC016687.2', 'RPRM', 'LRP1B', 'EPHA3', 'SEMA3E', 'SORCS3', 'SLC9A9', 'MAMDC2', 'CDH10', 'FSTL4', 'SORCS1')
#genes2 <- c('AL356295.1', 'TSHZ2', 'AC073578.2', 'USP9Y', 'UTY', 'RORA', 'DCC', 'GABRG1', 'AC243829.1', 'NLGN4Y', 'AL450352.1', 'MERTK', 'ALCAM')
#genes10 <- c('GALR1', 'KCNC2', 'AC124254.2', 'ATP10A', 'PDE1C', 'AC060765.2', 'CADPS2', 'GRIK2', 'LINC01680', 'PCDH7', 'XIST', 'ARHGAP29')
#genes11 <- c('LINC02196', 'AC079793.1', 'ADGRL4', 'ARHGAP15', 'MEGF10', 'CTNNA2', 'RUNX2', 'PDZRN4', 'PTPRG', 'SULF2', 'AC107208.1', 'TRIM22', 'TOX', 'AC016687.2')
#genes12 <- c('CCBE1', 'LINGO2', 'SORCS1', 'DPP10', 'AC021594.2', 'AL391117.1', 'AC099517.1', 'MIR99AHG', 'CTNNA2', 'ROBO2', 'OPCML', 'AC109466.1', 'ADTRP', 'PRUNE2', 'CNTNAP5')
#genes13 <- c('SLC22A10', 'TSHZ2', 'PCLO', 'DCC', 'LINC01725', 'IL1RAPL2', 'RARB', 'MEG3', 'ALCAM', 'COBLL1')
#genes14 <- c('CMTM8', 'GRIN3A', 'AC021134.1', 'CACNG5', 'PHACTR1', 'RBM20', 'KIRREL3', 'CELF4', 'CDH12', 'PRKY', 'AC109466.1', 'USP9Y', 'KCNH7', 'COL21A1', 'TTTY14')
#genes15 <- c('AC109466.1', 'AC008415.1', 'DPP10', 'CASC15', 'AC060765.2', 'CELF4', 'LRRK1', 'PRUNE2', 'AC073091.3', 'CAMK2D', 'CNTN3')
#genes16 <- c('NPFFR2', 'LINC02196', 'AC079380.1', 'AC099517.1', 'ADGRL4', 'OXR1', 'ADAMTS17', 'TOX', 'DACH1', 'TNNT2', 'CTNNA2', 'TENM2', 'VAV3')
#genes17 <- c('MT-ND3', 'RPRM', 'AP003066.1', 'AC016687.2', 'KCNC2', 'CADM2', 'MT-ATP6', 'NTRK3', 'MDGA2', 'MT-ND1', 'NELL1', 'SYT1', 'MT-ND4', 'MT-CO2', 'SEMA3E')
#genes18 <- c('MT-ND3', 'CACNG5', 'AL445207.1', 'AL589693.1', 'PRUNE2', 'DNAH6', 'MT-ATP6', 'CAMK2D', 'MALAT1', 'COL21A1', 'GRIN3A', 'MT-ND1', 'RAB27B')
#genes19 <- c('AL359636.2', 'CHST8', 'AC016687.2', 'CNGB3', 'GALR1', 'ASTN2', 'GRIK2', 'AC016687.3', 'PDE1C', 'GALNT17', 'SORCS3', 'LINC01378', 'SEMA5A', 'KCNC2', 'NRSN1')
#genes20 <- c('AC112770.1', 'RORA', 'BACH2', 'AC243829.1', 'LDB2', 'KCNH8', 'VAV3', 'MAGI2', 'TENM2', 'LAMA4', 'PCDH11Y', 'LOXHD1', 'WLS', 'TRPC3')
#genes21 <- c('SOX6', 'IL1RAPL2', 'CASC15', 'PARD3', 'SLC24A3', 'ERBB4', 'AC073578.2', 'SLC1A2')
#genes22 <- c('ST18', 'RNF220', 'DOCK5', 'FRMD4B', 'CERCAM', 'BCAS1', 'CASC15', 'LINC01608', 'IL1RAPL2', 'LPAR1', 'NCKAP5', 'SHROOM4', 'PLP1', 'UGT8', 'FA2H')
#
#
#genes <- c(genes0,genes1,genes2,
#           genes10,genes11,genes12,genes13,genes14,genes15,genes16,genes17,genes18,genes19,
#           genes20,genes21,genes22)
##RORB_marker.pdf 10x60
#DotPlot(seu, features=unique(genes), group.by = 'seurat_clusters') + RotatedAxis()
#
#pdf('RORB_cluster22_marker.pdf',width = 8,height = 6)
#DotPlot(seu, features=unique(genes22), group.by = 'seurat_clusters') + RotatedAxis()
#dev.off()
#
#
#
#
#
################# COMET ########################
#seu <- readRDS('seu.harmony.RORB.rds')
#marker <- as.matrix(seu@assays$RNA@data)
#write.table(marker, 'COMET/RORB_marker.txt', col.names=T, row.names=T, sep='\t',quote=F)
#
#vis <- seu@reductions$umap@cell.embeddings
#write.table(vis, 'COMET/RORB_vis.txt', col.names=F, row.names=T, sep='\t',quote=F)
#
#cluster <- cbind(colnames(seu),paste0('RORB_',seu$seurat_clusters))
#write.table(cluster, 'COMET/RORB_cluster.txt', col.names=F, row.names=F, sep='\t',quote=F)
#
### run COMET
### 命令行下输入
### Comet marker.txt vis.txt cluster.txt output/
#
#
#
############### dend Markers ##########################
#seu <- readRDS('seu.harmony.RORB.rds')
#
#norm.dat <- seu@assays$RNA@data
#
#cl.clean <- paste0('C',seu$seurat_clusters)
#names(cl.clean) <- colnames(seu)
#
#cl.clean[cl.clean=='C17'] <- 'C1'
#
#select.markers <- select_markers(norm.dat, 
#                                 cl.clean, 
#                                 #de.genes = de.genes,
#                                 n.markers = 20)
#
#marker.genes <- select.markers$markers
#de.genes <- select.markers$de.genes
#
#cl.med <- get_cl_medians(norm.dat[marker.genes,], 
#                          cl.clean)
#
## Build the dendrogram
#dend.result <- build_dend(cl.med,
#                          #l.rank, 
#                          #l.color,
#                          nboot = 100)
#saveRDS(dend.result, 'seu.harmony.RORB.dend.rds')
#
#dendNew <- select_dend_markers(dend=dend.result$dend, cl=cl.clean, de.genes = de.genes, norm.dat = norm.dat)
#
#
#
#
#select_dend_markers <- function(dend, cl, de.genes,norm.dat=NULL,up.gene.score=NULL, down.gene.score=NULL,...)
#{
#  require(dendextend)
#  require(randomForest)
#  print(dend)
#  if(length(dend)>1){
#    cl.g = sapply(dend, labels,simplify=F)
#    markers = c()
#    for(i in 1:(length(cl.g)-1)){
#      for(j in (i+1):length(cl.g)){
#        g = select_markers_pair_group(cl=cl, cl.g[[i]],cl.g[[j]],de.genes=de.genes,up.gene.score=up.gene.score, down.gene.score=down.gene.score,...)$genes
#        markers=union(markers, g)          
#      }
#    }
#    tmp.cl = setNames(rep(1:length(cl.g),sapply(cl.g,length)),unlist(cl.g))
#    select.cells = names(cl)[cl %in% names(tmp.cl)]
#    if(length(markers)==0){
#      next
#    }
#    if(!is.null(norm.dat)){
#      select.cells = unlist(tapply(select.cells, cl[select.cells],function(x)sample(x,min(length(x), 50))))
#      tmp.cl = setNames(tmp.cl[as.character(cl[select.cells])],select.cells)
#      rf = randomForest(t(as.matrix(norm.dat[markers,select.cells])),factor(tmp.cl))
#      w = importance(rf)
#      attr(dend, "markers")=w[,1]
#    }
#    else{
#      attr(dend, "markers")=setNames(rep(1,length(markers)),markers)
#    }
#    for(i in 1:length(dend)){
#      dend[[i]] = select_dend_markers(dend[[i]], cl=cl, norm.dat=norm.dat, de.genes=de.genes,up.gene.score=up.gene.score, down.gene.score=down.gene.score,...)
#    }
#  }
#  return(dend)
#}
#
#
#
#
#
#
################ cluster merge ####################################################
#seu <- readRDS('seu.harmony.RORB.rds')
#
#seu$seurat_clusters_merge <- as.vector(seu$seurat_clusters)
#seu$seurat_clusters_merge[seu$seurat_clusters %in% c('1','17')] <- '1,17'
#seu$seurat_clusters_merge[seu$seurat_clusters %in% c('10','19')] <- '10,19'
#seu$seurat_clusters_merge[seu$seurat_clusters %in% c('4','8','14','18')] <- '4,8,14,18'
#seu$seurat_clusters_merge[seu$seurat_clusters %in% c('2','9','20')] <- '2,20'
#seu$seurat_clusters_merge[seu$seurat_clusters %in% c('3','6','7')] <- '3,6,7'
#
#seu$seurat_clusters_merge <- factor(seu$seurat_clusters_merge)
#
#DimPlot(seu, group.by = 'seurat_clusters', label = T)
#DimPlot(seu, group.by = 'seurat_clusters_merge', label = T)
#
#seu@active.ident <- seu$seurat_clusters
#markers <- FindAllMarkers(seu)
#saveRDS(markers, 'seu.harmony.RORB.markers.rds')
#
#seu@active.ident <- seu$seurat_clusters_merge
#markers <- FindAllMarkers(seu)
#saveRDS(markers, 'seu.harmony.RORB.merge.markers.rds')
#
#
#markers1 <- readRDS('seu.harmony.RORB.markers.rds')
#markers2 <- readRDS('seu.harmony.RORB.merge.markers.rds')
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



markers <- readRDS('seu.harmony.RORB.markers.rds')
seu <- readRDS('seu.harmony.RORB.rds')

DimPlot(seu, group.by = 'hicat_cluster_merge', label = T)

top <- top_n(group_by(markers,cluster), 5, avg_log2FC)
genes <- top$gene

DotPlot(seu, features=unique(c('LINC00507',genes)), group.by = 'hicat_cluster_merge') + RotatedAxis()


markers2 <- markers[,c("cluster", "gene","avg_log2FC", "pct.1", "pct.2","p_val",  "p_val_adj")]
write.csv(markers2, 'RORB.markers.csv',row.names = FALSE, quote = F)

