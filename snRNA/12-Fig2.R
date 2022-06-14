library(Seurat)
setwd("/home/luomeng/brain/final_v1/data/")
################ GABA CGE ######################################################

seu <- readRDS('../data/seu.harmony.GABA.CGE.rds')

GABA <- c("86", "54", "33", "32", "109",
          "83", "19", "77", "121","46", "15", "68", "52", "78", "29", 
          "132","92", "51", "49", "44", "20", "99", "94", "87", "101",
          "75", "73", "122","120","90", "130","40", "126","45", "118")
GABA_color <- c("#C36AE6","#C263CC","#D97C80","#FF8995","#B36C76",  
                "#F2859E","#F48C9E","#BE61D4","#9E56A6","#B65FBF","#B363CB","#BD76DC","#B967D9","#DD6DF2","#A45FBF",
                "#9A5AB3","#AF6FCC","#C667D9","#D26AE6","#AC63BE","#9256A6","#7F4C8C","#DF70FF","#D270FF","#CD26FF",
                "#CD6DF2","#DA81F7","#EA86FF","#E96DF2","#CC7E09","#B98327","#C38326","#B67B16","#A66D0D","#B57014")
names(GABA_color) <- GABA

region_color <- c('#3F4587','#8562AA','#EC8561','#B97CB5','#D43046','#F0592B','#ED4A96','#593C97','#A54486','#FBDE13','#299FAE','#75CCE3','#0C6939','#0D9547')
region <- c('FPPFC','DLPFC','VLPFC','M1','S1','S1E','PoCG','SPL','SMG','AG','VISP','ITG','STG','ACC')
regionID <- c('A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','A11','A12','A13','A14')
names(region_color) <- region
names(region) <- regionID

############## UMAP ############################
seu$hicat_cluster_merge_newID <- factor(seu$hicat_cluster_merge_newID,
                                        levels = GABA)
#UMAP_GABA_CGE_cluster.pdf 5x7
DimPlot(seu,group.by = 'hicat_cluster_merge_newID',cols = GABA_color)

############## UMAP ############################
regionID <- c('A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','A11','A12','A13','A14')
regionName <- c('FPPFC','DLPFC','VLPFC','M1','S1','S1E','PoCG','SPL','SMG','AG','VISP','ITG','STG','ACC')
names(regionName) <- regionID

seu$regionName <- regionName[seu$region]

seu$regionName <- factor(seu$regionName,
                         levels=regionName)
#UMAP_GABA_CGE_region.pdf 5x6.5
DimPlot(seu,group.by = 'regionName',cols = region_color)


############## cell composition ###################
cluster_dis<-table(seu$regionName,seu$hicat_cluster_merge_newID)
cluster_dis<-cluster_dis[,GABA]
cluster_dis_r <- cluster_dis/rowSums(cluster_dis)

#GABA_CGE_region_proportion_tree.pdf 4x10
pheatmap(cluster_dis_r,cluster_cols = F)

p <- pheatmap(cluster_dis_r,cluster_cols = F)

cluster_dis_r <- as.data.frame(cluster_dis_r)
colnames(cluster_dis_r) <- c('region','cluster','ratio')
cluster_dis_r$region <- factor(cluster_dis_r$region,
                               levels = rev(p$tree_row$labels[p$tree_row$order]))

# GABA_CGE_region_proportion.pdf 4x7.05
ggplot(cluster_dis_r,aes(cluster,region),showCategory=8)+
  geom_point(aes(size=ratio,color=region))+
  scale_color_manual(values=region_color[rev(p$tree_row$labels[p$tree_row$order])])+
  scale_size(limits = c(0.0005506608,0.4))+
  theme_classic() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))

aa <- rep(1,14)
names(aa) <- region[rev(p$tree_row$labels[p$tree_row$order])]
# GABA_CGE_region_color.pdf 4x12
barplot(aa,col = region_color[rev(p$tree_row$labels[p$tree_row$order])])



####### DotPlot ######################
seu$hicat_cluster_merge_newID <- factor(seu$hicat_cluster_merge_newID,
                                        levels = GABA)
seu@active.ident <- seu$hicat_cluster_merge_newID

#genes <- c('GAD1','ADARB2',
#           'LAMP5','SV2C','SFTA3','MDGA1','ROR2','TBL1Y','TRPC3', #LAMP5
#           'PAX6','NMBR','AC091576.1','SYT10','SLCO5A1','ANGPT1','DDR2','TGFBR2','CILP2','LIX1','TNFAIP8L3','RXFP1','CHRNB3','HCRTR2', #PAX6
#           'VIP','PENK','ABI3BP','ADAM12','ENPEP','DACH2','HTR2C','SCML4','SCTR','NOX4','LINC01630','CHRNA2','MSR1' #VIP
#)


ID <- cbind(as.vector(seu$hicat_cluster_merge),as.vector(seu$hicat_cluster_merge_newID),as.vector(seu$hicat_cluster_supertypes))
ID <- unique(ID)
colnames(ID) <- c('oldID','newID','anno')
rownames(ID) <- ID[,2]
ID <- ID[GABA,]


markers <- read.csv('markers_global_v1_luomeng.csv')

library(dplyr)
top2 <- top_n(group_by(markers,clusterName), 2, X0_scaled)
top2 <- top2[top2$clusterName %in% ID[,1],]

top2_order <- c()
for(i in 1:dim(ID)[1]){
  top2_order <- rbind(top2_order, top2[top2$clusterName == ID[i,1],])
}
genes <- unique(top2_order$marker)


seu_markers_m<-seu@assays$RNA@data[rownames(seu@assays$RNA@data) %in% genes,]
seu_markers_m<-t(as.matrix(seu_markers_m))
seu_markers_m<-as.data.frame(seu_markers_m)
seu_markers_m<-cbind(seu_markers_m,as.vector(seu$hicat_cluster_merge_newID))
colnames(seu_markers_m)[dim(seu_markers_m)[2]]<-'cluster'
#seu_markers_m$cluster<-as.vector(seu_markers_m$cluster)
seu_markers_m2<-aggregate(seu_markers_m[,1:(dim(seu_markers_m)[2]-1)],by=list(cluster=seu_markers_m$cluster),FUN=mean)
seu_markers_m3<-seu_markers_m2[,-1]

seu_markers_m4<-matrix(unlist(seu_markers_m3),dim(seu_markers_m3)[1],dim(seu_markers_m3)[2])
colnames(seu_markers_m4)<-colnames(seu_markers_m3)
rownames(seu_markers_m4)<-as.vector(seu_markers_m2[,1])
seu_markers_m4<-log2(seu_markers_m4+1)
seu_markers_m4 <- seu_markers_m4[rev(GABA),]

##GABA_markers_heatmap_CGE.pdf 7.1x10
pheatmap(seu_markers_m4,
         cluster_rows=F,
         cluster_cols=T,
         angle_col = 90)

pp <- pheatmap(seu_markers_m4,
               cluster_rows=F,
               cluster_cols=T,
               angle_col = 90)

# GABA_markers_dotplot_CGE.pdf 9x15
DotPlot(seu,features = unique(pp$tree_col$labels[pp$tree_col$order]), assay = 'RNA', cols = c("yellow", "red")) + theme(axis.text.x = element_text(angle = 90, hjust = 1))


############## Chi-Test ########################
a <- table(seu$hicat_cluster_merge_newID,seu$regionName)

chisq <- chisq.test(a)
chisq
#X-squared = 1147.2, df = 364, p-value < 2.2e-16
#
chisq$observed
round(chisq$expected,2)
round(chisq$residuals, 3)

library(corrplot)
library(grDevices)
col2 <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                               "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                               "#4393C3", "#2166AC", "#053061")))
#CGE_cluster_region_chiTest_corrplot.pdf 8x5
corrplot(chisq$residuals, is.cor = FALSE, col =  col2(200))



############# geom_jitter #####################
tmp  <- data.frame(cells = colnames(seu),
                   cluster = seu$hicat_cluster_merge_newID,
                   region = seu$regionName)
tmp$cluster <- factor(tmp$cluster, levels = GABA)
tmp$region <- factor(tmp$region, levels = rev(p$tree_row$labels[p$tree_row$order]))
# CGE_geom_point.pdf 4x8
ggplot(tmp, aes(x=factor(cluster), y=region, fill=factor(cluster)) ) + 
  geom_jitter(position=position_jitter(width=0.3, height=0.4), size = 0.5, aes(colour=factor(cluster)), alpha=0.9)+
  scale_color_manual(values=GABA_color)+
  theme_classic()



############ subclasses re-clustering CGE ########################
seu <- readRDS('seu.harmony.GABA.CGE.rds')

table(seu$hicat_cluster_subclasses)

seu$hicat_cluster_subclasses <- factor(seu$hicat_cluster_subclasses,
                                       levels = c('LAMP5/SV2C','PAX6','LAMP5/SNCG','VIP'))

subclasss_color <- c('#f2798d','#E96DF2','#D96C07','#9C28CC')

## Inh_CGE_subclasses.pdf 6x8
DimPlot(seu, group.by = 'hicat_cluster_subclasses', cols = subclasss_color)


region_color <- c('#3F4587','#8562AA','#EC8561','#B97CB5','#D43046','#F0592B','#ED4A96','#593C97','#A54486','#FBDE13','#299FAE','#75CCE3','#0C6939','#0D9547')
region <- c('FPPFC','DLPFC','VLPFC','M1','S1','S1E','PoCG','SPL','SMG','AG','VISP','ITG','STG','ACC')
regionID <- c('A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','A11','A12','A13','A14')
names(region_color) <- region
names(region) <- regionID

seu$regionName <- region[seu$region]
seu$regionName <- factor(seu$regionName, levels = region)

## Inh_CGE_region.pdf 6x7
DimPlot(seu, group.by = 'regionName', cols = region_color)

library(harmony)
seu.list <- SplitObject(seu, split.by = "hicat_cluster_subclasses")
for (i in 1:length(seu.list)) {
  seu.list[[i]] <- NormalizeData(object = seu.list[[i]], normalization.method = "LogNormalize", scale.factor = 1e4)
  seu.list[[i]] <- FindVariableFeatures(seu.list[[i]], selection.method = "vst", nfeatures = 2000)
  
  all.genes <- rownames(seu.list[[i]])
  seu.list[[i]] <- ScaleData(seu.list[[i]], features = all.genes)
  seu.list[[i]] <- RunPCA(seu.list[[i]], features = VariableFeatures(object = seu.list[[i]]))
  
  seu.list[[i]] <- RunHarmony(seu.list[[i]], "sample")
  
  seu.list[[i]] <- FindNeighbors(seu.list[[i]], reduction = "harmony", dims = 1:20)
  seu.list[[i]] <- FindClusters(seu.list[[i]], resolution = 1)
  seu.list[[i]] <- RunUMAP(seu.list[[i]], reduction = "harmony", dims = 1:20)
}

saveRDS(seu.list, 'seu.harmony.GABA.CGE.list.rds')



GABA <- c("86", "54", "33", "32", "109",
          "83", "19", "77", "121","46", "15", "68", "52", "78", "29", 
          "132","92", "51", "49", "44", "20", "99", "94", "87", "101",
          "75", "73", "122","120","90", "130","40", "126","45", "118")
GABA_color <- c("#C36AE6","#C263CC","#D97C80","#FF8995","#B36C76",  
                "#F2859E","#F48C9E","#BE61D4","#9E56A6","#B65FBF","#B363CB","#BD76DC","#B967D9","#DD6DF2","#A45FBF",
                "#9A5AB3","#AF6FCC","#C667D9","#D26AE6","#AC63BE","#9256A6","#7F4C8C","#DF70FF","#D270FF","#CD26FF",
                "#CD6DF2","#DA81F7","#EA86FF","#E96DF2","#CC7E09","#B98327","#C38326","#B67B16","#A66D0D","#B57014")
names(GABA_color) <- GABA

seu.list <- readRDS('seu.harmony.GABA.CGE.list.rds')

subclasses <- c('LAMP5/SV2C','PAX6','LAMP5/SNCG','VIP')

for(i in 1:length(subclasses)){
  
  #this_col <- GABA_color[levels(seu.list[[subclasses[i]]]$hicat_cluster_merge_newID)]
  DimPlot(seu.list[[subclasses[i]]],group.by = 'hicat_cluster_merge_newID') + ggtitle(subclasses[i])
  ggsave(path = "UMAP", filename=paste0('Inh_CGE_clusters_',i,'.pdf'), width=5, height=4)
  
  x1 <- min(seu.list[[subclasses[i]]]@reductions$umap@cell.embeddings[,1])
  x2 <- max(seu.list[[subclasses[i]]]@reductions$umap@cell.embeddings[,1])
  y1 <- min(seu.list[[subclasses[i]]]@reductions$umap@cell.embeddings[,2])
  y2 <- max(seu.list[[subclasses[i]]]@reductions$umap@cell.embeddings[,2])
  
  x1 <- x1 - (x2-x1)/20
  x2 <- x2 + (x2-x1)/20
  y1 <- y1 - (y2-y1)/20
  y2 <- y2 + (y2-y1)/20
  
  DimPlot(seu.list[[subclasses[i]]], group.by = 'regionName', cols = region_color) + ggtitle(subclasses[i]) + xlim(x1,x2) + ylim(y1,y2)
  ggsave(path = "UMAP", filename=paste0('Inh_CGE_subclasses_',i,'.pdf'), width=5, height=4)
  
  mm <- data.frame(cells = colnames(seu.list[[subclasses[i]]]),
                   UMAP_1 = seu.list[[subclasses[i]]]@reductions$umap@cell.embeddings[,1],
                   UMAP_2 = seu.list[[subclasses[i]]]@reductions$umap@cell.embeddings[,2],
                   cluster = seu.list[[subclasses[i]]]$hicat_cluster_merge_newID,
                   regionName = seu.list[[subclasses[i]]]$regionName)
  
  cc <- unique(seu.list[[subclasses[i]]]$hicat_cluster_merge_newID)
  for(j in 1:length(cc)){
    
    #cells_to_highlight <- list()
    #for(k in 1:length(names(region_color))){
    #  cells_to_highlight[[names(region_color)[k]]] <- colnames(seu.list[[subclasses[i]]])[seu.list[[subclasses[i]]]$hicat_cluster_merge_newID == cc[j] &
    #                                                                                        seu.list[[subclasses[i]]]$regionName == names(region_color)[k]]
    #}
    #
    #DimPlot(seu.list[[subclasses[i]]], cells.highlight = cells_to_highlight, group.by = 'regionName', cols.highlight = region_color) + ggtitle(paste0('cl: ',cc[j])) + xlim(x1,x2) + ylim(y1,y2)
    
    mm$color <- as.vector(mm$regionName)
    mm$color[!(mm$cluster == cc[j])] <- 'Unselected'
    
    mm$color <- factor(mm$color, levels = c(names(region_color),'Unselected'))
    mm_color <- c(region_color, 'grey70')
    names(mm_color) <- c(names(region_color),'Unselected')
    
    mm$size <- rep(2,dim(mm)[1])
    mm$size[!(mm$cluster == cc[j])] <- 1
    
    mm <- rbind(mm[mm$color != 'Unselected',],
                mm[mm$color == 'Unselected',])
    
    ggplot(mm, aes(x = UMAP_1, y = UMAP_2, colour = color, size = size)) +
      geom_point() + #scale_color_tableau() + 
      scale_color_manual(values = mm_color)+
      scale_size_continuous( range= c(0.1,1) )+
      theme_classic() + 
      ggtitle(paste0('cl: ',cc[j])) + 
      theme(plot.title = element_text(hjust = 0.5)) +
      xlim(x1,x2) + 
      ylim(y1,y2)
    ggsave(path = "UMAP", filename=paste0('Inh_CGE_subclasses_',i,'_',j,'.pdf'), width=4.9, height=3.65) 
  }
}











#### To XuChang ####################
devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")
remotes::install_github("mojaveazure/seurat-disk")
library(Seurat)
library(hdf5r)
library(loomR)
library(scater)
library(SeuratDisk)



seu0 <- readRDS('seu.harmony.GABA.CGE.rds')
seu <- CreateSeuratObject(
  seu0@assays$RNA@counts,
  min.cells = 0,
  min.features = 0
)
seu$hicat_cluster_merge_newID <- seu0$hicat_cluster_merge_newID[colnames(seu)]
seu <- NormalizeData(object = seu, normalization.method = "LogNormalize", scale.factor = 1e4)
pan.loom <- as.loom(x = seu, filename = "CGE.loom", verbose = FALSE)
pan.loom$close_all()







#####################################################################################





#################### GABA MGE #######################################################

seu <- readRDS('seu.harmony.GABA.MGE.rds')

GABA <- c("67", "100","97", "102","25", "24", "110","55", "112","58",     ##MGE
          "21", "28", "37", "135","141","23", "64", "105","36", "50", 
          "88", "113","8")
GABA_color <- c("#F27503","#E67105","#EC820D","#E69A05","#D99207","#FF7A00","#CE8024","#BF670B","#C36F1E","#DC7C16",     ##MGE
                "#CC7209","#E89420","#FFA900","#AB5B42","#D0494E","#B1614A","#D92D43","#EC2D46","#8C2A27","#873C46",
                "#B9342C","#D93031","#B92C3E")
names(GABA_color) <- GABA

region_color <- c('#3F4587','#8562AA','#EC8561','#B97CB5','#D43046','#F0592B','#ED4A96','#593C97','#A54486','#FBDE13','#299FAE','#75CCE3','#0C6939','#0D9547')
region <- c('FPPFC','DLPFC','VLPFC','M1','S1','S1E','PoCG','SPL','SMG','AG','VISP','ITG','STG','ACC')
regionID <- c('A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','A11','A12','A13','A14')
names(region_color) <- region
names(region) <- regionID

############## UMAP ############################
seu$hicat_cluster_merge_newID <- factor(seu$hicat_cluster_merge_newID,
                                        levels = GABA)
#UMAP_GABA_MGE_cluster.pdf 5x7
DimPlot(seu,group.by = 'hicat_cluster_merge_newID',cols = GABA_color)

############## UMAP ############################
regionID <- c('A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','A11','A12','A13','A14')
regionName <- c('FPPFC','DLPFC','VLPFC','M1','S1','S1E','PoCG','SPL','SMG','AG','VISP','ITG','STG','ACC')
names(regionName) <- regionID

seu$regionName <- regionName[seu$region]

seu$regionName <- factor(seu$regionName,
                         levels=regionName)
#UMAP_GABA_MGE_region.pdf 5x6.5
DimPlot(seu,group.by = 'regionName',cols = region_color)


############## cell composition ###################
cluster_dis<-table(seu$regionName,seu$hicat_cluster_merge_newID)
cluster_dis<-cluster_dis[,GABA]
cluster_dis_r <- cluster_dis/rowSums(cluster_dis)

#GABA_MGE_region_proportion_tree.pdf 4x10
pheatmap(cluster_dis_r,cluster_cols = F)

p <- pheatmap(cluster_dis_r,cluster_cols = F)

cluster_dis_r <- as.data.frame(cluster_dis_r)
colnames(cluster_dis_r) <- c('region','cluster','ratio')
cluster_dis_r$region <- factor(cluster_dis_r$region,
                               levels = rev(p$tree_row$labels[p$tree_row$order]))

# GABA_MGE_region_proportion.pdf 4x6
ggplot(cluster_dis_r,aes(cluster,region),showCategory=8)+
  geom_point(aes(size=ratio,color=region))+
  scale_color_manual(values=region_color[rev(p$tree_row$labels[p$tree_row$order])])+
  scale_size(limits = c(0.0005202914,0.4))+
  theme_classic() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))

aa <- rep(1,14)
names(aa) <- region[rev(p$tree_row$labels[p$tree_row$order])]
# GABA_MGE_region_color.pdf 4x12
barplot(aa,col = region_color[rev(p$tree_row$labels[p$tree_row$order])])



####### DotPlot ######################
seu$hicat_cluster_merge_newID <- factor(seu$hicat_cluster_merge_newID,
                                        levels = GABA)
seu@active.ident <- seu$hicat_cluster_merge_newID

#genes <- c('GAD1','ADARB2',
#           'LAMP5','SV2C','SFTA3','MDGA1','ROR2','TBL1Y','TRPC3', #LAMP5
#           'PAX6','NMBR','AC091576.1','SYT10','SLCO5A1','ANGPT1','DDR2','TGFBR2','CILP2','LIX1','TNFAIP8L3','RXFP1','CHRNB3','HCRTR2', #PAX6
#           'VIP','PENK','ABI3BP','ADAM12','ENPEP','DACH2','HTR2C','SCML4','SCTR','NOX4','LINC01630','CHRNA2','MSR1' #VIP
#)


ID <- cbind(as.vector(seu$hicat_cluster_merge),as.vector(seu$hicat_cluster_merge_newID),as.vector(seu$hicat_cluster_supertypes))
ID <- unique(ID)
colnames(ID) <- c('oldID','newID','anno')
rownames(ID) <- ID[,2]
ID <- ID[GABA,]


markers <- read.csv('markers_global_v1_luomeng.csv')

library(dplyr)
top2 <- top_n(group_by(markers,clusterName), 2, X0_scaled)
top2 <- top2[top2$clusterName %in% ID[,1],]

top2_order <- c()
for(i in 1:dim(ID)[1]){
  top2_order <- rbind(top2_order, top2[top2$clusterName == ID[i,1],])
}
genes <- unique(top2_order$marker)


seu_markers_m<-seu@assays$RNA@data[rownames(seu@assays$RNA@data) %in% genes,]
seu_markers_m<-t(as.matrix(seu_markers_m))
seu_markers_m<-as.data.frame(seu_markers_m)
seu_markers_m<-cbind(seu_markers_m,as.vector(seu$hicat_cluster_merge_newID))
colnames(seu_markers_m)[dim(seu_markers_m)[2]]<-'cluster'
#seu_markers_m$cluster<-as.vector(seu_markers_m$cluster)
seu_markers_m2<-aggregate(seu_markers_m[,1:(dim(seu_markers_m)[2]-1)],by=list(cluster=seu_markers_m$cluster),FUN=mean)
seu_markers_m3<-seu_markers_m2[,-1]

seu_markers_m4<-matrix(unlist(seu_markers_m3),dim(seu_markers_m3)[1],dim(seu_markers_m3)[2])
colnames(seu_markers_m4)<-colnames(seu_markers_m3)
rownames(seu_markers_m4)<-as.vector(seu_markers_m2[,1])
seu_markers_m4<-log2(seu_markers_m4+1)
seu_markers_m4 <- seu_markers_m4[rev(GABA),]

##GABA_markers_heatmap_MGE.pdf 7x7
pheatmap(seu_markers_m4,
         cluster_rows=F,
         cluster_cols=T,
         angle_col = 90)

pp <- pheatmap(seu_markers_m4,
               cluster_rows=F,
               cluster_cols=T,
               angle_col = 90)

# GABA_markers_dotplot_MGE.pdf 7x10
DotPlot(seu,features = unique(pp$tree_col$labels[pp$tree_col$order]), assay = 'RNA', cols = c("yellow", "red")) + theme(axis.text.x = element_text(angle = 90, hjust = 1))


#### Vln plot 

VlnPlot(seu, features = genes, group.by = 'hicat_cluster_merge_newID')

## GABA_MGE_vlnplot.pdf 40x10
StackedVlnPlot(obj = seu, 
               features = unique(genes),
               group.by = 'hicat_cluster_merge_newID')

modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}




############## Chi-Test ########################
a <- table(seu$hicat_cluster_merge_newID,seu$regionName)

chisq <- chisq.test(a)
chisq
#X-squared = 1147.2, df = 364, p-value < 2.2e-16
#
chisq$observed
round(chisq$expected,2)
round(chisq$residuals, 3)

library(corrplot)
library(grDevices)
col2 <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                               "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                               "#4393C3", "#2166AC", "#053061")))
#MGE_cluster_region_chiTest_corrplot.pdf 8x5
corrplot(chisq$residuals, is.cor = FALSE, col =  col2(200))



############# geom_jitter #####################
tmp  <- data.frame(cells = colnames(seu),
                   cluster = seu$hicat_cluster_merge_newID,
                   region = seu$regionName)
tmp$cluster <- factor(tmp$cluster, levels = GABA)
tmp$region <- factor(tmp$region, levels = rev(p$tree_row$labels[p$tree_row$order]))
# MGE_geom_point.pdf 3x6
ggplot(tmp, aes(x=factor(cluster), y=region, fill=factor(cluster)) ) + 
  geom_jitter(position=position_jitter(width=0.3, height=0.4), size = 0.5, aes(colour=factor(cluster)), alpha=0.9)+
  scale_color_manual(values=GABA_color)+
  theme_classic()




############ subclasses re-clustering MGE ########################
seu <- readRDS('seu.harmony.GABA.MGE.rds')

table(seu$hicat_cluster_subclasses)

seu$hicat_cluster_subclasses <- factor(seu$hicat_cluster_subclasses,
                                       levels = c('SST','PVALB'))

subclasss_color <- c('#F2B003','#FF2D4E')

## Inh_MGE_subclasses.pdf 6x7
DimPlot(seu, group.by = 'hicat_cluster_subclasses', cols = subclasss_color)


region_color <- c('#3F4587','#8562AA','#EC8561','#B97CB5','#D43046','#F0592B','#ED4A96','#593C97','#A54486','#FBDE13','#299FAE','#75CCE3','#0C6939','#0D9547')
region <- c('FPPFC','DLPFC','VLPFC','M1','S1','S1E','PoCG','SPL','SMG','AG','VISP','ITG','STG','ACC')
regionID <- c('A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','A11','A12','A13','A14')
names(region_color) <- region
names(region) <- regionID

seu$regionName <- region[seu$region]
seu$regionName <- factor(seu$regionName, levels = region)

## Inh_MGE_region.pdf 6x7
DimPlot(seu, group.by = 'regionName', cols = region_color)

library(harmony)
seu.list <- SplitObject(seu, split.by = "hicat_cluster_subclasses")
for (i in 1:length(seu.list)) {
  seu.list[[i]] <- NormalizeData(object = seu.list[[i]], normalization.method = "LogNormalize", scale.factor = 1e4)
  seu.list[[i]] <- FindVariableFeatures(seu.list[[i]], selection.method = "vst", nfeatures = 2000)
  
  all.genes <- rownames(seu.list[[i]])
  seu.list[[i]] <- ScaleData(seu.list[[i]], features = all.genes)
  seu.list[[i]] <- RunPCA(seu.list[[i]], features = VariableFeatures(object = seu.list[[i]]))
  
  seu.list[[i]] <- RunHarmony(seu.list[[i]], "sample")
  
  seu.list[[i]] <- FindNeighbors(seu.list[[i]], reduction = "harmony", dims = 1:20)
  seu.list[[i]] <- FindClusters(seu.list[[i]], resolution = 1)
  seu.list[[i]] <- RunUMAP(seu.list[[i]], reduction = "harmony", dims = 1:20)
}

saveRDS(seu.list, 'seu.harmony.GABA.MGE.list.rds')


seu.list <- readRDS('seu.harmony.GABA.MGE.list.rds')

subclasses <- c('SST','PVALB')

for(i in 1:length(subclasses)){
  
  DimPlot(seu.list[[subclasses[i]]],group.by = 'hicat_cluster_merge_newID') + ggtitle(subclasses[i])
  ggsave(path = "UMAP", filename=paste0('Inh_MGE_clusters_',i,'.pdf'), width=5, height=4)
  
  x1 <- min(seu.list[[subclasses[i]]]@reductions$umap@cell.embeddings[,1])
  x2 <- max(seu.list[[subclasses[i]]]@reductions$umap@cell.embeddings[,1])
  y1 <- min(seu.list[[subclasses[i]]]@reductions$umap@cell.embeddings[,2])
  y2 <- max(seu.list[[subclasses[i]]]@reductions$umap@cell.embeddings[,2])
  
  x1 <- x1 - (x2-x1)/20
  x2 <- x2 + (x2-x1)/20
  y1 <- y1 - (y2-y1)/20
  y2 <- y2 + (y2-y1)/20
  
  DimPlot(seu.list[[subclasses[i]]], group.by = 'regionName', cols = region_color) + ggtitle(subclasses[i]) + xlim(x1,x2) + ylim(y1,y2)
  ggsave(path = "UMAP", filename=paste0('Inh_MGE_subclasses_',i,'.pdf'), width=5, height=4)

  mm <- data.frame(cells = colnames(seu.list[[subclasses[i]]]),
                   UMAP_1 = seu.list[[subclasses[i]]]@reductions$umap@cell.embeddings[,1],
                   UMAP_2 = seu.list[[subclasses[i]]]@reductions$umap@cell.embeddings[,2],
                   cluster = seu.list[[subclasses[i]]]$hicat_cluster_merge_newID,
                   regionName = seu.list[[subclasses[i]]]$regionName)
  
  cc <- unique(seu.list[[subclasses[i]]]$hicat_cluster_merge_newID)
  for(j in 1:length(cc)){
    
    #cells_to_highlight <- list()
    #for(k in 1:length(names(region_color))){
    #  cells_to_highlight[[names(region_color)[k]]] <- colnames(seu.list[[subclasses[i]]])[seu.list[[subclasses[i]]]$hicat_cluster_merge_newID == cc[j] &
    #                                                                                        seu.list[[subclasses[i]]]$regionName == names(region_color)[k]]
    #}
    #DimPlot(seu.list[[subclasses[i]]], cells.highlight = cells_to_highlight, cols.highlight = region_color) + ggtitle(paste0('cl: ',cc[j])) + xlim(x1,x2) + ylim(y1,y2)
    
    mm$color <- as.vector(mm$regionName)
    mm$color[!(mm$cluster == cc[j])] <- 'Unselected'
    
    mm$color <- factor(mm$color, levels = c(names(region_color),'Unselected'))
    mm_color <- c(region_color, 'grey70')
    names(mm_color) <- c(names(region_color),'Unselected')
    
    mm$size <- rep(2,dim(mm)[1])
    mm$size[!(mm$cluster == cc[j])] <- 1
    
    ggplot(mm, aes(x = UMAP_1, y = UMAP_2, colour = color, size = size)) +
      geom_point() + #scale_color_tableau() + 
      scale_color_manual(values = mm_color)+
      scale_size_continuous( range= c(0.1,1) )+
      theme_classic() + 
      ggtitle(paste0('cl: ',cc[j])) + 
      theme(plot.title = element_text(hjust = 0.5)) +
      xlim(x1,x2) + 
      ylim(y1,y2)
    ggsave(path = "UMAP", filename=paste0('Inh_MGE_subclasses_',i,'_',j,'.pdf'), width=4.9, height=3.65) 
  }
}









