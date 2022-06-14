library(Seurat)

#seu0 <- readRDS('seu.harmony.anno.rds')
seu <- readRDS('seu.harmony.Glu.rds')


#################### Glu #######################################################

seu <- readRDS('seu.harmony.Glu.rds')

Glu <- c("111","119","140","116","104","60", "136",                      ##Glu
         "128","139","123","85", "134","76", "138","93", "91", "131",
         "133","98", "80", "107","142","72", "114","63", "81", "53", 
         "124","115","137","61", "106","74", "22", "47", "43", "71", 
         "16", "11", "10", "5",  "38", "3",  "9",  "7",  "17", "26", 
         "34", "18", "30", "108","13", "41", "48", "39", "65", "62", 
         "12", "14")
Glu_color <- c("#93C43B","#AFE32A","#B1E67B","#C0F27F","#A4D374","#88A662","#8D8C20",                                   ##Glu
               "#898325","#A7A322","#B0CE1F","#B1B10C","#DABE23","#CEC823","#C5A93D","#91910E","#C8B323","#0E8B8B",
               "#09EEDC","#19E3BE","#2B8AA5","#7944AA","#2FBCE5","#2E84BE","#2E9FCB","#296B98","#2C93B2","#2B77A5",
               "#69419D","#7E47B6","#5E419D","#71419D","#7643A4","#6444AA","#0BBEBE","#52ADC0","#4A9E9E","#18C5C8",
               "#02E2F9","#09CCC6","#07C6D9","#2E9FBB","#378695","#138C98","#0DA3A3","#0E9797","#07D8D8","#5BDFDC",
               "#637384","#5B7893","#61867A","#5991A4","#4EA8AC","#41AEAE","#4E9F8C","#99FFCC","#52B8AA","#4E9E9E",
               "#55C5B5","#4BBCAD")
names(Glu_color) <- Glu

region_color <- c('#3F4587','#8562AA','#EC8561','#B97CB5','#D43046','#F0592B','#ED4A96','#593C97','#A54486','#FBDE13','#299FAE','#75CCE3','#0C6939','#0D9547')
region <- c('FPPFC','DLPFC','VLPFC','M1','S1','S1E','PoCG','SPL','SMG','AG','VISP','ITG','STG','ACC')
regionID <- c('A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','A11','A12','A13','A14')
names(region_color) <- region
names(region) <- regionID

############## UMAP ############################
seu$hicat_cluster_merge_newID <- factor(seu$hicat_cluster_merge_newID,
                                        levels = Glu)
#UMAP_Glu_cluster.pdf 5x8
DimPlot(seu,group.by = 'hicat_cluster_merge_newID',cols = Glu_color)

############## UMAP ############################
regionID <- c('A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','A11','A12','A13','A14')
regionName <- c('FPPFC','DLPFC','VLPFC','M1','S1','S1E','PoCG','SPL','SMG','AG','VISP','ITG','STG','ACC')
names(regionName) <- regionID

seu$regionName <- regionName[seu$region]

seu$regionName <- factor(seu$regionName,
                         levels=regionName)
#UMAP_Glu_region.pdf 5x6.5
DimPlot(seu,group.by = 'regionName',cols = region_color)


############## cell composition ###################
cluster_dis<-table(seu$regionName,seu$hicat_cluster_merge_newID)
cluster_dis<-cluster_dis[,Glu]
cluster_dis_r <- cluster_dis/rowSums(cluster_dis)

#Glu_region_proportion_tree.pdf 4x15
pheatmap(cluster_dis_r,cluster_cols = F)

p <- pheatmap(cluster_dis_r,cluster_cols = F)

cluster_dis_r <- as.data.frame(cluster_dis_r)
colnames(cluster_dis_r) <- c('region','cluster','ratio')
cluster_dis_r$region <- factor(cluster_dis_r$region,
                               levels = rev(p$tree_row$labels[p$tree_row$order]))

# Glu_region_proportion.pdf 3.5x10
ggplot(cluster_dis_r,aes(cluster,region),showCategory=8)+
  geom_point(aes(size=ratio,color=region))+
  scale_color_manual(values=region_color[rev(p$tree_row$labels[p$tree_row$order])])+
  scale_size(limits = c(0.000107032,0.3))+
  theme_classic() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))

aa <- rep(1,14)
names(aa) <- region[rev(p$tree_row$labels[p$tree_row$order])]
# Glu_region_color.pdf 4x12
barplot(aa,col = region_color[rev(p$tree_row$labels[p$tree_row$order])])



####### DotPlot ######################
seu$hicat_cluster_merge_newID <- factor(seu$hicat_cluster_merge_newID,
                                        levels = Glu)
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
ID <- ID[Glu,]


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
seu_markers_m4 <- seu_markers_m4[rev(Glu),]

##Glu_markers_heatmap.pdf 7x7
pheatmap(seu_markers_m4,
         cluster_rows=F,
         cluster_cols=T,
         angle_col = 90)

pp <- pheatmap(seu_markers_m4,
               cluster_rows=F,
               cluster_cols=T,
               angle_col = 90)

# Glu_markers_dotplot.pdf 14x24
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
#Glu_cluster_region_chiTest_corrplot.pdf 16x10
corrplot(chisq$residuals, is.cor = FALSE, col =  col2(200))



############# geom_jitter #####################
tmp  <- data.frame(cells = colnames(seu),
                   cluster = seu$hicat_cluster_merge_newID,
                   region = seu$regionName)
tmp$cluster <- factor(tmp$cluster, levels = Glu)
tmp$region <- factor(tmp$region, levels = rev(p$tree_row$labels[p$tree_row$order]))
# Glu_geom_point.pdf 3x12
ggplot(tmp, aes(x=factor(cluster), y=region, fill=factor(cluster)) ) + 
  geom_jitter(position=position_jitter(width=0.3, height=0.4), size = 0.5, aes(colour=factor(cluster)), alpha=0.9)+
  scale_color_manual(values=Glu_color)+
  theme_classic()








#######################################################



Glu <- c('44','46','92','67','30',
         '35','23','36','65','61','12','1','18','6','21',
         '7','11','63','19')
Glu_color <- c('#299373','#2FD3B0','#2CAD7A','#2FCA96','#8AA3CC',
               '#90A7D9','#30E6BA','#2BA087','#268063','#9EE61D','#BCFF1A','#8DD31E','#02F970','#0E993B','#0BBF45',
               '#13A23E','#00FF34','#07D945','#09CC4F')
names(Glu_color) <- Glu

a<-rep(1,19)
names(a) <- Glu

#Glu_cluster_color.pdf 3x10
barplot(a,col = Glu_color)





################################################################
seu <- readRDS('seu.harmony.anno.de80.rds')

#result <- readRDS('seu.hicat.run_consensus_clust.de20.niter100.de80Exc.rds')
#dend.result <- readRDS('seu.hicat.run_consensus_clust.de20.niter100.de80Exc.dend.rds')

result <- readRDS('seu.hicat.run_consensus_clust.de40.niter100.de80Exc.rds')
dend.result <- readRDS('seu.hicat.run_consensus_clust.de40.niter100.de80Exc.dend.rds')

#dend_Exc_de20.pdf 10x3
#dend_Exc_de40.pdf 10x3
plot(dend.result$dend,horiz = T)

seu$tmp <- result$cl.result$cl[colnames(seu)]
seu$tmp[is.na(seu$tmp)]<-'unknown'
seu$tmp <- factor(seu$tmp,
                  levels = c(labels(dend.result$dend),'unknown'))

cluster_dis<-cbind(as.vector(seu$tmp),as.vector(seu$individual))
colnames(cluster_dis)<-c('cluster','individual')
cluster_dis<-as.data.frame(cluster_dis)
cluster_dis<-cluster_dis[!cluster_dis[,1]=='unknown',]
cluster_dis$cluster<-factor(cluster_dis$cluster,
                            levels=rev(labels(dend.result$dend)))

library(ggplot2)
library(RColorBrewer)
##cluster_individual_composition_Exc_de20.pdf
##cluster_individual_composition_Exc_de40.pdf
##pdf 3x20
ggplot(cluster_dis, aes(cluster)) + 
  geom_bar(aes(fill=individual), position="fill", color="gray20",alpha = 0.7) +#
  scale_fill_manual(values=brewer.pal(5,"Dark2"))+ #values=colors[c(1:4,6)]
  ylab('Percentage')+
  theme(axis.line = element_line(colour = "gray10"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + RotatedAxis()

genes <- c('SLC17A7',
           'FEZF2','ATP6V1C2','LINC00922','SGCG','SYT6','TSPAN18','SCUBE1','MDFIC','TNIP3','SEMA3D',
           'LINC00507','ACVR1C','THEMIS',
           'RORB','TRABD2A','GALR1','CRISPLD2','COL22A1','COL22A1','COL5A2','LINC02388','AC079793.1','AL355672.1','AL356295.1'
)

#marker_Exc_de20.pdf 10x10
#marker_Exc_de40.pdf 10x10
DotPlot(seu,features = unique(genes),group.by = 'tmp') + RotatedAxis()










######################## UMAP #############################################

table(seu$hicat_cluster_subclasses)
seu$hicat_cluster_subclasses[seu$hicat_cluster_merge_newID=='107'] <- 'L6B/FEZF2'
seu$hicat_cluster_subclasses[seu$hicat_cluster_merge_newID=='80'] <- 'L6CT/FEZF2'

table(seu$hicat_cluster_subclasses)
# ET L1-L3 IT LINC00507         L5 IT/RORB          L6B/FEZF2            L6 CAR3         L6CT/FEZF2              L6/IT                 NP 
#196              59481              24244                923                 90               1563                802                739



seu$hicat_cluster_subclasses <- factor(seu$hicat_cluster_subclasses,
                                       levels = c('NP', 'L6CT/FEZF2', 'L6B/FEZF2', 'L5 IT/RORB', 'L1-L3 IT LINC00507', 'ET', 'L6 CAR3', 'L6/IT'))

subclasss_color <- c('#93C43B', '#2D8CB8', '#7944AA', '#09B2B2', '#07D8D8', '#19E3BE', '#898325', '#CEC823')

## Exc_subclasses.pdf 6x8
DimPlot(seu, group.by = 'hicat_cluster_subclasses', cols = subclasss_color)




region_color <- c('#3F4587','#8562AA','#EC8561','#B97CB5','#D43046','#F0592B','#ED4A96','#593C97','#A54486','#FBDE13','#299FAE','#75CCE3','#0C6939','#0D9547')
region <- c('FPPFC','DLPFC','VLPFC','M1','S1','S1E','PoCG','SPL','SMG','AG','VISP','ITG','STG','ACC')
regionID <- c('A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','A11','A12','A13','A14')
names(region_color) <- region
names(region) <- regionID

seu$regionName <- region[seu$region]
seu$regionName <- factor(seu$regionName, levels = region)

## Exc_region.pdf 6x7
DimPlot(seu, group.by = 'regionName', cols = region_color)

saveRDS(seu, 'seu.harmony.Glu.rds')


############ subclasses re-clustering ########################
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

saveRDS(seu.list, 'seu.harmony.Glu.list.rds')

seu.list <- readRDS('seu.harmony.Glu.list.rds')

subclasses <- c('NP', 'L6CT/FEZF2', 'L6B/FEZF2', 'L5 IT/RORB', 'L1-L3 IT LINC00507', 'ET', 'L6 CAR3', 'L6/IT')

for(i in 1:length(subclasses)){
  
  DimPlot(seu.list[[subclasses[i]]],group.by = 'hicat_cluster_merge_newID') + ggtitle(subclasses[i])
  ggsave(path = "UMAP", filename=paste0('Exc_clusters_',i,'.pdf'), width=5, height=4)
  
  x1 <- min(seu.list[[subclasses[i]]]@reductions$umap@cell.embeddings[,1])
  x2 <- max(seu.list[[subclasses[i]]]@reductions$umap@cell.embeddings[,1])
  y1 <- min(seu.list[[subclasses[i]]]@reductions$umap@cell.embeddings[,2])
  y2 <- max(seu.list[[subclasses[i]]]@reductions$umap@cell.embeddings[,2])
  
  x1 <- x1 - (x2-x1)/20
  x2 <- x2 + (x2-x1)/20
  y1 <- y1 - (y2-y1)/20
  y2 <- y2 + (y2-y1)/20
  
  DimPlot(seu.list[[subclasses[i]]], group.by = 'regionName', cols = region_color) + ggtitle(subclasses[i]) + xlim(x1,x2) + ylim(y1,y2)
  ggsave(path = "UMAP", filename=paste0('Exc_subclasses_',i,'.pdf'), width=5, height=4)
  
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
    ggsave(path = "UMAP", filename=paste0('Exc_subclasses_',i,'_',j,'.pdf'), width=4.9, height=3.65) 
  }
}





