library(Seurat)


#################### NonNeuron #######################################################

seu <- readRDS('seu.harmony.nonNeuron.rds')

NonNeuron <- c("84", "82", "6",  "129","117","79", "35", "125","59", "96",    ##Non-Neuron
               "31", "103","70", "69", "27", "2",  "42", "127","1",  "56", 
               "95", "57", "89", "66", "4")
NonNeuron_color <- c("#382650","#2F2144","#513577","#4D4639","#6E7C6F","#535944","#697255","#3B3F32","#A77C70","#A6906F",      ##Non-Neuron
                     "#8D6041","#604B47","#2E3E39","#4D483C","#332F26","#5A503F","#72B09C","#6CA491","#47645C","#4F7165",
                     "#65978A","#3F574E","#5E8A79","#577D70","#374A45")
names(NonNeuron_color) <- NonNeuron

region_color <- c('#3F4587','#8562AA','#EC8561','#B97CB5','#D43046','#F0592B','#ED4A96','#593C97','#A54486','#FBDE13','#299FAE','#75CCE3','#0C6939','#0D9547')
region <- c('FPPFC','DLPFC','VLPFC','M1','S1','S1E','PoCG','SPL','SMG','AG','VISP','ITG','STG','ACC')
regionID <- c('A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','A11','A12','A13','A14')
names(region_color) <- region
names(region) <- regionID

############## UMAP ############################
seu$hicat_cluster_merge_newID <- factor(seu$hicat_cluster_merge_newID,
                                        levels = NonNeuron)
#UMAP_NonNeuron_cluster.pdf 5x7
DimPlot(seu,group.by = 'hicat_cluster_merge_newID',cols = NonNeuron_color)

############## UMAP ############################
regionID <- c('A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','A11','A12','A13','A14')
regionName <- c('FPPFC','DLPFC','VLPFC','M1','S1','S1E','PoCG','SPL','SMG','AG','VISP','ITG','STG','ACC')
names(regionName) <- regionID

seu$regionName <- regionName[seu$region]

seu$regionName <- factor(seu$regionName,
                         levels=regionName)
#UMAP_NonNeuron_region.pdf 5x6.5
DimPlot(seu,group.by = 'regionName',cols = region_color)


############## cell composition ###################
cluster_dis<-table(seu$regionName,seu$hicat_cluster_merge_newID)
cluster_dis<-cluster_dis[,NonNeuron]
cluster_dis_r <- cluster_dis/rowSums(cluster_dis)

#NonNeuron_region_proportion_tree.pdf 4x7
pheatmap(cluster_dis_r,cluster_cols = F)

p <- pheatmap(cluster_dis_r,cluster_cols = F)

cluster_dis_r <- as.data.frame(cluster_dis_r)
colnames(cluster_dis_r) <- c('region','cluster','ratio')
cluster_dis_r$region <- factor(cluster_dis_r$region,
                               levels = rev(p$tree_row$labels[p$tree_row$order]))

# NonNeuron_region_proportion.pdf 3.5x6
ggplot(cluster_dis_r,aes(cluster,region),showCategory=8)+
  geom_point(aes(size=ratio,color=region))+
  scale_color_manual(values=region_color[rev(p$tree_row$labels[p$tree_row$order])])+
  scale_size(limits = c(0.0001494098,0.65))+
  theme_classic() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))

aa <- rep(1,14)
names(aa) <- region[rev(p$tree_row$labels[p$tree_row$order])]
# NonNeuron_region_color.pdf 4x12
barplot(aa,col = region_color[rev(p$tree_row$labels[p$tree_row$order])])



####### DotPlot ######################
seu$hicat_cluster_merge_newID <- factor(seu$hicat_cluster_merge_newID,
                                        levels = NonNeuron)
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
ID <- ID[NonNeuron,]


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
seu_markers_m4 <- seu_markers_m4[rev(NonNeuron),]

##NonNeuron_markers_heatmap.pdf 5x9
pheatmap(seu_markers_m4,
         cluster_rows=F,
         cluster_cols=T,
         angle_col = 90)

pp <- pheatmap(seu_markers_m4,
               cluster_rows=F,
               cluster_cols=T,
               angle_col = 90)

# NonNeuron_markers_dotplot.pdf 6x12
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
#NonNeuron_cluster_region_chiTest_corrplot.pdf 16x10
corrplot(chisq$residuals, is.cor = FALSE, col =  col2(200))



############# geom_jitter #####################
tmp  <- data.frame(cells = colnames(seu),
                   cluster = seu$hicat_cluster_merge_newID,
                   region = seu$regionName)
tmp$cluster <- factor(tmp$cluster, levels = NonNeuron)
tmp$region <- factor(tmp$region, levels = rev(p$tree_row$labels[p$tree_row$order]))
# NonNeuron_geom_point.pdf 3x7
ggplot(tmp, aes(x=factor(cluster), y=region, fill=factor(cluster)) ) + 
  geom_jitter(position=position_jitter(width=0.3, height=0.4), size = 0.5, aes(colour=factor(cluster)), alpha=0.9)+
  scale_color_manual(values=NonNeuron_color)+
  theme_classic()








#