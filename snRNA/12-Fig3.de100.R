library(Seurat)

#seu0 <- readRDS('seu.harmony.anno.rds')
seu <- readRDS('seu.harmony.anno.Glu.rds')

Glu <- c('44','46','92','67','30',
         '35','23','36','65','61','12','1','18','6','21',
         '7','11','63','19')
Glu_color <- c('#299373','#2FD3B0','#2CAD7A','#2FCA96','#8AA3CC',
               '#90A7D9','#30E6BA','#2BA087','#268063','#9EE61D','#BCFF1A','#8DD31E','#02F970','#0E993B','#0BBF45',
               '#13A23E','#00FF34','#07D945','#09CC4F')

region_color <- c('#3F4587','#8562AA','#EC8561','#B97CB5','#D43046','#F0592B','#ED4A96','#593C97','#A54486','#FBDE13','#299FAE','#75CCE3','#0C6939','#0D9547')
region <- c('FPPFC','DLPFC','VLPFC','M1','S1','S1E','PoCG','SPL','SMG','AG','VISP','ITG','STG','ACC')
regionID <- c('A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','A11','A12','A13','A14')
names(region_color) <- regionID
names(region) <- regionID

############## UMAP ############################
seu$hicat_cluster_newID <- factor(seu$hicat_cluster_newID,
                            levels = Glu)
#UMAP_Glu_cluster.pdf 5x7
DimPlot(seu,group.by = 'hicat_cluster_newID',cols = Glu_color)

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
cluster_dis<-table(seu$region,seu$hicat_cluster_newID)
cluster_dis<-cluster_dis[,Glu]
cluster_dis_r <- cluster_dis/rowSums(cluster_dis)

#Glu_region_proportion_tree.pdf 4x10
pheatmap(cluster_dis_r,cluster_cols = F)

cluster_dis_r <- as.data.frame(cluster_dis_r)
colnames(cluster_dis_r) <- c('region','cluster','ratio')
cluster_dis_r$region <- factor(cluster_dis_r$region,
                               levels = rev(c('A9','A7','A8','A4','A13','A12','A1','A2','A10','A11','A5','A14','A3','A6')))

# Glu_region_proportion.pdf 4x5.2
ggplot(cluster_dis_r,aes(cluster,region),showCategory=8)+
  geom_point(aes(size=ratio,color=region))+
  scale_color_manual(values=region_color[c('A9','A7','A8','A4','A13','A12','A1','A2','A10','A11','A5','A14','A3','A6')])+
  scale_size(limits = c(0.0001322052,0.8))+
  theme_classic() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))



aa <- rep(1,14)
names(aa) <- region[c('A9','A7','A8','A4','A13','A12','A1','A2','A10','A11','A5','A14','A3','A6')]
# Glu_region_color.pdf 4x12
barplot(aa,col = region_color[c('A9','A7','A8','A4','A13','A12','A1','A2','A10','A11','A5','A14','A3','A6')])



####### DotPlot ######################
seu$hicat_cluster_newID <- factor(seu$hicat_cluster_newID,
                                  levels = c('44','46','92','67','30',
                                             '35','23','36','65','61','12','1','18','6','21',
                                             '7','11','63','19'))
seu@active.ident <- seu$hicat_cluster_newID

genes <- c('SLC17A7',
           'FEZF2','ATP6V1C2','LINC00922','SGCG','SYT6','TSPAN18','SCUBE1','MDFIC','TNIP3','SEMA3D',
           'LINC00507','ACVR1C','THEMIS',
           'RORB','TRABD2A','GALR1','CRISPLD2','COL22A1','COL22A1','COL5A2','LINC02388','AC079793.1','AL355672.1','AL356295.1'
           )

seu_markers_m<-seu@assays$RNA@data[rownames(seu@assays$RNA@data) %in% genes,]
seu_markers_m<-t(as.matrix(seu_markers_m))
seu_markers_m<-as.data.frame(seu_markers_m)
seu_markers_m<-cbind(seu_markers_m,as.vector(seu$hicat_cluster_newID))
colnames(seu_markers_m)[dim(seu_markers_m)[2]]<-'cluster'
#seu_markers_m$cluster<-as.vector(seu_markers_m$cluster)
seu_markers_m2<-aggregate(seu_markers_m[,1:(dim(seu_markers_m)[2]-1)],by=list(cluster=seu_markers_m$cluster),FUN=mean)
seu_markers_m3<-seu_markers_m2[,-1]

seu_markers_m4<-matrix(unlist(seu_markers_m3),dim(seu_markers_m3)[1],dim(seu_markers_m3)[2])
colnames(seu_markers_m4)<-colnames(seu_markers_m3)
rownames(seu_markers_m4)<-as.vector(seu_markers_m2[,1])
seu_markers_m4<-log2(seu_markers_m4+1)
seu_markers_m4 <- seu_markers_m4[rev(c('44','46','92','67','30',
                                       '35','23','36','65','61','12','1','18','6','21',
                                       '7','11','63','19')),]

##Glu_markers_heatmap.pdf 7.1x8
pheatmap(seu_markers_m4,
         cluster_rows=F,
         cluster_cols=T,
         angle_col = 90)

pp <- pheatmap(seu_markers_m4,
               cluster_rows=F,
               cluster_cols=T,
               angle_col = 90)

# Glu_markers_dotplot.pdf 6x9.5
DotPlot(seu,features = unique(pp$tree_col$labels[pp$tree_col$order]), assay = 'RNA', cols = c("yellow", "red")) + theme(axis.text.x = element_text(angle = 90, hjust = 1))




############## Chi-Test ########################
a <- table(seu$hicat_cluster_newID,seu$regionName)

chisq <- chisq.test(a)
chisq
#X-squared = 1995.7, df = 299, p-value < 2.2e-16
#
chisq$observed
round(chisq$expected,2)
round(chisq$residuals, 3)

library(corrplot)
library(grDevices)
col2 <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                               "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                               "#4393C3", "#2166AC", "#053061")))
#Glu_cluster_region_chiTest_corrplot.pdf 8x5
corrplot(chisq$residuals, is.cor = FALSE, col =  col2(200))


############# geom_jitter #####################
tmp  <- data.frame(cells = colnames(seu),
                   cluster = seu$hicat_cluster_newID,
                   region = seu$regionName)
tmp$cluster <- factor(tmp$cluster, levels = Glu)
# Glu_geom_point.pdf 8x7
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
seu <- readRDS('seu.harmony.anno.Glu.rds')
result <- readRDS('seu.hicat.run_consensus_clust.de40.niter100.Exc.rds')
dend.result <- readRDS('seu.hicat.run_consensus_clust.de40.niter100.Exc.dend.rds')
#dend_Exc.pdf 10x3
plot(dend.result$dend,horiz = T)

seu$tmp <- result$cl.result$cl[colnames(seu)]
seu$tmp <- factor(seu$tmp,
                  levels = c(labels(dend.result$dend)))

cluster_dis<-cbind(as.vector(seu$tmp),as.vector(seu$individual))
colnames(cluster_dis)<-c('cluster','individual')
cluster_dis<-as.data.frame(cluster_dis)
cluster_dis$cluster<-factor(cluster_dis$cluster,
                            levels=rev(labels(dend.result$dend)))

library(ggplot2)
library(RColorBrewer)
##cluster_individual_composition_Exc.pdf
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
#marker_Exc.pdf 10x10
DotPlot(seu,features = unique(genes),group.by = 'tmp') + RotatedAxis()

