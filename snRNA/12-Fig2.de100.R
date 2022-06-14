library(Seurat)

################ GABA CGE ######################################################

#seu0 <- readRDS('seu.harmony.anno.de100.rds')
seu <- readRDS('seu.harmony.anno.de100.GABA.CGE.rds')

GABA <- c('24','79',
          '15','95','88','58','53','37','34','83','31','33',
          '14','93','10','75','41','66','26','90','32','72',
          '29','69','74','86','85','60','57')
GABA_color <- c('#F2859E','#D97C80',
                '#F48C9E','#FF8995','#A45FBF','#BD76DC','#9E56A6','#AF6FCC','#B967D9','#BE61D4','#B65FBF','#DD6DF2',
                '#9A5AB3','#D270FF','#D26AE6','#D92D43','#E69A05','#B98327','#C39317','#B363CB','#A6770D','#C38326',
                '#996A0E','#AC7913','#A66D0D','#8C620E','#D99207','#B57014','#B67B16')

region_color <- c('#3F4587','#8562AA','#EC8561','#B97CB5','#D43046','#F0592B','#ED4A96','#593C97','#A54486','#FBDE13','#299FAE','#75CCE3','#0C6939','#0D9547')
region <- c('FPPFC','DLPFC','VLPFC','M1','S1','S1E','PoCG','SPL','SMG','AG','VISP','ITG','STG','ACC')
regionID <- c('A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','A11','A12','A13','A14')
names(region_color) <- regionID
names(region) <- regionID

############## UMAP ############################
seu$hicat_cluster_newID <- factor(seu$hicat_cluster_newID,
                            levels = GABA)
#UMAP_GABA_CGE_cluster.pdf 5x7
DimPlot(seu,group.by = 'hicat_cluster_newID',cols = GABA_color)

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
cluster_dis<-table(seu$region,seu$hicat_cluster_newID)
cluster_dis<-cluster_dis[,GABA]
cluster_dis_r <- cluster_dis/rowSums(cluster_dis)

#GABA_CGE_region_proportion_tree.pdf 4x10
pheatmap(cluster_dis_r,cluster_cols = F)

cluster_dis_r <- as.data.frame(cluster_dis_r)
colnames(cluster_dis_r) <- c('region','cluster','ratio')
cluster_dis_r$region <- factor(cluster_dis_r$region,
                               levels = rev(c('A5','A11','A8','A14','A3','A6','A4','A9','A12','A1','A2','A13','A7','A10')))

# GABA_CGE_region_proportion.pdf 4x7.05
ggplot(cluster_dis_r,aes(cluster,region),showCategory=8)+
  geom_point(aes(size=ratio,color=region))+
  scale_color_manual(values=region_color[c('A5','A11','A8','A14','A3','A6','A4','A9','A12','A1','A2','A13','A7','A10')])+
  scale_size(limits = c(0.0006397953,0.4))+
  theme_classic() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))

aa <- rep(1,14)
names(aa) <- region[c('A5','A11','A8','A14','A3','A6','A4','A9','A12','A1','A2','A13','A7','A10')]
# GABA_CGE_region_color.pdf 4x12
barplot(aa,col = region_color[c('A5','A11','A8','A14','A3','A6','A4','A9','A12','A1','A2','A13','A7','A10')])



####### DotPlot ######################
seu$hicat_cluster_newID <- factor(seu$hicat_cluster_newID,
                            levels = c('24','79',
                                       '15','95','88','58','53','37','34','83','31','33',
                                       '14','93','10','75','41','66','26','90','32','72',
                                       '29','69','74','86','85','60','57'))
seu@active.ident <- seu$hicat_cluster_newID

genes <- c('GAD1','ADARB2',
           'LAMP5','SV2C','SFTA3','MDGA1','ROR2','TBL1Y','TRPC3', #LAMP5
           'PAX6','NMBR','AC091576.1','SYT10','SLCO5A1','ANGPT1','DDR2','TGFBR2','CILP2','LIX1','TNFAIP8L3','RXFP1','CHRNB3','HCRTR2', #PAX6
           'VIP','PENK','ABI3BP','ADAM12','ENPEP','DACH2','HTR2C','SCML4','SCTR','NOX4','LINC01630','CHRNA2','MSR1' #VIP
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
seu_markers_m4 <- seu_markers_m4[rev(c('24','79',
                                       '15','95','88','58','53','37','34','83','31','33',
                                       '14','93','10','75','41','66','26','90','32','72',
                                       '29','69','74','86','85','60','57')),]

##GABA_markers_heatmap_CGE.pdf 7.1x8
pheatmap(seu_markers_m4,
         cluster_rows=F,
         cluster_cols=T,
         angle_col = 90)

pp <- pheatmap(seu_markers_m4,
               cluster_rows=F,
               cluster_cols=T,
               angle_col = 90)

# GABA_markers_dotplot_CGE.pdf 7.1x10
DotPlot(seu,features = unique(pp$tree_col$labels[pp$tree_col$order]), assay = 'RNA', cols = c("yellow", "red")) + theme(axis.text.x = element_text(angle = 90, hjust = 1))


############## Chi-Test ########################
a <- table(seu$hicat_cluster_newID,seu$regionName)

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
                   cluster = seu$hicat_cluster_newID,
                   region = seu$regionName)
tmp$cluster <- factor(tmp$cluster, levels = GABA)
# CGE_geom_point.pdf 8x8
ggplot(tmp, aes(x=factor(cluster), y=region, fill=factor(cluster)) ) + 
  geom_jitter(position=position_jitter(width=0.3, height=0.4), size = 0.5, aes(colour=factor(cluster)), alpha=0.9)+
  scale_color_manual(values=GABA_color)+
  theme_classic()
              

#####################################################################################





#################### GABA MGE #######################################################

#seu0 <- readRDS('seu.harmony.anno.de100.rds')
seu <- readRDS('seu.harmony.anno.de100.GABA.MGE.rds')

GABA <- c('97','16','80','96','8','76',
          '49','50','59','81','84','42','43','77','28','13',
          '38','91','87','20','17','51','27','39')
GABA_color <- c('#F2B003','#D62228','#A62A2F','#EC2D46','#8C2734','#B92C3E',
                '#B9342C','#D93031','#EC820D','#F28403','#D96C07','#FFB800','#CC6709','#E67105','#F2841A','#E89420',
                '#FFA900','#FF7A00','#F27503','#D97807','#FF9900','#EC3031','#FF2D4E','#C60C0F')

region_color <- c('#3F4587','#8562AA','#EC8561','#B97CB5','#D43046','#F0592B','#ED4A96','#593C97','#A54486','#FBDE13','#299FAE','#75CCE3','#0C6939','#0D9547')
region <- c('FPPFC','DLPFC','VLPFC','M1','S1','S1E','PoCG','SPL','SMG','AG','VISP','ITG','STG','ACC')
regionID <- c('A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','A11','A12','A13','A14')
names(region_color) <- regionID
names(region) <- regionID

############## UMAP ############################
seu$hicat_cluster_newID <- factor(seu$hicat_cluster_newID,
                                  levels = GABA)
#UMAP_GABA_MGE_cluster.pdf 5x7
DimPlot(seu,group.by = 'hicat_cluster_newID',cols = GABA_color)

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
cluster_dis<-table(seu$region,seu$hicat_cluster_newID)
cluster_dis<-cluster_dis[,GABA]
cluster_dis_r <- cluster_dis/rowSums(cluster_dis)

#GABA_MGE_region_proportion_tree.pdf 4x10
pheatmap(cluster_dis_r,cluster_cols = F)

cluster_dis_r <- as.data.frame(cluster_dis_r)
colnames(cluster_dis_r) <- c('region','cluster','ratio')
cluster_dis_r$region <- factor(cluster_dis_r$region,
                               levels = rev(c('A2','A13','A8','A7','A9','A3','A5','A6','A14','A4','A10','A1','A11','A12')))

# GABA_MGE_region_proportion.pdf 4x6
ggplot(cluster_dis_r,aes(cluster,region),showCategory=8)+
  geom_point(aes(size=ratio,color=region))+
  scale_color_manual(values=region_color[c('A2','A13','A8','A7','A9','A3','A5','A6','A14','A4','A10','A1','A11','A12')])+
  scale_size(limits = c(0.0006064281,0.5))+
  theme_classic() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))


aa <- rep(1,14)
names(aa) <- region[c('A2','A13','A8','A7','A9','A3','A5','A6','A14','A4','A10','A1','A11','A12')]
# GABA_MGE_region_color.pdf 4x12
barplot(aa,col = region_color[c('A2','A13','A8','A7','A9','A3','A5','A6','A14','A4','A10','A1','A11','A12')])

####### DotPlot ######################
seu$hicat_cluster_newID <- factor(seu$hicat_cluster_newID,
                                  levels = c('97','16','80','96','8','76',
                                             '49','50','59','81','84','42','43','77','28','13',
                                             '38','91','87','20','17','51','27','39'))
seu@active.ident <- seu$hicat_cluster_newID

genes <- c('GAD1','LHX6',
           'PVALB','COL15A1','AC046195.2','MEPE','TAC1','PIEZO2','CALB1','LHX6','AGBL1','FOLH1','SULF1','BLNK', #PVALB
           'SST','GRIK3','NPY','CALB1','SLC27A6','AC087379.2','KLHL14','CNBD1','ST18','TRPC6','TAC1','STK32B','STK32A','ADGRG6','SLC9A2' #SST
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
seu_markers_m4 <- seu_markers_m4[rev(c('97','16','80','96','8','76',
                                       '49','50','59','81','84','42','43','77','28','13',
                                       '38','91','87','20','17','51','27','39')),]

##GABA_markers_heatmap_MGE.pdf 6x6
pheatmap(seu_markers_m4,
         cluster_rows=F,
         cluster_cols=T,
         angle_col = 90)

pp <- pheatmap(seu_markers_m4,
               cluster_rows=F,
               cluster_cols=T,
               angle_col = 90)

# GABA_markers_dotplot_MGE.pdf 6x8
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
#MGE_cluster_region_chiTest_corrplot.pdf 8x5
corrplot(chisq$residuals, is.cor = FALSE, col =  col2(200))


############# geom_jitter #####################
tmp  <- data.frame(cells = colnames(seu),
                   cluster = seu$hicat_cluster_newID,
                   region = seu$regionName)
tmp$cluster <- factor(tmp$cluster, levels = GABA)
# MGE_geom_point.pdf 8x8
ggplot(tmp, aes(x=factor(cluster), y=region, fill=factor(cluster)) ) + 
  geom_jitter(position=position_jitter(width=0.3, height=0.4), size = 0.5, aes(colour=factor(cluster)), alpha=0.9)+
  scale_color_manual(values=GABA_color)+
  theme_classic()


######################################################################################


genes <- c('ADARB2','PROX1','LHX6','RXFP3','NTF3','LAMP5','PDLIM5','NDNF',
           'RXFP1','DOCK5','LSP1','SLC35D3','JAM2','EGLN3','TAFA1','NPY2R',
           'PAX6','KRT73','SNCG','SERPINF1','SLC17A8','CALCB','NPFFR1','NTNG1',
           'VIP','PTHLH','PCDH11X','CP','MYBPC1','GPC3','SLC5A7','CBLN4','CHAT',
           'RSPO1','LMO1','TMEM176A','QRFPR','IGFBP6')






supertype_color1 <- c('#F2B003','#D62228','#A62A2F','#EC2D46','#8C2734','#B92C3E',
                      '#B9342C','#D93031','#EC820D','#F28403','#D96C07','#FFB800','#CC6709','#E67105','#F2841A','#E89420',
                      '#FFA900','#FF7A00','#F27503','#D97807','#FF9900','#EC3031','#FF2D4E','#C60C0F','#F2859E','#D97C80',
                      '#F48C9E','#FF8995','#A45FBF','#BD76DC','#9E56A6','#AF6FCC','#B967D9','#BE61D4','#B65FBF','#DD6DF2',
                      '#9A5AB3','#D270FF','#D26AE6','#D92D43','#E69A05','#B98327','#C39317','#B363CB','#A6770D','#C38326',
                      '#996A0E','#AC7913','#A66D0D','#8C620E','#D99207','#B57014','#B67B16')
names(supertype_color1) <- c('97','16','80','96','8','76',
                             '49','50','59','81','84','42','43','77','28','13',
                             '38','91','87','20','17','51','27','39','24','79',
                             '15','95','88','58','53','37','34','83','31','33',
                             '14','93','10','75','41','66','26','90','32','72',
                             '29','69','74','86','85','60','57')

#Inh_cluster_color.pdf 3x10
barplot(rep(1,53),col = supertype_color1)













