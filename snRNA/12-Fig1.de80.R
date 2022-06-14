library(scrattch.hicat)
library(dendextend)
library(Seurat)

#seu <- readRDS('seu.harmony.anno.de80.rds')
seu <- readRDS('seu.harmony.rds')

result<-readRDS('seu.hicat.de80.niter100.rds')
dend.result<-readRDS('seu.hicat.de80.niter100.dend.rds')

seu$de80_hicat_cluster_newID <- result$cl.result$cl[colnames(seu)]
seu$de80_hicat_cluster_newID <- factor(seu$de80_hicat_cluster_newID, 
                                       levels = rev(labels(dend.result$dend)))

#UMAP_cluster.pdf 5x9
DimPlot(seu, group.by = 'de80_hicat_cluster_newID', cols = rev(supertype_color))

seu$de80_hicat_cluster_classes <- factor(seu$de80_icat_cluster_classes,
                                         levels = c("Glutamatergic neuronal class",
                                                    "GABAergic neuronal class",
                                                    "Astrocyte/oligodendrocyte non-neuronal class",
                                                    "Immune/vascular non-neuronal class",
                                                    "Unassigned" ))

#UMAP_classes.pdf 5x10
DimPlot(seu, group.by = 'hicat_cluster_classes', cols = c('#00ADEE','#F05A28','#873C46','#808080','#101010'))


seu$individual <- factor(seu$individual,
                         levels=c('S0206','S0406','S0426','S0531'))
#UMAP_individual.pdf 5x6
DimPlot(seu, group.by = 'individual', cols = c('#1B9E77','#D95F02','#7570B3','#E7298A'))

regionID <- c('A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','A11','A12','A13','A14')
regionName <- c('FPPFC','DLPFC','VLPFC','M1','S1','S1E','PoCG','SPL','SMG','AG','VISP','ITG','STG','ACC')
names(regionName) <- regionID

seu$regionName <- regionName[seu$region]

seu$regionName <- factor(seu$regionName,
                         levels=c('FPPFC','DLPFC','VLPFC','M1','S1','S1E','PoCG','SPL','SMG','AG','VISP','ITG','STG','ACC'))
#UMAP_region.pdf 5x6
DimPlot(seu, group.by = 'regionName', cols = c('#3F4587','#8562AA','#EC8561','#B97CB5','#D43046','#F0592B',
                                               '#ED4A96','#593C97','#A54486','#FBDE13','#299FAE','#75CCE3',
                                               '#0C6939','#0D9547'))

DotPlot(seu,features = c('FEZF2', 'SGCG', 'SYT6', 'TSPAN18'),group.by = 'hicat_cluster_newID')



dend.result<-readRDS('seu.hicat.de80.niter100.dend.rds')
dend <- dend.result$dend
#dend_de80.pdf 25x4
plot(dend, horiz = T)
plot(remove_branches_edgePar(dend))


str(dend)
write.table(dend,'dend.txt')

dend.reorder <- reorder(dend,10:1)
plot(dend.reorder,horiz = T)

dend <- dend.result$dend
anno <- unique(cbind(as.vector(seu$hicat_cluster_newID),seu$hicat_cluster_supertypes))
rownames(anno) <- anno[,1]
anno[,2] <- paste0(anno[,1],':',anno[,2])

dend.label <- dend
labels(dend.label) <- anno[labels(dend),2]
#dend.pdf 20x4
plot(dend,horiz = T)
plot(dend.label,horiz = T)

# pdf dend_tree.pdf 4x20
plot(remove_branches_edgePar(dend.reorder))


library(dendsort)
plot(dendsort(dend))



############## cell composition ###################
cluster_dis<-cbind(as.vector(seu$de80_hicat_cluster_newID),as.vector(seu$individual))
colnames(cluster_dis)<-c('cluster','individual')
cluster_dis<-as.data.frame(cluster_dis)
cluster_dis$cluster<-factor(cluster_dis$cluster,
                            levels=rev(labels(dend.result$dend)))

library(ggplot2)
library(RColorBrewer)
##cluster_individual_composition_de80.pdf
##pdf 3x22
ggplot(cluster_dis, aes(cluster)) + 
  geom_bar(aes(fill=individual), position="fill", color="gray20",alpha = 0.7) +#
  scale_fill_manual(values=brewer.pal(6,"Dark2"))+ #values=colors[c(1:4,6)]
  ylab('Percentage')+
  theme(axis.line = element_line(colour = "gray10"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + RotatedAxis()

cluster_dis<-cbind(as.vector(seu$de80_hicat_cluster_newID),as.vector(seu$region))
colnames(cluster_dis)<-c('cluster','region')
cluster_dis<-as.data.frame(cluster_dis)
cluster_dis$cluster<-factor(cluster_dis$cluster,
                            levels=rev(labels(dend.result$dend)))
cluster_dis$region<-factor(cluster_dis$region,
                           levels=c('A1','A2','A3','A4','A5',
                                    'A6','A7','A8','A9','A10',
                                    'A11','A12','A13','A14'))
library(ggplot2)
library(RColorBrewer)
##cluster_region_composition_de80.pdf
##pdf 3x22
ggplot(cluster_dis, aes(cluster)) + 
  geom_bar(aes(fill=region), position="fill", color="gray20") +#,alpha = 0.7
  scale_fill_manual(values=c(brewer.pal(12,"Set3"),brewer.pal(2,"Set2")))+ #values=colors[c(1:4,6)]
  ylab('Percentage')+
  theme(axis.line = element_line(colour = "gray10"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + RotatedAxis()


#### % percent in each sample ##############
cluster_dis<-table(seu$de80_hicat_cluster_newID,seu$region)
cluster_dis_r <- cluster_dis/rowSums(cluster_dis)
cluster_dis_r <- cluster_dis_r[rev(labels(dend.result$dend)),]

cluster_dis_r[cluster_dis_r>0.3]<-0.3
myColors=brewer.pal(8,"Reds")[1:8]

cluster_dis_r <- cluster_dis_r[,c('A8','A9','A3','A5','A14','A12','A4','A7','A10','A13','A11','A1','A2','A6')]

library(pheatmap)
library(RColorBrewer)
#pdf cluster_region_composition_ratio_de80.pdf 22x3.5
#pdf cluster_region_composition_ratio_tree.pdf 20x3.5
pheatmap(cluster_dis_r,
         cluster_rows = F,
         cluster_cols = F, # T for Tree
         color = colorRampPalette(colors = myColors)(100))


cluster_dis_r_samples <- cluster_dis_r
for(i in 1:dim(cluster_dis_r)[1]){
  for(j in 1:dim(cluster_dis_r)[2]){
    cluster_dis_r_samples[i,j]<-length(unique(seu$individual[seu$de80_hicat_cluster_newID==rownames(cluster_dis_r)[i] & seu$region==colnames(cluster_dis_r)[j]]))
  }
}
cluster_dis_r_samples <- cluster_dis_r_samples[,c('A8','A9','A3','A5','A14','A12','A4','A7','A10','A13','A11','A1','A2','A6')]

#pdf cluster_region_composition_samples_de80.pdf 22x3.5
pheatmap(cluster_dis_r_samples,
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = T,
         number_format = '%.f')


cluster_dis<-cluster_dis[rev(labels(dend.result$dend)),
                         c('A8','A9','A3','A5','A14','A12','A4','A7','A10','A13','A11','A1','A2','A6')]
#pdf cluster_region_composition_cells_de80.pdf 22x3.5
pheatmap(log2(cluster_dis+1),
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = T,
         number_format = '%.f')


# cell_number_log2_de80.pdf 3x22
barplot(log2(table(seu$de80_hicat_cluster_newID)[rev(labels(dend.result$dend))]),
        col = supertype_color)

#cluster_color.pdf 3x20
barplot(rep(1,97),col = supertype_color)


################# pheatmap #############################
genes1 <- c('PTRPC','NOSTRIN','PDGRFB','MBP','PDGFRA','FEZF2','LINC00507','COL5A2','LAMP5','SV2C','PAX6','VIP','PVALB','SST')

genes2 <- c('TYROBP','NOSTRIN','OPALIN','FGFR3','PDGFRA','SLC1A3','FEZF2','THEMIS','RORB','LINC00507','SLC17A7','PVALB','SST','LHX6','VIP','LAMP5','PAX6','ADARB2','GAD1')

genes1[!genes1 %in% genes2]
genes2[!genes2 %in% genes1]

genes <- c('PTPRC','TYROBP','NOSTRIN','PDGFRB','OPALIN','FGFR3','PDGFRA','SLC1A3','FEZF2','THEMIS','LINC00507','SLC17A7','RORB','LAMP5','SV2C','PAX6','VIP','ADARB2','PVALB','SST','LHX6','GAD1')
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
seu_markers_m4 <- seu_markers_m4[rev(c('9','94','64','62','5','78','55','48','47','40','89','71','2',
                                       '52','73','22','82','45','25','56','70','3','68','54','4','44',
                                       '46','92','67','30','35','23','36','65','61','12','1','18','6',
                                       '21','7','11','63','19','97','16','80','96','8','76','49','50',
                                       '59','81','84','42','43','77','28','13','38','91','87','20','17',
                                       '51','27','39','24','79','15','95','88','58','53','37','34','83',
                                       '31','33','14','93','10','75','41','66','26','90','32','72','29',
                                       '69','74','86','85','60','57')),]
seu_markers_m4 <- seu_markers_m4[,genes]

anno <- unique(cbind(as.vector(seu$hicat_cluster_newID),seu$hicat_cluster_supertypes))
rownames(anno)<-anno[,1]

#rownames(seu_markers_m4) <- paste0(rownames(seu_markers_m4),':',anno[rownames(seu_markers_m4),2])
rownames(seu_markers_m4) <- anno[rownames(seu_markers_m4),2]

##marker_heatmap.pdf 19.4x7
pheatmap(seu_markers_m4,
         cluster_rows=F,
         cluster_cols=F,
         angle_col = 45)

seu_markers_m4_z <- t(seu_markers_m4)
for(i in 1:dim(seu_markers_m4_z)[1]){
  tmp <- scale(seu_markers_m4_z[i,],center=T,scale=T)
  seu_markers_m4_z[i,]<-tmp[,1]
}

#seu_markers_m4_z2 <- t(seu_markers_m4)
#for(i in 1:dim(seu_markers_m4_z2)[1]){
#  zscore = (seu_markers_m4_z2[i,] - mean(seu_markers_m4_z2[i,]))/sd(seu_markers_m4_z2[i,])
#  seu_markers_m4_z2[i,]<-zscore
#}

seu_markers_m4_z[seu_markers_m4_z>5] = 5

#marker_zscore_heatmap.pdf 22x7
pheatmap(t(seu_markers_m4_z),
         cluster_cols = F,
         cluster_rows = F,
         angle_col = c("45"))



######### DotPlot ############################
genes <- c('PTPRC','TYROBP','NOSTRIN','PDGFRB','OPALIN','FGFR3','PDGFRA','SLC1A3','FEZF2','THEMIS','LINC00507','SLC17A7','RORB','LAMP5','SV2C','PAX6','VIP','ADARB2','PVALB','SST','LHX6','GAD1')
seu$de80_hicat_cluster_newID <- factor(seu$de80_hicat_cluster_newID,
                                       levels=labels(dend.result$dend))
seu@active.ident <- seu$de80_hicat_cluster_newID
##marker_dotplot_de80.pdf 22x9
pdf('marker_dotplot_de80.pdf',width = 9, height = 22)
DotPlot(seu,features = genes, group.by = 'de80_hicat_cluster_newID') + RotatedAxis()
dev.off()




##################################################
#chi-test in R
#http://www.sthda.com/english/wiki/chi-square-test-of-independence-in-r

a <- table(seu$hicat_cluster_newID,seu$regionName)

chisq <- chisq.test(a)
chisq
#X-squared = 26865, df = 1248, p-value < 2.2e-16
#
chisq$observed
round(chisq$expected,2)
round(chisq$residuals, 3)


library(corrplot)
library(grDevices)
col2 <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                               "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                               "#4393C3", "#2166AC", "#053061")))
#cluster_region_chiTest_corrplot.pdf 20x5
corrplot(chisq$residuals, is.cor = FALSE, col =  col2(200))







