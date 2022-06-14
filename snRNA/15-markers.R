library(scrattch.hicat)
library(dendextend)
library(Seurat)

seu <- readRDS('seu.v2.rds')
meta <- readRDS('seu.harmony.anno.meta.v2.rds')
dend.result <- readRDS('Brain_subClusteringID_seurat_merge_dend_renames.rds')

seu$hicat_cluster_newID <- meta[colnames(seu),'hicat_cluster_merge']
seu$hicat_cluster_newID <- factor(seu$hicat_cluster_newID,
                                  levels = rev(labels(dend.result$dend)))
seu@active.ident <- seu$hicat_cluster_newID

clusters <- levels(seu$hicat_cluster_newID)
length(clusters)
#142

markers <- c()
for(i in 1:length(clusters)){
  tmp <- FindMarkers(seu, ident.1 = clusters[i])
  tmp$gene <- rownames(tmp)
  tmp$cluster <- rep(clusters[i],dim(tmp)[1])
  markers <- rbind(markers,tmp)
}
saveRDS(markers,'markers.rds')


m1 <- readRDS('markers1.rds')
m2 <- readRDS('markers2.rds')
m3 <- readRDS('markers3.rds')
m4 <- readRDS('markers4.rds')
m5 <- readRDS('markers5.rds')
m6 <- readRDS('markers6.rds')
m7 <- readRDS('markers7.rds')
m8 <- readRDS('markers8.rds')
m9 <- readRDS('markers9.rds')
m10 <- readRDS('markers10.rds')
m11 <- readRDS('markers11.rds')
m12 <- readRDS('markers12.rds')
m13 <- readRDS('markers13.rds')
m14 <- readRDS('markers14.rds')

m14_1 <- readRDS('markers446.rds')
m14_2 <- readRDS('markers453.rds')
m14_3 <- readRDS('markers454.rds')
m14_4 <- readRDS('markers456.rds')
m14_5 <- readRDS('markers463.rds')
m14_6 <- readRDS('markers471.rds')

colnames(m8) <- colnames(m1)
colnames(m9) <- colnames(m1)
colnames(m11) <- colnames(m1)
colnames(m13) <- colnames(m1)
colnames(m14) <- colnames(m1)

m <- rbind(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m14_1,m14_2,m14_3,m14_4,m14_5,m14_6)
saveRDS(m, 'hicat_cluster_merge_markers.rds')

t1 <- readRDS('markers_celltype1.rds')
t2 <- readRDS('markers_celltype2.rds')
t3 <- readRDS('markers_celltype3.rds')
t4 <- readRDS('markers_celltype4.rds')
t5 <- readRDS('markers_celltype5.rds')

t <- rbind(t1,t2,t3,t4,t5)
saveRDS(t, 'hicat_cluster_celltype_markers.rds')

c <- readRDS('markers_classes.rds')



############# DotPlot ####################################
library(dplyr)
c_top5 <- top_n(group_by(c,cluster), 5, avg_logFC)
c_top5_genes <- unique(c_top5$gene)

t_top5 <- top_n(group_by(t,cluster), 5, avg_logFC)
t_top5_genes <- unique(t_top5$gene)

m_top5 <- top_n(group_by(m,cluster), 5, avg_logFC)
m_top5_genes <- unique(m_top5$gene)


seu$hicat_cluster_subclass <- meta[colnames(seu),'hicat_cluster_anno']
seu$hicat_cluster_subclass <- factor(seu$hicat_cluster_subclass,
                                  levels = c("OPC","Oligo L4-L6 MOBP","Oligo L4-L6 OPALIN",
                                             "Astrocyte FGFR3","NP","X","L6CT/FEZF2","L6B/FEZF2",
                                             "L5 IT/RORB","L1-L3 IT LINC00507","ET","L6 CAR3",
                                             "L6/IT","LAMP5/SV2C","PAX6","LAMP5/SNCG","VIP",
                                             "SST","PVALB","OPC2","Microglia","VLMC L1-L3","ENDO L2-L5"))

bigclasses <- unique(c$cluster)
classes <- levels(seu$hicat_cluster_subclass)
clusters <- rev(labels(dend.result$dend))


genes0 <- c()
g0 <- table(c$gene)
for(i in 1:length(bigclasses)){
  tmp <- c[c$cluster==bigclasses[i],]
  if(dim(tmp)[1]<1){
    print(i)
    next
  }
  j=0
  for(k in 1:dim(tmp)[1]){
    if(!(tmp$gene[k] %in% genes1)){
      if(g0[tmp$gene[k]]>1){
        next
      }
      genes0 <- c(genes0, tmp$gene[k])
      j=j+1
      if(j==5){
        break
      }
    }
  }
}

g<-table(m$gene)
genes1 <- genes0
for(i in 1:length(classes)){
  tmp <- t[t$cluster==classes[i],]
  if(dim(tmp)[1]<1){
    print(i)
    next
  }
  j=0
  for(k in 1:dim(tmp)[1]){
    if(!(tmp$gene[k] %in% genes1)){
      if(g[tmp$gene[k]]>10){
        next
      }
      genes1 <- c(genes1, tmp$gene[k])
      j=j+1
      if(j==5){
        break
      }
    }
  }
}

genes2 <- genes1
for(i in 1:length(clusters)){
  tmp <- m[m$cluster==clusters[i],]
  if(dim(tmp)[1]<1){
    print(clusters[i])
    next
  }
  j=0
  for(k in 1:dim(tmp)[1]){
    if(!(tmp$gene[k] %in% genes2)){
      if(g[tmp$gene[k]]>10){
        next
      }
      genes2 <- c(genes2, tmp$gene[k])
      j=j+1
      if(j==5){
        break
      }
    }
  }
}

length(genes2)
#755

exclude <- c('MALAT1')
genes3 <- genes2[!(genes2 %in% exclude)]

seu$hicat_cluster_newID <- factor(seu$hicat_cluster_newID,
                                  levels = labels(dend.result$dend))
## markers_dotplot.pdf 25x200
pdf('markers_dotplot.pdf',height = 25,width = 200)
DotPlot(seu, features = genes3, group.by = 'hicat_cluster_newID') + RotatedAxis()
dev.off()



############ pheatmap #########################################

seu_markers_m<-seu@assays$RNA@data[rownames(seu@assays$RNA@data) %in% genes3,]
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
seu_markers_m4 <- seu_markers_m4[clusters,genes2]

seu_markers_m5<-log2(seu_markers_m4+1)

##markers_heatmap.pdf 25x200
pdf('markers_heatmap.pdf',height = 25,width = 200)
pheatmap(seu_markers_m5,
         cluster_rows=F,
         cluster_cols=F,
         angle_col = 45) 
dev.off()


seu_markers_m6 <- t(seu_markers_m4)
for(i in 1:dim(seu_markers_m6)[1]){
  tmp <- scale(seu_markers_m6[i,],center=T,scale=T)
  seu_markers_m6[i,]<-tmp[,1]
}

seu_markers_m7 <- seu_markers_m6
seu_markers_m7[seu_markers_m7<0]<-0
seu_markers_m7[seu_markers_m7>10]<-10

pdf('markers_heatmap_zscore.pdf',height = 25,width = 200)
pheatmap(t(seu_markers_m7),
         cluster_cols = F,
         cluster_rows = F,
         angle_col = c("45"))
dev.off()

