library(scrattch.hicat)
library(Seurat)
source('KNN.graph.R')


seu <- readRDS('seu.harmony.GABA.CGE.rds')

#cl.center.df <- read.csv('cl.center.df.csv')
#knn.cl.df <- read.csv('knn.cl.df.csv')

norm.dat<-seu@assays$RNA@data

DE <- readRDS('Brain_subClusteringID_seurat_merge_DE20.rds')

norm.dat2 <- norm.dat[DE$markers,]

cl <- seu$hicat_cluster_merge_newID
cl.df <-  data.frame(cluster_id = unique(cl[colnames(norm.dat2)]),
                     cluster_label = paste0("cluster_",unique(cl[colnames(norm.dat2)])),
                     cluster_color = rainbow(length(unique(cl[colnames(norm.dat2)]))))

rd.dat <- t(as.matrix(norm.dat2))

a <- get_knn_graph(rd.dat, cl, cl.df)
saveRDS(a,'get_knn_graph.CGE.rds')




a<-readRDS('get_knn_graph.CGE.rds')
knn.cl.df2 <- a$knn.cl.df

clusters <- c("86", "54", "33", "32", "109",
              "83", "19", "77", "121","46", "15", "68", "52", "78", "29", 
              "132","92", "51", "49", "44", "20", "99", "94", "87", "101",
              "75", "73", "122","120","90", "130","40", "126","45", "118")
cluster_color <- c("#C36AE6","#C263CC","#D97C80","#FF8995","#B36C76",  
                   "#F2859E","#F48C9E","#BE61D4","#9E56A6","#B65FBF","#B363CB","#BD76DC","#B967D9","#DD6DF2","#A45FBF",
                   "#9A5AB3","#AF6FCC","#C667D9","#D26AE6","#AC63BE","#9256A6","#7F4C8C","#DF70FF","#D270FF","#CD26FF",
                   "#CD6DF2","#DA81F7","#EA86FF","#E96DF2","#CC7E09","#B98327","#C38326","#B67B16","#A66D0D","#B57014")
names(cluster_color)<-clusters

cl.center.df2 <- data.frame(x = rep(0,length(clusters)),
                            y = rep(0,length(clusters)),
                            cl = clusters,
                            cluster_color = cluster_color,
                            cluster_size = as.vector(table(seu$hicat_cluster_merge_newID)[clusters]))

for(i in 1:length(clusters)){
  x1 <- mean(seu@reductions$umap@cell.embeddings[seu$hicat_cluster_merge_newID==clusters[i],1])
  y1 <- mean(seu@reductions$umap@cell.embeddings[seu$hicat_cluster_merge_newID==clusters[i],2])
  cl.center.df2[i,1:2] <- c(x1,y1)
}


anno <- unique(cbind(as.vector(seu$hicat_cluster_merge_newID),seu$hicat_cluster_subclasses))
colnames(anno) <- c('cluster','anno')
rownames(anno) <- anno[,1]

anno_color <- cbind(c('LAMP5/SNCG', 'LAMP5/SV2C', 'PAX6', 'VIP'),
                    c('#D96C07','#f2798d','#E96DF2','#9C28CC'))
rownames(anno_color) <- anno_color[,1]
colnames(anno_color) <- c('anno','color')
anno_color <- as.data.frame(anno_color)
anno_color$id <- c(1:4)

cl.center.df2$clade_label <- anno[cl.center.df2$cl,2]
cl.center.df2$clade_id <- anno_color[cl.center.df2$clade_label,3]
cl.center.df2$clade_color <- anno_color[cl.center.df2$clade_label,2]

cl.center.df2$cluster_id <- cl.center.df2$cl
cl.center.df2$cluster_label <- cl.center.df2$cl

plot_constellation(knn.cl.df = knn.cl.df2, cl.center.df = cl.center.df2, out.dir = "plot_constellation", node.dodge=TRUE, plot.hull=c(1,2), label_repel=TRUE) 





