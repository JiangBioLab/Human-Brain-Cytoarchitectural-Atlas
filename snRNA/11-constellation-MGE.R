library(scrattch.hicat)
library(Seurat)
source('KNN.graph.R')


seu <- readRDS('seu.harmony.GABA.MGE.rds')

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
saveRDS(a,'get_knn_graph.MGE.rds')




a<-readRDS('get_knn_graph.MGE.rds')
knn.cl.df2 <- a$knn.cl.df

clusters <- c("67", "100","97", "102","25", "24", "110","55", "112","58",     ##MGE
              "21", "28", "37", "135","141","23", "64", "105","36", "50", 
              "88", "113","8")
cluster_color <- c("#F27503","#E67105","#EC820D","#E69A05","#D99207","#FF7A00","#CE8024","#BF670B","#C36F1E","#DC7C16",     ##MGE
                   "#CC7209","#E89420","#FFA900","#AB5B42","#D0494E","#B1614A","#D92D43","#EC2D46","#8C2A27","#873C46",
                   "#B9342C","#D93031","#B92C3E")
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

anno_color <- cbind(c("PVALB", "SST"),
                    c('#FF2D4E','#F2B003'))
rownames(anno_color) <- anno_color[,1]
colnames(anno_color) <- c('anno','color')
anno_color <- as.data.frame(anno_color)
anno_color$id <- c(1:2)

cl.center.df2$clade_label <- anno[cl.center.df2$cl,2]
cl.center.df2$clade_id <- anno_color[cl.center.df2$clade_label,3]
cl.center.df2$clade_color <- anno_color[cl.center.df2$clade_label,2]

cl.center.df2$cluster_id <- cl.center.df2$cl
cl.center.df2$cluster_label <- cl.center.df2$cl

plot_constellation(knn.cl.df = knn.cl.df2, cl.center.df = cl.center.df2, out.dir = "plot_constellation", node.dodge=TRUE, plot.hull=c(1,2), label_repel=TRUE) 





