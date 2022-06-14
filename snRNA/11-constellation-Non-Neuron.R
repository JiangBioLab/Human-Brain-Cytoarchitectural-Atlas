library(scrattch.hicat)
library(Seurat)
source('KNN.graph.R')


seu <- readRDS('seu.harmony.nonNeuron.rds')

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
saveRDS(a,'get_knn_graph.NonNeuron.rds')




a<-readRDS('get_knn_graph.NonNeuron.rds')
knn.cl.df2 <- a$knn.cl.df

clusters <- c("84", "82", "6",  "129","117","79", "35", "125","59", "96",    ##Non-Neuron
              "31", "103","70", "69", "27", "2",  "42", "127","1",  "56", 
              "95", "57", "89", "66", "4")
cluster_color <- c("#382650","#2F2144","#513577","#4D4639","#6E7C6F","#535944","#697255","#3B3F32","#A77C70","#A6906F",      ##Non-Neuron
                   "#8D6041","#604B47","#2E3E39","#4D483C","#332F26","#5A503F","#72B09C","#6CA491","#47645C","#4F7165",
                   "#65978A","#3F574E","#5E8A79","#577D70","#374A45")
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

seu$hicat_cluster_subclasses[seu$hicat_cluster_merge=='OPC2']="OPC"

anno <- unique(cbind(as.vector(seu$hicat_cluster_merge_newID),seu$hicat_cluster_subclasses))
colnames(anno) <- c('cluster','anno')
rownames(anno) <- anno[,1]

anno_color <- cbind(c("Microglia", "VLMC L1-L3", "ENDO L2-L5", "OPC", "Oligo L4-L6 MOBP", "Oligo L4-L6 OPALIN", "Astrocyte FGFR3"),
                    c("#513577", "#697255", "#604B47", "#2E3E39", "#5E8A79", "#72B09C", "#665C47"))
rownames(anno_color) <- anno_color[,1]
colnames(anno_color) <- c('anno','color')
anno_color <- as.data.frame(anno_color)
anno_color$id <- c(1:7)

cl.center.df2$clade_label <- anno[cl.center.df2$cl,2]
cl.center.df2$clade_id <- anno_color[cl.center.df2$clade_label,3]
cl.center.df2$clade_color <- anno_color[cl.center.df2$clade_label,2]

cl.center.df2$cluster_id <- cl.center.df2$cl
cl.center.df2$cluster_label <- cl.center.df2$cl

plot_constellation(knn.cl.df = knn.cl.df2, cl.center.df = cl.center.df2, out.dir = "plot_constellation", node.dodge=TRUE, plot.hull=c(1,2), label_repel=TRUE) 





