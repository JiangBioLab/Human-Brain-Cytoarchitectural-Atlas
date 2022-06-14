library(scrattch.hicat)
library(Seurat)
source('KNN.graph.R')


seu <- readRDS('seu.harmony.Glu.rds')

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
saveRDS(a,'get_knn_graph.Glu.rds')




a<-readRDS('get_knn_graph.Glu.rds')
knn.cl.df2 <- a$knn.cl.df

clusters <- c("111","119","140","116","104","60", "136",                      ##Glu
              "128","139","123","85", "134","76", "138","93", "91", "131",
              "133","98", "80", "107","142","72", "114","63", "81", "53", 
              "124","115","137","61", "106","74", "22", "47", "43", "71", 
              "16", "11", "10", "5",  "38", "3",  "9",  "7",  "17", "26", 
              "34", "18", "30", "108","13", "41", "48", "39", "65", "62", 
              "12", "14")
cluster_color <- c("#93C43B","#AFE32A","#B1E67B","#C0F27F","#A4D374","#88A662","#8D8C20",                                   ##Glu
                   "#898325","#A7A322","#B0CE1F","#B1B10C","#DABE23","#CEC823","#C5A93D","#91910E","#C8B323","#0E8B8B",
                   "#09EEDC","#19E3BE","#2B8AA5","#7944AA","#2FBCE5","#2E84BE","#2E9FCB","#296B98","#2C93B2","#2B77A5",
                   "#69419D","#7E47B6","#5E419D","#71419D","#7643A4","#6444AA","#0BBEBE","#52ADC0","#4A9E9E","#18C5C8",
                   "#02E2F9","#09CCC6","#07C6D9","#2E9FBB","#378695","#138C98","#0DA3A3","#0E9797","#07D8D8","#5BDFDC",
                   "#637384","#5B7893","#61867A","#5991A4","#4EA8AC","#41AEAE","#4E9F8C","#99FFCC","#52B8AA","#4E9E9E",
                   "#55C5B5","#4BBCAD")
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

seu$hicat_cluster_subclasses[seu$hicat_cluster_merge=='FEZF2.5_1']="L6CT/FEZF2"
seu$hicat_cluster_subclasses[seu$hicat_cluster_merge=='FEZF2.4_3']="L6B/FEZF2"

anno <- unique(cbind(as.vector(seu$hicat_cluster_merge_newID),seu$hicat_cluster_subclasses))
colnames(anno) <- c('cluster','anno')
rownames(anno) <- anno[,1]

anno_color <- cbind(c("NP", "L6CT/FEZF2", "L6B/FEZF2", "L5 IT/RORB", "L1-L3 IT LINC00507", "ET", "L6 CAR3", "L6/IT"),
                    c("#93C43B", "#2D8CB8", "#7944AA", "#09B2B2", "#07D8D8", "#19E3BE", "#898325", "#CEC823"))
rownames(anno_color) <- anno_color[,1]
colnames(anno_color) <- c('anno','color')
anno_color <- as.data.frame(anno_color)
anno_color$id <- c(1:8)

cl.center.df2$clade_label <- anno[cl.center.df2$cl,2]
cl.center.df2$clade_id <- anno_color[cl.center.df2$clade_label,3]
cl.center.df2$clade_color <- anno_color[cl.center.df2$clade_label,2]

cl.center.df2$cluster_id <- cl.center.df2$cl
cl.center.df2$cluster_label <- cl.center.df2$cl

plot_constellation(knn.cl.df = knn.cl.df2, cl.center.df = cl.center.df2, out.dir = "plot_constellation", node.dodge=TRUE, plot.hull=c(1,2), label_repel=TRUE) 





