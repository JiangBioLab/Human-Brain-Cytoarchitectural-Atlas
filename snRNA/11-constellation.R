library(scrattch.hicat)
source('KNN.graph.R')

seu <- readRDS('seu.harmony.anno.v2.rds')
#meta <- readRDS('seu.harmony.anno.meta.v2.rds')
#
#
################## Color ################################
#supertype_color <- c("#577D70", "#374A45", "#5E8A79", "#3F574E", "#65978A", "#47645C", "#4F7165", "#6CA491", "#72B09C", "#332F26", 
#                     "#5A503F", "#4D483C", "#A4D374", "#88A662", "#C0F27F", "#B1E67B", "#AFE32A", "#93C43B", "#2B8AA5", "#7944AA", 
#                     "#2C93B2", "#2B77A5", "#296B98", "#2E9FCB", "#2E84BE", "#2FBCE5", "#69419D", "#7E47B6", "#5E419D", "#71419D", 
#                     "#7643A4", "#6444AA", "#52ADC0", "#4A9E9E", "#0BBEBE", "#4EA8AC", "#41AEAE", "#5991A4", "#5B7893", "#61867A", 
#                     "#637384", "#52B8AA", "#4E9E9E", "#55C5B5", "#4BBCAD", "#4E9F8C", "#99FFCC", "#5BDFDC", "#0E9797", "#07D8D8", 
#                     "#0DA3A3", "#138C98", "#378695", "#07C6D9", "#2E9FBB", "#09CCC6", "#02E2F9", "#18C5C8", "#09EEDC", "#19E3BE", 
#                     "#0E8B8B", "#8D8C20", "#898325", "#B0CE1F", "#B1B10C", "#A7A322", "#91910E", "#C8B323", "#C5A93D", "#DABE23", 
#                     "#CEC823", "#F2859E", "#F48C9E", "#B36C76", "#FF8995", "#D97C80", "#C36AE6", "#C263CC", "#CD6DF2", "#DA81F7", 
#                     "#EA86FF", "#E96DF2", "#D270FF", "#CD26FF", "#A66D0D", "#B57014", "#B67B16", "#B98327", "#C38326", "#CC7E09", 
#                     "#7F4C8C", "#DF70FF", "#9256A6", "#AC63BE", "#C667D9", "#D26AE6", "#9A5AB3", "#AF6FCC", "#DD6DF2", "#A45FBF", 
#                     "#B967D9", "#BD76DC", "#B65FBF", "#B363CB", "#BE61D4", "#9E56A6", "#E67105", "#EC820D", "#F27503", "#D99207", 
#                     "#FF7A00", "#E69A05", "#BF670B", "#C36F1E", "#CE8024", "#E89420", "#FFA900", "#CC7209", "#DC7C16", "#D93031", 
#                     "#B92C3E", "#B9342C", "#8C2A27", "#873C46", "#D92D43", "#EC2D46", "#D0494E", "#B1614A", "#AB5B42", "#2E3E39", 
#                     "#2F2144", "#513577", "#382650", "#697255", "#3B3F32", "#535944", "#6E7C6F", "#4D4639", "#8D6041", "#604B47", 
#                     "#A6906F", "#A77C70")
#names(supertype_color) <- c("66", "4", "89", "57", "95", "1", "56", "127", "42", "27", 
#                            "2", "69", "104", "60", "116", "140", "119", "111", "80", "107", 
#                            "81", "53", "63", "114", "72", "142", "124", "115", "137", "61", 
#                            "106", "74", "47", "43", "22", "13", "41", "108", "18", "30", 
#                            "34", "65", "62", "12", "14", "48", "39", "26", "7", "17", 
#                            "9", "3", "38", "10", "5", "11", "16", "71", "133", "98", 
#                            "131", "136", "128", "123", "85", "139", "93", "91", "138", "134", 
#                            "76", "83", "19", "109", "32", "33", "86", "54", "75", "73", 
#                            "122", "120", "87", "101", "45", "118", "126", "130", "40", "90", 
#                            "99", "94", "20", "44", "51", "49", "132", "92", "78", "29", 
#                            "52", "68", "46", "15", "77", "121", "100", "97", "67", "25", 
#                            "24", "102", "55", "112", "110", "28", "37", "21", "58", "113", 
#                            "8", "88", "36", "50", "64", "105", "141", "23", "135", "70", 
#                            "82", "6", "84", "35", "125", "79", "117", "129", "31", "103", 
#                            "96", "59")
#
#
##cl.center.df <- read.csv('cl.center.df.csv')
##knn.cl.df <- read.csv('knn.cl.df.csv')
#
#norm.dat<-readRDS('norm.dat.v2.rds')
#
#DE <- readRDS('Brain_subClusteringID_seurat_merge_DE20.rds')
#
#norm.dat2 <- norm.dat[as.vector(DE$markers),]
#
#cl <- meta$hicat_cluster_merge_newID
#names(cl) <- rownames(meta)
#cl.df <-  data.frame(cluster_id = unique(as.vector(cl[colnames(norm.dat2)])),
#                     cluster_label = unique(as.vector(cl[colnames(norm.dat2)])),
#                     cluster_color = supertype_color[unique(as.vector(cl[colnames(norm.dat2)]))])
#
#rd.dat <- t(as.matrix(norm.dat2))
#
#saveRDS(rd.dat,'rd.dat.rds')
#saveRDS(cl,'cl.rds')
#saveRDS(cl.df,'cl.df.rds')


rd.dat<-readRDS('rd.dat.rds')
cl<-readRDS('cl.rds')
cl.df<-readRDS('cl.df.rds')


a <- get_knn_graph(rd.dat, cl, cl.df)
saveRDS(a,'get_knn_graph.rds')




a<-readRDS('get_knn_graph.rds')
knn.cl.df2 <- a$knn.cl.df


cl.center.df2 <- data.frame(x = rep(0,length(supertype_color)),
                            y = rep(0,length(supertype_color)),
                            cl = names(supertype_color),
                            cluster_color = supertype_color,
                            cluster_size = as.vector(table(seu$hicat_cluster_merge_newID)[names(supertype_color)]))

for(i in 1:length(supertype_color)){
  x1 <- mean(seu@reductions$umap@cell.embeddings[seu$hicat_cluster_merge_newID==names(supertype_color)[i],1])
  y1 <- mean(seu@reductions$umap@cell.embeddings[seu$hicat_cluster_merge_newID==names(supertype_color)[i],2])
  cl.center.df2[i,1:2] <- c(x1,y1)
}

anno <- unique(cbind(as.vector(seu$hicat_cluster_merge_newID),seu$hicat_cluster_subclasses))
colnames(anno) <- c('cluster','anno')
rownames(anno) <- anno[,1]

anno_color <- cbind(c("OPC", "Oligo L4-L6 MOBP", "Oligo L4-L6 OPALIN", "Astrocyte FGFR3", "NP", "L6CT/FEZF2", "L6B/FEZF2", "L5 IT/RORB", "L1-L3 IT LINC00507", "ET", 
                      "L6 CAR3", "L6/IT", "LAMP5/SV2C", "PAX6", "LAMP5/SNCG", "VIP", "SST", "PVALB", "Microglia", "VLMC L1-L3", 
                      "ENDO L2-L5"),
                    c("#2E3E39", "#5E8A79", "#72B09C", "#665C47", "#93C43B", "#2D8CB8", "#7944AA", "#09B2B2", "#07D8D8", "#19E3BE", 
                      "#898325", "#CEC823", "#f2798d", "#E96DF2", "#D96C07", "#9C28CC", "#F2B003", "#FF2D4E", "#513577", "#697255", 
                      "#604B47"))

rownames(anno_color) <- anno_color[,1]
colnames(anno_color) <- c('anno','color')
anno_color <- as.data.frame(anno_color)
anno_color$id <- c(1:21)

cl.center.df2$clade_label <- anno[cl.center.df2$cl,2]
cl.center.df2$clade_id <- anno_color[cl.center.df2$clade_label,3]
cl.center.df2$clade_color <- anno_color[cl.center.df2$clade_label,2]

cl.center.df2$cluster_id <- cl.center.df2$cl
cl.center.df2$cluster_label <- cl.center.df2$cl

plot_constellation(knn.cl.df = knn.cl.df2, cl.center.df = cl.center.df2, out.dir = "plot_constellation", node.dodge=TRUE, plot.hull=c(1,2), label_repel=TRUE) 





