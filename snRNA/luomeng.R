

a <- read.csv('luomeng_Exc.csv', header = T)
rownames(a) <- a$cells
table(a$clusterID,a$type)
library(pheatmap)
#cluster_data.pdf 30x3
pheatmap(table(a$clusterID,a$type),display_numbers=T,number_format='%.0f')



b <- read.csv('metadata.csv',header = T)

c <- cbind(b$sample_name,b$cluster_label)
colnames(c) <- c('cells','label')
rownames(c) <- c[,1]

aa1 <- a[a$type=='allen',]
aa1$annotation <- c[aa1$cells,2]
aaa1 <- table(aa1$clusterID,aa1$annotation)
#cluster_allen_annotation.pdf 25x12
pheatmap(log2(aaa1+1),display_numbers=T,number_format='%.0f')



result<-readRDS('seu.hicat.de80.niter100.rds')

aa2 <- a[a$type=='us',]
aa2$cluster_de80 <- result$cl.result$cl[aa2$cells]
aaa2 <- table(aa2$clusterID,aa2$cluster_de80)
#cluster_us_annotation.pdf 25x12
pheatmap(log2(aaa2+1),display_numbers=T,number_format='%.0f')


table(rownames(aaa2) %in% rownames(aaa1))

common <- rownames(aaa1)[rownames(aaa1) %in% rownames(aaa2)]
cc1 <- rownames(aaa1)[!rownames(aaa1) %in% rownames(aaa2)]
cc2 <- rownames(aaa2)[!rownames(aaa2) %in% rownames(aaa1)]

ccc1 <- c(common,cc1)
ccc2 <- c(common,cc2)

aaa1 <- table(aa1$clusterID,aa1$annotation)
aaa2 <- table(aa2$clusterID,aa2$cluster_de80)

aaa1 <- aaa1[ccc1,]
aaa2 <- aaa2[ccc2,]
#cluster_allen_annotation_Exc.pdf 8x12
pheatmap(log2(aaa1+1),cluster_rows = F,display_numbers=T,number_format='%.0f')

clusterOrder <- rev(c("25","30","29","31","58","425","68","70","40","36","39","51","42",
                      "48","49","72","93","125","111","104","109","98","115",
                      "123","124","88","65","77","86","89"))
aaa2 <- aaa2[,clusterOrder]
#cluster_us_annotation_Exc.pdf 8x8
pheatmap(log2(aaa2+1),cluster_rows = F,cluster_cols = F,display_numbers=T,number_format='%.0f')




########### Inh ###################################################
a <- read.csv('luomeng_Inh.csv', header = T)
rownames(a) <- a$cells
table(a$clusterID,a$type)
library(pheatmap)
#cluster_data.pdf 30x3
pheatmap(table(a$clusterID,a$type),display_numbers=T,number_format='%.0f')



b <- read.csv('metadata.csv',header = T)

c <- cbind(b$sample_name,b$cluster_label)
colnames(c) <- c('cells','label')
rownames(c) <- c[,1]

aa1 <- a[a$type=='allen',]
aa1$annotation <- c[aa1$cells,2]
aaa1 <- table(aa1$clusterID,aa1$annotation)
#cluster_allen_annotation.pdf 25x12
pheatmap(log2(aaa1+1),display_numbers=T,number_format='%.0f')



result<-readRDS('seu.hicat.de80.niter100.rds')
dend.result<-readRDS('seu.hicat.de80.niter100.dend.rds')

aa2 <- a[a$type=='us',]
aa2$cluster_de80 <- result$cl.result$cl[aa2$cells]
aaa2 <- table(aa2$clusterID,aa2$cluster_de80)
#cluster_us_annotation.pdf 25x12
pheatmap(log2(aaa2+1),display_numbers=T,number_format='%.0f')


table(rownames(aaa2) %in% rownames(aaa1))

common <- rownames(aaa1)[rownames(aaa1) %in% rownames(aaa2)]
cc1 <- rownames(aaa1)[!rownames(aaa1) %in% rownames(aaa2)]
cc2 <- rownames(aaa2)[!rownames(aaa2) %in% rownames(aaa1)]

ccc1 <- c(common,cc1)
ccc2 <- c(common,cc2)

aaa1 <- table(aa1$clusterID,aa1$annotation)
aaa2 <- table(aa2$clusterID,aa2$cluster_de80)

aaa1 <- aaa1[ccc1,]
aaa2 <- aaa2[ccc2,]
#cluster_allen_annotation_Inh.pdf 12x10
pheatmap(log2(aaa1+1), cluster_rows = F, display_numbers=T, number_format='%.0f')

clusterOrder <- rev(c("297","284","380","342","368","334","373","339","267","384","327","329","333","299","301","309",
                  "303","311","383","377","375","378","268","336","127","160","162","148","144","156","146","157",
                  "172","173","190","175","180","185","206",
                  "204","211","201","202","198","199","192","226","214","217","220","224","270","256",
                  "245","235","243","264","282","222","279","280","229",
                  "247","232","233"))
aaa2 <- aaa2[,clusterOrder]
#cluster_us_annotation_Inh.pdf 12x12
pheatmap(log2(aaa2+1), cluster_rows = F, cluster_cols = F, display_numbers=T, number_format='%.0f')



#############################################################
