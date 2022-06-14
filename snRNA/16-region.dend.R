#library(tasic2016data)
library(Seurat)
library(scrattch.hicat)
library(dendextend)
library(dplyr)
library(matrixStats)
library(Matrix)

##############################################
################### Exc ######################
##############################################
seu <- readRDS('seu.harmony.Glu.rds')

seu$hicat_cluster_subclasses[seu$hicat_cluster_merge_newID=='107'] <- 'L6B/FEZF2'
seu$hicat_cluster_subclasses[seu$hicat_cluster_merge_newID=='80'] <- 'L6CT/FEZF2'

region_color <- c('#3F4587','#8562AA','#EC8561','#B97CB5','#D43046','#F0592B','#ED4A96','#593C97','#A54486','#FBDE13','#299FAE','#75CCE3','#0C6939','#0D9547')
region <- c('FPPFC','DLPFC','VLPFC','M1','S1','S1E','PoCG','SPL','SMG','AG','VISP','ITG','STG','ACC')
regionID <- c('A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','A11','A12','A13','A14')
names(region_color) <- region
names(region) <- regionID

seu$regionName <- region[seu$region]
seu$regionName <- factor(seu$regionName, levels = region)

region.clean <- seu$regionName


### dend ###
cl.med <- get_cl_medians(seu@assays$RNA@data, 
                         region.clean)
dend.result <- build_dend(cl.med,
                          #l.rank, 
                          #l.color,
                          nboot = 100)
saveRDS(dend.result, 'Exc.dend.result.rds')

# Exc_dend.pdf 3x5
plot(dend.result$dend)

dend.result <- readRDS('Exc.dend.result.rds')



### cor heatmap ###
library(pheatmap)
library(RColorBrewer)
subclasses <- c('NP', 'L6CT/FEZF2', 'L6B/FEZF2', 'L5 IT/RORB', 'L1-L3 IT LINC00507', 'ET', 'L6 CAR3', 'L6/IT')
for(i in 1:length(subclasses)){
  norm.dat <- seu@assays$RNA@data[,seu$hicat_cluster_subclasses==subclasses[i]]
  cl.clean <- seu$regionName[seu$hicat_cluster_subclasses==subclasses[i]]
  tmp.med <- get_cl_medians(norm.dat, 
                            cl.clean)
  cl.cor <- cor(tmp.med)
  
  cl.cor <- cl.cor[rev(labels(dend.result$dend)),labels(dend.result$dend)]
  
  pdf(paste0('Exc_cor_',i,'.pdf'), width = 4.5, height = 4)
  pheatmap(cl.cor, 
           cluster_rows = F,
           cluster_cols = F,
           main = subclasses[i],
           color = brewer.pal(9,'Blues')[3:9])
  dev.off()
}

aa <- rep(1,14)
names(aa) <- labels(dend.result$dend)
# Exc_dend_region_color.pdf 4x5
barplot(aa,col = region_color[labels(dend.result$dend)])


### confusionMatrix ###
library(caret)
subclasses <- c('NP', 'L6CT/FEZF2', 'L6B/FEZF2', 'L5 IT/RORB', 'L1-L3 IT LINC00507', 'ET', 'L6 CAR3', 'L6/IT')
#allcells <- seu$regionName
#for(i in 1:length(subclasses)){
#  this_cells <- seu$regionName[seu$hicat_cluster_subclasses==subclasses[i]]
#  
#  rr0 <- table(allcells)*(length(this_cells)/length(allcells))
#  rr1 <- round(table(allcells)*(length(this_cells)/length(allcells)))
#  
#  ref_cells <- c()
#  if(round(sum(rr0)-sum(rr1))>0){
#    nn <- sum(rr0)-sum(rr1)
#    dd <- rr0-rr1
#    ss <- sort(dd,decreasing = T)
#    
#    rr2 <- rr1
#    rr2[names(ss)[1:nn]] <- rr2[names(ss)[1:nn]] + 1
#    
#    for(j in 1:length(rr2)){
#      ref_cells <- c(ref_cells, rep(names(rr2)[j],rr2[j]))
#    }
#  }
#  else if(round(sum(rr0)-sum(rr1))<0){
#    nn <- sum(rr1)-sum(rr0)
#    dd <- rr1-rr0
#    ss <- sort(dd,decreasing = T)
#    
#    rr2 <- rr1
#    rr2[names(ss)[1:nn]] <- rr2[names(ss)[1:nn]] - 1
#    
#    for(j in 1:length(rr2)){
#      ref_cells <- c(ref_cells, rep(names(rr2)[j],rr2[j]))
#    }
#  }
#  else{
#    for(j in 1:length(rr1)){
#      ref_cells <- c(ref_cells, rep(names(rr1)[j],rr1[j]))
#    }
#  }
#  
#  ref_cells <- factor(ref_cells, levels = levels(this_cells))
#  
#  CC <- confusionMatrix(this_cells,ref_cells)  
#  
#  MM0 <- CC$table[rev(labels(dend.result$dend)),labels(dend.result$dend)]
#  
#  MM1 <- MM0/rowSums(MM0)
#  
#  pdf(paste0('Exc_confusion_',i,'.pdf'), width = 4.5, height = 4)
#  pheatmap(MM1, 
#           cluster_rows = F,
#           cluster_cols = F,
#           main = subclasses[i],
#           color = brewer.pal(9,'Purples')[3:9])
#  dev.off()
#}

region <- c('FPPFC','DLPFC','VLPFC','M1','S1','S1E','PoCG','SPL','SMG','AG','VISP','ITG','STG','ACC')
regionID <- c('A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','A11','A12','A13','A14')
names(region) <- regionID

subclasses <- c('NP', 'L6CT_FEZF2', 'L6B_FEZF2', 'L5_IT_RORB', 'L1-L3_IT_LINC00507', 'ET', 'L6_CAR3', 'L6_IT')
for(i in 1:length(subclasses)){
  file <- paste0('Exc.',subclasses[i],'.cm.rds')
  cm <- readRDS(file)
  
  colnames(cm) <- region[colnames(cm)]
  rownames(cm) <- region[rownames(cm)]
  
  cm1 <- cm[rev(labels(dend.result$dend)),labels(dend.result$dend)]
  cm2 <- cm1/rowSums(cm1)
  
  pdf(paste0('Exc_confusion_',i,'.pdf'), width = 4.5, height = 4)
  pheatmap(cm2, 
           cluster_rows = F,
           cluster_cols = F,
           main = subclasses[i],
           color = brewer.pal(9,'Purples')[3:9])
  dev.off()
}





#### DE genes ################
subclasses <- c('NP', 'L6CT/FEZF2', 'L6B/FEZF2', 'L5 IT/RORB', 'L1-L3 IT LINC00507', 
                'ET', 'L6 CAR3', 'L6/IT')
region <- c('FPPFC','DLPFC','VLPFC','M1','S1','S1E','PoCG','SPL','SMG','AG','VISP','ITG','STG','ACC')
DE <- c()
for(i in 1:length(subclasses)){
  for(j in 1:length(region)){
    cells.1 <- colnames(seu)[seu$hicat_cluster_subclasses == subclasses[i] & seu$regionName == region[j]]
    cells.2 <- colnames(seu)[seu$hicat_cluster_subclasses == subclasses[i] & seu$regionName != region[j]]
    if(length(cells.1) < 3){
      next
    }
    tmp <- FindMarkers(seu,
                       ident.1=cells.1,
                       ident.2=cells.2,
                       min.pct = 0.2)
    tmp$subclasses <- rep(subclasses[i],dim(tmp)[1])
    tmp$region <- rep(region[j],dim(tmp)[1])
    tmp$gene <- rownames(tmp)
    DE <- rbind(DE,tmp)
  }
}
saveRDS(DE, 'Exc_subclass_region_DE.rds')

DE <-readRDS('Exc_subclass_region_DE.rds')
genes <- c()
for(i in 1:length(subclasses)){
  tmp <- DE[DE$subclasses==subclasses[i],]
  top1 <- top_n(group_by(tmp,region), 1, avg_log2FC)
  genes <- top1$gene
  
  seu_markers_m<-seu@assays$RNA@data[rownames(seu@assays$RNA@data) %in% genes,seu$hicat_cluster_subclasses == subclasses[i]]
  seu_markers_m<-t(as.matrix(seu_markers_m))
  seu_markers_m<-as.data.frame(seu_markers_m)
  seu_markers_m<-cbind(seu_markers_m,as.vector(seu$regionName[seu$hicat_cluster_subclasses == subclasses[i]]))
  colnames(seu_markers_m)[dim(seu_markers_m)[2]]<-'cluster'
  #seu_markers_m$cluster<-as.vector(seu_markers_m$cluster)
  seu_markers_m2<-aggregate(seu_markers_m[,1:(dim(seu_markers_m)[2]-1)],by=list(cluster=seu_markers_m$cluster),FUN=mean)
  seu_markers_m3<-seu_markers_m2[,-1]
  
  seu_markers_m4<-matrix(unlist(seu_markers_m3),dim(seu_markers_m3)[1],dim(seu_markers_m3)[2])
  colnames(seu_markers_m4)<-colnames(seu_markers_m3)
  rownames(seu_markers_m4)<-as.vector(seu_markers_m2[,1])
  seu_markers_m4<-log2(seu_markers_m4+1)
  seu_markers_m4 <- seu_markers_m4[labels(dend.result$dend),]
  pdf(paste0('Exc_heatmap_',i,'.pdf'), width = 6, height = 4.5)
  pheatmap(t(seu_markers_m4),
           cluster_rows=T,
           cluster_cols=F,
           main = subclasses[i],
           color = brewer.pal(9,'Oranges')[1:8],
           angle_col = 90)
  dev.off()
}





##############################################
############# Inh CGE ########################
##############################################
seu <- readRDS('seu.harmony.GABA.CGE.rds')

region_color <- c('#3F4587','#8562AA','#EC8561','#B97CB5','#D43046','#F0592B','#ED4A96','#593C97','#A54486','#FBDE13','#299FAE','#75CCE3','#0C6939','#0D9547')
region <- c('FPPFC','DLPFC','VLPFC','M1','S1','S1E','PoCG','SPL','SMG','AG','VISP','ITG','STG','ACC')
regionID <- c('A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','A11','A12','A13','A14')
names(region_color) <- region
names(region) <- regionID

seu$regionName <- region[seu$region]
seu$regionName <- factor(seu$regionName, levels = region)

region.clean <- seu$regionName


### dend ###
cl.med <- get_cl_medians(seu@assays$RNA@data, 
                         region.clean)
dend.result <- build_dend(cl.med,
                          #l.rank, 
                          #l.color,
                          nboot = 100)
saveRDS(dend.result, 'Inh.CGE.dend.result.rds')


# Inh_CGE_dend.pdf 3x5
plot(dend.result$dend)


dend.result<-readRDS('Inh.CGE.dend.result.rds')


### cor heatmap ###
library(pheatmap)
library(RColorBrewer)
subclasses <- c('LAMP5/SNCG', 'LAMP5/SV2C', 'PAX6', 'VIP')
for(i in 1:length(subclasses)){
  norm.dat <- seu@assays$RNA@data[,seu$hicat_cluster_subclasses==subclasses[i]]
  cl.clean <- seu$regionName[seu$hicat_cluster_subclasses==subclasses[i]]
  tmp.med <- get_cl_medians(norm.dat, 
                            cl.clean)
  cl.cor <- cor(tmp.med)
  
  cl.cor <- cl.cor[rev(labels(dend.result$dend)),labels(dend.result$dend)]
  
  pdf(paste0('Inh_CGE_cor_',i,'.pdf'), width = 4.5, height = 4)
  pheatmap(cl.cor, 
           cluster_rows = F,
           cluster_cols = F,
           main = subclasses[i],
           color = brewer.pal(9,'Blues')[3:9])
  dev.off()
}

aa <- rep(1,14)
names(aa) <- labels(dend.result$dend)
# Inh_CGE_dend_region_color.pdf 4x5
barplot(aa,col = region_color[labels(dend.result$dend)])


### confusionMatrix ###
library(caret)
subclasses <- c('LAMP5/SNCG', 'LAMP5/SV2C', 'PAX6', 'VIP')
#allcells <- seu$regionName
#for(i in 1:length(subclasses)){
#  this_cells <- seu$regionName[seu$hicat_cluster_subclasses==subclasses[i]]
#  
#  rr0 <- table(allcells)*(length(this_cells)/length(allcells))
#  rr1 <- round(table(allcells)*(length(this_cells)/length(allcells)))
#  
#  ref_cells <- c()
#  if(round(sum(rr0)-sum(rr1))>0){
#    nn <- round(sum(rr0)-sum(rr1))
#    dd <- rr0-rr1
#    ss <- sort(dd,decreasing = T)
#    
#    rr2 <- rr1
#    rr2[names(ss)[1:nn]] <- rr2[names(ss)[1:nn]] + 1
#    
#    for(j in 1:length(rr2)){
#      ref_cells <- c(ref_cells, rep(names(rr2)[j],rr2[j]))
#    }
#  } else if(round(sum(rr0)-sum(rr1))<0){
#    nn <- sum(rr1)-sum(rr0)
#    dd <- rr1-rr0
#    ss <- sort(dd,decreasing = T)
#    
#    rr2 <- rr1
#    rr2[names(ss)[1:nn]] <- rr2[names(ss)[1:nn]] - 1
#    
#    for(j in 1:length(rr2)){
#      ref_cells <- c(ref_cells, rep(names(rr2)[j],rr2[j]))
#    }
#  } else{
#    for(j in 1:length(rr1)){
#      ref_cells <- c(ref_cells, rep(names(rr1)[j],rr1[j]))
#    }
#  }
#    
#  
#  ref_cells <- factor(ref_cells, levels = levels(this_cells))
#  
#  CC <- confusionMatrix(this_cells,ref_cells)  
#  
#  MM0 <- CC$table[rev(labels(dend.result$dend)),labels(dend.result$dend)]
#  
#  MM1 <- MM0/rowSums(MM0)
#  
#  pdf(paste0('Inh_CGE_confusion_',i,'.pdf'), width = 4.5, height = 4)
#  pheatmap(MM1, 
#           cluster_rows = F,
#           cluster_cols = F,
#           main = subclasses[i],
#           color = brewer.pal(9,'Purples')[3:9])
#  dev.off()
#}

region <- c('FPPFC','DLPFC','VLPFC','M1','S1','S1E','PoCG','SPL','SMG','AG','VISP','ITG','STG','ACC')
regionID <- c('A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','A11','A12','A13','A14')
names(region) <- regionID

subclasses <- c('LAMP5_SNCG', 'LAMP5_SV2C', 'PAX6', 'VIP')
for(i in 1:length(subclasses)){
  file <- paste0('Inh.CGE.',subclasses[i],'.cm.rds')
  cm <- readRDS(file)
  
  colnames(cm) <- region[colnames(cm)]
  rownames(cm) <- region[rownames(cm)]
  
  cm1 <- cm[rev(labels(dend.result$dend)),labels(dend.result$dend)]
  cm2 <- cm1/rowSums(cm1)
  
  pdf(paste0('Inh_CGE_confusion_',i,'.pdf'), width = 4.5, height = 4)
  pheatmap(cm2, 
           cluster_rows = F,
           cluster_cols = F,
           main = subclasses[i],
           color = brewer.pal(9,'Purples')[3:9])
  dev.off()
}


#### DE genes ################
subclasses <- c('LAMP5/SNCG', 'LAMP5/SV2C', 'PAX6', 'VIP')
region <- c('FPPFC','DLPFC','VLPFC','M1','S1','S1E','PoCG','SPL','SMG','AG','VISP','ITG','STG','ACC')
DE <- c()
for(i in 1:length(subclasses)){
  for(j in 1:length(region)){
    cells.1 <- colnames(seu)[seu$hicat_cluster_subclasses == subclasses[i] & seu$regionName == region[j]]
    cells.2 <- colnames(seu)[seu$hicat_cluster_subclasses == subclasses[i] & seu$regionName != region[j]]
    if(length(cells.1) < 3){
      next
    }
    tmp <- FindMarkers(seu,
                       ident.1=cells.1,
                       ident.2=cells.2,
                       min.pct = 0.2)
    tmp$subclasses <- rep(subclasses[i],dim(tmp)[1])
    tmp$region <- rep(region[j],dim(tmp)[1])
    tmp$gene <- rownames(tmp)
    DE <- rbind(DE,tmp)
  }
}
saveRDS(DE, 'Inh_CGE_subclass_region_DE.rds')

genes <- c()
for(i in 1:length(subclasses)){
  tmp <- DE[DE$subclasses==subclasses[i],]
  top1 <- top_n(group_by(tmp,region), 1, avg_log2FC)
  genes <- top1$gene
  
  seu_markers_m<-seu@assays$RNA@data[rownames(seu@assays$RNA@data) %in% genes,seu$hicat_cluster_subclasses == subclasses[i]]
  seu_markers_m<-t(as.matrix(seu_markers_m))
  seu_markers_m<-as.data.frame(seu_markers_m)
  seu_markers_m<-cbind(seu_markers_m,as.vector(seu$regionName[seu$hicat_cluster_subclasses == subclasses[i]]))
  colnames(seu_markers_m)[dim(seu_markers_m)[2]]<-'cluster'
  #seu_markers_m$cluster<-as.vector(seu_markers_m$cluster)
  seu_markers_m2<-aggregate(seu_markers_m[,1:(dim(seu_markers_m)[2]-1)],by=list(cluster=seu_markers_m$cluster),FUN=mean)
  seu_markers_m3<-seu_markers_m2[,-1]
  
  seu_markers_m4<-matrix(unlist(seu_markers_m3),dim(seu_markers_m3)[1],dim(seu_markers_m3)[2])
  colnames(seu_markers_m4)<-colnames(seu_markers_m3)
  rownames(seu_markers_m4)<-as.vector(seu_markers_m2[,1])
  seu_markers_m4<-log2(seu_markers_m4+1)
  seu_markers_m4 <- seu_markers_m4[labels(dend.result$dend),]
  pdf(paste0('Inh_CGE_heatmap_',i,'.pdf'), width = 6, height = 4.5)
  pheatmap(t(seu_markers_m4),
           cluster_rows=T,
           cluster_cols=F,
           main = subclasses[i],
           color = brewer.pal(9,'Oranges')[1:8],
           angle_col = 90)
  dev.off()
}






##############################################
############# Inh MGE ########################
##############################################
seu <- readRDS('seu.harmony.GABA.MGE.rds')

region_color <- c('#3F4587','#8562AA','#EC8561','#B97CB5','#D43046','#F0592B','#ED4A96','#593C97','#A54486','#FBDE13','#299FAE','#75CCE3','#0C6939','#0D9547')
region <- c('FPPFC','DLPFC','VLPFC','M1','S1','S1E','PoCG','SPL','SMG','AG','VISP','ITG','STG','ACC')
regionID <- c('A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','A11','A12','A13','A14')
names(region_color) <- region
names(region) <- regionID

seu$regionName <- region[seu$region]
seu$regionName <- factor(seu$regionName, levels = region)

region.clean <- seu$regionName


### dend ###
cl.med <- get_cl_medians(seu@assays$RNA@data, 
                         region.clean)
dend.result <- build_dend(cl.med,
                          #l.rank, 
                          #l.color,
                          nboot = 100)
saveRDS(dend.result, 'Inh.MGE.dend.result.rds')


# Inh_MGE_dend.pdf 3x5
plot(dend.result$dend)



dend.result <- readRDS('Inh.MGE.dend.result.rds')



### cor heatmap ###
library(pheatmap)
library(RColorBrewer)
subclasses <- c('PVALB', 'SST')
for(i in 1:length(subclasses)){
  norm.dat <- seu@assays$RNA@data[,seu$hicat_cluster_subclasses==subclasses[i]]
  cl.clean <- seu$regionName[seu$hicat_cluster_subclasses==subclasses[i]]
  tmp.med <- get_cl_medians(norm.dat, 
                            cl.clean)
  cl.cor <- cor(tmp.med)
  
  cl.cor <- cl.cor[rev(labels(dend.result$dend)),labels(dend.result$dend)]
  
  pdf(paste0('Inh_MGE_cor_',i,'.pdf'), width = 4.5, height = 4)
  pheatmap(cl.cor, 
           cluster_rows = F,
           cluster_cols = F,
           main = subclasses[i],
           color = brewer.pal(9,'Blues')[3:9])
  dev.off()
}

aa <- rep(1,14)
names(aa) <- labels(dend.result$dend)
# Inh_MGE_dend_region_color.pdf 4x5
barplot(aa,col = region_color[labels(dend.result$dend)])


### confusionMatrix ###
library(caret)
subclasses <- c('PVALB', 'SST')
#allcells <- seu$regionName
#for(i in 1:length(subclasses)){
#  this_cells <- seu$regionName[seu$hicat_cluster_subclasses==subclasses[i]]
#  
#  rr0 <- table(allcells)*(length(this_cells)/length(allcells))
#  rr1 <- round(table(allcells)*(length(this_cells)/length(allcells)))
#  
#  ref_cells <- c()
#  if(round(sum(rr0)-sum(rr1))>0){
#    nn <- sum(rr0)-sum(rr1)
#    dd <- rr0-rr1
#    ss <- sort(dd,decreasing = T)
#    
#    rr2 <- rr1
#    rr2[names(ss)[1:nn]] <- rr2[names(ss)[1:nn]] + 1
#    
#    for(j in 1:length(rr2)){
#      ref_cells <- c(ref_cells, rep(names(rr2)[j],rr2[j]))
#    }
#  }
#  else if(round(sum(rr0)-sum(rr1))<0){
#    nn <- sum(rr1)-sum(rr0)
#    dd <- rr1-rr0
#    ss <- sort(dd,decreasing = T)
#    
#    rr2 <- rr1
#    rr2[names(ss)[1:nn]] <- rr2[names(ss)[1:nn]] - 1
#    
#    for(j in 1:length(rr2)){
#      ref_cells <- c(ref_cells, rep(names(rr2)[j],rr2[j]))
#    }
#  }
#  else{
#    for(j in 1:length(rr1)){
#      ref_cells <- c(ref_cells, rep(names(rr1)[j],rr1[j]))
#    }
#  }
#  
#  ref_cells <- factor(ref_cells, levels = levels(this_cells))
#  
#  CC <- confusionMatrix(this_cells,ref_cells)  
#  
#  MM0 <- CC$table[rev(labels(dend.result$dend)),labels(dend.result$dend)]
#  
#  MM1 <- MM0/rowSums(MM0)
#  
#  pdf(paste0('Inh_MGE_confusion_',i,'.pdf'), width = 4.5, height = 4)
#  pheatmap(MM1, 
#           cluster_rows = F,
#           cluster_cols = F,
#           main = subclasses[i],
#           color = brewer.pal(9,'Purples')[3:9])
#  dev.off()
#}

subclasses <- c('PVALB', 'SST')
for(i in 1:length(subclasses)){
  file <- paste0('Inh.MGE.',subclasses[i],'.cm.rds')
  cm <- readRDS(file)
  
  colnames(cm) <- region[colnames(cm)]
  rownames(cm) <- region[rownames(cm)]
  
  cm1 <- cm[rev(labels(dend.result$dend)),labels(dend.result$dend)]
  cm2 <- cm1/rowSums(cm1)
  
  pdf(paste0('Inh_MGE_confusion_',i,'.pdf'), width = 4.5, height = 4)
  pheatmap(cm2, 
           cluster_rows = F,
           cluster_cols = F,
           main = subclasses[i],
           color = brewer.pal(9,'Purples')[3:9])
  dev.off()
}







#### DE genes ################
subclasses <- c('PVALB', 'SST')
region <- c('FPPFC','DLPFC','VLPFC','M1','S1','S1E','PoCG','SPL','SMG','AG','VISP','ITG','STG','ACC')
DE <- c()
for(i in 1:length(subclasses)){
  for(j in 1:length(region)){
    cells.1 <- colnames(seu)[seu$hicat_cluster_subclasses == subclasses[i] & seu$regionName == region[j]]
    cells.2 <- colnames(seu)[seu$hicat_cluster_subclasses == subclasses[i] & seu$regionName != region[j]]
    if(length(cells.1) < 3){
      next
    }
    tmp <- FindMarkers(seu,
                       ident.1=cells.1,
                       ident.2=cells.2,
                       min.pct = 0.2)
    tmp$subclasses <- rep(subclasses[i],dim(tmp)[1])
    tmp$region <- rep(region[j],dim(tmp)[1])
    tmp$gene <- rownames(tmp)
    DE <- rbind(DE,tmp)
  }
}
saveRDS(DE, 'Inh_MGE_subclass_region_DE.rds')

genes <- c()
for(i in 1:length(subclasses)){
  tmp <- DE[DE$subclasses==subclasses[i],]
  top1 <- top_n(group_by(tmp,region), 1, avg_log2FC)
  genes <- top1$gene
  
  seu_markers_m<-seu@assays$RNA@data[rownames(seu@assays$RNA@data) %in% genes,seu$hicat_cluster_subclasses == subclasses[i]]
  seu_markers_m<-t(as.matrix(seu_markers_m))
  seu_markers_m<-as.data.frame(seu_markers_m)
  seu_markers_m<-cbind(seu_markers_m,as.vector(seu$regionName[seu$hicat_cluster_subclasses == subclasses[i]]))
  colnames(seu_markers_m)[dim(seu_markers_m)[2]]<-'cluster'
  #seu_markers_m$cluster<-as.vector(seu_markers_m$cluster)
  seu_markers_m2<-aggregate(seu_markers_m[,1:(dim(seu_markers_m)[2]-1)],by=list(cluster=seu_markers_m$cluster),FUN=mean)
  seu_markers_m3<-seu_markers_m2[,-1]
  
  seu_markers_m4<-matrix(unlist(seu_markers_m3),dim(seu_markers_m3)[1],dim(seu_markers_m3)[2])
  colnames(seu_markers_m4)<-colnames(seu_markers_m3)
  rownames(seu_markers_m4)<-as.vector(seu_markers_m2[,1])
  seu_markers_m4<-log2(seu_markers_m4+1)
  seu_markers_m4 <- seu_markers_m4[labels(dend.result$dend),]
  pdf(paste0('Inh_MGE_heatmap_',i,'.pdf'), width = 6, height = 4.5)
  pheatmap(t(seu_markers_m4),
           cluster_rows=T,
           cluster_cols=F,
           main = subclasses[i],
           color = brewer.pal(9,'Oranges')[1:8],
           angle_col = 90)
  dev.off()
}



