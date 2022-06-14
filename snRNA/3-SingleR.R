library(Seurat)
library(dplyr)
library(pheatmap)
library(SingleCellExperiment)
library(ggplot2)

library(SingleR)
library(celldex)

sample <- c('S0128-A2','S0128-A6','S0128-A7','S0128-A8','S0128-A9','S0128-A10','S0128-A11','S0128-A12','S0128-A13','S0128-A14'
            ,'S0206-A1','S0206-A2','S0206-A3','S0206-A4','S0206-A5','S0206-A6','S0206-A7','S0206-A11','S0206-A12'
            ,'S0406-A1','S0406-A3','S0406-A4','S0406-A5','S0406-A6', 'S0406-A7', 'S0406-A8', 'S0406-A9', 'S0406-A10','S0406-A14'
            ,'S0426-A5', 'S0426-A8', 'S0426-A9', 'S0426-A10','S0426-A11','S0426-A12','S0426-A13','S0426-A14'
            ,'S0531-A1', 'S0531-A2', 'S0531-A3','S0531-A4', 'S0531-A13',
            'S0916-A2','S0916-A6','S0916-A7','S0916-A8','S0916-A9','S0916-A10','S0916-A11','S0916-A12','S0916-A13','S0916-A14')

#ref <- BlueprintEncodeData()
#for(i in 1:length(sample)){
#  file <- paste0('counts/',sample[i],'.counts.rds')
#  seu <- readRDS(file)
#  seu.singleR<-SingleR(test=seu, 
#                       ref=ref,
#                       assay.type.test=1,
#                       labels = ref$label.fine)
#  outfile <- paste0(sample[i],'.singleR.blue.rds')
#  saveRDS(seu.singleR, outfile)
#}

#ref <- DatabaseImmuneCellExpressionData()
#for(i in 1:length(sample)){
#  file <- paste0('counts/',sample[i],'.counts.rds')
#  seu <- readRDS(file)
#  seu.singleR<-SingleR(test=seu, 
#                       ref=ref,
#                       assay.type.test=1,
#                       labels = ref$label.fine)
#  outfile <- paste0(sample[i],'.singleR.data.rds')
#  saveRDS(seu.singleR, outfile)
#}

ref <- HumanPrimaryCellAtlasData()
for(i in 1:length(sample)){
  file <- paste0('counts/',sample[i],'.counts.rds')
  seu <- readRDS(file)
  seu.singleR<-SingleR(test=seu, 
                       ref=ref,
                       assay.type.test=1,
                       labels = ref$label.fine)
  outfile <- paste0(sample[i],'.singleR.human.rds')
  saveRDS(seu.singleR, outfile)
}

ref <- ImmGenData()
for(i in 1:length(sample)){
  file <- paste0('counts/',sample[i],'.counts.rds')
  seu <- readRDS(file)
  seu.singleR<-SingleR(test=seu, 
                       ref=ref,
                       assay.type.test=1,
                       labels = ref$label.fine)
  outfile <- paste0(sample[i],'.singleR.imm.rds')
  saveRDS(seu.singleR, outfile)
}

ref <- MonacoImmuneData()
for(i in 1:length(sample)){
  file <- paste0('counts/',sample[i],'.counts.rds')
  seu <- readRDS(file)
  seu.singleR<-SingleR(test=seu, 
                       ref=ref,
                       assay.type.test=1,
                       labels = ref$label.fine)
  outfile <- paste0(sample[i],'.singleR.monaco.rds')
  saveRDS(seu.singleR, outfile)
}

ref <- NovershternHematopoieticData()
for(i in 1:length(sample)){
  file <- paste0('counts/',sample[i],'.counts.rds')
  seu <- readRDS(file)
  seu.singleR<-SingleR(test=seu, 
                       ref=ref,
                       assay.type.test=1,
                       labels = ref$label.fine)
  outfile <- paste0(sample[i],'.singleR.nover.rds')
  saveRDS(seu.singleR, outfile)
}


anno <- data.frame()
for(i in 1:length(sample)){
  file1 <- paste0('singleR/',sample[i],'.singleR.blue.rds')
  file2 <- paste0('singleR/',sample[i],'.singleR.data.rds')
  file3 <- paste0('singleR/',sample[i],'.singleR.human.rds')
  file4 <- paste0('singleR/',sample[i],'.singleR.imm.rds')
  file5 <- paste0('singleR/',sample[i],'.singleR.monaco.rds')
  file6 <- paste0('singleR/',sample[i],'.singleR.nover.rds')
 
  anno1 <- readRDS(file1)
  anno2 <- readRDS(file2)
  anno3 <- readRDS(file3)
  anno4 <- readRDS(file4)
  anno5 <- readRDS(file5)
  anno6 <- readRDS(file6)
 
  cells <- rownames(anno1)
  anno <- rbind(anno, cbind(cells, 
                            anno1[cells,'labels'],
                            anno2[cells,'labels'],
                            anno3[cells,'labels'],
                            anno4[cells,'labels'],
                            anno5[cells,'labels'],
                            anno6[cells,'labels']))
}

colnames(anno) <- c('cells','blue','data','human','imm','monaco','nover')
rownames(anno) <- anno$cells

saveRDS(anno,'singleR/merge_singleR.rds')

dim(anno)
254701      7

dim(sc)
#33538 254701


sc$singleR.blue <- anno[colnames(sc),'blue']
sc$singleR.data <- anno[colnames(sc),'data']
sc$singleR.human <- anno[colnames(sc),'human']
sc$singleR.imm <- anno[colnames(sc),'imm']
sc$singleR.monaco <- anno[colnames(sc),'monaco']
sc$singleR.nover <- anno[colnames(sc),'nover']

saveRDS(sc,'sc.rds')



table(sc$singleR.blue,sc$sample)










