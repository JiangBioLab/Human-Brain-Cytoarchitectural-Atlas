#library(tasic2016data)
library(Seurat)
library(scrattch.hicat)
library(dendextend)
library(dplyr)
library(matrixStats)
library(Matrix)

seu <- readRDS('seu.harmony.anno.de80.rds')

############# Alter Tree wpp #########################################
result <- readRDS('hicat_old/seu.hicat.run_consensus_clust.de80.niter100.rds')
norm.dat0 <- readRDS('norm.dat.rds')

seu <- FindVariableFeatures(seu, nfeatures = 10000)

norm.dat <- norm.dat0[seu@assays$RNA@var.features,]

cl.clean <- result$cl.result$cl
#de.genes <- result$cl.result$de.genes
de.genes <- seu@assays$RNA@var.features

cl.clean[cl.clean %in% c('146','149','152')] <- '146_149_152'
cl.clean[cl.clean %in% c('305','313')] <- '305_313'
cl.clean[cl.clean %in% c('164','165')] <- '164_165'
cl.clean[cl.clean %in% c('131','135')] <- '131_135'
cl.clean[cl.clean %in% c('118','129')] <- '118_129'
cl.clean[cl.clean %in% c('355','357')] <- '355_357'
cl.clean[cl.clean %in% c('47','52','56')] <- '47_52_56'
cl.clean[cl.clean %in% c('311','328')] <- '311_328'
cl.clean[cl.clean %in% c('180','186')] <- '180_186'
cl.clean[cl.clean %in% c('42','43')] <- '42_43'
cl.clean[cl.clean %in% c('360','361')] <- '360_361'
cl.clean[cl.clean %in% c('238','249')] <- '238_249'
cl.clean[cl.clean %in% c('269','271')] <- '269_271'
cl.clean[cl.clean %in% c('60','64')] <- '60_64'

cl.clean.count <- table(cl.clean)
cl.clean.count <- cl.clean.count[order(cl.clean.count)]
cl.clean.count <- rev(cl.clean.count)
cl.newID <- data.frame(old=names(cl.clean.count),
                       new=c(1:length(cl.clean.count)),
                       cells=cl.clean.count)
rownames(cl.newID) <- cl.newID[,1]
write.table(cl.newID,'cl.newID.txt',col.names = F,row.names = F,sep='\t',quote = F)

cl.clean2 <- cl.newID[cl.clean,2]
names(cl.clean2) <- names(cl.clean)

select.markers <- select_markers(norm.dat, 
                                 cl.clean2, 
                                 #de.genes = de.genes,
                                 n.markers = 50)

marker.genes <- select.markers$markers
de.genes <- select.markers$de.genes

#cl.med <- get_cl_medians(norm.dat[marker.genes,], 
#                         cl.clean)
cl.med2 <- get_cl_medians(norm.dat[marker.genes,], 
                          cl.clean2)

##The prefered order for the leaf nodes.
#l.rank <- setNames(1:nrow(consensus.cl.df), 
#                   row.names(consensus.cl.df))

##Color of the leaf nodes.
#l.color <- setNames(as.character(consensus.cl.df$cluster_color), row.names(consensus.cl.df))

# Build the dendrogram
dend.result <- build_dend(cl.med2,
                          #l.rank, 
                          #l.color,
                          nboot = 100)
saveRDS(dend.result, 'seu.hicat.run_consensus_clust.de80.niter100.dend.alter3.rds')


dend.result<-readRDS('seu.hicat.run_consensus_clust.de80.niter100.dend.alter3.rds')
plot(dend.result$dend)

dend <- dend.result$dend
anno <- unique(cbind(as.vector(seu$de80_hicat_cluster_newID),seu$hicat_cluster_subclasses))
rownames(anno) <- anno[,1]
anno[,2] <- paste0(anno[,1],':',anno[,2])

dend.label <- dend
labels(dend.label) <- anno[labels(dend),2]
#dend.pdf 20x5
plot(dend,horiz = T)
plot(dend.label,horiz = T)



a<-table(seu$hicat_cluster_newID,seu$de80_hicat_cluster_newID)
write.table(a,'cmp.txt',col.names = T,row.names = T,sep='\t',quote = F)



dend.result1 <- readRDS('seu.hicat.run_consensus_clust.de100.niter100.dend.alter.rds')
dend.result2 <- readRDS('seu.hicat.run_consensus_clust.de80.niter100.dend.alter.rds')
a <- table(seu$hicat_cluster_newID,seu$de80_hicat_cluster_newID)
a <- a[labels(dend.result1$dend),labels(dend.result2$dend)]
b <- unique(cbind(as.vector(seu$hicat_cluster_newID),seu$hicat_cluster_supertypes))
rownames(b) <- b[,1]
b <- b[labels(dend.result1$dend),]
b <- cbind(b,paste0(b[,1],': ',b[,2]))
rownames(a) <- b[,3]
myColors=brewer.pal(8,"Reds")[1:8]
# de100_vs_de80_tree.pdf 20x25
pheatmap(log2(a+1),
         cluster_rows = F,
         cluster_cols = F,
         color = colorRampPalette(colors = myColors)(100))
# de100_tree.pdf 20x5
plot(dend.result1$dend,horiz = T)
# de80_tree.pdf 5x20
plot(dend.result2$dend)

