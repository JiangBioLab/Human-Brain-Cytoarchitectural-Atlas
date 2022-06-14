#library(tasic2016data)
library(Seurat)
library(scrattch.hicat)
library(dendextend)
library(dplyr)
library(matrixStats)
library(Matrix)

#seu <- readRDS('seu.harmony.anno.rds')
#norm.dat <- seu@assays$RNA@data
#saveRDS(norm.dat,'norm.dat.rds')


#data <- seu@data$RNA@counts
#norm.dat <- logCPM(data)
#saveRDS(norm.dat,'norm.dat.logCPM.rds')
  

norm.dat<-readRDS('norm.dat.rds')

de.param <- de_param(padj.th     = 0.05, 
                     lfc.th      = 1, 
                     low.th      = 1, 
                     q1.th       = 0.3,
                     q2.th       = NULL,
                     q.diff.th   = 0.7, 
                     de.score.th = 150,
                     min.cells = 20)

#gene.counts <- colSums(norm.dat > 0)
#rm.eigen <- matrix(log2(gene.counts), ncol = 1)
#row.names(rm.eigen) <- names(gene.counts)
#colnames(rm.eigen) <- "log2GeneCounts"

#onestep.result <- onestep_clust(norm.dat, 
#                                #select.cells = select.cells, 
#                                dim.method = "pca", 
#                                de.param = de.param)  #strict.param
#                                #rm.eigen = rm.eigen)
#
#saveRDS(onestep.result,'seu.hicat.onestep.result.rds')

#pheatmap(log(table(seu$clusterAnno, onestep.result$cl)+1))
#pheatmap(log(table(seu$individual, onestep.result$cl)+1))

#iter.result <- iter_clust(norm.dat, 
#                          #select.cells = select.cells, 
#                          dim.method = "pca", 
#                          de.param = de.param#, 
#                          #rm.eigen = rm.eigen
#                          )
#saveRDS(iter.result,'seu.hicat.iter.result.rds')



#set.seed(123)

#result <- run_consensus_clust(norm.dat, 
#                              niter = 20, 
#                              #select.cells = select.cells,
#                              de.param = de.param, 
#                              #rm.eigen = rm.eigen, 
#                              dim.method = "pca", 
#                              output_dir = "subsample_PCA_de150_niter20", 
#                              mc.cores = 10)
#saveRDS(result,'seu.hicat.run_consensus_clust.de150.niter20.rds')



############# Tree #########################################
#cl.clean <- result$cl.result$cl
#de.genes <- result$cl.result$de.genes
#
#select.markers <- select_markers(norm.dat, 
#                                 cl.clean, 
#                                 de.genes = de.genes,
#                                 n.markers = 50)
#
#marker.genes <- select.markers$markers
#de.genes <- select.markers$de.genes
#
#cl.med <- get_cl_medians(norm.dat[marker.genes,], 
#                         cl.clean)
#
###The prefered order for the leaf nodes.
##l.rank <- setNames(1:nrow(consensus.cl.df), 
##                   row.names(consensus.cl.df))
#
###Color of the leaf nodes.
##l.color <- setNames(as.character(consensus.cl.df$cluster_color), row.names(consensus.cl.df))
#
## Build the dendrogram
#dend.result <- build_dend(cl.med,
#                          #l.rank, 
#                          #l.color,
##                          nboot = 100)
#saveRDS(dend.result, 'seu.hicat.run_consensus_clust.de150.niter20.dend.rds')
##dend <- dend.result$dend
#####attach cluster labels to the leaves of the tree 
##dend.labeled <- dend
###labels(dend.labeled) <- consensus.cl.df[labels(dend), "cluster_label"]
##plot(dend.labeled)
#
#
#cl <- result$cl.result$cl
#cl.df <-  data.frame(cluster_id = unique(cl),
#                     cluster_label = paste0("cluster_",unique(cl)),
#                     cluster_color = rainbow(length(unique(cl))))
#
#tsne.result <- plot_tsne_cl(norm.dat, marker.genes, cl, cl.df, fn.size=5, cex=1)
#
#saveRDS(tsne.result, 'seu.hicat.run_consensus_clust.de150.niter20.tsne.rds')
##tsne.result$g





set.seed(123)

result <- run_consensus_clust(norm.dat, 
                              niter = 100, 
                              #select.cells = select.cells,
                              de.param = de.param, 
                              #rm.eigen = rm.eigen, 
                              dim.method = "pca", 
                              output_dir = "subsample_PCA_de150_niter100", 
                              mc.cores = 10)
saveRDS(result,'seu.hicat.run_consensus_clust.de150.niter100.rds')



#set.seed(123)
#result <- run_consensus_clust_wpp(norm.dat, 
#                                  niter = 100, 
#                                  #select.cells = select.cells,
#                                  de.param = de.param, 
#                                  #rm.eigen = rm.eigen, 
#                                  dim.method = "pca", 
#                                  output_dir = "subsample_PCA_de150_niter100", 
#                                  mc.cores = 10)
#saveRDS(result,'seu.hicat.run_consensus_clust.de150.niter100.rds')


############# Tree #########################################
cl.clean <- result$cl.result$cl
de.genes <- result$cl.result$de.genes

select.markers <- select_markers(norm.dat, 
                                 cl.clean, 
                                 de.genes = de.genes,
                                 n.markers = 50)

marker.genes <- select.markers$markers
de.genes <- select.markers$de.genes

cl.med <- get_cl_medians(norm.dat[marker.genes,], 
                         cl.clean)

##The prefered order for the leaf nodes.
#l.rank <- setNames(1:nrow(consensus.cl.df), 
#                   row.names(consensus.cl.df))

##Color of the leaf nodes.
#l.color <- setNames(as.character(consensus.cl.df$cluster_color), row.names(consensus.cl.df))

# Build the dendrogram
dend.result <- build_dend(cl.med,
                          #l.rank, 
                          #l.color,
                          nboot = 100)
saveRDS(dend.result, 'seu.hicat.run_consensus_clust.niter100.dend.rds')
#dend <- dend.result$dend
####attach cluster labels to the leaves of the tree 
#dend.labeled <- dend
##labels(dend.labeled) <- consensus.cl.df[labels(dend), "cluster_label"]
#plot(dend.labeled)


cl <- result$cl.result$cl
cl.df <-  data.frame(cluster_id = unique(cl),
                     cluster_label = paste0("cluster_",unique(cl)),
                     cluster_color = rainbow(length(unique(cl))))

tsne.result <- plot_tsne_cl(norm.dat, marker.genes, cl, cl.df, fn.size=5, cex=1)

saveRDS(tsne.result, 'seu.hicat.run_consensus_clust.niter100.tsne.rds')
#tsne.result$g







############# Alter Tree #########################################
result <- readRDS('seu.hicat.run_consensus_clust.de100.niter100.rds')
norm.dat <- readRDS('norm.dat.rds')

cl.clean <- result$cl.result$cl
de.genes <- result$cl.result$de.genes

cl.clean[cl.clean %in% c('261','266')] <- '261_266'
cl.clean[cl.clean %in% c('191','192')] <- '191_192'
cl.clean[cl.clean %in% c('117','119')] <- '117_119'
cl.clean[cl.clean %in% c('57','58')] <- '57_58'
cl.clean[cl.clean %in% c('44','45')] <- '44_45'
cl.clean[cl.clean %in% c('59','61')] <- '59_61'
cl.clean[cl.clean %in% c('25','30','32')] <- '25_30_32'
cl.clean[cl.clean %in% c('38','41','42')] <- '38_41_42'

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
                                 de.genes = de.genes,
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
saveRDS(dend.result, 'seu.hicat.run_consensus_clust.de100.niter100.dend.alter.rds')
plot(dend.result$dend)
#dend <- dend.result$dend
####attach cluster labels to the leaves of the tree 
#dend.labeled <- dend
##labels(dend.labeled) <- consensus.cl.df[labels(dend), "cluster_label"]
#plot(dend.labeled)












load('subsample_PCA_niter100/cells.16.rda')
result <- scrattch.hicat::iter_clust(norm.dat=norm.dat, dim.method = "pca", select.cells=select.cells, de.param = de.param)
save(result, file='subsample_PCA_niter100/result.16.rda')




