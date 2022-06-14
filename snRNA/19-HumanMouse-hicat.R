#library(tasic2016data)
library(Seurat)
library(scrattch.hicat)
library(dendextend)
library(dplyr)
library(matrixStats)
library(Matrix)
source('consensusCluster.R')



####### Exc ########

seu <- readRDS('EXC_h_m_SCT.rds')
norm.dat <- seu@assays$SCT@data

de.param <- de_param(padj.th     = 0.05, 
                     lfc.th      = 1, 
                     low.th      = 1, 
                     q1.th       = 0.3,
                     q2.th       = NULL,
                     q.diff.th   = 0.7, 
                     de.score.th = 80,
                     min.cells = 20)

set.seed(12345)
#result <- run_consensus_clust_wpp(norm.dat, 
#                                  niter = 100, 
#                                  #select.cells = select.cells,
#                                  de.param = de.param, 
#                                  #rm.eigen = rm.eigen, 
#                                  dim.method = "pca", 
#                                  output_dir = "EXC_h_m_SCT.hicat.de80", 
#                                  mc.cores = 1)

#for(i in 1:10){
#  infile <- paste0('EXC_h_m_SCT.hicat.de80/cells.',i,'.rda')
#  outfile <- paste0('EXC_h_m_SCT.hicat.de80/result.',i,'.rda')
#  load(infile)
#  result <- scrattch.hicat::iter_clust(norm.dat=norm.dat, dim.method = "pca", select.cells=select.cells, de.param = de.param)
#  save(result, file=outfile)
#}

result <- run_consensus_clust_wpp(norm.dat, 
                                  niter = 100, 
                                  #select.cells = select.cells,
                                  de.param = de.param, 
                                  #rm.eigen = rm.eigen, 
                                  dim.method = "pca", 
                                  output_dir = "EXC_h_m_SCT.hicat.de80", 
                                  mc.cores = 1)

saveRDS(result,'EXC_h_m_SCT.hicat.de80.rds')



#### dend
result <- readRDS('EXC_h_m_SCT.hicat.de80.rds')
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

# Build the dendrogram
dend.result <- build_dend(cl.med,
                          #l.rank, 
                          #l.color,
                          nboot = 100)

# EXC_h_m_SCT_dend.pdf 4x20
plot(dend.result$dend)
saveRDS(dend.result, 'EXC_h_m_SCT.hicat.de80.dend.rds')


result <- readRDS('EXC_h_m_SCT.hicat.de80.rds')
dend.result <- readRDS('EXC_h_m_SCT.hicat.de80.dend.rds')

cl <- result$cl.result$cl
cells <- names(cl)
M <- cells[grep('F',cells)]
H <- cells[!(cells %in% M)]

cc <- data.frame(cl = cl,
                 sample = cl)
rownames(cc) <- names(cl)
cc[M,'sample'] <- 'M'
cc[H,'sample'] <- 'H'

dd <- table(cc$cl,cc$sample)
ll <- dd[,1] > dd[,2]
ll2 <- ll
ll2[ll] <- paste0(names(ll2)[ll],'(H)')
ll2[!ll] <- paste0(names(ll2)[!ll],'(M)')


labels(dend.result$dend) <- ll2[labels(dend.result$dend)]
# EXC_h_m_SCT_dend.pdf 20x4
plot(dend.result$dend,
     horiz = T)





####### Inh ########
seu <- readRDS('INH_h_m_SCT.rds')
norm.dat <- seu@assays$SCT@data

de.param <- de_param(padj.th     = 0.05, 
                     lfc.th      = 1, 
                     low.th      = 1, 
                     q1.th       = 0.3,
                     q2.th       = NULL,
                     q.diff.th   = 0.7, 
                     de.score.th = 80,
                     min.cells = 20)

set.seed(12345)

result <- run_consensus_clust_wpp(norm.dat,
                              niter = 100, 
                              #select.cells = select.cells,
                              de.param = de.param, 
                              #rm.eigen = rm.eigen, 
                              dim.method = "pca", 
                              output_dir = "INH_h_m_SCT.hicat.de80", 
                              mc.cores = 1)
saveRDS(result,'INH_h_m_SCT.hicat.de80.rds')


#### dend
result <- readRDS('INH_h_m_SCT.hicat.de80.rds')
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

# Build the dendrogram
dend.result <- build_dend(cl.med,
                          #l.rank, 
                          #l.color,
                          nboot = 100)

# INH_h_m_SCT_dend.pdf 4x20
plot(dend.result$dend)

saveRDS(dend.result, 'INH_h_m_SCT.hicat.de80.dend.rds')




result <- readRDS('INH_h_m_SCT.hicat.de80.rds')
dend.result <- readRDS('INH_h_m_SCT.hicat.de80.dend.rds')

cl <- result$cl.result$cl
cells <- names(cl)
M <- cells[grep('F',cells)]
H <- cells[!(cells %in% M)]

cc <- data.frame(cl = cl,
                 sample = cl)
rownames(cc) <- names(cl)
cc[M,'sample'] <- 'M'
cc[H,'sample'] <- 'H'

dd <- table(cc$cl,cc$sample)
ll <- dd[,1] > dd[,2]
ll2 <- ll
ll2[ll] <- paste0(names(ll2)[ll],'(H)')
ll2[!ll] <- paste0(names(ll2)[!ll],'(M)')


labels(dend.result$dend) <- ll2[labels(dend.result$dend)]
# INH_h_m_SCT_dend.pdf 20x4
plot(dend.result$dend,
     horiz = T)



