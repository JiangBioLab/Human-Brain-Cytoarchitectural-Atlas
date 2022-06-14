library(DoubletFinder)
library(Seurat)


sample <- c(#'S0128-A2','S0128-A6','S0128-A7','S0128-A8','S0128-A9','S0128-A10','S0128-A11','S0128-A13','S0128-A12','S0128-A14',
            #'S0206-A1','S0206-A2','S0206-A3','S0206-A4','S0206-A5','S0206-A6','S0206-A7','S0206-A11','S0206-A12',
            #'S0406-A1','S0406-A3','S0406-A4','S0406-A5','S0406-A6', 'S0406-A7', 'S0406-A8', 'S0406-A9', 'S0406-A10','S0406-A14',
            #'S0426-A5', 'S0426-A8', 'S0426-A9', 'S0426-A10','S0426-A11','S0426-A12','S0426-A13','S0426-A14',
            #'S0531-A1', 'S0531-A2', 'S0531-A3','S0531-A4', 'S0531-A13',
	    'S0916-A2','S0916-A6','S0916-A7','S0916-A8','S0916-A9','S0916-A10','S0916-A11','S0916-A13','S0916-A12','S0916-A14')

for(i in 1:length(sample)){
  file <- paste0('counts/',sample[i],'.counts.rds')
  count <- readRDS(file)
  seu <- CreateSeuratObject(counts = count,  min.cells = 0, min.features = 0)
  
  seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(seu)
  seu <- ScaleData(seu, features = all.genes)
  seu <- RunPCA(seu, features = VariableFeatures(object = seu))
  seu <- FindNeighbors(seu, dims = 1:15)
  seu <- FindClusters(seu, resolution = 1)
  
  fileSeu <- paste0(sample[i],'_seu.rds')
  saveRDS(seu,fileSeu)
  
  ################### Doublet ###########################
  sweep.res.list <- paramSweep_v3(seu, PCs = 1:15, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  mpK<-as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  
  # 找最佳 nExp
  annotations <- seu$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(0.075*ncol(seu@assays$RNA@data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  # 找 Doublet
  seurat_filterDouble <- doubletFinder_v3(seu, PCs = 1:15, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = F)
  fileout <- paste0(sample[i],'_doublet.rds')
  saveRDS(seurat_filterDouble,fileout)
}

#remove.packages("Matrix")
#packageurl <- "https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.3-2.tar.gz"
#install.packages(packageurl, repos=NULL, type="source")



S0128_A2  <- readRDS('doublet/S0128-A2_doublet.rds')
S0128_A6  <- readRDS('doublet/S0128-A6_doublet.rds')
S0128_A7  <- readRDS('doublet/S0128-A7_doublet.rds')
S0128_A8  <- readRDS('doublet/S0128-A8_doublet.rds')
S0128_A9  <- readRDS('doublet/S0128-A9_doublet.rds')
S0128_A10 <- readRDS('doublet/S0128-A10_doublet.rds')
S0128_A11 <- readRDS('doublet/S0128-A11_doublet.rds')
S0128_A12 <- readRDS('doublet/S0128-A12_doublet.rds')
S0128_A13 <- readRDS('doublet/S0128-A13_doublet.rds')
S0128_A14 <- readRDS('doublet/S0128-A14_doublet.rds')
S0206_A1  <- readRDS('doublet/S0206-A1_doublet.rds')
S0206_A2  <- readRDS('doublet/S0206-A2_doublet.rds')
S0206_A3  <- readRDS('doublet/S0206-A3_doublet.rds')
S0206_A4  <- readRDS('doublet/S0206-A4_doublet.rds')
S0206_A5  <- readRDS('doublet/S0206-A5_doublet.rds')
S0206_A6  <- readRDS('doublet/S0206-A6_doublet.rds')
S0206_A7  <- readRDS('doublet/S0206-A7_doublet.rds')
S0206_A11 <- readRDS('doublet/S0206-A11_doublet.rds')
S0206_A12 <- readRDS('doublet/S0206-A12_doublet.rds')
S0406_A1  <- readRDS('doublet/S0406-A1_doublet.rds')
S0406_A3  <- readRDS('doublet/S0406-A3_doublet.rds')
S0406_A4  <- readRDS('doublet/S0406-A4_doublet.rds')
S0406_A5  <- readRDS('doublet/S0406-A5_doublet.rds')
S0406_A6  <- readRDS('doublet/S0406-A6_doublet.rds')
S0406_A7  <- readRDS('doublet/S0406-A7_doublet.rds')
S0406_A8  <- readRDS('doublet/S0406-A8_doublet.rds')
S0406_A9  <- readRDS('doublet/S0406-A9_doublet.rds')
S0406_A10 <- readRDS('doublet/S0406-A10_doublet.rds')
S0406_A14 <- readRDS('doublet/S0406-A14_doublet.rds')
S0426_A5  <- readRDS('doublet/S0426-A5_doublet.rds')
S0426_A8  <- readRDS('doublet/S0426-A8_doublet.rds')
S0426_A9  <- readRDS('doublet/S0426-A9_doublet.rds')
S0426_A10 <- readRDS('doublet/S0426-A10_doublet.rds')
S0426_A11 <- readRDS('doublet/S0426-A11_doublet.rds')
S0426_A12 <- readRDS('doublet/S0426-A12_doublet.rds')
S0426_A13 <- readRDS('doublet/S0426-A13_doublet.rds')
S0426_A14 <- readRDS('doublet/S0426-A14_doublet.rds')
S0531_A1  <- readRDS('doublet/S0531-A1_doublet.rds')
S0531_A2  <- readRDS('doublet/S0531-A2_doublet.rds')
S0531_A3  <- readRDS('doublet/S0531-A3_doublet.rds')
S0531_A4  <- readRDS('doublet/S0531-A4_doublet.rds')
S0531_A13 <- readRDS('doublet/S0531-A13_doublet.rds')
S0916_A2  <- readRDS('doublet/S0916-A2_doublet.rds')
S0916_A6  <- readRDS('doublet/S0916-A6_doublet.rds')
S0916_A7  <- readRDS('doublet/S0916-A7_doublet.rds')
S0916_A8  <- readRDS('doublet/S0916-A8_doublet.rds')
S0916_A9  <- readRDS('doublet/S0916-A9_doublet.rds')
S0916_A10 <- readRDS('doublet/S0916-A10_doublet.rds')
S0916_A11 <- readRDS('doublet/S0916-A11_doublet.rds')
S0916_A12 <- readRDS('doublet/S0916-A12_doublet.rds')
S0916_A13 <- readRDS('doublet/S0916-A13_doublet.rds')
S0916_A14 <- readRDS('doublet/S0916-A14_doublet.rds')
#
doublet <- c()
doublet <- rbind(cbind(names(S0128_A2$orig.ident), S0128_A2$DF.classifications_0.25_0.03_190),
                 cbind(names(S0128_A6$orig.ident), S0128_A6$DF.classifications_0.25_0.25_266),
                 cbind(names(S0128_A7$orig.ident), S0128_A7$DF.classifications_0.25_0.005_282),
                 cbind(names(S0128_A8$orig.ident), S0128_A8$DF.classifications_0.25_0.01_205),
                 cbind(names(S0128_A9$orig.ident), S0128_A9$DF.classifications_0.25_0.02_98),
                 cbind(names(S0128_A10$orig.ident),S0128_A10$DF.classifications_0.25_0.02_125),
                 cbind(names(S0128_A11$orig.ident),S0128_A11$DF.classifications_0.25_0.005_174),
                 cbind(names(S0128_A12$orig.ident),S0128_A12$DF.classifications_0.25_0.14_94),
                 cbind(names(S0128_A13$orig.ident),S0128_A13$DF.classifications_0.25_0.28_146),
                 cbind(names(S0128_A14$orig.ident),S0128_A14$DF.classifications_0.25_0.02_175),
                 cbind(names(S0206_A1$orig.ident), S0206_A1$DF.classifications_0.25_0.02_410),
                 cbind(names(S0206_A2$orig.ident), S0206_A2$DF.classifications_0.25_0.29_398),
                 cbind(names(S0206_A3$orig.ident), S0206_A3$DF.classifications_0.25_0.07_433),
                 cbind(names(S0206_A4$orig.ident), S0206_A4$DF.classifications_0.25_0.02_476),
                 cbind(names(S0206_A5$orig.ident), S0206_A5$DF.classifications_0.25_0.27_389),
                 cbind(names(S0206_A6$orig.ident), S0206_A6$DF.classifications_0.25_0.01_263),
                 cbind(names(S0206_A7$orig.ident), S0206_A7$DF.classifications_0.25_0.23_508),
                 cbind(names(S0206_A11$orig.ident),S0206_A11$DF.classifications_0.25_0.22_595),
                 cbind(names(S0206_A12$orig.ident),S0206_A12$DF.classifications_0.25_0.01_523),
                 cbind(names(S0406_A1$orig.ident), S0406_A1$DF.classifications_0.25_0.09_843),
                 cbind(names(S0406_A3$orig.ident), S0406_A3$DF.classifications_0.25_0.005_677),
                 cbind(names(S0406_A4$orig.ident), S0406_A4$DF.classifications_0.25_0.09_577),
                 cbind(names(S0406_A5$orig.ident), S0406_A5$DF.classifications_0.25_0.005_562),
                 cbind(names(S0406_A6$orig.ident), S0406_A6$DF.classifications_0.25_0.27_801),
                 cbind(names(S0406_A7$orig.ident), S0406_A7$DF.classifications_0.25_0.3_663),
                 cbind(names(S0406_A8$orig.ident), S0406_A8$DF.classifications_0.25_0.3_1104),
                 cbind(names(S0406_A9$orig.ident), S0406_A9$DF.classifications_0.25_0.21_801),
                 cbind(names(S0406_A10$orig.ident),S0406_A10$DF.classifications_0.25_0.3_733),
                 cbind(names(S0406_A14$orig.ident),S0406_A14$DF.classifications_0.25_0.15_705),
                 cbind(names(S0426_A5$orig.ident), S0426_A5$DF.classifications_0.25_0.05_427),
                 cbind(names(S0426_A8$orig.ident), S0426_A8$DF.classifications_0.25_0.3_523),
                 cbind(names(S0426_A9$orig.ident), S0426_A9$DF.classifications_0.25_0.25_430),
                 cbind(names(S0426_A10$orig.ident),S0426_A10$DF.classifications_0.25_0.17_443),
                 cbind(names(S0426_A11$orig.ident),S0426_A11$DF.classifications_0.25_0.03_420),
                 cbind(names(S0426_A12$orig.ident),S0426_A12$DF.classifications_0.25_0.03_586),
                 cbind(names(S0426_A13$orig.ident),S0426_A13$DF.classifications_0.25_0.3_496),
                 cbind(names(S0426_A14$orig.ident),S0426_A14$DF.classifications_0.25_0.18_437),
                 cbind(names(S0531_A1$orig.ident), S0531_A1$DF.classifications_0.25_0.15_606),
                 cbind(names(S0531_A2$orig.ident), S0531_A2$DF.classifications_0.25_0.13_655),
                 cbind(names(S0531_A3$orig.ident), S0531_A3$DF.classifications_0.25_0.12_594),
                 cbind(names(S0531_A4$orig.ident), S0531_A4$DF.classifications_0.25_0.06_610),
                 cbind(names(S0531_A13$orig.ident),S0531_A13$DF.classifications_0.25_0.17_516),
                 cbind(names(S0916_A2$orig.ident), S0916_A2$DF.classifications_0.25_0.06_499),
                 cbind(names(S0916_A6$orig.ident), S0916_A6$DF.classifications_0.25_0.3_575),
                 cbind(names(S0916_A7$orig.ident), S0916_A7$DF.classifications_0.25_0.23_430),
                 cbind(names(S0916_A8$orig.ident), S0916_A8$DF.classifications_0.25_0.04_475),
                 cbind(names(S0916_A9$orig.ident), S0916_A9$DF.classifications_0.25_0.16_499),
                 cbind(names(S0916_A10$orig.ident),S0916_A10$DF.classifications_0.25_0.2_506),
                 cbind(names(S0916_A11$orig.ident),S0916_A11$DF.classifications_0.25_0.3_456),
                 cbind(names(S0916_A12$orig.ident),S0916_A12$DF.classifications_0.25_0.005_472),
                 cbind(names(S0916_A13$orig.ident),S0916_A13$DF.classifications_0.25_0.12_529),
                 cbind(names(S0916_A14$orig.ident),S0916_A14$DF.classifications_0.25_0.09_437))


colnames(doublet) <- c('cells','doublet')
rownames(doublet) <- doublet[,1]
saveRDS(doublet,'doublet/merge_doublet.rds')

#sc$doublet <- doublet[colnames(sc),2]

#saveRDS(sc,'sc.rds')#
