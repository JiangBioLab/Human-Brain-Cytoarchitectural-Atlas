library(Seurat)
library(ggplot2)
library(harmony)

sc <- readRDS('sc.qc2.rds')
dim(sc)
#36601 252296

#seu <- CreateSeuratObject(
#  counts(sc),
#  min.cells = 10,
#  min.features = 200
#)
## min.cells = 10, min.features =200


seu <- CreateSeuratObject(
  counts(sc)[,sc$individual %in% c("S0206", "S0406", "S0426", "S0531", "S0916")],
  min.cells = 10,
  min.features = 200
)
# min.cells = 10, min.features =200
# 33261 234645

meta <- data.frame(cells = colnames(sc),
                   sample = sc$sample,
                   individual = sc$individual,
                   region = sc$region,
                   singleR.blue = sc$singleR.blue,
                   singleR.data = sc$singleR.data,
                   singleR.human = sc$singleR.human,
                   singleR.imm = sc$singleR.imm,
                   singleR.monaco = sc$singleR.monaco,
                   singleR.nover = sc$singleR.nover)
rownames(meta) <- meta$cells

seu$sample <- meta[colnames(seu),'sample']
seu$individual <- meta[colnames(seu),'individual']
seu$region <- meta[colnames(seu),'region']
seu$singleR.blue <- meta[colnames(seu),'singleR.blue']
seu$singleR.data <- meta[colnames(seu),'singleR.data']
seu$singleR.human <- meta[colnames(seu),'singleR.human']
seu$singleR.imm <- meta[colnames(seu),'singleR.imm']
seu$singleR.monaco <- meta[colnames(seu),'singleR.monaco']
seu$singleR.nover <- meta[colnames(seu),'singleR.nover']


seu$sample <- factor(seu$sample, levels = c(#'S0128-A2','S0128-A6','S0128-A7','S0128-A8','S0128-A9','S0128-A10','S0128-A11','S0128-A12','S0128-A13','S0128-A14',
                                            'S0206-A1','S0206-A2','S0206-A3','S0206-A4','S0206-A5','S0206-A6','S0206-A7','S0206-A11','S0206-A12',
                                            'S0406-A1','S0406-A3','S0406-A4','S0406-A5','S0406-A6', 'S0406-A7', 'S0406-A8', 'S0406-A9', 'S0406-A10','S0406-A14',
                                            'S0426-A5', 'S0426-A8', 'S0426-A9', 'S0426-A10','S0426-A11','S0426-A12','S0426-A13','S0426-A14',
                                            'S0531-A1', 'S0531-A2', 'S0531-A3','S0531-A4', 'S0531-A13',
                                            'S0916-A2','S0916-A6','S0916-A7','S0916-A8','S0916-A9','S0916-A10','S0916-A11','S0916-A13','S0916-A12','S0916-A14'))
seu$individual <- factor(seu$individual, levels = c(#"S0128", 
                                                    "S0206", "S0406", "S0426","S0531", "S0916"))
seu$region <- factor(seu$region, levels = c("A1",  "A2",  "A3",  "A4",  "A5",  "A6", "A7", "A8", "A9", "A10", "A11",  "A12",  "A13",  "A14"))



seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

#data <- data.frame(sample = seu$sample,
#                   genes = seu$nFeature_RNA,
#                   counts = seu$nCount_RNA,
#                   percent.mt = seu$percent.mt)
#
#ggplot(data, aes(x = sample, y = genes)) +
#  geom_violin(aes(fill = sample)) +
#  ylim(0, 3000) + 
#  theme_classic() + RotatedAxis()
#
#ggplot(data, aes(x = sample, y = counts)) +
#  geom_violin(aes(fill = sample)) +
#  ylim(0, 10000) + 
#  theme_classic() + RotatedAxis()
#
#ggplot(data, aes(x = sample, y = percent.mt)) +
#  geom_violin(aes(fill = sample)) +
#  ylim(0, 20) + 
#  theme_classic() + RotatedAxis()


#plot1 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = 'sample')
#plot2 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = 'sample')
#plot1 + plot2


seu <- NormalizeData(object = seu, normalization.method = "LogNormalize", scale.factor = 1e4)

#norm.dat <- seu@assays$RNA@data
#saveRDS(norm.dat,'norm.dat.rds')

seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

seu <- CellCycleScoring(seu, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
seu$CC.Difference <- seu$S.Score - seu$G2M.Score

seu <- ScaleData(seu, vars.to.regress = "CC.Difference", features = rownames(seu), verbose = FALSE)
#seu <- ScaleData(seu, features = rownames(seu), verbose = FALSE)

saveRDS(seu,'seu.scaleMT.harmony.rds')

seu <- RunPCA(seu, features = VariableFeatures(object = seu))

ElbowPlot(seu)

seu <- RunHarmony(seu, "sample")
#harmony_embeddings <- Embeddings(seu, 'harmony')

seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:50) #1:20
seu <- FindClusters(seu, resolution = 2.0)  #1.0
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:50)
seu <- RunTSNE(seu, reduction = "harmony", dims = 1:50)

seu <- RunUMAP(seu, reduction = "harmony", dims = 1:40)

#DimPlot(seu)
#DimPlot(seu, reduction = 'tsne')


saveRDS(seu,'seu.scaleMT.harmony.rds')
#saveRDS(seu,'seu.harmony.rds')





