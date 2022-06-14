library(Seurat)
library(ggplot2)
library(harmony)


sc <- readRDS('sc.qc2.rds')
dim(sc)
#36601 252296

seu <- CreateSeuratObject(
  counts(sc)[,sc$individual %in% c("S0206", "S0406", "S0426", "S0531", "S0916")],
  min.cells = 10,
  min.features = 200
)

# min.cells = 10, min.features =200
# 33407 252294

meta <- data.frame(cells = colnames(sc),
                   Sample = sc$sample,
                   Individual = sc$individual,
                   Region = sc$region,
                   singleR.blue = sc$singleR.blue,
                   singleR.data = sc$singleR.data,
                   singleR.human = sc$singleR.human,
                   singleR.imm = sc$singleR.imm,
                   singleR.monaco = sc$singleR.monaco,
                   singleR.nover = sc$singleR.nover)
rownames(meta) <- meta$cells

seu$sample <- meta[colnames(seu),'Sample']
seu$individual <- meta[colnames(seu),'Individual']
seu$region <- meta[colnames(seu),'Region']
seu$singleR.blue <- meta[colnames(seu),'singleR.blue']
seu$singleR.data <- meta[colnames(seu),'singleR.data']
seu$singleR.human <- meta[colnames(seu),'singleR.human']
seu$singleR.imm <- meta[colnames(seu),'singleR.imm']
seu$singleR.monaco <- meta[colnames(seu),'singleR.monaco']
seu$singleR.nover <- meta[colnames(seu),'singleR.nover']


seu$sample <- factor(seu$sample, levels = c(#'S0128_A2','S0128_A6','S0128_A7','S0128_A8','S0128_A9','S0128_A10','S0128_A11','S0128_A12','S0128_A13','S0128_A14',
                                            'S0206_A1','S0206_A2','S0206_A3','S0206_A4','S0206_A5','S0206_A6','S0206_A7','S0206_A11','S0206_A12',
                                            'S0406_A1','S0406_A3','S0406_A4','S0406_A5','S0406_A6', 'S0406_A7', 'S0406_A8', 'S0406_A9', 'S0406_A10','S0406_A14',
                                            'S0426_A5', 'S0426_A8', 'S0426_A9', 'S0426_A10','S0426_A11','S0426_A12','S0426_A13','S0426_A14',
                                            'S0531_A1', 'S0531_A2', 'S0531_A3','S0531_A4', 'S0531_A13',
                                            'S0916-A2','S0916-A6','S0916-A7','S0916-A8','S0916-A9','S0916-A10','S0916-A11','S0916-A13','S0916-A12','S0916-A14'))
seu$individual <- factor(seu$individual, levels = c(#"S0128", 
                                                    "S0206", "S0406", "S0426","S0531","S0916"))
seu$region <- factor(seu$region, levels = c("A1",  "A2",  "A3",  "A4",  "A5",  "A6", "A7", "A8", "A9", "A10", "A11",  "A12",  "A13",  "A14"))


#saveRDS(seu,'seu.rds')


seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

seu <- NormalizeData(object = seu, normalization.method = "LogNormalize", scale.factor = 1e4)

seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(seu)
seu <- ScaleData(seu, features = all.genes)

seu <- RunPCA(seu, features = VariableFeatures(object = seu))

ElbowPlot(seu)

#seu <- FindNeighbors(seu, dims = 1:20)
#seu <- FindClusters(seu, resolution = 1)

#seu <- RunUMAP(seu, dims = 1:20)
#DimPlot(seu, reduction = "umap")
#DimPlot(seu, reduction = "umap", group.by = 'tissue')

#saveRDS(seu,'seu.pre.rds')


seu <- RunHarmony(seu, "sample")
#harmony_embeddings <- Embeddings(seu, 'harmony')

seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:20)
seu <- FindClusters(seu, resolution = 1)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:20)

saveRDS(seu,'seu.harmony.rds')





