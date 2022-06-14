library(Seurat)
library(harmony)

############# GABA CGE###############
counts0 <- readRDS('counts.v2.rds')
meta <- readRDS('seu.harmony.anno.meta.v2.rds')

clusters <- c("118", "45", "126", "40", "130", "90", "120", "122", "73", "75", "101", "87", "94", "99", 
              "20", "44", "49", "51", "92", "132", "29", "78", "52", "68", "15", "46", "121", "77", "19",
              "83", "109", "32", "33", "54", "86")
cells <- rownames(meta)[meta$hicat_cluster_merge_newID %in% clusters]

counts <- counts0[,cells]

seu <- CreateSeuratObject(
  counts,
  min.cells = 0,
  min.features = 0
)

# 33014 36013

seu$sample <- meta[colnames(seu),'sample']
seu$individual <- meta[colnames(seu),'individual']
seu$region <- meta[colnames(seu),'region']
seu$hicat_cluster_merge <- meta[colnames(seu),'hicat_cluster_merge']
seu$hicat_cluster_merge_newID <- meta[colnames(seu),'hicat_cluster_merge_newID']
seu$hicat_cluster_classes <- meta[colnames(seu),'hicat_cluster_classes']
seu$hicat_cluster_subclasses <- meta[colnames(seu),'hicat_cluster_subclasses']
seu$hicat_cluster_supertypes <- meta[colnames(seu),'hicat_cluster_supertypes']

seu$hicat_cluster_merge_newID <- factor(seu$hicat_cluster_merge_newID, levels = clusters)

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

seu <- NormalizeData(object = seu, normalization.method = "LogNormalize", scale.factor = 1e4)

seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(seu)
seu <- ScaleData(seu, features = all.genes)

seu <- RunPCA(seu, features = VariableFeatures(object = seu))

seu <- RunHarmony(seu, "sample")

ElbowPlot(seu)

seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:20)
seu <- FindClusters(seu, resolution = 1)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:20)

saveRDS(seu,'seu.harmony.GABA.CGE.rds')




############## GABA MEG ###############
counts0 <- readRDS('counts.v2.rds')
meta <- readRDS('seu.harmony.anno.meta.v2.rds')

clusters <- c("8", "113", "88", "50", "36", "105", "64", "23", "141", "135", "37", "28", "21", 
              "58", "112", "55", "110", "24", "25", "102", "97", "100", "67")
cells <- rownames(meta)[meta$hicat_cluster_merge_newID %in% clusters]

counts <- counts0[,cells]

seu <- CreateSeuratObject(
  counts,
  min.cells = 0,
  min.features = 0
)

# 33014 36013

seu$sample <- meta[colnames(seu),'sample']
seu$individual <- meta[colnames(seu),'individual']
seu$region <- meta[colnames(seu),'region']
seu$hicat_cluster_merge <- meta[colnames(seu),'hicat_cluster_merge']
seu$hicat_cluster_merge_newID <- meta[colnames(seu),'hicat_cluster_merge_newID']
seu$hicat_cluster_classes <- meta[colnames(seu),'hicat_cluster_classes']
seu$hicat_cluster_subclasses <- meta[colnames(seu),'hicat_cluster_subclasses']
seu$hicat_cluster_supertypes <- meta[colnames(seu),'hicat_cluster_supertypes']

seu$hicat_cluster_merge_newID <- factor(seu$hicat_cluster_merge_newID, levels = clusters)

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

seu <- NormalizeData(object = seu, normalization.method = "LogNormalize", scale.factor = 1e4)

seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(seu)
seu <- ScaleData(seu, features = all.genes)

seu <- RunPCA(seu, features = VariableFeatures(object = seu))

seu <- RunHarmony(seu, "sample")

ElbowPlot(seu)

seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:20)
seu <- FindClusters(seu, resolution = 1)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:20)

saveRDS(seu,'seu.harmony.GABA.MGE.rds')



############## Glu ###############

counts0 <- readRDS('counts.v2.rds')
meta <- readRDS('seu.harmony.anno.meta.v2.rds')

clusters <- c("14", "12", "62", "65", "39", "48", "41", "13", "108", "30", "18", "34", "26", 
              "17", "7", "9", "3", "38", "5", "10", "11", "16", "71", "43", 
              "47", "22", "74", "106", "61", "137", "115", "124", "53", "81", "63", "114", 
              "72", "142", "107", "80", "98", "133", "131", "91", "93", "138", "76", "134",
              "85", "123", "139", "128", "136", "60", "104", "116", "140", "119", "111")
cells <- rownames(meta)[meta$hicat_cluster_merge_newID %in% clusters]

counts <- counts0[,cells]

seu <- CreateSeuratObject(
  counts,
  min.cells = 0,
  min.features = 0
)

# 33014 36013

seu$sample <- meta[colnames(seu),'sample']
seu$individual <- meta[colnames(seu),'individual']
seu$region <- meta[colnames(seu),'region']
seu$hicat_cluster_merge <- meta[colnames(seu),'hicat_cluster_merge']
seu$hicat_cluster_merge_newID <- meta[colnames(seu),'hicat_cluster_merge_newID']
seu$hicat_cluster_classes <- meta[colnames(seu),'hicat_cluster_classes']
seu$hicat_cluster_subclasses <- meta[colnames(seu),'hicat_cluster_subclasses']
seu$hicat_cluster_supertypes <- meta[colnames(seu),'hicat_cluster_supertypes']

seu$hicat_cluster_merge_newID <- factor(seu$hicat_cluster_merge_newID, levels = clusters)

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

seu <- NormalizeData(object = seu, normalization.method = "LogNormalize", scale.factor = 1e4)

seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(seu)
seu <- ScaleData(seu, features = all.genes)

seu <- RunPCA(seu, features = VariableFeatures(object = seu))

seu <- RunHarmony(seu, "sample")

ElbowPlot(seu)

seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:20)
seu <- FindClusters(seu, resolution = 1)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:20)

saveRDS(seu,'seu.harmony.Glu.rds')


############### non-neuron ###############
counts0 <- readRDS('counts.v2.rds')
meta <- readRDS('seu.harmony.anno.meta.v2.rds')

clusters <- c("4", "66", "89", "57", "95", "56", "1", "127", "42", "2", "27", 
              "69", "70", "103", "31", "96", "59", "125", "35", "79", "117", 
              "129", "6", "82", "84")
cells <- rownames(meta)[meta$hicat_cluster_merge_newID %in% clusters]

counts <- counts0[,cells]

seu <- CreateSeuratObject(
  counts,
  min.cells = 0,
  min.features = 0
)

# 33014 36013

seu$sample <- meta[colnames(seu),'sample']
seu$individual <- meta[colnames(seu),'individual']
seu$region <- meta[colnames(seu),'region']
seu$hicat_cluster_merge <- meta[colnames(seu),'hicat_cluster_merge']
seu$hicat_cluster_merge_newID <- meta[colnames(seu),'hicat_cluster_merge_newID']
seu$hicat_cluster_classes <- meta[colnames(seu),'hicat_cluster_classes']
seu$hicat_cluster_subclasses <- meta[colnames(seu),'hicat_cluster_subclasses']
seu$hicat_cluster_supertypes <- meta[colnames(seu),'hicat_cluster_supertypes']

seu$hicat_cluster_merge_newID <- factor(seu$hicat_cluster_merge_newID, levels = clusters)

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

seu <- NormalizeData(object = seu, normalization.method = "LogNormalize", scale.factor = 1e4)

seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(seu)
seu <- ScaleData(seu, features = all.genes)

seu <- RunPCA(seu, features = VariableFeatures(object = seu))

seu <- RunHarmony(seu, "sample")

ElbowPlot(seu)

seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:20)
seu <- FindClusters(seu, resolution = 1)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:20)

saveRDS(seu,'seu.harmony.nonNeuron.rds')


