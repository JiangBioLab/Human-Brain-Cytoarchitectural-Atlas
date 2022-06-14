## https://theislab.github.io/scanpy-in-R/#setting-up


install.packages("renv")

renv::init()

renv::install("reticulate")

renv::use_python()

pkgs <- c(
  "renv",
  "reticulate",
  "png",
  "ggplot2",
  "BiocManager",
  "Seurat"
)

bioc_pkgs <- c(
  "SingleCellExperiment",
  "scater",
  "multtest"
)

# If you are using an {renv} environment
renv::install(pkgs)

# Otherwise do it the normal way
#install.packages(pkgs)

# Install Bioconductor packages
BiocManager::install(bioc_pkgs, update = FALSE)


py_pkgs <- c(
  "scanpy",
  "python-igraph",
  "louvain"
)

reticulate::py_install(py_pkgs)


renv::snapshot()


suppressPackageStartupMessages({
  library("reticulate")
  library("ggplot2")
  library("SingleCellExperiment")
  library("scater")
  library("Seurat")
})

renv::snapshot()

seurat <- readRDS('/home/wangpp/DATA/Liver/LiverAnalysis/seu.harmony.anno.End.v6.rds')
exprs <- GetAssayData(seurat)
meta <- seurat[[]]
feature_meta <- GetAssay(seurat)[[]]
embedding <- Embeddings(seurat, "umap")

sc <- import("scanpy")

adata_seurat = sc$AnnData(X = t(exprs), obs = meta, var = feature_meta)
adata_seurat$obsm['umap'] = embedding


plt <- import("matplotlib.pyplot")




ax = sc.pl.stacked_violin(pbmc, marker_genes_dict, groupby='clusters', swap_axes=False, dendrogram=True)

sc$pl$stacked_violin(adata_seurat, c('PDPN', 'PROX1', 'VEGFC', 'BMX', 'SEMA3G'),groupby='seurat_clusters')


