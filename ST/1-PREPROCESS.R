setwd("/home/luomeng/DATA/BRAIN/ST/")
options(bitmapType = "cairo")

# LOAD PACKAGES -----------------------------------------------------------
{
  library(Seurat)
  library(ggsignif)
  library(purrr)
  library(ggplot2)
  library(tidyverse)
  library(RColorBrewer)
  library(ggtern)
  library(scatterpie)
  library(magrittr)
  library(ggnewscale)
  library(cowplot)
  library(ggridges)
  library(DelayedMatrixStats)
  library(grid)
  library(readbitmap)
  library(foreach)
  library(parallel)
  library(doParallel)
  library(patchwork)
  library(ggpubr)
}

# SAMPLES TO LOAD ---------------------------------------------------------
samples <-
  c(
    "ITG",
    "V1_1",
    "V1_2",
    "V1_3",
    "M1",
    "FPPFC_2",
    "FPPFC_1",
    "S1_1",
    "S1_2",
    "DLPFC_1",
    "DLPFC_2"
  )

sample_rep <-
  c("ITG",
    "V1_1",
    "FPPFC_1",
    "M1",
    "S1_1",
    "DLPFC_2")

sample_rep2 <-
  c("ITG",
    "V1_1",
    "FPPFC_1",
    "S1_1",
    "DLPFC_2")

# LOAD FUNCTIONS ----------------------------------------
source("./1-ANALY/func.R", echo = T)

# Seurat Reading H5 data --------------------------------
ST_paths <- str_c("./0-DATA/", samples)
ST_list <- ReadData(samples = samples, ST_paths = ST_paths)

# SCT Integration ---------------------------------------
st.combined.sct <- integratedData(ST_list)
st.combined.sct$globalKey <- rownames(st.combined.sct@meta.data)

# Add Layer Annotation ----------------------------------
Anno_paths <- str_c("1-ANALY/Anno/", samples, ".csv")

# ST_list %<>% purrr::map( ~ addAnno(., Anno_paths, samples))
# 
# ST_list %<>% purrr::map(function(x) {
#   x$Layer2 <- x$Layer %>% str_remove("[ab]")
#   x$globalKey <- rownames(x@meta.data)
#   x
# })
# ST_list <- purrr::map(ST_list,
#                       function(x) {
#                         x$globalKey <- 1
#                         x$globalKey <- NULL
#                         x
#                       })

# Seurat standard pipeline ------------------------------------------------
st.combined.sct %<>% RunPCA(verbose = FALSE)
st.combined.sct %<>% RunUMAP(
  n.neighbors = 50,
  min.dist = 0.05,
  reduction = "pca",
  dims = 1:20,
  verbose = FALSE
)

# Add Layer Annotation ----------------------------------
Anno_paths <- str_c("1-ANALY/Anno/", samples, ".csv")

st.combined.sct <- addAnno(st.combined.sct, Anno_paths, samples)

st.combined.sct$Layer2 <-
  st.combined.sct$Layer %>% str_remove("[ab]")

Layers <-  st.combined.sct$Layer2 %>% unique %>% sort

# Add depth information ---------------------------------------------------

annotations <- st.combined.sct$Layer %>% unique %>% sort
annotations_wt <- lapply(str_extract_all(annotations, "[1-6]"),
                         function(x) {
                           mean(as.numeric(x))
                         }) %>% unlist

names(annotations_wt) <- annotations
annotations_wt["WM"] <- 7 #

st.combined.sct$depth <- annotations_wt[st.combined.sct$Layer]

# Add deconvolution information -------------------------------------------

Decon_paths <- str_c("1-ANALY/Decon/", samples, ".csv")

decon <- 
  purrr::map2(
  Decon_paths,
  samples,
  ~ read_csv(.x) %>%
    dplyr::rename(Barcode = 1) %>%
    mutate(sample = .y,
           globalKey = str_c(.y, Barcode, sep = "@"))
) %>% purrr::reduce(rbind)  %>% column_to_rownames("globalKey")

st.combined.sct@meta.data %<>% cbind(decon[rownames(st.combined.sct@meta.data),])

# Add cell type information -----------------------------------------------
cellType.meta <- read_delim("0-DATA/CellType_meta.txt")
rownames(cellType.meta) <- cellType.meta$anno
cellTypes <- cellType.meta$anno

# Add color configuration -------------------------------------------------
color_layer <-
  c(ggsci::pal_aaas()(10), ggsci::pal_npg()(10))[1:length(annotations)]

names(color_layer) <- annotations

# color for different levels of celltypes

color_cluster <-
  read_delim("0-DATA/cluster_color.txt") %$%
  {
    color %>%
      set_names(cellType.meta$anno %>%
                  set_names(cellType.meta$newID) %>%
                  `[`(as.character(cluster)))
  }

color_subclass <-
  color_subclass <-
  read_delim("0-DATA/subclass_color.txt") %$% {
    color %>% set_names(subclass)
  }

others <- c(EXC = "#70A3BE", INH = "#C638A9")

color_ct <- c(color_cluster, color_subclass, others)

color_sample <-
  c(ggsci::pal_npg()(10), ggsci::pal_aaas()(10))[1:length(samples)]
names(color_sample) <- samples

ColorBar <- colorRampPalette(rev(c(brewer.pal(11, "Spectral"),"#808080","#FFFFFF")))
ColorBar2 <- colorRampPalette(c("#4E0707",brewer.pal(5, "Spectral")))

ct.order <- names(color_cluster)

#!!! DeepST Umap -------------------------------------------------------------
## process original umap
# {umap_sample <- read.csv("./1-ANALY/Umap_DeepST/batch_name.csv")
# # #umap坐标转换46948行xy坐标
# umap_xy <-
#   read.csv("./1-ANALY/Umap_DeepST/X_umap.csv") %>% left_join(umap_sample, by = c("X" = "X"))
# cell.id <-
#   stringr::str_split(umap_xy$X, "-", n = 3, simplify = T)[, c(1, 2)]
# cell.id <- stringr::str_c(cell.id[, 1], cell.id[, 2], sep = "-")
# umap_xy$cell.id <- cell.id
# names(umap_xy) <-
#   c("orig.id", "umap_x", "umap_y", "sample", "cell.id")
# trans_name <-
#   c(
#     'DLPFC_1',
#     'DLPFC_2',
#     'FPPFC_1',
#     'FPPFC_2',
#     'ITG',
#     'M1',
#     'S1_1',
#     'S1_2',
#     'V1_1',
#     'V1_2',
#     'V1_3'
#   )
# names(trans_name) <-
#   c('DLPFC',
#     'DLPFC2',
#     'FPPFC',
#     'FPPFC2',
#     'ITG',
#     'M1',
#     'S1',
#     'S1_2',
#     'V1',
#     'V1_2',
#     'V1_3')
# umap_xy$sample <- trans_name[umap_xy$sample]
# umap_xy$globalKey <-
#   str_c(umap_xy$sample, umap_xy$cell.id, sep = "@")
# 
# write.csv(umap_xy, "1-ANALY/Umap_DeepST/umap_DeepST.csv")
# }

# replace the original umap from Seurat
umap_xy <- read.csv( "1-ANALY/Umap_DeepST/umap_DeepST.csv")

embedding <- umap_xy[c("umap_x", "umap_y", "globalKey")] %>%
  column_to_rownames("globalKey") %>%
  as.matrix() %>%
  `[`(rownames(st.combined.sct@meta.data), )

spotAll <- merge(spotAll,embedding, by = 0)
# attributes(embedding) <-
#   attributes(st.combined.sct@reductions$umap@cell.embeddings)
# st.combined.sct@reductions$umap@cell.embeddings <- embedding

# DeepST reduction --------------------------------------------------------

#
# reduction_DpST <-
#   read_csv(file = "1-ANALY/Umap_DeepST/all_11_brain_deepst_data.csv") %>%
#   left_join(umap_sample, by = c("...1" = "X"))
# cell.id <-
#   stringr::str_split(reduction_DpST$...1, "-", n = 3, simplify = T)[, c(1, 2)]
# cell.id <- stringr::str_c(cell.id[, 1], cell.id[, 2], sep = "-")
# reduction_DpST$cell.id <- cell.id
# names(reduction_DpST)[1] <- c("orig.id")
# reduction_DpST$sample <- trans_name[reduction_DpST$batch_name]
# reduction_DpST$globalKey <-
#   str_c(reduction_DpST$sample, reduction_DpST$cell.id, sep = "@")
#

# Extract data for convenient plotting --------------------------------------------------------------
coords <-
  purrr::map(
    samples,
    ~ GetTissueCoordinates(st.combined.sct, image = .) 
  ) %>%
  purrr::reduce(rbind)

spotAll <-
  Seurat::FetchData(st.combined.sct,
            vars = c(colnames(st.combined.sct@meta.data))) %>%
  cbind(coords[rownames(.),]) %>% rownames_to_column("globalKey") %>% 
  left_join(umap_xy,by = c("globalKey" = "globalKey"))

spotAll_long <- spotAll %>%
  pivot_longer(cols = cellTypes,
               names_to = "cellType",
               values_to = "probability") %>%
  left_join(cellType.meta, by = c("cellType" = "anno"))

# spotLayerSample <- spotAll %>% group_split(Layer, sample)
#
# # keep only top 90%
# spotLayer <-
#   spotAll %>% group_split(Layer2) %>% purrr::map( ~ slice_max(., order_by = nCount_Spatial, prop = .9))
# # spotLayer2 <- spotAll %>% group_split(Layer2) %>% purrr::map(~ slice_max(., order_by = nFeature_SCT, prop = .9))
#
# filtered_spots <-
#   spotLayer %>% purrr::map( ~ .$globalKey)  %>% unlist()
# # filtered_spots2 <- spotLayer2 %>% purrr::map(~.$globalKey)  %>% unlist()
# st.combined.sct2 <-
#   subset(st.combined.sct, globalKey %in% filtered_spots)

# Caculating the sd and mean of depth per sample -------------------------------------

cellType_stat <- spotAll_long %>%
  group_by(level1, level2, level3, level4, cellType, sample) %>%
  dplyr::summarise(# sum_p = sum(probability.norm),
    sd = weightedSd(depth, w = probability),
    m =  weighted.mean(depth, w = probability)) %>%
  mutate(min_ = m - sd,
         max_ = m + sd)


cellType_stat$cellType %<>% factor(
  levels = cellType_stat %>% group_by(cellType) %>% summarise(median = median(sd)) %>% arrange(median) %$% cellType)

# a complicated method to get the mode of Layer for each cell type
# stat based on downsample
cellType_stat2 <- {
  spotAll_long %>%
    group_by(sample, Layer2, cellType) %>%
    slice_sample(
      n =  spotAll %>%
        group_by(sample, Layer2) %>%
        summarise(Layer_n = n(), ) %>%
        ungroup() %$% max(Layer_n) * 2,
      replace = T
    )
} %>%
  group_by(sample, cellType) %>%
  slice_max(prop = 0.05, order_by = probability.norm)  %>%
  dplyr::summarise(
    sd = weightedSd(depth, w = 1000 * probability.norm),
    # num is too small so expand it and this will not affect the result while using weightSD
    m =  weighted.mean(depth, w = probability.norm),
    m_ = mean(depth),
    median = median(depth),
    most = which.max(table(depth))
  )  %>% left_join(cellType.meta, by = c("cellType" = "anno")) %>%
  mutate(min_ = m - sd,
         max_ = m + sd)

cellType_depth <- cellType_stat2 %>% 
  group_by(cellType, level4) %>% 
  summarise(depth_m = mean(m)) %>%
  column_to_rownames("cellType")

# Read Outline ------------------------------------------------------------

line_paths <-
  str_c("1-ANALY/outline/",
        samples, ".png")

outline <- list()

for (i in 1:length(samples)) {
  # print(line_paths[[i]])
  outline[[samples[[i]]]] = LoadOutline(st.combined.sct, samples[[i]], line_paths[[i]])
  print(outline[[i]] %>%
          environment()  %>%
          as.list %$% line_paths)
}

# Read image information --------------------------------------------------------------

ST_image_paths <- file.path(ST_paths, "spatial")
images_info <- list()
for (i in 1:length(samples)) {
  print(samples[i])
  images_info[[samples[i]]] <-
    Read10X_Image(ST_image_paths[i], filter.matrix = F)
}

# scale_f <- 0.05
images_lim <- purrr::map(
  images_info,
  ~ .@scale.factors$lowres * data.frame(
    y_max = max(.@coordinates["imagerow"]),
    y_min = min(.@coordinates["imagerow"]),
    x_max = max(.@coordinates["imagecol"]),
    x_min = min(.@coordinates["imagecol"])
  )
)%>% purrr::reduce(rbind) %>% mutate(
  sample = samples,
  dis_y = (y_max - y_min) / 78,
  dis_x = (x_max - x_min) / 64,
  down = y_min - dis_y * 3,
  up = y_max + dis_y * 3,
  left = x_min - dis_x * 3,
  right = x_max + dis_x * 3
) %>%
  column_to_rownames("sample")

# Read image --------------------------------------------------------------

images <- list()
for (sample in samples) {
  # print(sample)
  images[[sample]] = alphaImg(st.combined.sct,
                              sample)
  print(images[[sample]] %>%
          environment()  %>%
          as.list %$% sample)
}

# gene list ---------------------------------------------------------------

gene_list <-
  c(
    "MBP",
    "SNAP25",
    "CCN2",
    # "PLP1",
    "THEMIS",
    # "NTNG2",
    "PCP4",
    "RORB",
    "COL5A2",
    # "ADCYAP1",
    "CUX2",
    "LAMP5",
    "PCDH8",
    "LINC00507",
    # "AQP4",
    "RELN"
  )



# Theme configuration -----------------------------------------------------

theme_umap <- function(gg) {
  gg + coord_cartesian(expand = F) +
    theme_cowplot() + theme(
      plot.margin = margin(
        t = 0,
        r = 0,
        b = 0,
        l = 0
      ),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.direction = "horizontal",
      legend.position = "top",
      legend.justification = "center",
      legend.margin = margin(
        t = 0,
        r = 0,
        b = 0,
        l = 0
      ),
      legend.box.spacing = unit(0, units = "pt")
    )
}


# Umap of ST --------------------------------------------------------------

labels.loc <-
  FetchData(st.combined.sct, vars = c("UMAP_1", "UMAP_2", "ident")) %>% group_by(ident) %>% summarise(x = mean(UMAP_1), y = mean(UMAP_2))


# Read alignment of each cells across all slides  -----------------------

cellAlign<- str_c("0-DATA/cellAlignment/",samples) %>% 
  {.[file.exists(.)]} %>%
  purrr::map(~readRDS(.)) 

cl <- parallel::makeCluster(10)
doParallel::registerDoParallel(cl)
cellAlign <-
  foreach(data = cellAlign,.packages = c("tidyverse","purrr")) %dopar% {
    tmp <-reduce(data,rbind)  
    }
parallel::stopCluster(cl)

cellAlign <- reduce(cellAlign,rbind)


# EXC tsne Data -----------------------------------------------------------
{
  EXC_6region <-
    readRDS("../../REMOTE_Broken/EXC_6region/EXC_6region_harmony.rds")
  
  EXC_6region@meta.data <- EXC_6region@meta.data[c("hicat_cluster_merge_anno","region")]
  
  EXC_6region@meta.data$sample <- 
    read.csv("0-DATA/Region.csv", row.names = "No") %>% 
    `[`(EXC_6region@meta.data$region,"Abbr")
  
  EXC_6region@meta.data$aver.depth <- cell_average_depth[rownames(EXC_6region@meta.data),"aver.depth"]

}


# Read SnRNA data ------------------

#cell alignment to spots and calculate the average depth across 11 slides
cell_to_spot <- list.files("/home/luomeng/DATA/BRAIN/SnRNA/",pattern = "cellAlignment.rds") %>% purrr::map(readRDS) %>% purrr::reduce(rbind)
cell_to_spot$globalKey <- with(cell_to_spot,paste(sample, spotID, sep = "@"))
cell_to_spot$depth <- spotAll["depth"][cell_to_spot$globalKey,]
cell_average_depth <- cell_to_spot  %>% drop_na(depth) %>%  group_by(cellID) %>% summarise(aver.depth = mean(depth)) %>% column_to_rownames("cellID")

umap_all <- readRDS("../SnRNA/seu.harmony.anno.meta.v2.rds") %>% rownames_to_column("cells") %>% select(cells,UMAP_1,UMAP_2)
umap_INHCGE <- readRDS("../SnRNA/INH_CGE_UMAP.rds")
umap_INHMGE <- readRDS("../SnRNA/INH_MGE_UMAP.rds")
umap_EXC <- readRDS("../SnRNA/EXC_UMAP.rds")
umap_NoN <- readRDS("../SnRNA/NonNeuron_UMAP.rds")

umap_all$aver.depth <- cell_average_depth[umap_all$cells,"aver.depth"]
umap_INHCGE$aver.depth <- cell_average_depth[umap_INHCGE$cells,"aver.depth"]
umap_INHMGE$aver.depth <- cell_average_depth[umap_INHMGE$cells,"aver.depth"]
umap_EXC$aver.depth <- cell_average_depth[umap_EXC$cells,"aver.depth"]
umap_NoN$aver.depth <- cell_average_depth[umap_NoN$cells,"aver.depth"]


