library(speciesTree)
library(conos)
library(pagoda2)


# read example files
# gene expression profile
expression.matrix0 <- readRDS("./example_data/integrated_matrix_all_cells.RDS")
# mete cell clusters
meta_clusters0 <- readRDS("./example_data/all_cells_cluster_membership.RDS")
# read cell annotations
upperlevelinfo0 <- readRDS("./example_data/upperlevelinfo.rds")
species_annot0 <- sapply(strsplit(names(meta_clusters0), "_"), `[`, 1)
names(species_annot0) <- names(meta_clusters0)
# read subsampled cell clusters
subsampledClusters0 <- readRDS("./example_data/subsampledClusters.rds")

# build tree
d0 <- cluster.matrix.expression.distances(t(expression.matrix0), groups=meta_clusters0, dist="cor", 
                                         useVariablegenes=FALSE,  use.scaled.data=TRUE)
dendr0 <- hclust(as.dist(d0), method='ward.D2')
dend0 <- as.dendrogram(dendr0)
tree0 <- buildSpeciesTree(dend=dend0, 
                         expMatrix=t(expression.matrix0), 
                         subsampleClusters=subsampledClusters0,
                         cls.groups=meta_clusters0, 
                         upperlevelannot=upperlevelinfo0, 
                         species=species_annot0,
                         n.cores=10)

# prune tree
tree0$dend <- TreeEntropy(tree0$dend, entropy.cutoff = 2.9)
dend_pruned0 <- pruneTreeEntropy(tree0$dend, cutoff=2.9)

s

# get homologous clusters
homo.clusters <- getClusters(dend_pruned, plotTree=TRUE)




### EXC #################

seu <- readRDS('/data/share/Brain/h_m_SCT/EXC_h_mouse10x_newCluster_SCT_100_2_marker.rds')

expression.matrix <- as(seu@assays$integrated@scale.data, "dgCMatrix")
meta_clusters <- paste0('c',as.vector(seu$seurat_clusters))
names(meta_clusters) <- colnames(seu)
#upperlevelinfo <- seu$level3
species_annot <- seu$species
#upperlevelinfo <- c(seu$level3[seu$species=='human'],
#                    seu$subclass_label[seu$species=='mouse'])
#upperlevelinfo[upperlevelinfo %in% c('L2/3 IT CTX-1','L2/3 IT CTX-2')] <- 'L2/3_IT'
#upperlevelinfo[upperlevelinfo %in% c('L5 IT CTX', 'L5_IT_RORB')] <- 'L5_IT'
#upperlevelinfo[upperlevelinfo %in% c('L6b CTX','L6B_FEZF2')] <- 'L6B'
#upperlevelinfo[upperlevelinfo %in% c('L6 CT CTX','L6CT_FEZF2')] <- 'L6_CT'
#upperlevelinfo[upperlevelinfo %in% c('L6_IT','L6 IT CTX')] <- 'L6_IT'
#upperlevelinfo[upperlevelinfo %in% c('L1-L3_IT_LINC00507')] <- 'L1-L3_IT'
#upperlevelinfo[upperlevelinfo %in% c('L5 PT CTX')] <- 'L5_PT'
#upperlevelinfo[upperlevelinfo %in% c('L5 NP CTX')] <- 'L5_NP'
#upperlevelinfo[upperlevelinfo %in% c('L4/5 IT CTX')] <- 'L4/5_IT'

upperlevelinfo <- as.vector(seu$seurat_clusters)

upperlevelinfo[upperlevelinfo %in% c("41","43","8", "18")] <- 'X1'
upperlevelinfo[upperlevelinfo %in% c("23","44","39","5", "6")] <- 'X2'
upperlevelinfo[upperlevelinfo %in% c("36","7", "24","19","27")] <- 'X3'
upperlevelinfo[upperlevelinfo %in% c("35","34","37")] <- 'X4'
upperlevelinfo[upperlevelinfo %in% c("47","15","4","20","0", "33")] <- 'X5'
upperlevelinfo[upperlevelinfo %in% c("21","46","17","31")] <- 'X6'
upperlevelinfo[upperlevelinfo %in% c("2", "42","11","3", "26")] <- 'X7'
upperlevelinfo[upperlevelinfo %in% c("10","25","16","45")] <- 'X8'
upperlevelinfo[upperlevelinfo %in% c("28","14","22","30")] <- 'X9'
upperlevelinfo[upperlevelinfo %in% c("9", "12","29","38","40")] <- 'X10'
upperlevelinfo[upperlevelinfo %in% c("48","1", "13","32")] <- 'X11'


subsampledClusters <- list()
for(i in 1:10){
  subsampledClusters[[i]] <- sample(meta_clusters,28224)
}

# build tree
d <- cluster.matrix.expression.distances(t(expression.matrix), groups=meta_clusters, dist="cor", 
                                         useVariablegenes=FALSE,  use.scaled.data=TRUE)
dendr <- hclust(as.dist(d), method='ward.D2')
dend_ <- as.dendrogram(dendr)

saveRDS(dend_, 'EXC_h_mouse10x_speciesTree_dend.rds')

# EXC_h_m_SCT_dend_speciesTree.pdf 4*15
plot(dend_)


tree <- buildSpeciesTree(dend=dend_, 
                          expMatrix=t(expression.matrix), 
                          subsampleClusters=subsampledClusters,
                          cls.groups=meta_clusters, 
                          upperlevelannot=upperlevelinfo, 
                          species=species_annot,
                          n.cores=10)

# prune tree
tree$dend <- TreeEntropy(tree$dend, entropy.cutoff = 2.9)
dend_pruned <- pruneTreeEntropy(tree$dend, cutoff=2.9)


## pie 
piec<-c()
aa <- table(seu$seurat_clusters,seu$species)
bb <- aa/rowSums(aa)
piec <- melt(bb)

colnames(piec)<-c('cluster','species','value')


piec$title<-factor(piec$cluster,
                   levels = labels(dend_))

#color<-c('#4E79A7','#F28E2B','#E15759','#76B7B2','#59A14F','#EDC948','#B07AA1')
# EXC_h_m_SCT_pie_speciesTree.pdf 3*30
#ggplot(piec, aes("", value, fill = Cluster)) + 
ggplot(piec, aes(x=1, y=value, fill=species)) + 
  geom_bar(stat = "identity", color = "gray30", size = 0.5) +
  #geom_text(aes(label = paste0(value * 100, "%")), 
  #          position = position_stack(vjust = 0.5), 
  #          color = "white", size = 3) +
  coord_polar(theta = "y") +
  facet_wrap(~ title, ncol = 49) +
  #scale_fill_manual(values = color) +
  theme_void()


tree <- buildSpeciesTree(dend=dend_, 
                         expMatrix=t(expression.matrix), 
                         subsampleClusters=subsampledClusters,
                         cls.groups=meta_clusters, 
                         upperlevelannot=upperlevelinfo, 
                         species=species_annot,
                         n.cores=10)

dend=dend_
expMatrix=t(expression.matrix)
subsampleClusters=subsampledClusters
cls.groups=meta_clusters
upperlevelannot=upperlevelinfo
species=species_annot
n.cores=10

mappingcells=FALSE
cellannot=NULL
#upperlevelannot=NULL
renameCluster=TRUE
plot=TRUE


buildSpeciesTree <- function(dend, expMatrix, subsampleClusters=NULL, cls.groups, mappingcells=FALSE, cellannot=NULL, species,
                             upperlevelannot=NULL, renameCluster=TRUE, plot=TRUE, n.cores=10){
  dendr <- TransferDend(dend, renameCluster=renameCluster, cls.groups = cls.groups)
  cls.groups <- dendr$new.groups
  
  dend <- dendr$dendrogram
  leafcontent <- dendr$leafcontent
  
  if(!is.null(subsampleClusters)){
    subsampled.dend <- subSampleTree(expMatrix, subsample.groups=subsampleClusters)
    stability.measurements <- TreeStabilityDend(dend, cls.groups=cls.groups, subsampled.dend, n.cores=n.cores)
    dend <- stability.measurements$dendrogram
  } else{
    stability.measurements = NULL
  }
  
  # add cluster attribute to dendrogram
  dend <- AddTreeAttribute(dend, species, leafcontent)
  dend <- dendSetWidthBysize(dend, scale=8)
  
  if(mappingcells){
    if(is.null(cellannot)){
      stop("Please provide cell annotations")
    }
    leaflabels <- mappingLabel(dend, leafcontent, cellannot, humanAnnot=T)
    # set labels
    dend <- set_labels(dend, paste(dend %>% labels(), leaflabels, sep=" "))
  }
  
  # add upperlevel info to each nodes
  if(!is.null(upperlevelannot[1])){
    dend <- UpperLevelInfo(dend, cellannot=upperlevelannot, leafcontent, propCutoff = 0.1)
    upperLevelnodes <- getUpperLevelNode(dend, cutoff=0.65)
    # normalize Tree
    dend <- NormTree(dend, upperLevelnodes, upperlevelannot, species)
    dend <- dendSetColorByNormMixing(dend)
  } else{
    # need to modify
    colorpallete <- colorRampPalette(c("blue", "grey", "grey", "grey", "red"))(101)
    dend <- dendSetColorByMixing(dend, species, leafContent, colorpallete)
    upperLevelnodes = NULL
  }
  if(plot){
    if(!is.null(upperlevelannot) & !is.null(subsampleClusters)){
      #par(cex=1, mar=c(20, 5, 0, 10))
      plot(dend)
      text(upperLevelnodes$xy, labels=upperLevelnodes$upperlabel, adj=c(0.5, -1.2), cex=0.8, col="red")
      text(get_nodes_xy(dend), labels=get_nodes_attr(dend, "stability"), adj=c(0.4,0.4), cex=0.5, col="red")
    } else if(!is.null(upperlevelannot)){
      #par(cex=1, mar=c(20, 5, 0, 10))
      plot(dend)
      text(get_nodes_xy(dend), labels=get_nodes_attr(dend, "stability"), adj=c(0.4,0.4), cex=0.5, col="red")
    } else if(!is.null(subsampleClusters)){
      #par(cex=1, mar=c(20, 5, 0, 10))
      plot(dend)
      text(stability.measurements$stability.loc, labels=stability.measurements$stability.labels,
           adj=c(0.4, 0.1), cex=0.35, col="red")
    }
    else{
      plot(dend)
    }
  }
  return(list(dend=dend, upperLevelnodes=upperLevelnodes, leafcontent=leafcontent))
}


dendSetColorByNormMixing <- function(d, withinFacNorm=FALSE){
  cc2col <- function(cc, rate=15, base=0.001) {
    if(length(cc)==2) { # 2-color
      cv <- cc
      cv <- dexp(c(cv[1], 0, cv[2]), rate) / base * (1-0.001)
      #cv <- c(cc[1],0,cc[2])+base; cv <- cv/max(cv) * (1-base)
      #rgb(base+cc[2],base,base+cc[3],1)
      adjustcolor(rgb(cv[1], cv[2], cv[3], 1), offset = c(0.5, 0.5, 0.5, 0.1))
    } else if(length(cc)==3) { # 3-color
      #cv <- 1 - cc
      cv <- cc
      cv <- dexp(cv, rate)
      #cv <- cv/max(cv) * (1-0.2)
      cv <- cv/rate * (1-base)
      adjustcolor(rgb(cv[1],cv[2],cv[3], 1), offset = c(0.5, 0.5, 0.5, 0.1))
    }
  }
  
  cbm <- function(d,fac){
    if(is.leaf(d)) {
      if(withinFacNorm){
        normpercent <- attr(d, "withinFacPercent")
      } else{
        normpercent <- attr(d, "normPercentage")
      }
      
      col <- cc2col(normpercent)
      attr(d,"edgePar") <- c(attr(d,"edgePar"),list(col=col))
      return(d);
    } else {
      oa <- attributes(d);
      d <- lapply(d,cbm);
      attributes(d) <- oa;
      if(withinFacNorm){
        normpercent <- attr(d, "withinFacPercent")
      } else{
        normpercent <- attr(d, "normPercentage")
      }
      
      col <- cc2col(normpercent)
      attr(d,"edgePar") <- c(attr(d,"edgePar"),list(col=col))
      return(d);
    }
  }
  cbm(d)
}


cc2col <- function(cc, rate=15, base=0.001){
  if(sum(cc)==0) {
    cc <- rep(1, length(cc))
  } else {
    # normalized by total number of cells
    if(length(cc) == 3){
      if(normTofac){
        cc <- round((cc[2:3]/totalCells)/sum(cc[2:3]/totalCells), 2)
      } else{
        cc <- round((cc[2:3])/sum(cc[2:3]), 2)
      }
      cv <- round(cc[1]*100, 0)
      colorpallete <- colorRampPalette(c("blue", "grey", "grey", "grey", "red"))(101)
      col <- colorpallete[cv + 1]
      
    } else{
      if(normTofac){
        cc <- round((cc[2:4]/totalCells)/sum(cc[2:4]/totalCells), 2)
      } else{
        cc <- round((cc[2:4])/sum(cc[2:4]), 2)
      }
      cv <- cc
      cv <- dexp(cv, rate)
      cv <- cv/rate * (1-base)
      col <- adjustcolor(rgb(cv[1],cv[2],cv[3], 1), offset = c(0.5, 0.5, 0.5, 0.1))
    }
  }
  return(col)
}


# prune tree
tree$dend <- TreeEntropy(tree$dend, entropy.cutoff = 2.9)
dend_pruned <- pruneTreeEntropy(tree$dend, cutoff=2.9)

# get homologous clusters
homo.clusters <- getClusters(dend_pruned, plotTree=TRUE)








### INH #################

seu <- readRDS('/data/share/Brain/h_m_SCT/INH_h_mouse10x_newCluster_SCT_100_2_marker.rds')

expression.matrix <- as(seu@assays$integrated@scale.data, "dgCMatrix")
meta_clusters <- seu$seurat_clusters
#upperlevelinfo <- seu$level3
species_annot <- seu$species
upperlevelinfo <- c(seu$level4[seu$species=='human'],
                    seu$subclass_label[seu$species=='mouse'])
upperlevelinfo[upperlevelinfo %in% c('Lamp5','LAMP5/SV2C')] <- 'LAMP5'
upperlevelinfo[upperlevelinfo %in% c('Pvalb')] <- 'PVALB'
upperlevelinfo[upperlevelinfo %in% c('Sncg','LAMP5/SNCG')] <- 'SNCG'
upperlevelinfo[upperlevelinfo %in% c('Sst','Sst Chodl')] <- 'SST'
upperlevelinfo[upperlevelinfo %in% c('Vip')] <- 'VIP'

subsampledClusters <- list()
for(i in 1:10){
  subsampledClusters[[i]] <- sample(meta_clusters,30631)
}

# build tree
d <- cluster.matrix.expression.distances(t(expression.matrix), groups=meta_clusters, dist="cor", 
                                         useVariablegenes=FALSE,  use.scaled.data=TRUE)
dendr <- hclust(as.dist(d), method='ward.D2')
dend_ <- as.dendrogram(dendr)

saveRDS(dend_, 'INH_h_mouse10x_speciesTree_dend.rds')

# INH_h_m_SCT_dend_speciesTree.pdf 4*15
plot(dend_)

## pie 
piec<-c()
aa <- table(seu$seurat_clusters,seu$species)
bb <- aa/rowSums(aa)
piec <- melt(bb)

colnames(piec)<-c('cluster','species','value')


piec$title<-factor(piec$cluster,
                   levels = labels(dend_))

#color<-c('#4E79A7','#F28E2B','#E15759','#76B7B2','#59A14F','#EDC948','#B07AA1')
# INH_h_m_SCT_pie_speciesTree.pdf 3*30
#ggplot(piec, aes("", value, fill = Cluster)) + 
ggplot(piec, aes(x=1, y=value, fill=species)) + 
  geom_bar(stat = "identity", color = "gray30", size = 0.5) +
  #geom_text(aes(label = paste0(value * 100, "%")), 
  #          position = position_stack(vjust = 0.5), 
  #          color = "white", size = 3) +
  coord_polar(theta = "y") +
  facet_wrap(~ title, ncol = 53) +
  #scale_fill_manual(values = color) +
  theme_void()


