library(Seurat)

seu<-readRDS('seu.harmony.rds')

result <- readRDS('seu.hicat.de80.niter100.rds')
dend.result <- readRDS('seu.hicat.de80.niter100.dend.rds')

seu$hicat_cluster <- result$cl.result$cl[colnames(seu)]
seu$hicat_cluster[seu$hicat_cluster %in% c(279,280)] <- '279_280'
seu$hicat_cluster[seu$hicat_cluster %in% c(235,243)] <- '235_243'
seu$hicat_cluster[seu$hicat_cluster %in% c(204,211)] <- '204_211'
seu$hicat_cluster[seu$hicat_cluster %in% c(146,157)] <- '146_157'
seu$hicat_cluster[seu$hicat_cluster %in% c(144,156)] <- '144_156'
seu$hicat_cluster[seu$hicat_cluster %in% c(160,162)] <- '160_162'
seu$hicat_cluster[seu$hicat_cluster %in% c(329,333)] <- '329_333'

seu$hicat_cluster <- factor(seu$hicat_cluster,
                            levels = labels(dend.result$dend))
seu@active.ident <- seu$hicat_cluster

markers <- FindAllMarkers(seu)
saveRDS(markers, 'seu.harmony.markers.rds')

#seu<-readRDS('seu.harmony.anno.Exc.rds')
#markers <- FindAllMarkers(seu)
#saveRDS(markers, 'seu.harmony.anno.Exc.markers.rds')
#
#seu<-readRDS('seu.harmony.anno.Inh_SV2C.rds')
#markers <- FindAllMarkers(seu)
#saveRDS(markers, 'seu.harmony.anno.Inh_SV2C.markers.rds')
#
#seu<-readRDS('seu.harmony.anno.Inh_unassigned1.rds')
#markers <- FindAllMarkers(seu)
#saveRDS(markers, 'seu.harmony.anno.Inh_unassigned1.markers.rds')
#
#seu<-readRDS('seu.harmony.anno.Inh_VIP.rds')
#markers <- FindAllMarkers(seu)
#saveRDS(markers, 'seu.harmony.anno.Inh_VIP.markers.rds')
#
#seu<-readRDS('seu.harmony.anno.Inh_SST.rds')
#markers <- FindAllMarkers(seu)
#saveRDS(markers, 'seu.harmony.anno.Inh_SST.markers.rds')
#
#seu<-readRDS('seu.harmony.anno.Inh_PVALB.rds')
#markers <- FindAllMarkers(seu)
#saveRDS(markers, 'seu.harmony.anno.Inh_PVALB.markers.rds')
#
#
#
#
## Inh_unassigned1 vs Inh_VIP
#cell1 <- colnames(seu)[seu$hicat_cluster %in% c('93','183','177','180','186','191','192',  # Inh_unassigned1
#                                                '203','204','200','201','139','194')]
#cell2 <- colnames(seu)[seu$hicat_cluster %in% c('129','131','162','174','126','127','133',  # Inh_VIP
#                                                '137','135','167','141','318')]
#de<-FindMarkers(seu,
#                ident.1 = c('93','183','177','180','186','191','192',  # Inh_unassigned1
#                            '203','204','200','201','139','194'))
#saveRDS(de,'Inh_unassigned1_de.rds')
#
#
#de <- readRDS('Inh_unassigned1_vs_VIP_de.rds')
#genes <- rownames(de)[de$avg_log2FC>0][1:40]
#DotPlot(seu, features = unique(genes), assay = 'RNA') + RotatedAxis()
#
#
#
#dend <- readRDS('seu.hicat.run_consensus_clust.de100.niter100.dend.rds')
#


#seu<-readRDS('seu.harmony.anno.rds')
##seu@active.ident <- seu$hicat_cluster_newID
#markers <- FindAllMarkers(seu)
#saveRDS(markers, 'seu.harmony.anno.newID.markers.rds')
#
#
#
#seu<-readRDS('seu.harmony.anno.Astro.rds')
##seu@active.ident <- seu$hicat_cluster_newID
#markers <- FindAllMarkers(seu)
#saveRDS(markers, 'seu.harmony.anno.Astro.markers.rds')
#
#seu<-readRDS('seu.harmony.anno.Micro.rds')
##seu@active.ident <- seu$hicat_cluster_newID
#markers <- FindAllMarkers(seu)
#saveRDS(markers, 'seu.harmony.anno.Micro.markers.rds')
#
#
##seu<-readRDS('seu.harmony.anno.Glu.rds')
##seu@active.ident <- seu$hicat_cluster
##markers <- FindAllMarkers(seu)
##saveRDS(markers, 'seu.harmony.anno.Glu.markers.rds')
#
#
#
#CGE <- readRDS('seu.harmony.anno.GABA.CGE.newID.markers.rds')
#MGE <- readRDS('seu.harmony.anno.GABA.MGE.newID.markers.rds')
#Glu <- readRDS('seu.harmony.anno.Glu.newID.markers.rds')
#Non <- readRDS('seu.harmony.anno.nonNeuron.newID.markers.rds')
#
#CGE_top5 <- CGE %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
#
##CGE_top5_dotplot 20x25
#DotPlot(seu,features = unique(CGE_top5$gene)) + RotatedAxis()
#
#
#MGE_top5 <- MGE %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
#
##MGE_top5_dotplot 20x25
#DotPlot(seu,features = unique(MGE_top5$gene)) + RotatedAxis()
#
#
#Glu_top5 <- Glu %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
#
##Glu_top5_dotplot 20x25
#DotPlot(seu,features = unique(Glu_top5$gene)) + RotatedAxis()
#
#
#Non_top5 <- Non %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
#
##NonNueron_top5_dotplot 20x25
#DotPlot(seu,features = unique(Non_top5$gene)) + RotatedAxis()
#
#
#
#
#
#
#