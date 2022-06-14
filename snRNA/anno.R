library(Seurat)

seu <- readRDS('seu.harmony.v2.rds')

result<-readRDS('seu.hicat.de80.niter100.rds')
dend.result<-readRDS('Brain_subClusteringID_seurat_merge_dend_renames.rds')

barcode_df2 <- readRDS('Brain_subClusteringID_seurat_merge_v2.rds')

seu$hicat_cluster <- result$cl.result$cl[colnames(seu)]
seu$hicat_cluster <- factor(seu$hicat_cluster,
                            levels = rev(labels(dend.result$dend)))

seu$hicat_cluster_merge <- barcode_df2[colnames(seu),'subClusteringID_merge']
seu$hicat_cluster_merge <- factor(seu$hicat_cluster_merge, 
                                  levels = rev(labels(dend.result$dend)))

seu@active.ident <- seu$hicat_cluster_merge




clusters <- c("22", "19", "17", "1", "21", "395", 
              "392", "389", "386", "477", "476", "478", 
              "FEZF2.3_1", "FEZF2.3_0,5,6,7", "FEZF2.3_3", "FEZF2.3_8", "FEZF2.3_4", "FEZF2.3_2", 
              "FEZF2.5_1", "FEZF2.4_3", "FEZF2.5_3", "FEZF2.5_0,4,5,11", "FEZF2.5_2,9,10", "FEZF2.5_8", 
              "FEZF2.5_6,7", "FEZF2.5_12", "FEZF2.4_7", "FEZF2.4_6", "FEZF2.4_8", "FEZF2.4_1,2,5", 
              "FEZF2.4_4", "FEZF2.4_0,9", "RORB_15", "RORB_13", "LINC00507_14", "RORB_4,8,14,18", 
              "RORB_12", "RORB_23", "RORB_1,17", "RORB_10_19", "RORB_5", "RORB_22", 
              "RORB_21", "RORB_3,6,7", "RORB_2,9,20", "RORB_16", "RORB_11", "RORB_0", 
              "LINC00507_4,7,12", "LINC00507_10", "LINC00507_5,15", "LINC00507_1,3,11", "LINC00507_17", "LINC00507_2", 
              "LINC00507_0,6,16", "LINC00507_8,13,19", "LINC00507_9", "LINC00507_18", "FEZF2.2_2", "FEZF2.2_0,3", 
              "FEZF2.2_1", "THEMIS_8", "THEMIS_6", "THEMIS_4", "THEMIS_0", "THEMIS_10", 
              "THEMIS_3", "THEMIS_2", "THEMIS_9", "THEMIS_7", "THEMIS_1,5", "146_157", 
              "144_156", "148", "160_162", "127", "206", "204_211", 
              "226", "192", "199", "198", "202", "201", 
              "185", "180", "175", "173", "172", "190", 
              "233", "232", "247", "229", "279_280", "222", 
              "282", "264", "245", "235_243", "256", "270", 
              "224", "220", "217", "214", "336", "268", 
              "383", "378", "375", "377", "311", "303", 
              "309", "329_333", "301", "327", "299", "384", 
              "267", "339", "373", "334", "368", "342", 
              "380", "284", "297", "413", "404", "403", 
              "401", "469", "466", "461", "471", "463", 
              "456", "454", "453", "446")

seu$hicat_cluster_anno <- as.vector(seu$hicat_cluster_merge)
seu$hicat_cluster_anno[seu$hicat_cluster_merge %in% c("22", "19", "17", "1")] <- 'OPC'
seu$hicat_cluster_anno[seu$hicat_cluster_merge %in% c("21")] <- 'Oligo L4-L6 MOBP'
seu$hicat_cluster_anno[seu$hicat_cluster_merge %in% c("395", "392", "389", "386")] <- 'Oligo L4-L6 OPALIN'
seu$hicat_cluster_anno[seu$hicat_cluster_merge %in% c("477", "476", "478")] <- 'Astrocyte FGFR3'
seu$hicat_cluster_anno[seu$hicat_cluster_merge %in% c("FEZF2.3_1", "FEZF2.3_0,5,6,7", "FEZF2.3_3", "FEZF2.3_8", "FEZF2.3_4", "FEZF2.3_2")] <- 'NP'
seu$hicat_cluster_anno[seu$hicat_cluster_merge %in% c("FEZF2.5_1", "FEZF2.4_3")] <- 'X'
seu$hicat_cluster_anno[seu$hicat_cluster_merge %in% c("FEZF2.5_3", "FEZF2.5_0,4,5,11", "FEZF2.5_2,9,10", 
                                                      "FEZF2.5_8", "FEZF2.5_6,7", "FEZF2.5_12")] <- 'L6CT/FEZF2'
seu$hicat_cluster_anno[seu$hicat_cluster_merge %in% c("FEZF2.4_7", "FEZF2.4_6", "FEZF2.4_8", 
                                                      "FEZF2.4_1,2,5", "FEZF2.4_4", "FEZF2.4_0,9")] <- 'L6B/FEZF2'
seu$hicat_cluster_anno[seu$hicat_cluster_merge %in% c("RORB_15", "RORB_13", "LINC00507_14", "RORB_4,8,14,18", 
                                                      "RORB_12", "RORB_23", "RORB_1,17", "RORB_10_19", "RORB_5", "RORB_22", 
                                                      "RORB_21", "RORB_3,6,7", "RORB_2,9,20", "RORB_16", "RORB_11", "RORB_0")] <- 'L5 IT/RORB'
seu$hicat_cluster_anno[seu$hicat_cluster_merge %in% c("LINC00507_4,7,12", "LINC00507_10", "LINC00507_5,15", 
                                                      "LINC00507_1,3,11", "LINC00507_17", "LINC00507_2", 
                                                      "LINC00507_0,6,16", "LINC00507_8,13,19", "LINC00507_9", 
                                                      "LINC00507_18")] <- 'L1-L3 IT LINC00507'
seu$hicat_cluster_anno[seu$hicat_cluster_merge %in% c("FEZF2.2_2", "FEZF2.2_0,3", "FEZF2.2_1")] <- 'ET'
seu$hicat_cluster_anno[seu$hicat_cluster_merge %in% c("THEMIS_8", "THEMIS_6")] <- 'L6 CAR3'
seu$hicat_cluster_anno[seu$hicat_cluster_merge %in% c("THEMIS_4", "THEMIS_0", "THEMIS_10", 
                                                      "THEMIS_3", "THEMIS_2", "THEMIS_9", 
                                                      "THEMIS_7", "THEMIS_1,5")] <- 'L6/IT'
seu$hicat_cluster_anno[seu$hicat_cluster_merge %in% c("146_157", "144_156", "148", "160_162", "127")] <- 'LAMP5/SV2C'
seu$hicat_cluster_anno[seu$hicat_cluster_merge %in% c("206", "204_211", "226", "192", 
                                                      "199", "198", "202", "201")] <- 'PAX6'
seu$hicat_cluster_anno[seu$hicat_cluster_merge %in% c("185", "180", "175", "173", "172", "190")] <- 'LAMP5/SNCG'
seu$hicat_cluster_anno[seu$hicat_cluster_merge %in% c("233", "232", "247", "229", "279_280", "222", 
                                                      "282", "264", "245", "235_243", "256", "270", 
                                                      "224", "220", "217", "214")] <- 'VIP'
seu$hicat_cluster_anno[seu$hicat_cluster_merge %in% c("336", "268", 
                                                      "383", "378", "375", "377", "311", "303", 
                                                      "309", "329_333", "301", "327", "299")] <- 'SST'
seu$hicat_cluster_anno[seu$hicat_cluster_merge %in% c("384", "267", "339", "373", "334", "368", "342", 
                                                      "380", "284", "297")] <- 'PVALB'
seu$hicat_cluster_anno[seu$hicat_cluster_merge %in% c("413")] <- 'OPC2'
seu$hicat_cluster_anno[seu$hicat_cluster_merge %in% c("404", "403", "401")] <- 'Microglia'
seu$hicat_cluster_anno[seu$hicat_cluster_merge %in% c("469", "466", "461", "471", "463")] <- 'VLMC L1-L3'
seu$hicat_cluster_anno[seu$hicat_cluster_merge %in% c("456", "454", "453", "446")] <- 'ENDO L2-L5'


DimPlot(seu,group.by = 'hicat_cluster_anno', label=T)

#saveRDS(seu, 'seu.harmony.anno.v2.rds')



##################################################################################
##################################################################################

######## classes #######################
seu$hicat_cluster_classes <- as.vector(seu$hicat_cluster_merge)
seu$hicat_cluster_classes[seu$hicat_cluster_merge %in% c("22", "19", "17", "1", "21", "395", 
                                                         "392", "389", "386", "477", "476", "478",
                                                         "413")] <- 'Astrocyte/oligodendrocyte non-neuronal class'
seu$hicat_cluster_classes[seu$hicat_cluster_merge %in% c("FEZF2.3_1", "FEZF2.3_0,5,6,7", "FEZF2.3_3", "FEZF2.3_8", "FEZF2.3_4", "FEZF2.3_2", 
                                                         "FEZF2.5_1", "FEZF2.4_3", "FEZF2.5_3", "FEZF2.5_0,4,5,11", "FEZF2.5_2,9,10", "FEZF2.5_8", 
                                                         "FEZF2.5_6,7", "FEZF2.5_12", "FEZF2.4_7", "FEZF2.4_6", "FEZF2.4_8", "FEZF2.4_1,2,5", 
                                                         "FEZF2.4_4", "FEZF2.4_0,9", "RORB_15", "RORB_13", "LINC00507_14", "RORB_4,8,14,18", 
                                                         "RORB_12", "RORB_23", "RORB_1,17", "RORB_10_19", "RORB_5", "RORB_22", 
                                                         "RORB_21", "RORB_3,6,7", "RORB_2,9,20", "RORB_16", "RORB_11", "RORB_0", 
                                                         "LINC00507_4,7,12", "LINC00507_10", "LINC00507_5,15", "LINC00507_1,3,11", "LINC00507_17", "LINC00507_2", 
                                                         "LINC00507_0,6,16", "LINC00507_8,13,19", "LINC00507_9", "LINC00507_18", "FEZF2.2_2", "FEZF2.2_0,3", 
                                                         "FEZF2.2_1", "THEMIS_8", "THEMIS_6", "THEMIS_4", "THEMIS_0", "THEMIS_10", 
                                                         "THEMIS_3", "THEMIS_2", "THEMIS_9", "THEMIS_7", "THEMIS_1,5")] <- 'Glutamatergic neuronal class'
seu$hicat_cluster_classes[seu$hicat_cluster_merge %in% c("146_157", "144_156", "148", "160_162", "127", "206", "204_211", 
                                                         "226", "192", "199", "198", "202", "201", 
                                                         "185", "180", "175", "173", "172", "190", 
                                                         "233", "232", "247", "229", "279_280", "222", 
                                                         "282", "264", "245", "235_243", "256", "270", 
                                                         "224", "220", "217", "214", "336", "268", 
                                                         "383", "378", "375", "377", "311", "303", 
                                                         "309", "329_333", "301", "327", "299", "384", 
                                                         "267", "339", "373", "334", "368", "342", 
                                                         "380", "284", "297")] <- 'GABAergic neuronal class'
seu$hicat_cluster_classes[seu$hicat_cluster_merge %in% c("404", "403", "401", "469", "466", "461", "471", "463", 
                                                         "456", "454", "453", "446")] <- 'Immune/vascular non-neuronal class'




######## subclasses #######################
seu$hicat_cluster_subclasses <- seu$hicat_cluster_anno




######## supertypes #######################
seu$hicat_cluster_supertypes <- as.vector(seu$hicat_cluster_merge)

a <- table(seu$hicat_cluster_merge)
b <- a[rev(order(a))]

oldID <- names(b)
newID <- c(1:142)
names(newID) <- oldID

seu$hicat_cluster_newID <- newID[as.character(seu$hicat_cluster_merge)]

#saveRDS(seu,'seu.harmony.anno.v2.rds')
saveRDS(seu@meta.data,'seu.harmony.anno.meta.v2.rds')

seu0 <- CreateSeuratObject(
  seu@assays$RNA@counts,
  min.cells = 0,
  min.features = 0
)
seu0 <- NormalizeData(object = seu0, normalization.method = "LogNormalize", scale.factor = 1e4)
saveRDS(seu0,'seu.v2.rds')


saveRDS(seu,'seu.harmony.anno.v2.rds')
