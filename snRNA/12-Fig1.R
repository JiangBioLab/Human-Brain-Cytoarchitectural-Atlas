library(scrattch.hicat)
library(dendextend)
library(Seurat)

#seu <- readRDS('seu.harmony.anno.de80.rds')
seu <- readRDS('seu.harmony.anno.v2.rds')

dend.result <- readRDS('Brain_subClusteringID_seurat_merge_dend_renames_newID.rds')
dend <- dend.result$dend
dend_sort <- readRDS('Brain_subClusteringID_seurat_merge_dend_renames_newID_sort.rds') 

ID <- readRDS('ID.rds')

seu$hicat_cluster_merge_newID <- factor(seu$hicat_cluster_merge_newID ,
                                        levels = rev(labels(dend_sort)))

################# Color ################################
supertype_color <- c("#382650","#2F2144","#513577","#4D4639","#6E7C6F","#535944","#697255","#3B3F32","#A77C70","#A6906F",      ##Non-Neuron
                     "#8D6041","#604B47","#2E3E39","#4D483C","#332F26","#5A503F","#72B09C","#6CA491","#47645C","#4F7165",
                     "#65978A","#3F574E","#5E8A79","#577D70","#374A45",
                     
                     "#C36AE6","#C263CC","#D97C80","#FF8995","#B36C76",                                                        ##CGE
                     "#F2859E","#F48C9E","#BE61D4","#9E56A6","#B65FBF","#B363CB","#BD76DC","#B967D9","#DD6DF2","#A45FBF",
                     "#9A5AB3","#AF6FCC","#C667D9","#D26AE6","#AC63BE","#9256A6","#7F4C8C","#DF70FF","#D270FF","#CD26FF",
                     "#CD6DF2","#DA81F7","#EA86FF","#E96DF2","#CC7E09","#B98327","#C38326","#B67B16","#A66D0D","#B57014",
                     
                     "#F27503","#E67105","#EC820D","#E69A05","#D99207","#FF7A00","#CE8024","#BF670B","#C36F1E","#DC7C16",     ##MGE
                     "#CC7209","#E89420","#FFA900","#AB5B42","#D0494E","#B1614A","#D92D43","#EC2D46","#8C2A27","#873C46",
                     "#B9342C","#D93031","#B92C3E",
                     
                     "#93C43B","#AFE32A","#B1E67B","#C0F27F","#A4D374","#88A662","#8D8C20",                                   ##Glu
                     "#898325","#A7A322","#B0CE1F","#B1B10C","#DABE23","#CEC823","#C5A93D","#91910E","#C8B323","#0E8B8B",
                     "#09EEDC","#19E3BE","#2B8AA5","#7944AA","#2FBCE5","#2E84BE","#2E9FCB","#296B98","#2C93B2","#2B77A5",
                     "#69419D","#7E47B6","#5E419D","#71419D","#7643A4","#6444AA","#0BBEBE","#52ADC0","#4A9E9E","#18C5C8",
                     "#02E2F9","#09CCC6","#07C6D9","#2E9FBB","#378695","#138C98","#0DA3A3","#0E9797","#07D8D8","#5BDFDC",
                     "#637384","#5B7893","#61867A","#5991A4","#4EA8AC","#41AEAE","#4E9F8C","#99FFCC","#52B8AA","#4E9E9E",
                     "#55C5B5","#4BBCAD")

names(supertype_color) <- c("84", "82", "6",  "129","117","79", "35", "125","59", "96",    ##Non-Neuron
                            "31", "103","70", "69", "27", "2",  "42", "127","1",  "56", 
                            "95", "57", "89", "66", "4",  
                            
                            "86", "54", "33", "32", "109",                                 ##CGE
                            "83", "19", "77", "121","46", "15", "68", "52", "78", "29", 
                            "132","92", "51", "49", "44", "20", "99", "94", "87", "101",
                            "75", "73", "122","120","90", "130","40", "126","45", "118",
                            
                            "67", "100","97", "102","25", "24", "110","55", "112","58",     ##MGE
                            "21", "28", "37", "135","141","23", "64", "105","36", "50", 
                            "88", "113","8",  
                            
                            "111","119","140","116","104","60", "136",                      ##Glu
                            "128","139","123","85", "134","76", "138","93", "91", "131",
                            "133","98", "80", "107","142","72", "114","63", "81", "53", 
                            "124","115","137","61", "106","74", "22", "47", "43", "71", 
                            "16", "11", "10", "5",  "38", "3",  "9",  "7",  "17", "26", 
                            "34", "18", "30", "108","13", "41", "48", "39", "65", "62", 
                            "12", "14")

subclass_color <- c("#2E3E39", "#5E8A79", "#72B09C", "#665C47", "#93C43B", "#2D8CB8", "#7944AA", "#09B2B2", "#07D8D8", "#19E3BE", 
                    "#898325", "#CEC823", "#f2798d", "#E96DF2", "#D96C07", "#9C28CC", "#F2B003", "#FF2D4E", "#513577", "#697255", 
                    "#604B47")
names(subclass_color) <- c("OPC", "Oligo L4-L6 MOBP", "Oligo L4-L6 OPALIN", "Astrocyte FGFR3", "NP", "L6CT/FEZF2", "L6B/FEZF2", "L5 IT/RORB", "L1-L3 IT LINC00507", "ET", 
                           "L6 CAR3", "L6/IT", "LAMP5/SV2C", "PAX6", "LAMP5/SNCG", "VIP", "SST", "PVALB", "Microglia", "VLMC L1-L3", 
                           "ENDO L2-L5")
class_color <- c('#00ADEE', '#F05A28', '#808080', '#808080')
names(class_color) <- c('Glutamatergic neuronal class', 'GABAergic neuronal class', 'Astrocyte/oligodendrocyte non-neuronal class', 'Immune/vascular non-neuronal class')

########################################################



################### dend #####################################
dend.result <- readRDS('Brain_subClusteringID_seurat_merge_dend_renames_newID.rds')
dend <- dend.result$dend

#dend.pdf 30x4
plot(dend, horiz = T)
plot(remove_branches_edgePar(dend))


library(dendsort)
#dend_sort.pdf 30x4
plot(dendsort(dend),
     horiz = T)

saveRDS(dendsort(dend), 'Brain_subClusteringID_seurat_merge_dend_renames_newID_sort.rds')


str(dend)
write.table(dend,'dend.txt')

dend.reorder <- reorder(dend,10:1)
plot(dend.reorder,horiz = T)

dend <- dend.result$dend
anno <- unique(cbind(as.vector(seu$hicat_cluster_newID),seu$hicat_cluster_supertypes))
rownames(anno) <- anno[,1]
anno[,2] <- paste0(anno[,1],':',anno[,2])

dend.label <- dend
labels(dend.label) <- anno[labels(dend),2]

plot(dend,horiz = T)
plot(dend.label,horiz = T)

plot(remove_branches_edgePar(dend.reorder))





############## Chi-Test ########################
a <- table(seu$hicat_cluster_merge_newID,seu$region)

chisq <- chisq.test(a)
chisq
#X-squared = 1147.2, df = 364, p-value < 2.2e-16
#
chisq$observed
round(chisq$expected,2)
round(chisq$residuals, 3)

library(corrplot)
library(grDevices)
col2 <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                               "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                               "#4393C3", "#2166AC", "#053061")))
#cluster_region_chiTest_corrplot.pdf 30x5
corrplot(chisq$residuals, is.cor = FALSE, col =  col2(200))



#UMAP_cluster.pdf 10x18
DimPlot(seu, group.by = 'hicat_cluster_merge_newID', cols = rev(supertype_color[levels(seu$hicat_cluster_merge_newID)]))

seu$hicat_cluster_classes <- factor(seu$hicat_cluster_classes,
                                    levels = c("Glutamatergic neuronal class",
                                               "GABAergic neuronal class",
                                               "Astrocyte/oligodendrocyte non-neuronal class",
                                               "Immune/vascular non-neuronal class"))

#UMAP_classes.pdf 5x10
DimPlot(seu, group.by = 'hicat_cluster_classes', cols = c('#00ADEE','#F05A28','#808080','#505050'))


seu$individual <- factor(seu$individual,
                         levels=c('S0206','S0406','S0426','S0531','S0916'))
#UMAP_individual.pdf 5x6
DimPlot(seu, group.by = 'individual', cols = c('#1B9E77','#D95F02','#7570B3','#E7298A','#66A61E'))


regionID <- c('A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','A11','A12','A13','A14')
regionName <- c('FPPFC','DLPFC','VLPFC','M1','S1','S1E','PoCG','SPL','SMG','AG','VISP','ITG','STG','ACC')
names(regionName) <- regionID

seu$regionName <- regionName[seu$region]

seu$regionName <- factor(seu$regionName,
                         levels=c('FPPFC','DLPFC','VLPFC','M1','S1','S1E','PoCG','SPL','SMG','AG','VISP','ITG','STG','ACC'))
#UMAP_region.pdf 5x6
DimPlot(seu, group.by = 'regionName', cols = c('#3F4587','#8562AA','#EC8561','#B97CB5','#D43046','#F0592B',
                                               '#ED4A96','#593C97','#A54486','#FBDE13','#299FAE','#75CCE3',
                                               '#0C6939','#0D9547'))



DotPlot(seu,features = c('FEZF2', 'SGCG', 'SYT6', 'TSPAN18'),group.by = 'hicat_cluster_newID')





############## cell composition ###################

cluster_dis<-cbind(as.vector(seu$hicat_cluster_merge_newID),as.vector(seu$individual))
colnames(cluster_dis)<-c('cluster','individual')
cluster_dis<-as.data.frame(cluster_dis)
cluster_dis$cluster<-factor(cluster_dis$cluster,
                            levels=rev(labels(dend_sort)))

library(ggplot2)
library(RColorBrewer)
##cluster_individual_composition_de80.pdf
##pdf 3x30
ggplot(cluster_dis, aes(cluster)) + 
  geom_bar(aes(fill=individual), position="fill", color="gray20",alpha = 0.7) +#
  scale_fill_manual(values=brewer.pal(6,"Dark2"))+ #values=colors[c(1:4,6)]
  ylab('Percentage')+
  theme(axis.line = element_line(colour = "gray10"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + RotatedAxis()

cluster_dis<-cbind(as.vector(seu$hicat_cluster_merge_newID),as.vector(seu$regionName))
colnames(cluster_dis)<-c('cluster','region')
cluster_dis<-as.data.frame(cluster_dis)
cluster_dis$cluster<-factor(cluster_dis$cluster,
                            levels=rev(labels(dend_sort)))
cluster_dis$region<-factor(cluster_dis$region,
                           levels=c('FPPFC','DLPFC','VLPFC','M1','S1','S1E','PoCG','SPL','SMG','AG','VISP','ITG','STG','ACC'))
library(ggplot2)
library(RColorBrewer)
##cluster_region_composition_de80.pdf
##pdf 3x30
ggplot(cluster_dis, aes(cluster)) + 
  geom_bar(aes(fill=region), position="fill", color="gray20") +#,alpha = 0.7
  scale_fill_manual(values=c(brewer.pal(12,"Set3"),brewer.pal(2,"Set2")))+ #values=colors[c(1:4,6)]
  ylab('Percentage')+
  theme(axis.line = element_line(colour = "gray10"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + RotatedAxis()


#### % percent in each sample ##############
cluster_dis<-table(seu$hicat_cluster_merge_newID,seu$regionName)
cluster_dis_r <- cluster_dis/rowSums(cluster_dis)
cluster_dis_r <- cluster_dis_r[rev(labels(dend_sort)),]

cluster_dis_r[cluster_dis_r>0.3]<-0.3
myColors=brewer.pal(8,"Reds")[1:8]

cluster_dis_r <- cluster_dis_r[,c('FPPFC','DLPFC','VLPFC','M1','S1','S1E','PoCG','SPL','SMG','AG','VISP','ITG','STG','ACC')]

library(pheatmap)
library(RColorBrewer)
#pdf cluster_region_composition_ratio_de80.pdf 30x4
#pdf cluster_region_composition_ratio_tree.pdf 20x3.5
pheatmap(cluster_dis_r,
         cluster_rows = F,
         cluster_cols = F, # T for Tree
         color = colorRampPalette(colors = myColors)(100))


cluster_dis_r_samples <- cluster_dis_r
for(i in 1:dim(cluster_dis_r)[1]){
  for(j in 1:dim(cluster_dis_r)[2]){
    cluster_dis_r_samples[i,j]<-length(unique(seu$individual[seu$hicat_cluster_merge_newID==rownames(cluster_dis_r)[i] & seu$regionName==colnames(cluster_dis_r)[j]]))
  }
}
cluster_dis_r_samples <- cluster_dis_r_samples[,c('FPPFC','DLPFC','VLPFC','M1','S1','S1E','PoCG','SPL','SMG','AG','VISP','ITG','STG','ACC')]

#pdf cluster_region_composition_samples_de80.pdf 30x5
pheatmap(cluster_dis_r_samples,
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = T,
         number_format = '%.f')


cluster_dis<-cluster_dis[rev(labels(dend_sort)),
                         c('FPPFC','DLPFC','VLPFC','M1','S1','S1E','PoCG','SPL','SMG','AG','VISP','ITG','STG','ACC')]
#pdf cluster_region_composition_cells_de80.pdf 30x5
pheatmap(log2(cluster_dis+1),
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = T,
         number_format = '%.f')


# cell_number_log2_de80.pdf 3x30
barplot(log2(table(seu$hicat_cluster_merge_newID)[rev(labels(dend_sort))]),
        col = supertype_color[rev(labels(dend_sort))])

ID <- readRDS('ID.rds')
rownames(ID) <- ID$hicat_cluster_merge_newID

#cluster_color.pdf 3x20
barplot(rep(1,142),
        col = supertype_color[rev(labels(dend_sort))],
        xlab = ID[rev(labels(dend_sort)),'ID'])


################# pheatmap #############################
genes1 <- c('PTRPC','NOSTRIN','PDGRFB','MBP','PDGFRA','FEZF2','LINC00507','COL5A2','LAMP5','SV2C','PAX6','VIP','PVALB','SST')

genes2 <- c('TYROBP','NOSTRIN','OPALIN','FGFR3','PDGFRA','SLC1A3','FEZF2','THEMIS','RORB','LINC00507','SLC17A7','PVALB','SST','LHX6','VIP','LAMP5','PAX6','ADARB2','GAD1')

genes1[!genes1 %in% genes2]
genes2[!genes2 %in% genes1]

genes <- c('PTPRC','TYROBP','NOSTRIN','PDGFRB','OPALIN','FGFR3','PDGFRA','SLC1A3','FEZF2','THEMIS','LINC00507','SLC17A7','RORB','LAMP5','SV2C','CTGF','SNCG','LINC01116','CCN2','PAX6','VIP','ADARB2','PVALB','SST','LHX6','GAD1','ADARB2','AQP4')
seu_markers_m<-seu@assays$RNA@data[rownames(seu@assays$RNA@data) %in% genes,]
seu_markers_m<-t(as.matrix(seu_markers_m))
seu_markers_m<-as.data.frame(seu_markers_m)
seu_markers_m<-cbind(seu_markers_m,as.vector(seu$hicat_cluster_newID))
colnames(seu_markers_m)[dim(seu_markers_m)[2]]<-'cluster'
#seu_markers_m$cluster<-as.vector(seu_markers_m$cluster)
seu_markers_m2<-aggregate(seu_markers_m[,1:(dim(seu_markers_m)[2]-1)],by=list(cluster=seu_markers_m$cluster),FUN=mean)
seu_markers_m3<-seu_markers_m2[,-1]

seu_markers_m4<-matrix(unlist(seu_markers_m3),dim(seu_markers_m3)[1],dim(seu_markers_m3)[2])
colnames(seu_markers_m4)<-colnames(seu_markers_m3)
rownames(seu_markers_m4)<-as.vector(seu_markers_m2[,1])
seu_markers_m4<-log2(seu_markers_m4+1)
seu_markers_m4 <- seu_markers_m4[rev(c('9','94','64','62','5','78','55','48','47','40','89','71','2',
                                       '52','73','22','82','45','25','56','70','3','68','54','4','44',
                                       '46','92','67','30','35','23','36','65','61','12','1','18','6',
                                       '21','7','11','63','19','97','16','80','96','8','76','49','50',
                                       '59','81','84','42','43','77','28','13','38','91','87','20','17',
                                       '51','27','39','24','79','15','95','88','58','53','37','34','83',
                                       '31','33','14','93','10','75','41','66','26','90','32','72','29',
                                       '69','74','86','85','60','57')),]
seu_markers_m4 <- seu_markers_m4[,genes]

anno <- unique(cbind(as.vector(seu$hicat_cluster_newID),seu$hicat_cluster_supertypes))
rownames(anno)<-anno[,1]

#rownames(seu_markers_m4) <- paste0(rownames(seu_markers_m4),':',anno[rownames(seu_markers_m4),2])
rownames(seu_markers_m4) <- anno[rownames(seu_markers_m4),2]

##marker_heatmap.pdf 19.4x7
pheatmap(seu_markers_m4,
         cluster_rows=F,
         cluster_cols=F,
         angle_col = 45)

seu_markers_m4_z <- t(seu_markers_m4)
for(i in 1:dim(seu_markers_m4_z)[1]){
  tmp <- scale(seu_markers_m4_z[i,],center=T,scale=T)
  seu_markers_m4_z[i,]<-tmp[,1]
}

#seu_markers_m4_z2 <- t(seu_markers_m4)
#for(i in 1:dim(seu_markers_m4_z2)[1]){
#  zscore = (seu_markers_m4_z2[i,] - mean(seu_markers_m4_z2[i,]))/sd(seu_markers_m4_z2[i,])
#  seu_markers_m4_z2[i,]<-zscore
#}

seu_markers_m4_z[seu_markers_m4_z>5] = 5

#marker_zscore_heatmap.pdf 22x7
pheatmap(t(seu_markers_m4_z),
         cluster_cols = F,
         cluster_rows = F,
         angle_col = c("45"))



######### DotPlot ############################
genes <- c('PTPRC','TYROBP','NOSTRIN','PDGFRB','OPALIN','FGFR3','PDGFRA','SLC1A3','FEZF2','THEMIS','LINC00507','SLC17A7','RORB','LAMP5','SV2C','PAX6','VIP','ADARB2','PVALB','SST','LHX6','GAD1')

genes <- c('PTPRC','TYROBP','NOSTRIN','PDGFRB','OPALIN','FGFR3','PDGFRA','SLC1A3','FEZF2','THEMIS','LINC00507','SLC17A7','RORB','LAMP5','SV2C','CTGF','SNCG','LINC01116','CCN2','PAX6','VIP','ADARB2','PVALB','SST','LHX6','GAD1','ADARB2','AQP4')

seu$hicat_cluster_merge_newID <- factor(seu$hicat_cluster_merge_newID,
                            levels=labels(dend_sort))
#seu@active.ident <- seu$hicat_cluster
##marker_dotplot_de80.pdf 30x12
pdf('marker_dotplot_de80.pdf',width = 12, height = 30)
DotPlot(seu,features = unique(genes), group.by = 'hicat_cluster_merge_newID') + RotatedAxis()
dev.off()




genes <- c('PTPRC','TYROBP','NOSTRIN','PDGFRB','OPALIN','FGFR3','PDGFRA','SLC1A3','FEZF2',
           'THEMIS','LINC00507','SLC17A7','RORB','LAMP5','SV2C','CTGF','SNCG','LINC01116',
           'CCN2','PAX6','VIP','ADARB2','PVALB','SST','LHX6','GAD1',
           'FRMD7','ADAMTS19','SCGN','TH','LHX1OS','LHX5','LHX1','EBF3')



##################################################
#chi-test in R
#http://www.sthda.com/english/wiki/chi-square-test-of-independence-in-r

a <- table(seu$hicat_cluster_merge_newID,seu$regionName)

chisq <- chisq.test(a)
chisq
#X-squared = 26865, df = 1248, p-value < 2.2e-16
#
chisq$observed
round(chisq$expected,2)
round(chisq$residuals, 3)


library(corrplot)
library(grDevices)
col2 <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                               "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                               "#4393C3", "#2166AC", "#053061")))
#cluster_region_chiTest_corrplot.pdf 20x5
corrplot(chisq$residuals, is.cor = FALSE, col =  col2(200))








############## check some clusters #####################################
### 384 339
DimPlot(seu,group.by = 'hicat_cluster_anno',label=T)

DimPlot(seu, 
        group.by="hicat_cluster_anno", 
        cells.highlight=findCells(seu, 'hicat_cluster_merge', '384'))

DimPlot(seu, cells = unlist(findCells(seu, 'hicat_cluster_merge', '384')))

FeaturePlot(seu, 
            cells = unlist(findCells(seu, 'hicat_cluster_merge', '384')),
            reduction = "umap", 
            features = c('OPALIN','SLC1A3','PVALB','LHX6'),
            ncol = 2,
            cols = c('grey','tomato','tomato4'))


findCells <- function(obj, column, values, name=NULL) {
  ## Given a Seurat OBJ, return a list with the names of the cells where
  ## the specified meta.data COLUMN equals any of the strings specified
  ## in VALUES (both must be characters or factors). Name of the list member
  ## must be specified using the NAME argument if length(values)>1
  stopifnot(is(obj, "Seurat"))
  stopifnot(is.character(column))
  stopifnot(column %in% names(obj@meta.data))
  col <- obj@meta.data[[column]]
  stopifnot(is.character(col) || is.factor(col))
  values <- unique(values)
  stopifnot(is.character(values) || is.factor(values))
  if (length(values)>1 && is.null(name))
    stop("findCells: specify a name to be used for the selection")
  if(is.null(name))
    name <- values
  stopifnot(is.character(name))
  rem <- setdiff(c(values), col)
  if(length(rem)>0)stop("findCells: requested value(s) never occurs in this column: ", rem)
  l <- list(colnames(obj)[ col %in% values ])
  names(l) <- name
  l
}                                       #findCells



########################################################################



