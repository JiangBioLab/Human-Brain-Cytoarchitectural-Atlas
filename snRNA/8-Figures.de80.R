library(Seurat)

seu <- readRDS('seu.harmony.anno.rds')


##### color #########################
individual_color <- c('S0206','S0406','S0426','S0531')
names(individual_color) <- c('#1B9E77','#D95F02','#7570B3','#E7298A')

region_color <- c('#3F4587','#8562AA','#EC8561','#B97CB5','#D43046','#F0592B','#ED4A96','#593C97','#A54486','#FBDE13','#299FAE','#75CCE3','#0C6939','#0D9547')
names(region_color) <- c('A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','A11','A12','A13','A14')

class_color <- c('#00ADEE','#F05A28','#808080','#101010')
names(class_color) <- c('Glutamatergic','GABAergic','non-neuronal','Unassigned')

subclass_color <- c('#808080','#2EBF5E','#58D2CF','#0EED60','#8D6C62','#94AF97','#53776C','#697255','#807059','#DA808C','#B1B10C','#D93137','#FF9900','#B864CC')
names(subclass_color) <- c('Unassigned','COL5A2','FEZF2','LINC00507','Endo','Mig','Olig','OPC','Peri','LAMP5/SV2C','PAX6','PVALB','SST','VIP')

supertype_color <- c('#86A7B9','#8BA38D','#6E7C6F','#535944','#3B3F32','#87A68B','#697255','#374A45','#6CA491','#577D70',
                     '#3F574E','#5E8A79','#2E3E39','#3B3F32','#A77C70','#A6906F','#535944','#604B47','#697255','#768BA6',
                     '#72B09C','#7BA6C4','#47645C','#4F7165','#65978A','#299373','#2FD3B0','#2CAD7A','#2FCA96','#8AA3CC',
                     '#90A7D9','#30E6BA','#2BA087','#268063','#9EE61D','#BCFF1A','#8DD31E','#02F970','#0E993B','#0BBF45',
                     '#13A23E','#00FF34','#07D945','#09CC4F','#F2B003','#D62228','#A62A2F','#EC2D46','#8C2734','#B92C3E',
                     '#B9342C','#D93031','#EC820D','#F28403','#D96C07','#FFB800','#CC6709','#E67105','#F2841A','#E89420',
                     '#FFA900','#FF7A00','#F27503','#D97807','#FF9900','#EC3031','#FF2D4E','#C60C0F','#F2859E','#D97C80',
                     '#F48C9E','#FF8995','#A45FBF','#BD76DC','#9E56A6','#AF6FCC','#B967D9','#BE61D4','#B65FBF','#DD6DF2',
                     '#9A5AB3','#D270FF','#D26AE6','#D92D43','#E69A05','#B98327','#C39317','#B363CB','#A6770D','#C38326',
                     '#996A0E','#AC7913','#A66D0D','#8C620E','#D99207','#B57014','#B67B16')
names(supertype_color) <- c('9','94','64','62','5','78','55','48','47','40',
                            '89','71','2','52','73','22','82','45','25','56',
                            '70','3','68','54','4','44','46','92','67','30',
                            '35','23','36','65','61','12','1','18','6','21',
                            '7','11','63','19','97','16','80','96','8','76',
                            '49','50','59','81','84','42','43','77','28','13',
                            '38','91','87','20','17','51','27','39','24','79',
                            '15','95','88','58','53','37','34','83','31','33',
                            '14','93','10','75','41','66','26','90','32','72',
                            '29','69','74','86','85','60','57')

#####################################



seu$hicat_cluster <- factor(seu$hicat_cluster,
                            levels = c('380','378','382','359','353','372','333','337','342','341','339','343','392','396','259','400',
                                       '386','43','76','78','38','41','42','35','25','30','32','87','65','73','63','59','61','44','45','69',
                                       '72','56','57','58','110','114','121','117','119','93','183','177','180','186','191','192',
                                       '203','204','200','201','139','194','129','131','162','174','126','127','133','137','135',
                                       '167','141','318','285','321','309','216','297','277','245','278','251','276','293','263',
                                       '269','261','266','212','248','249','222','224','306','295','271','211','229','322','329',
                                       '331','327','323','325','348','20','1','16','21','22'))
seu@active.ident <- seu$hicat_cluster


DimPlot(seu,group.by = 'clusterAnno', label=T)
DimPlot(seu, cells.highlight = colnames(seu)[seu$hicat_cluster=='127'])

DimPlot(seu, label = T)
DimPlot(seu,label = T,reduction='tsne')
DimPlot(seu, group.by = 'individual')
DimPlot(seu, group.by = 'region')

neuron_marker <- c('SNAP25','SYT1','GRIN1', #neuron
                   'GAD1','GAD2','SLC6A1','SLC32A1', #inhibitory neuron
                   'SLC17A7','STAB2','NRGN','CAMK2A', #excitatory neuron
                   'RYR1','RASGRF2','CA8','PVALB','GRIA1', #purkinje neuron
                   'TIAM1','PPFIA2','HTR2C','TLE4','TSHZ2','FOXP2','GRM4','PCP4', #granule cell
                   'SST', #SST+
                   'VIP', #VIP+
                   'SV2C', #SV2C+
                   'PVALB', #PVALB+
                   'SLC17A6', #Glutamatergic neuron
                   'CUX2', #upper-layer
                   'RORB', #layer 4
                   'TLE4') #deep-layer

neuron_marker <- c('RYR1','RASGRF2','CA8','PVALB','GRIA1', #purkinje neuron
                   'TIAM1','PPFIA2','HTR2C','TLE4','TSHZ2','FOXP2','GRM4','PCP4', #granule cell
                   'SST', #SST+
                   'VIP', #VIP+
                   'SV2C', #SV2C+
                   'PVALB', #PVALB+
                   'SLC17A6', #Glutamatergic neuron
                   'CUX2', #upper-layer
                   'RORB', #layer 4
                   'TLE4')

#neuron_marker_dotplot.pdf 15x12
DotPlot(seu, features = unique(neuron_marker)) + RotatedAxis()

astrocytes_marker <- c('AQP4','SLC1A2','GJA1','GFAP','SLC1A3')
ODC_marker <- c('MBP','SOX10', 'PLP1', 'MOBP','MOG','SLC1A3','SLC1A3')
OPC_marker <- c('PDGFRA','PCDH15')
microglia_marker <- c('HLA-DRA', 'CX3CR1', 'CSF1R','SLC1A3','PTPRC', 'TMEM119', 'P2RY12')
endothelial_marker <- c('FLT1', 'CLDN5','SLC1A3','TEK', 'PDGFB', 'NOS3', 'ELTD1', 'PECAM1')
pericytes_marker <- c('AMBP','PDGFRB','COBLL1','KCNJ8', 'ABCC9', 'ATP13A5', 'ART3', 'PLA1A', 'ACE2')

markers <- c('GAD1','ADARB2','PAX6','LAMP5','VIP','LHX6','SST','PVALB','SLC17A','LINC00507','RORB','THEMIS','FEZF2','SLC1A3','PDGFRA','FGFR3','SLC14A1','GFAP','OPALIN','NOSTRIN','TYROBP')

markers <- c('SNAP25','SYT1','GRIN1', #neuron
             'AQP4','SLC1A2','GJA1','GFAP','SLC1A3', #astrocyte
             'MBP','SOX10', 'PLP1', 'MOBP','MOG','SLC1A3','SLC1A3', #ODC
             'PDGFRA','PCDH15', #OPC
             'HLA-DRA', 'CX3CR1', 'CSF1R','SLC1A3','PTPRC', 'TMEM119', 'P2RY12', #microglia
             'FLT1', 'CLDN5','SLC1A3','TEK', 'PDGFB', 'NOS3', 'ADGRL4', 'PECAM1', #endothelial
             'AMBP','PDGFRB','COBLL1','KCNJ8', 'ABCC9', 'ATP13A5', 'ART3', 'PLA1A', 'ACE2' #pericyte
             )


markers <- c('SNAP25','SYT1','GRIN1', #neuron
             'GAD1','GAD2','SLC6A1','SLC32A1', #inhibitory neuron
             'SLC17A7','STAB2','NRGN','CAMK2A', #excitatory neuron
             'RYR1','RASGRF2','CA8','PVALB','GRIA1', #purkinje neuron
             'TIAM1','PPFIA2','HTR2C','TLE4','TSHZ2','FOXP2','GRM4','PCP4', #granule cell
             'MEIS2',
             'SNCG',
             
             'SST', #SST+
             'VIP', #VIP+
             'SV2C', #SV2C+
             'LAMP5',
             'PVALB', #PVALB+
             'CUX2', #upper-layer
             'RORB', #layer 4
             'TLE4',
             'AQP4','SLC1A2','GJA1','GFAP','SLC1A3', #astrocyte
             'MBP','SOX10', 'PLP1', 'MOBP','MOG','SLC1A3','SLC1A3', #ODC
             'PDGFRA','PCDH15', #OPC
             'HLA-DRA', 'CX3CR1', 'CSF1R','SLC1A3','PTPRC', 'TMEM119', 'P2RY12', #microglia
             'FLT1', 'CLDN5','SLC1A3','TEK', 'PDGFB', 'NOS3', 'ADGRL4', 'PECAM1', #endothelial
             'AMBP','PDGFRB','COBLL1','KCNJ8', 'ABCC9', 'ATP13A5', 'ART3', 'PLA1A', 'ACE2' #pericyte
)

#celltype_dotplot.pdf 25x20
DotPlot(seu, features = unique(markers), assay = 'RNA') + theme(axis.text.x = element_text(angle = 90, hjust = 1))


markers <- c('MEIS2','HPF','SNCG','SERPINF1', 'KRT73', 'NTNG1', 'JAM2', 'NPY2R', 
             'LAMP5', 'LHX6', 'PDLIM5', 'EGLN3', 'PAX6',
             'VIP', 'PCDH11X', 'MYBPC1', 'CP', 'RSPO1', 'LMO1', 'CBLN4', 'HPF', 'IGFBP6',
             'SST', 'SYNDIG1L', 'CRH', 'NTS', 'LMO1', 'HPF', 'MYH8', 'ETV1', 'NMBR', 'HPSE', 'CALB2', 'MME', 'CTSC',
             'PVALB', 'TH', 'LPL', 'VIPR2')
markers <- c('LAMP5','PAX6','ADARB2','LINC00507','THEMIS','RORB','FEZF2','THEMIS','FEZF2')
DotPlot(seu, features = unique(markers), assay = 'RNA') + RotatedAxis()





FeaturePlot(seu, 
            #slot = 'IGLL5',
            reduction = "umap", 
            features = c('FGFR3','SLC14A1','GFAP'),
            #ncol = 3,
            cols = c('grey','tomato','tomato4'))



############## cell composition ###################
#bigger.clusters1
cluster_dis<-cbind(as.vector(seu$clusterAnno),as.vector(seu$sample))
colnames(cluster_dis)<-c('cluster','sample')
cluster_dis<-as.data.frame(cluster_dis)
cluster_dis<-cluster_dis[cluster_dis$cluster != 'Unassigned',]
#cluster_dis$cluster <- factor(cluster_dis$cluster,
#                              levels = c(0,1,2,3,4,5,6,7,8))
cluster_dis$sample <- factor(cluster_dis$sample,
                             levels = c('S0128-A2','S0128-A6','S0128-A7','S0128-A8','S0128-A9','S0128-A10','S0128-A11','S0128-A13','S0128-A12','S0128-A14',
                                        'S0206-A1','S0206-A2','S0206-A3','S0206-A4','S0206-A5','S0206-A6','S0206-A7','S0206-A11','S0206-A12',
                                        'S0406-A1','S0406-A3','S0406-A4','S0406-A5','S0406-A6', 'S0406-A7', 'S0406-A8', 'S0406-A9', 'S0406-A10','S0406-A14',
                                        'S0426-A5', 'S0426-A8', 'S0426-A9', 'S0426-A10','S0426-A11','S0426-A12','S0426-A13','S0426-A14',
                                        'S0531-A1', 'S0531-A2', 'S0531-A3','S0531-A4', 'S0531-A13'))

library(ggplot2)
library(RColorBrewer)
##Celltype_sample_composition.pdf
##pdf 8x3
ggplot(cluster_dis, aes(sample)) + 
  geom_bar(aes(fill=cluster), position="fill", color="gray20") +#,alpha = 0.7
  scale_fill_manual(values=c(brewer.pal(12,"Set3"),brewer.pal(3,"Set2")))+ #values=colors[c(1:4,6)]
  ylab('Percentage')+
  theme(axis.line = element_line(colour = "gray10"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + RotatedAxis()



################## chi-test #####################
observe<-table(seu$sample,seu$seurat_clusters)

chisq <- chisq.test(observe)
#X-squared = 2598.9, df = 64, p-value < 2.2e-16

round(chisq$expected,2)

ratio <- chisq$observed/chisq$expected

annotate <- data.frame(Individual=sub('.*-','',rownames(ratio)),
                       Region=sub('-.*','',rownames(ratio))
)
rownames(annotate) <- rownames(ratio)
colnames(annotate) <- c('Region','Individual')

# Roe.heatmap.pdf 8x10
pheatmap(t(ratio),
         cluster_cols = F,
         cluster_rows = F,
         annotation_col = annotate)

library(corrplot)
library(grDevices)
col2 <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                               "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                               "#4393C3", "#2166AC", "#053061")))
#residual_dotplot.pdf 8x10
corrplot(chisq$residuals, is.cor = FALSE, col =  col2(200))



########
observe<-table(seu$region,seu$seurat_clusters)

chisq <- chisq.test(observe)
#X-squared = 2598.9, df = 64, p-value < 2.2e-16

round(chisq$expected,2)

ratio <- chisq$observed/chisq$expected

annotate <- data.frame(Individual=sub('.*-','',rownames(ratio)),
                       Region=sub('-.*','',rownames(ratio))
)
rownames(annotate) <- rownames(ratio)
colnames(annotate) <- c('Region','Individual')

# Roe.heatmap.pdf 8x10
pheatmap(t(ratio),
         #annotation_col = annotate,
         cluster_cols = F,
         cluster_rows = F)

library(corrplot)
library(grDevices)
col2 <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                               "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                               "#4393C3", "#2166AC", "#053061")))
#residual_dotplot_region.pdf 5x10
corrplot(chisq$residuals, is.cor = FALSE, col =  col2(200))






