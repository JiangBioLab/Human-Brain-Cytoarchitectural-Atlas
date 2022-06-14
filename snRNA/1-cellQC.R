#BiocManager::install('Seurat',lib='/home/user/wangpingping/R/x86_64-pc-linux-gnu-library/4.0')


library(Seurat)

S0128_A2  <- Read10X('../1-matrix/S0128-A2')
S0128_A6  <- Read10X('../1-matrix/S0128-A6')
S0128_A7  <- Read10X('../1-matrix/S0128-A7')
S0128_A8  <- Read10X('../1-matrix/S0128-A8')
S0128_A9  <- Read10X('../1-matrix/S0128-A9')
S0128_A10 <- Read10X('../1-matrix/S0128-A10')
S0128_A11 <- Read10X('../1-matrix/S0128-A11')
S0128_A12 <- Read10X('../1-matrix/S0128-A12')
S0128_A13 <- Read10X('../1-matrix/S0128-A13')
S0128_A14 <- Read10X('../1-matrix/S0128-A14')
S0206_A1  <- Read10X('../1-matrix/S0206-A1')
S0206_A2  <- Read10X('../1-matrix/S0206-A2')
S0206_A3  <- Read10X('../1-matrix/S0206-A3')
S0206_A4  <- Read10X('../1-matrix/S0206-A4')
S0206_A5  <- Read10X('../1-matrix/S0206-A5')
S0206_A6  <- Read10X('../1-matrix/S0206-A6')
S0206_A7  <- Read10X('../1-matrix/S0206-A7')
S0206_A11 <- Read10X('../1-matrix/S0206-A11')
S0206_A12 <- Read10X('../1-matrix/S0206-A12')
S0406_A1  <- Read10X('../1-matrix/S0406-A1')
S0406_A3  <- Read10X('../1-matrix/S0406-A3')
S0406_A4  <- Read10X('../1-matrix/S0406-A4')
S0406_A5  <- Read10X('../1-matrix/S0406-A5')
S0406_A6  <- Read10X('../1-matrix/S0406-A6')
S0406_A7  <- Read10X('../1-matrix/S0406-A7')
S0406_A8  <- Read10X('../1-matrix/S0406-A8')
S0406_A9  <- Read10X('../1-matrix/S0406-A9')
S0406_A10 <- Read10X('../1-matrix/S0406-A10')
S0406_A14 <- Read10X('../1-matrix/S0406-A14')
S0426_A5  <- Read10X('../1-matrix/S0426-A5')
S0426_A8  <- Read10X('../1-matrix/S0426-A8')
S0426_A9  <- Read10X('../1-matrix/S0426-A9')
S0426_A10 <- Read10X('../1-matrix/S0426-A10')
S0426_A11 <- Read10X('../1-matrix/S0426-A11')
S0426_A12 <- Read10X('../1-matrix/S0426-A12')
S0426_A13 <- Read10X('../1-matrix/S0426-A13')
S0426_A14 <- Read10X('../1-matrix/S0426-A14')
S0531_A1  <- Read10X('../1-matrix/S0531-A1')
S0531_A2  <- Read10X('../1-matrix/S0531-A2')
S0531_A3  <- Read10X('../1-matrix/S0531-A3')
S0531_A4  <- Read10X('../1-matrix/S0531-A4')
S0531_A13 <- Read10X('../1-matrix/S0531-A13')
S0916_A2  <- Read10X('../1-matrix/S0916-A2')
S0916_A6  <- Read10X('../1-matrix/S0916-A6')
S0916_A7  <- Read10X('../1-matrix/S0916-A7')
S0916_A8  <- Read10X('../1-matrix/S0916-A8')
S0916_A9  <- Read10X('../1-matrix/S0916-A9')
S0916_A10 <- Read10X('../1-matrix/S0916-A10')
S0916_A11 <- Read10X('../1-matrix/S0916-A11')
S0916_A12 <- Read10X('../1-matrix/S0916-A12')
S0916_A13 <- Read10X('../1-matrix/S0916-A13')
S0916_A14 <- Read10X('../1-matrix/S0916-A14')

colnames(S0128_A2) <- paste0('S0128-A2_' ,colnames(S0128_A2))
colnames(S0128_A6) <- paste0('S0128-A6_' ,colnames(S0128_A6))
colnames(S0128_A7) <- paste0('S0128-A7_' ,colnames(S0128_A7))
colnames(S0128_A8) <- paste0('S0128-A8_' ,colnames(S0128_A8))
colnames(S0128_A9) <- paste0('S0128-A9_' ,colnames(S0128_A9))
colnames(S0128_A10) <- paste0('S0128-A10_',colnames(S0128_A10))
colnames(S0128_A11) <- paste0('S0128-A11_',colnames(S0128_A11))
colnames(S0128_A12) <- paste0('S0128-A12_',colnames(S0128_A12))
colnames(S0128_A13) <- paste0('S0128-A13_',colnames(S0128_A13))
colnames(S0128_A14) <- paste0('S0128-A14_',colnames(S0128_A14))
colnames(S0206_A1) <- paste0('S0206-A1_' ,colnames(S0206_A1))
colnames(S0206_A2) <- paste0('S0206-A2_' ,colnames(S0206_A2))
colnames(S0206_A3) <- paste0('S0206-A3_' ,colnames(S0206_A3))
colnames(S0206_A4) <- paste0('S0206-A4_' ,colnames(S0206_A4))
colnames(S0206_A5) <- paste0('S0206-A5_' ,colnames(S0206_A5))
colnames(S0206_A6) <- paste0('S0206-A6_' ,colnames(S0206_A6))
colnames(S0206_A7) <- paste0('S0206-A7_' ,colnames(S0206_A7))
colnames(S0206_A11) <- paste0('S0206-A11_',colnames(S0206_A11))
colnames(S0206_A12) <- paste0('S0206-A12_',colnames(S0206_A12))
colnames(S0406_A1) <- paste0('S0406-A1_' ,colnames(S0406_A1))
colnames(S0406_A3) <- paste0('S0406-A3_' ,colnames(S0406_A3))
colnames(S0406_A4) <- paste0('S0406-A4_' ,colnames(S0406_A4))
colnames(S0406_A5) <- paste0('S0406-A5_' ,colnames(S0406_A5))
colnames(S0406_A6) <- paste0('S0406-A6_' ,colnames(S0406_A6))
colnames(S0406_A7) <- paste0('S0406-A7_' ,colnames(S0406_A7))
colnames(S0406_A8) <- paste0('S0406-A8_' ,colnames(S0406_A8))
colnames(S0406_A9) <- paste0('S0406-A9_' ,colnames(S0406_A9))
colnames(S0406_A10) <- paste0('S0406-A10_',colnames(S0406_A10))
colnames(S0406_A14) <- paste0('S0406-A14_',colnames(S0406_A14))
colnames(S0426_A5) <- paste0('S0426-A5_' ,colnames(S0426_A5))
colnames(S0426_A8) <- paste0('S0426-A8_' ,colnames(S0426_A8))
colnames(S0426_A9) <- paste0('S0426-A9_' ,colnames(S0426_A9))
colnames(S0426_A10) <- paste0('S0426-A10_',colnames(S0426_A10))
colnames(S0426_A11) <- paste0('S0426-A11_',colnames(S0426_A11))
colnames(S0426_A12) <- paste0('S0426-A12_',colnames(S0426_A12))
colnames(S0426_A13) <- paste0('S0426-A13_',colnames(S0426_A13))
colnames(S0426_A14) <- paste0('S0426-A14_',colnames(S0426_A14))
colnames(S0531_A1) <- paste0('S0531-A1_' ,colnames(S0531_A1))
colnames(S0531_A2) <- paste0('S0531-A2_' ,colnames(S0531_A2))
colnames(S0531_A3) <- paste0('S0531-A3_' ,colnames(S0531_A3))
colnames(S0531_A4) <- paste0('S0531-A4_' ,colnames(S0531_A4))
colnames(S0531_A13) <- paste0('S0531-A13_',colnames(S0531_A13))
colnames(S0916_A2) <- paste0('S0916-A2_' ,colnames(S0916_A2))
colnames(S0916_A6) <- paste0('S0916-A6_' ,colnames(S0916_A6))
colnames(S0916_A7) <- paste0('S0916-A7_' ,colnames(S0916_A7))
colnames(S0916_A8) <- paste0('S0916-A8_' ,colnames(S0916_A8))
colnames(S0916_A9) <- paste0('S0916-A9_' ,colnames(S0916_A9))
colnames(S0916_A10) <- paste0('S0916-A10_',colnames(S0916_A10))
colnames(S0916_A11) <- paste0('S0916-A11_',colnames(S0916_A11))
colnames(S0916_A12) <- paste0('S0916-A12_',colnames(S0916_A12))
colnames(S0916_A13) <- paste0('S0916-A13_',colnames(S0916_A13))
colnames(S0916_A14) <- paste0('S0916-A14_',colnames(S0916_A14))


saveRDS(S0128_A2 ,'counts/S0128_A2.counts.rds') 
saveRDS(S0128_A6 ,'counts/S0128_A6.counts.rds') 
saveRDS(S0128_A7 ,'counts/S0128_A7.counts.rds') 
saveRDS(S0128_A8 ,'counts/S0128_A8.counts.rds') 
saveRDS(S0128_A9 ,'counts/S0128_A9.counts.rds') 
saveRDS(S0128_A10,'counts/S0128_A10.counts.rds')
saveRDS(S0128_A11,'counts/S0128_A11.counts.rds')
saveRDS(S0128_A12,'counts/S0128_A12.counts.rds')
saveRDS(S0128_A13,'counts/S0128_A13.counts.rds')
saveRDS(S0128_A14,'counts/S0128_A14.counts.rds')
saveRDS(S0206_A1 ,'counts/S0206_A1.counts.rds') 
saveRDS(S0206_A2 ,'counts/S0206_A2.counts.rds') 
saveRDS(S0206_A3 ,'counts/S0206_A3.counts.rds') 
saveRDS(S0206_A4 ,'counts/S0206_A4.counts.rds') 
saveRDS(S0206_A5 ,'counts/S0206_A5.counts.rds') 
saveRDS(S0206_A6 ,'counts/S0206_A6.counts.rds') 
saveRDS(S0206_A7 ,'counts/S0206_A7.counts.rds') 
saveRDS(S0206_A11,'counts/S0206_A11.counts.rds')
saveRDS(S0206_A12,'counts/S0206_A12.counts.rds')
saveRDS(S0406_A1 ,'counts/S0406_A1.counts.rds') 
saveRDS(S0406_A3 ,'counts/S0406_A3.counts.rds') 
saveRDS(S0406_A4 ,'counts/S0406_A4.counts.rds') 
saveRDS(S0406_A5 ,'counts/S0406_A5.counts.rds') 
saveRDS(S0406_A6 ,'counts/S0406_A6.counts.rds') 
saveRDS(S0406_A7 ,'counts/S0406_A7.counts.rds') 
saveRDS(S0406_A8 ,'counts/S0406_A8.counts.rds') 
saveRDS(S0406_A9 ,'counts/S0406_A9.counts.rds') 
saveRDS(S0406_A10,'counts/S0406_A10.counts.rds')
saveRDS(S0406_A14,'counts/S0406_A14.counts.rds')
saveRDS(S0426_A5 ,'counts/S0426_A5.counts.rds') 
saveRDS(S0426_A8 ,'counts/S0426_A8.counts.rds') 
saveRDS(S0426_A9 ,'counts/S0426_A9.counts.rds') 
saveRDS(S0426_A10,'counts/S0426_A10.counts.rds')
saveRDS(S0426_A11,'counts/S0426_A11.counts.rds')
saveRDS(S0426_A12,'counts/S0426_A12.counts.rds')
saveRDS(S0426_A13,'counts/S0426_A13.counts.rds')
saveRDS(S0426_A14,'counts/S0426_A14.counts.rds')
saveRDS(S0531_A1 ,'counts/S0531_A1.counts.rds') 
saveRDS(S0531_A2 ,'counts/S0531_A2.counts.rds') 
saveRDS(S0531_A3 ,'counts/S0531_A3.counts.rds') 
saveRDS(S0531_A4 ,'counts/S0531_A4.counts.rds') 
saveRDS(S0531_A13,'counts/S0531_A13.counts.rds')
saveRDS(S0916_A2 ,'counts/S0916_A2.counts.rds') 
saveRDS(S0916_A6 ,'counts/S0916_A6.counts.rds') 
saveRDS(S0916_A7 ,'counts/S0916_A7.counts.rds') 
saveRDS(S0916_A8 ,'counts/S0916_A8.counts.rds') 
saveRDS(S0916_A9 ,'counts/S0916_A9.counts.rds') 
saveRDS(S0916_A10,'counts/S0916_A10.counts.rds')
saveRDS(S0916_A11,'counts/S0916_A11.counts.rds')
saveRDS(S0916_A12,'counts/S0916_A12.counts.rds')
saveRDS(S0916_A13,'counts/S0916_A13.counts.rds')
saveRDS(S0916_A14,'counts/S0916_A14.counts.rds')


genes <- rownames(S0128_A2)

mtx <- cbind(S0128_A2[genes,],
             S0128_A6[genes,],
             S0128_A7[genes,],
             S0128_A8[genes,],
             S0128_A9[genes,],
             S0128_A10[genes,],
             S0128_A11[genes,],
             S0128_A12[genes,],
             S0128_A13[genes,],
             S0128_A14[genes,],
             S0206_A1[genes,],
             S0206_A2[genes,],
             S0206_A3[genes,],
             S0206_A4[genes,],
             S0206_A5[genes,],
             S0206_A6[genes,],
             S0206_A7[genes,],
             S0206_A11[genes,],
             S0206_A12[genes,],
             S0406_A1[genes,],
             S0406_A3[genes,],
             S0406_A4[genes,],
             S0406_A5[genes,],
             S0406_A6[genes,],
             S0406_A7[genes,],
             S0406_A8[genes,],
             S0406_A9[genes,],
             S0406_A10[genes,],
             S0406_A14[genes,],
             S0426_A5[genes,],
             S0426_A8[genes,],
             S0426_A9[genes,],
             S0426_A10[genes,],
             S0426_A11[genes,],
             S0426_A12[genes,],
             S0426_A13[genes,],
             S0426_A14[genes,],
             S0531_A1[genes,],
             S0531_A2[genes,],
             S0531_A3[genes,],
             S0531_A4[genes,],
             S0531_A13[genes,],
             S0916_A2[genes,],
             S0916_A6[genes,],
             S0916_A7[genes,],
             S0916_A8[genes,],
             S0916_A9[genes,],
             S0916_A10[genes,],
             S0916_A11[genes,],
             S0916_A12[genes,],
             S0916_A13[genes,],
             S0916_A14[genes,])
dim(mtx)
#36601 331146

library(SingleCellExperiment)
sc <- SingleCellExperiment(assays = list(counts = mtx)
                           #,colData <- data.frame(cell_names=cell[,2])
)

colData(sc)[1]<-sub('_.*','',colnames(mtx))
names(colData(sc))[1]<-'sample'

colData(sc)[2]<-sub('-.*','',colnames(mtx))
names(colData(sc))[2]<-'individual'

colData(sc)[3]<-sub('.*-','',sub('_.*','',colnames(mtx)))
names(colData(sc))[3]<-'region'


sample <- c('S0128-A2','S0128-A6','S0128-A7','S0128-A8','S0128-A9','S0128-A10','S0128-A11','S0128-A12','S0128-A13','S0128-A14',
            'S0206-A1','S0206-A2','S0206-A3','S0206-A4','S0206-A5','S0206-A6','S0206-A7','S0206-A11','S0206-A12',
            'S0406-A1','S0406-A3','S0406-A4','S0406-A5','S0406-A6', 'S0406-A7', 'S0406-A8', 'S0406-A9', 'S0406-A10','S0406-A14',
            'S0426-A5', 'S0426-A8', 'S0426-A9', 'S0426-A10','S0426-A11','S0426-A12','S0426-A13','S0426-A14',
            'S0531-A1', 'S0531-A2', 'S0531-A3','S0531-A4', 'S0531-A13',
            'S0916-A2','S0916-A6','S0916-A7','S0916-A8','S0916-A9','S0916-A10','S0916-A11','S0916-A12','S0916-A13','S0916-A14')


saveRDS(sc,'sc.rds')

doublet<-readRDS('doublet/merge_doublet.rds')
sc$doublet <- doublet[colnames(sc),2]

anno<-readRDS('singleR/merge_singleR.rds')
sc$singleR.blue <- anno[colnames(sc),'blue']
sc$singleR.data <- anno[colnames(sc),'data']
sc$singleR.human <- anno[colnames(sc),'human']
sc$singleR.imm <- anno[colnames(sc),'imm']
sc$singleR.monaco <- anno[colnames(sc),'monaco']
sc$singleR.nover <- anno[colnames(sc),'nover']



library(scater)

sc <- addPerCellQC(sc, 
                   subsets=list(Mito=grep("MT-", rownames(sc))))

plotColData(sc, x = "sum", y="detected", colour_by="sample") 

plotColData(sc, x = "sum", y="subsets_Mito_percent",colour_by="sample") 


#不要使用“一刀切”的模式，因为每种质控都有自己的范围，它们应该根据自己的特点去进行过滤。
#例如，不管细胞质量低不低，只要测序深，reads就多、有表达的feature数量也就多。
#还有，如果实验本身使用的spike-in RNA就多，那么最后得到的spike-in比例就高。
#我们过滤时就要先假设大部分细胞都是合格的，我们只需要剔除那些在QC结果中的离群点（Outliers）。
#这些离群点的定义是根据绝对中位差（median absolute deviation，MAD），具体方法是：
#Remove cells with log-library sizes that are more than 3 MADs below the median log-library size
#看到使用了log归一化处理，它可以让非常大的值和非常小的值能够放在可以比较的范围之内。
#下面过滤文库大小和有表达的feature数量（以低于中位数3倍MAD做为过滤条件）
# 一般也就是去掉低文库大小、低表达基因、高spike-in
libsize.drop1 <- isOutlier(sc$sum, nmads=3, type="lower", 
                           log=TRUE, batch=sc$sample)

feature.drop1 <- isOutlier(sc$detected, nmads=3, type="lower", 
                           log=TRUE, batch=sc$sample)

libsize.drop2 <- isOutlier(sc$sum, nmads=3, type="both", 
                           log=TRUE, batch=sc$sample)

feature.drop2 <- isOutlier(sc$detected, nmads=3, type="both", 
                           log=TRUE, batch=sc$sample)

spike.drop <- isOutlier(sc$subsets_Mito_percent, nmads=3, type="higher",
                        batch=sc$sample)

attr(libsize.drop1, "thresholds")
#       S0128-A10 S0128-A11 S0128-A12 S0128-A13 S0128-A14 S0128-A2 S0128-A6 S0128-A7 S0128-A8 S0128-A9 S0206-A1 S0206-A11 S0206-A12 S0206-A2
#lower   113.0991  72.30725  58.87049  37.48075  36.52515 58.02702 55.42414 59.58252  181.996 94.66247 320.2307  93.66148  173.1115 83.69372
#higher       Inf       Inf       Inf       Inf       Inf      Inf      Inf      Inf      Inf      Inf      Inf       Inf       Inf      Inf
#       S0206-A3 S0206-A4 S0206-A5 S0206-A6 S0206-A7 S0406-A1 S0406-A10 S0406-A14 S0406-A3 S0406-A4 S0406-A5 S0406-A6 S0406-A7 S0406-A8
#lower  30.39854 63.49327 75.97038  59.8857 98.58884 640.9131  149.7715  220.3431 151.2569 61.71395 123.8424 42.13712 51.76468  82.8142
#higher      Inf      Inf      Inf      Inf      Inf      Inf       Inf       Inf      Inf      Inf      Inf      Inf      Inf      Inf
       #S0406-A9 S0426-A10 S0426-A11 S0426-A12 S0426-A13 S0426-A14 S0426-A5 S0426-A8 S0426-A9 S0531-A1 S0531-A13 S0531-A2 S0531-A3 S0531-A4
#lower  165.7444   41.9713  13.19449  17.44055   59.9789  148.8561 20.12522 72.24173 22.65291 69.16277  64.95556 130.8256 84.70487 67.53269
#higher      Inf       Inf       Inf       Inf       Inf       Inf      Inf      Inf      Inf      Inf       Inf      Inf      Inf      Inf
#       S0916-A10 S0916-A11 S0916-A12 S0916-A13 S0916-A14 S0916-A2 S0916-A6 S0916-A7 S0916-A8 S0916-A9
#lower   53.67105  106.9311  118.1616   70.0582  85.61559 97.83657 82.75046 59.56572 39.48553 33.24996
#higher       Inf       Inf       Inf       Inf       Inf      Inf      Inf      Inf      Inf      Inf
attr(feature.drop1, "thresholds")
#       S0128-A10 S0128-A11 S0128-A12 S0128-A13 S0128-A14 S0128-A2 S0128-A6 S0128-A7 S0128-A8 S0128-A9 S0206-A1 S0206-A11 S0206-A12 S0206-A2
#lower   121.9811  114.5309  88.30656  79.10356  96.35006 84.96323 89.91536 94.15257 215.1607 110.8784 417.0106  171.5042  276.6429 161.1942
#higher       Inf       Inf       Inf       Inf       Inf      Inf      Inf      Inf      Inf      Inf      Inf       Inf       Inf      Inf
#       S0206-A3 S0206-A4 S0206-A5 S0206-A6 S0206-A7 S0406-A1 S0406-A10 S0406-A14 S0406-A3 S0406-A4 S0406-A5 S0406-A6 S0406-A7 S0406-A8
#lower  86.68063 128.8016 158.4796 151.5983 263.1201 949.9573  368.2864   451.289 366.7466 205.9443 371.1235 124.5338 165.1987 211.6575
#higher      Inf      Inf      Inf      Inf      Inf      Inf       Inf       Inf      Inf      Inf      Inf      Inf      Inf      Inf
#       S0406-A9 S0426-A10 S0426-A11 S0426-A12 S0426-A13 S0426-A14 S0426-A5 S0426-A8 S0426-A9 S0531-A1 S0531-A13 S0531-A2 S0531-A3 S0531-A4
#lower   372.729  187.1464  49.30719  62.04754  228.4928  242.0625 43.57721  248.665 84.38731 190.6164  171.5985 277.5694 223.0629 153.4458
#higher      Inf       Inf       Inf       Inf       Inf       Inf      Inf      Inf      Inf      Inf       Inf      Inf      Inf      Inf
#       S0916-A10 S0916-A11 S0916-A12 S0916-A13 S0916-A14 S0916-A2 S0916-A6 S0916-A7 S0916-A8 S0916-A9
#lower   146.7519  200.8503  291.9383  177.4455  168.7913 198.6731 214.9737 159.4675 137.2106 115.6947
#higher       Inf       Inf       Inf       Inf       Inf      Inf      Inf      Inf      Inf      Inf
attr(libsize.drop2, "thresholds")
#        S0128-A10    S0128-A11    S0128-A12    S0128-A13    S0128-A14    S0128-A2     S0128-A6     S0128-A7  S0128-A8    S0128-A9   S0206-A1
#lower    113.0991     72.30725     58.87049     37.48075     36.52515    58.02702     55.42414     59.58252   181.996    94.66247   320.2307
#higher 38363.6036 107575.95587 104810.61329 242369.62295 575178.21178 87787.52945 103752.69363 140759.66121 45069.657 44307.99346 75007.8120
#          S0206-A11  S0206-A12     S0206-A2     S0206-A3     S0206-A4     S0206-A5    S0206-A6     S0206-A7    S0406-A1   S0406-A10
#lower      93.66148   173.1115     83.69372     30.39854     63.49327     75.97038     59.8857     98.58884    640.9131    149.7715
#higher 239172.90885 90449.5264 329200.32144 929994.93866 232722.85553 276112.86745 784563.9485 761742.41090 197419.5145 451801.9100
#         S0406-A14    S0406-A3     S0406-A4    S0406-A5     S0406-A6     S0406-A7    S0406-A8    S0406-A9    S0426-A10    S0426-A11
#lower     220.3431    151.2569 6.171395e+01    123.8424     42.13712 5.176468e+01     82.8142    165.7444      41.9713 1.319449e+01
#higher 365406.1394 632421.9924 1.261246e+06 976870.8077 724173.17567 1.022427e+06 483993.8132 377497.5511 1927318.0471 1.761494e+06
#          S0426-A12    S0426-A13   S0426-A14     S0426-A5     S0426-A8     S0426-A9     S0531-A1    S0531-A13    S0531-A2     S0531-A3
#lower  1.744055e+01      59.9789    148.8561     20.12522 7.224173e+01 2.265291e+01 6.916277e+01     64.95556    130.8256     84.70487
#higher 1.552199e+06 1441696.6558 126186.0048 828965.19233 1.094610e+06 1.593173e+06 1.001584e+06 987877.13911 575236.1679 970552.02273
#           S0531-A4    S0916-A10   S0916-A11   S0916-A12   S0916-A13    S0916-A14     S0916-A2     S0916-A6     S0916-A7     S0916-A8
#lower      67.53269 5.367105e+01    106.9311    118.1616     70.0582     85.61559     97.83657     82.75046 5.956572e+01 3.948553e+01
#higher 698061.31405 1.128056e+06 469434.8193 828452.5289 982377.7346 381753.94645 514956.78063 857665.34246 1.051866e+06 1.956992e+06
#           S0916-A9
#lower  3.324996e+01
#higher 1.581701e+06

attr(spike.drop, "thresholds")
#       S0128-A10 S0128-A11 S0128-A12 S0128-A13 S0128-A14 S0128-A2 S0128-A6 S0128-A7 S0128-A8 S0128-A9  S0206-A1 S0206-A11 S0206-A12  S0206-A2
#lower       -Inf      -Inf      -Inf      -Inf      -Inf     -Inf     -Inf     -Inf     -Inf     -Inf      -Inf      -Inf      -Inf      -Inf
#higher  19.25042  5.053502    6.0889  6.242009  5.648428 6.459528 10.24398 7.785563 5.694916 10.99277 0.8440413 0.4017636 0.7290572 0.7268691
#        S0206-A3  S0206-A4  S0206-A5  S0206-A6  S0206-A7  S0406-A1 S0406-A10 S0406-A14  S0406-A3 S0406-A4  S0406-A5  S0406-A6 S0406-A7
#lower       -Inf      -Inf      -Inf      -Inf      -Inf      -Inf      -Inf      -Inf      -Inf     -Inf      -Inf      -Inf     -Inf
#higher 0.6719029 0.8338352 0.7904421 0.9108463 0.8134889 0.5619244  1.045953  1.261159 0.8112225 1.538097 0.6190881 0.7930609  1.25304
#       S0406-A8 S0406-A9 S0426-A10 S0426-A11 S0426-A12 S0426-A13 S0426-A14 S0426-A5 S0426-A8 S0426-A9 S0531-A1 S0531-A13 S0531-A2 S0531-A3
#lower      -Inf     -Inf      -Inf      -Inf      -Inf      -Inf      -Inf     -Inf     -Inf     -Inf     -Inf      -Inf     -Inf     -Inf
#higher 2.576055 1.279361  1.842093  1.411724  1.122135  1.301772   2.05095 5.771097 1.798565 2.108922 1.932489  2.066531  1.86022 1.675229
#       S0531-A4 S0916-A10 S0916-A11 S0916-A12 S0916-A13 S0916-A14 S0916-A2 S0916-A6 S0916-A7 S0916-A8 S0916-A9
#lower      -Inf      -Inf      -Inf      -Inf      -Inf      -Inf     -Inf     -Inf     -Inf     -Inf     -Inf
#higher 1.805154  2.990023  3.811565  1.886796  2.237937  2.559268 2.855018 2.786915 5.302234 4.811237 2.866671


sc$mito_3mads<- spike.drop
sc$libsize.drop1 <- libsize.drop1
sc$libsize.drop2 <- libsize.drop2
sc$feature.drop1 <- feature.drop1
sc$feature.drop2 <- feature.drop2

stats <- data.frame(cells = table(sc$sample)[sample]#,
                    ,doublet = table(sc$doublet=='Doublet',sc$sample)[2,][sample]#,
                    ,libsize_low = table(sc$libsize.drop1,sc$sample)[2,][sample]#,
                    ,libsize_both = table(sc$libsize.drop2,sc$sample)[2,][sample]#,
                    ,feature_low = table(sc$feature.drop1,sc$sample)[2,][sample]#,
                    ,feature_both = table(sc$feature.drop2,sc$sample)[2,][sample]#,
                    ,mito_3MADs = table(sc$mito_3mads,sc$sample)[2,][sample]#,
                    ,mito_10percent = table(sc$sample[sc$subsets_Mito_percent>=10])[sample]#,
                    ,mito_20percent = table(sc$sample[sc$subsets_Mito_percent>=20])[sample]#,
                    #,mito_30percent = table(sc$Sample[sc$subsets_Mito_percent>=30])
                    )

stats

keep1 <- !(sc$libsize.drop1 | sc$feature.drop1 | sc$mito_3mads | sc$doublet=='Doublet')
keep2 <- !(sc$libsize.drop2 | sc$feature.drop2 | sc$mito_3mads | sc$doublet=='Doublet')
data.frame(ByLibSize=sum(sc$libsize.drop1), 
           ByFeature=sum(sc$feature.drop1), 
           ByMito=sum(sc$mito_3mads), 
           ByDoublet=sum(sc$doublet=='Doublet'), 
           Remaining=sum(keep1))
#ByLibSize ByFeature ByMito ByDoublet Remaining
#      283      1238  55975     24837    252459

data.frame(ByLibSize=sum(sc$libsize.drop2), 
           ByFeature=sum(sc$feature.drop2), 
           ByMito=sum(sc$mito_3mads), 
           ByDoublet=sum(sc$doublet=='Doublet'), 
           Remaining=sum(keep2))
#ByLibSize ByFeature ByMito ByDoublet Remaining
#      529      1238  55975     24837    252296

sc$use1<-keep1
sc$use2<-keep2



#Gene QC
keep_genes1 <-apply(
  counts(sc[,colData(sc)$use1]),
  1,  # 1 by row, 2 by column
  function(x) length(x[x>1])>=2  #We define a gene as detectable if at least two cells contain more than 1 transcript from the gene.
)
keep_genes2 <-apply(
  counts(sc[,colData(sc)$use2]),
  1,  # 1 by row, 2 by column
  function(x) length(x[x>1])>=2  #We define a gene as detectable if at least two cells contain more than 1 transcript from the gene.
)
table(keep_genes1)
# FALSE  TRUE 
# 20734 18526 
table(keep_genes2)
#FALSE  TRUE 
#21075 18185 


rowData(sc)$use1 <- keep_genes1
rowData(sc)$use2 <- keep_genes2

sc.qc1<-sc[,colData(sc)$use1]
sc.qc2<-sc[,colData(sc)$use2]
dim(sc)
# 36601 331146
dim(sc.qc1)
# 36601 252459
dim(sc.qc2)
# 36601 252296


saveRDS(sc,'sc.rds')
saveRDS(sc.qc1,'sc.qc1.rds')
saveRDS(sc.qc2,'sc.qc2.rds')




logcounts(sc.qc) <- log2(counts(sc.qc) + 1)

#spike-in and MTsfinal_qc_stats.txt
png(file="F4_QC_totalGenes_PctMT_Frequency.png", width = 800,height = 600)
plotColData(
  sc.qc2,
  x = "total_features_by_counts",
  y = "pct_counts_MT",
  #shape_by ='sample',  #用户可以选择用哪个特征来绘制图像
  colour_by ='sample'   #用户可以选择用哪个特征来绘制图像
)
dev.off()

png(file="F4_QC_totalCounts_TotalGenes.png", width = 800,height = 600)
plotColData(
  sc.qc2,
  x = "total_counts",
  y = "total_features_by_counts",
  #shape_by ='sample',  #用户可以选择用哪个特征来绘制图像
  colour_by ='sample'   #用户可以选择用哪个特征来绘制图像
)
dev.off()





