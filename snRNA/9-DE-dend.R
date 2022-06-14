library(Seurat)

seu <- readRDS('seu.harmony.anno.rds')

compare <- rbind(c('380,378,382,359,353,372,333,337,342,341,339,343','392,396,259,400,386,43,76,78,38,41,42,35,25,30,32,87,65,73,63,59,61,44,45,69,72,56,57,58,110,114,121,117,119,93,183,177,180,186,191,192,203,204,200,201,139,194,129,131,162,174,126,127,133,137,135,167,141,318,285,321,309,216,297,277,245,278,251,276,293,263,269,261,266,212,248,249,222,224,306,295,271,211,229,322,329,331,327,323,325,348,20,1,16,21,22'),
                 c('380,378,382,359,353,372','333,337,342,341,339,343'),
                 c('380,378,382','359,353,372'),
                 c('380','378,382'),
                 c('378','382'),
                 c('359','353,372'),
                 c('353','372'),
                 c('333,337','342,341,339,343'),
                 c('333','337'),
                 c('342','341,339,343'),
                 c('341','339,343'),
                 c('339','343'),
                 
                 c('392,396','259,400,386,43,76,78,38,41,42,35,25,30,32,87,65,73,63,59,61,44,45,69,72,56,57,58,110,114,121,117,119,93,183,177,180,186,191,192,203,204,200,201,139,194,129,131,162,174,126,127,133,137,135,167,141,318,285,321,309,216,297,277,245,278,251,276,293,263,269,261,266,212,248,249,222,224,306,295,271,211,229,322,329,331,327,323,325,348,20,1,16,21,22'),
                 c('392','396'),
                 c('259,400,386,43,76,78,38,41,42,35,25,30,32,87,65,73,63,59,61,44,45,69,72,56,57,58,110,114,121,117,119,93,183,177,180,186,191,192,203,204,200,201,139,194,129,131,162,174,126,127,133,137,135,167,141,318,285,321,309,216,297,277,245,278,251,276,293,263,269,261,266,212,248,249,222,224,306,295,271,211,229,322','329,331,327,323,325,348,20,1,16,21,22'),
                 
                 c('329,331,327,323,325','348,20,1,16,21,22'),
                 c('329','331,327,323,325'),
                 c('331','327,323,325'),
                 c('327','323,325'),
                 c('323','325'),
                 c('348','20,1,16,21,22'),
                 c('20','1,16,21,22'),
                 c('1','16,21,22'),
                 c('16','21,22'),
                 c('21','22'),
                 
                 c('259,400','386,43,76,78,38,41,42,35,25,30,32,87,65,73,63,59,61,44,45,69,72,56,57,58,110,114,121,117,119,93,183,177,180,186,191,192,203,204,200,201,139,194,129,131,162,174,126,127,133,137,135,167,141,318,285,321,309,216,297,277,245,278,251,276,293,263,269,261,266,212,248,249,222,224,306,295,271,211,229,322'),
                 c('259','400'),
                 c('386,43,76,78,38,41,42,35,25,30,32,87,65,73,63,59,61,44,45,69,72,56,57,58','110,114,121,117,119,93,183,177,180,186,191,192,203,204,200,201,139,194,129,131,162,174,126,127,133,137,135,167,141,318,285,321,309,216,297,277,245,278,251,276,293,263,269,261,266,212,248,249,222,224,306,295,271,211,229,322'),
                 c('386','43,76,78,38,41,42,35,25,30,32,87,65,73,63,59,61,44,45,69,72,56,57,58'),
                 c('43','76,78,38,41,42,35,25,30,32,87,65,73,63,59,61,44,45,69,72,56,57,58'),
                 c('76,78','38,41,42,35,25,30,32,87,65,73,63,59,61,44,45,69,72,56,57,58'),
                 c('76','78'),
                 c('38,41,42,35,25,30,32','87,65,73,63,59,61,44,45,69,72,56,57,58'),
                 c('38,41,42','35,25,30,32'),
                 c('38','41,42'),
                 c('41','42'),
                 c('35','25,30,32'),
                 c('25','30,32'),
                 c('30','32'),
                 c('87,65,73','63,59,61,44,45,69,72,56,57,58'),
                 c('87','65,73'),
                 c('65','73'),
                 c('63,59,61','44,45,69,72,56,57,58'),
                 c('63','59,61'),
                 c('59','61'),
                 c('44,45','69,72,56,57,58'),
                 c('44','45'),
                 c('69,72','56,57,58'),
                 c('69','72'),
                 c('56','57,58'),
                 c('57','58'),
                 
                 c('110,114,121,117,119,93,183,177,180,186,191,192,203,204,200,201,139,194,129,131,162,174,126,127,133,137,135,167,141,318','285,321,309,216,297,277,245,278,251,276,293,263,269,261,266,212,248,249,222,224,306,295,271,211,229,322'),
                 c('110,114,121,117,119','93,183,177,180,186,191,192,203,204,200,201,139,194,129,131,162,174,126,127,133,137,135,167,141,318'),
                 c('110','114,121,117,119'),
                 c('114','121,117,119'),
                 c('121','117,119'),
                 c('117','119'),
                 c('93,183,177,180,186,191,192,203,204,200,201,139,194','129,131,162,174,126,127,133,137,135,167,141,318'),
                 c('93,183,177,180','186,191,192,203,204,200,201,139,194'),
                 c('93','183,177,180'),
                 c('183','177,180'),
                 c('177','180'),
                 c('186,191,192','203,204,200,201,139,194'),
                 c('186','191,192'),
                 c('191','192'),
                 c('203,204','200,201,139,194'),
                 c('203','204'),
                 c('200,201','139,194'),
                 c('200','201'),
                 c('139','194'),
                 c('129,131','162,174,126,127,133,137,135,167,141,318'),
                 c('129','131'),
                 c('162','174,126,127,133,137,135,167,141,318'),
                 c('174,126,127','133,137,135,167,141,318'),
                 c('174','126,127'),
                 c('126','127'),
                 c('133,137','135,167,141,318'),
                 c('133','137'),
                 c('135,167','141,318'),
                 c('135','167'),
                 c('141','318'),
                 
                 c('285','321,309,216,297,277,245,278,251,276,293,263,269,261,266,212,248,249,222,224,306,295,271,211,229,322'),
                 c('321,309,216,297,277,245,278,251,276,293,263,269,261,266','212,248,249,222,224,306,295,271,211,229,322'),
                 c('321,309,216,297','277,245,278,251,276,293,263,269,261,266'),
                 c('321','309,216,297'),
                 c('309','216,297'),
                 c('216','297'),
                 c('277,245,278','251,276,293,263,269,261,266'),
                 c('277','245,278'),
                 c('245','278'),
                 c('251','276,293,263,269,261,266'),
                 c('276,293','263,269,261,266'),
                 c('276','293'),
                 c('263','269,261,266'),
                 c('269','261,266'),
                 c('261','266'),
                 c('212','248,249,222,224,306,295,271,211,229,322'),
                 c('248,249','222,224,306,295,271,211,229,322'),
                 c('248','249'),
                 c('222,224,306','295,271,211,229,322'),
                 c('222','224,306'),
                 c('224','306'),
                 c('295','271,211,229,322'),
                 c('271','211,229,322'),
                 c('211','229,322'),
                 c('229','322'))

colnames(compare) <- c('clusters1','clusters2')

DE <- c()
for(i in 1:106){
  ident1 <- strsplit(compare[i,1],',')
  ident2 <- strsplit(compare[i,2],',')
  de<-FindMarkers(seu,
                  ident.1 = ident1$clusters1,
                  ident.2 = ident2$clusters2)
  de <- cbind(rep(compare[i,1],dim(de)[1]),rep(compare[i,2],dim(de)[1]),rownames(de),de)
  colnames(de)[1:3] <- c('clusters1','clusters2','genes')
  DE <- rbind(DE, de)
}
saveRDS(DE,'DE_dend.rds')



#dend1 <- readRDS('DE_dend1.rds')
#dend2 <- readRDS('DE_dend2.rds')
#dend3 <- readRDS('DE_dend3.rds')
#dend4 <- readRDS('DE_dend4.rds')
#dend5 <- readRDS('DE_dend5.rds')
#dend6 <- readRDS('DE_dend6.rds')
#dend7 <- readRDS('DE_dend7.rds')
#dend8 <- readRDS('DE_dend8.rds')
#dend9 <- readRDS('DE_dend9.rds')
#dend10 <- readRDS('DE_dend10.rds')
#
#dend <- rbind(dend1,
#              dend2,
#              dend3,
#              dend4,
#              dend5,
#              dend6,
#              dend7,
#              dend8,
#              dend9,
#              dend10)

#ii <- unique(paste0(dend$clusters1,';',dend$clusters2))
#jj <- paste0('DE',c(1:106))
#names(jj) <- ii
#
#dend$ID <- jj[paste0(dend$clusters1,';',dend$clusters2)]

#saveRDS(dend,'DE_dend.rds')
#write.table(dend,'DE_tree.txt',col.names = T,row.names = F,sep = '\t',quote = F)




library(Seurat)
seu <- readRDS('seu.harmony.anno.rds')

dend<-readRDS('DE_dend.rds')

id <- unique(dend$ID)
for(i in 2:10){
  tmp <- dend[dend$ID==id[i],]
  tmp_up <- tmp[tmp$avg_log2FC>0,]
  tmp_down <- tmp[tmp$avg_log2FC<0,]
  
  up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
  down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]
  
  up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
  down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()
  
  up_file <- paste0(id[i],'_up_top40.pdf')
  pdf(up_file, width=20, height=25)
  up_plot
  dev.off()
  
  down_file <- paste0(id[i],'_down_top40.pdf')
  pdf(down_file, width=20, height=25)
  down_plot
  dev.off()
}










i=10
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()








i=11
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()








i=12
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()








i=13
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()








i=14
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()








i=15
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()









i=16
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()









i=17
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()








i=18
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()








i=19
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()








i=20
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()








i=21
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()








i=22
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()








i=23
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()








i=24
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()








i=25
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()








i=26
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()








i=27
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()








i=28
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()








i=29
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()








i=30
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()








i=31
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()








i=32
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()







i=33
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()







i=34
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()







i=35
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()







i=36
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()







i=37
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()







i=38
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()







i=39
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()






i=40
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()








i=41
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()








i=42
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()







i=43
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()







i=44
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()







i=45
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()







i=46
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()







i=47
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()







i=48
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()







i=49
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()



i=50
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()








i=51
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()








i=52
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()







i=53
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()







i=54
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()







i=55
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()







i=56
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()







i=57
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()







i=58
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()







i=59
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()


i=60
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()








i=61
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()








i=62
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()







i=63
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()







i=64
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()







i=65
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()







i=66
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()







i=67
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()







i=68
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()







i=69
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()



i=70
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()








i=71
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()








i=72
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()







i=73
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()







i=74
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()







i=75
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()







i=76
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()







i=77
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()







i=78
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()







i=79
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()

i=80
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()








i=81
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()








i=82
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()







i=83
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()







i=84
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()







i=85
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()







i=86
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()







i=87
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()







i=88
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()







i=89
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()

i=90
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()








i=91
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()








i=92
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()







i=93
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()







i=94
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()







i=95
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()







i=96
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()







i=97
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()







i=98
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()







i=99
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()




i=100
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()


i=101
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()


i=102
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()


i=103
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()


i=104
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()


i=105
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()


i=106
tmp <- dend[dend$ID==id[i],]
tmp_up <- tmp[tmp$avg_log2FC>0,]
tmp_down <- tmp[tmp$avg_log2FC<0,]

up_genes <- tmp_up$genes[1:min(40,dim(tmp_up)[1])]
down_genes <- tmp_down$genes[1:min(40,dim(tmp_down)[1])]

up_plot <- DotPlot(seu, features = unique(up_genes), assay = 'RNA') + RotatedAxis()
down_plot <- DotPlot(seu, features = unique(down_genes), assay = 'RNA') + RotatedAxis()

up_file <- paste0(id[i],'_up_top40.pdf')
pdf(up_file, width=20, height=25)
up_plot
dev.off()

down_file <- paste0(id[i],'_down_top40.pdf')
pdf(down_file, width=20, height=25)
down_plot
dev.off()
