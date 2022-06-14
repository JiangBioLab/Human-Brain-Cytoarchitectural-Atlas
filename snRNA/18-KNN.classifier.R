# Installing Packages
install.packages("e1071")
install.packages("caTools")
install.packages("class")

# Loading package
library(e1071)
library(caTools)
library(class)

# Loading data
data(iris)
head(iris)

# Splitting data into train
# and test data
split <- sample.split(iris, SplitRatio = 0.7)
train_cl <- subset(iris, split == "TRUE")
test_cl <- subset(iris, split == "FALSE")

# Feature Scaling
train_scale <- scale(train_cl[, 1:4])
test_scale <- scale(test_cl[, 1:4])

# Fitting KNN Model 
# to training dataset
classifier_knn <- knn(train = train_scale,
                      test = test_scale,
                      cl = train_cl$Species,
                      k = 1)
classifier_knn

# Confusiin Matrix
cm <- table(test_cl$Species, classifier_knn)
cm

# Model Evaluation - Choosing K
# Calculate out of Sample error
misClassError <- mean(classifier_knn != test_cl$Species)
print(paste('Accuracy =', 1-misClassError))

# K = 3
classifier_knn <- knn(train = train_scale,
                      test = test_scale,
                      cl = train_cl$Species,
                      k = 3)
misClassError <- mean(classifier_knn != test_cl$Species)
print(paste('Accuracy =', 1-misClassError))

# K = 5
classifier_knn <- knn(train = train_scale,
                      test = test_scale,
                      cl = train_cl$Species,
                      k = 5)
misClassError <- mean(classifier_knn != test_cl$Species)
print(paste('Accuracy =', 1-misClassError))

# K = 7
classifier_knn <- knn(train = train_scale,
                      test = test_scale,
                      cl = train_cl$Species,
                      k = 7)
misClassError <- mean(classifier_knn != test_cl$Species)
print(paste('Accuracy =', 1-misClassError))

# K = 15
classifier_knn <- knn(train = train_scale,
                      test = test_scale,
                      cl = train_cl$Species,
                      k = 15)
misClassError <- mean(classifier_knn != test_cl$Species)
print(paste('Accuracy =', 1-misClassError))

# K = 19
classifier_knn <- knn(train = train_scale,
                      test = test_scale,
                      cl = train_cl$Species,
                      k = 19)
misClassError <- mean(classifier_knn != test_cl$Species)
print(paste('Accuracy =', 1-misClassError))


### output 
classifier_knn

### confusion matrix
cm




########################
seu <- readRDS('seu.harmony.GABA.CGE.rds')
DE <- readRDS('Inh_CGE_subclass_region_DE.rds')

subclasses <- unique(seu$hicat_cluster_subclasses)

for(cc in 1:length(subclasses)){
  iris <- seu@assays$RNA@data[,seu$hicat_cluster_subclasses==subclasses[cc]]
  #select <- rowSums(iris>0)>20
  
  DE0 <- DE[DE$subclasses==subclasses[cc],]
  select <- unique(DE0$gene)
  
  iris <- iris[select,]
  iris <- t(as.matrix(iris))
  iris <- as.data.frame(iris)
  
  iris$region <- seu$region[seu$hicat_cluster_subclasses==subclasses[cc]]
  
  # Splitting data into train
  # and test data
  aa <- names(table(iris$region))
  samples <- c()
  for(i in 1:length(aa)){
    tmp <- rownames(iris)[iris$region==aa[i]]
    tmp2 <- sample(tmp,length(tmp)/2)
    samples <- c(samples,tmp2)
  }
  
  train_cl <- iris[samples,]
  test_cl <- iris[!(rownames(iris) %in% samples),]
  
  # Feature Scaling
  train_scale <- scale(train_cl[, 1:length(select)])
  test_scale <- scale(test_cl[, 1:length(select)])
  
  # Fitting KNN Model 
  # to training dataset
  classifier_knn <- knn(train = train_scale,
                        test = test_scale,
                        cl = train_cl$region,
                        k = 15)
  classifier_knn
  
  # Confusiin Matrix
  cm <- table(test_cl$region, classifier_knn)
  cm
  
  subclasses[cc] <- sub('/','_',subclasses[cc])
  
  train_cl_file <- paste0('Inh.CGE.',subclasses[cc],'.train.rds')
  test_cl_file <- paste0('Inh.CGE.',subclasses[cc],'.test.rds')
  classifier_knn_file <- paste0('Inh.CGE.',subclasses[cc],'.knn.rds')
  cm_file <- paste0('Inh.CGE.',subclasses[cc],'.cm.rds')
  
  saveRDS(train_cl, train_cl_file)
  saveRDS(test_cl, test_cl_file)
  saveRDS(classifier_knn, classifier_knn_file)
  saveRDS(cm, cm_file)
}




######### MGE ##############################################
seu <- readRDS('seu.harmony.GABA.MGE.rds')
DE <- readRDS('Inh_MGE_subclass_region_DE.rds')

subclasses <- unique(seu$hicat_cluster_subclasses)

for(cc in 1:length(subclasses)){
  iris <- seu@assays$RNA@data[,seu$hicat_cluster_subclasses==subclasses[cc]]
  #select <- rowSums(iris>0)>20
  
  DE0 <- DE[DE$subclasses==subclasses[cc],]
  select <- unique(DE0$gene)
  
  iris <- iris[select,]
  iris <- t(as.matrix(iris))
  iris <- as.data.frame(iris)
  
  iris$region <- seu$region[seu$hicat_cluster_subclasses==subclasses[cc]]
  
  # Splitting data into train
  # and test data
  aa <- names(table(iris$region))
  samples <- c()
  for(i in 1:length(aa)){
    tmp <- rownames(iris)[iris$region==aa[i]]
    tmp2 <- sample(tmp,length(tmp)/2)
    samples <- c(samples,tmp2)
  }
  
  train_cl <- iris[samples,]
  test_cl <- iris[!(rownames(iris) %in% samples),]
  
  # Feature Scaling
  train_scale <- scale(train_cl[, 1:length(select)])
  test_scale <- scale(test_cl[, 1:length(select)])
  
  # Fitting KNN Model 
  # to training dataset
  classifier_knn <- knn(train = train_scale,
                        test = test_scale,
                        cl = train_cl$region,
                        k = 15)
  classifier_knn
  
  # Confusiin Matrix
  cm <- table(test_cl$region, classifier_knn)
  cm
  
  subclasses[cc] <- sub('/','_',subclasses[cc])
  
  train_cl_file <- paste0('Inh.MGE.',subclasses[cc],'.train.rds')
  test_cl_file <- paste0('Inh.MGE.',subclasses[cc],'.test.rds')
  classifier_knn_file <- paste0('Inh.MGE.',subclasses[cc],'.knn.rds')
  cm_file <- paste0('Inh.MGE.',subclasses[cc],'.cm.rds')
  
  saveRDS(train_cl, train_cl_file)
  saveRDS(test_cl, test_cl_file)
  saveRDS(classifier_knn, classifier_knn_file)
  saveRDS(cm, cm_file)
}




######### Glu ##############################################
seu <- readRDS('seu.harmony.Glu.rds')
DE <- readRDS('Exc_subclass_region_DE.rds')


seu$hicat_cluster_subclasses[seu$hicat_cluster_merge_newID=='107'] <- 'L6B/FEZF2'
seu$hicat_cluster_subclasses[seu$hicat_cluster_merge_newID=='80'] <- 'L6CT/FEZF2'

subclasses <- unique(as.vector(seu$hicat_cluster_subclasses))

subclasses <- c("ET","L6 CAR3")

for(cc in 1:length(subclasses)){
  iris <- seu@assays$RNA@data[,seu$hicat_cluster_subclasses==subclasses[cc]]
  #select <- rowSums(iris>0)>20
  
  DE0 <- DE[DE$subclasses==subclasses[cc],]
  select <- unique(DE0$gene)
  
  iris <- iris[select,]
  #iris[is.na(iris)] <- 0
  iris <- t(as.matrix(iris))
  iris <- as.data.frame(iris)
  
  iris$region <- as.vector(seu$region[seu$hicat_cluster_subclasses==subclasses[cc]])
  
  # Splitting data into train
  # and test data
  aa <- names(table(iris$region))
  samples <- c()
  for(i in 1:length(aa)){
    tmp <- rownames(iris)[iris$region==aa[i]]
    tmp2 <- sample(tmp,length(tmp)/2)
    samples <- c(samples,tmp2)
  }
  
  train_cl <- iris[samples,]
  test_cl <- iris[!(rownames(iris) %in% samples),]
  
  
  # Feature Scaling
  train_scale <- scale(train_cl[, 1:length(select)])
  test_scale <- scale(test_cl[, 1:length(select)])
  
  # Fitting KNN Model 
  # to training dataset
  classifier_knn <- knn(train = train_scale,
                        test = test_scale,
                        cl = train_cl$region,
                        k = 15)
  classifier_knn
  
  # Confusiin Matrix
  cm <- table(test_cl$region, classifier_knn)
  cm
  
  subclasses[cc] <- gsub('/','_',subclasses[cc])
  subclasses[cc] <- gsub(' ','_',subclasses[cc])
  
  train_cl_file <- paste0('Exc.',subclasses[cc],'.train.rds')
  test_cl_file <- paste0('Exc.',subclasses[cc],'.test.rds')
  classifier_knn_file <- paste0('Exc.',subclasses[cc],'.knn.rds')
  cm_file <- paste0('Exc.',subclasses[cc],'.cm.rds')
  
  saveRDS(train_cl, train_cl_file)
  saveRDS(test_cl, test_cl_file)
  saveRDS(classifier_knn, classifier_knn_file)
  saveRDS(cm, cm_file)
}



