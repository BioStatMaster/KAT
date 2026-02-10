library(energy)
library(Ball)
library(hierNet)
library(sprintr)
library(mlbench)
library(randomForest)
library(randomForestExplainer)
library(reticulate)
library(GSelection)
library(dplyr)
library(dHSIC)

#########################
####### Load data #######
#########################

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rawdata <- read.csv("Total_features.csv", header = T, na.strings = "")
rawdata2 <-t(rawdata)
colnames(rawdata2) <- rawdata2[1, ]
rawdata2 <- rawdata2[-c(1, 105:117),]
rawdata2[is.na(rawdata2)] <- 0
rawdata2<-as.data.frame(rawdata2)

rawdata2 <- lapply(rawdata2, function(x) {
  if (all(grepl("^[0-9.]+$", x))) as.numeric(x) else x
})

rawdata2 <- as.data.frame(rawdata2)
pf<-rowSums(rawdata2[,2:5])
pf<-log(pf)
rawdata2 <- rawdata2[,-(2:5)]
rawdata2 <- as.data.frame(cbind(pf,rawdata2))

rawdata2 <- rbind(rawdata2,rawdata2[85,])
rawdata2 <- rawdata2[order(rawdata2$group),]

ctldata<-log2(rawdata2[rep(1:8,12),-(1:2)])
trtdata<-log2(rawdata2[-c(1:8),-(1:2)])

data1<-(trtdata-ctldata)
data1<-as.data.frame(cbind(rawdata2[-c(1:8),1:2],data1))
names<-colnames(data1)
data1 <- data1 %>%
  group_by(group) %>%
  mutate(across(where(is.numeric), ~ as.numeric(scale(.x)))) %>%
  ungroup()

View(data1)

#######################
###### Screening ######
#######################



#############################
########## Fitting ##########
#############################


###### ALL GROUP ######

covariates<-data1[,3:ncol(data1)]
y<-data1$pf

hsic_scores <- sapply(1:ncol(covariates), function(j) {
  dhsic.test(covariates[, j], y)$statistic
})

ranked_vars <- order(hsic_scores, decreasing = TRUE)
screened<-head(ranked_vars, 25)

covariates1 = as.data.frame(covariates[,screened])

fit_hsic <- feature.selection(covariates1, y, d = min(5, ncol(covariates1)))
sel_idx  <- fit_hsic$hsic_selected_feature_index

covariates1 = covariates1[,sel_idx]
idx<-colnames(covariates1)
colnames(covariates1)=c("X1","X2","X3","X4","X5")
heatmap(cor(covariates1))
lmdata<-data.frame(y,covariates1)
group<-data1$group

### rf_interaction ### X1:X4

set.seed(5)
df <- data.frame(y = y, covariates1)
rf <- randomForest(y ~ ., data = df, ntree = 500)
rf$importance[order(rf$importance,decreasing = T),]

interactions <- min_depth_interactions(rf, colnames(covariates1))
interactions_sorted <- interactions[order(-interactions$occurrences), ]
interactions_sorted_depth <- interactions[order(interactions$mean_min_depth), ]
head(interactions_sorted, 3)
head(interactions_sorted_depth, 3)

### hierNet: NO QUADRATIC ### NULL

set.seed(199)
cvfit <- hierNet.path(covariates1, y, diagonal=F, strong = T)
cvfit <-hierNet.cv(cvfit,covariates1, nfolds=6, y)
bestlam <- cvfit$lamhat
fit <-hierNet(as.matrix(covariates1), y, lam=bestlam, strong = T, diagonal=F)
fit

lmfit<-lm(y ~ X1+X2+X3+X5+X1:X3, data=lmdata)
summary(lmfit)

### sprinter: QUADRATIC ### NULL
set.seed(1)
fit.cv <- cv.sprinter(
  x = covariates1,
  y = y,
  square = T,
  nfold=6,
  nlam = 100,
  lam_min_ratio = ifelse(nrow(covariates1) < ncol(covariates1), 0.01, 1e-04)
)

coefficients_table <- fit.cv$compact
print(coefficients_table)

lmfit<-lm(y ~ X1+X2+X3+X5+I(X3^2)+X1:X3, data=lmdata)
summary(lmfit)

# 결론: 인터랙션 없음

######  GROUP  1 ######

data2=data1[group=="le1" | group=="le2" | group=="le3" | group=="le4",]
covariates2<-data2[,3:ncol(data2)]
y<-data2$pf

hsic_scores <- sapply(1:ncol(covariates2), function(j) {
  dhsic.test(covariates2[, j], y)$statistic
})

ranked_vars <- order(hsic_scores, decreasing = TRUE)
screened<-head(ranked_vars, 25)
covariates2 = as.data.frame(covariates2[,screened])
fit_hsic <- feature.selection(covariates2, y, d = min(5, ncol(covariates2)))
sel_idx  <- fit_hsic$hsic_selected_feature_index

covariates2 = covariates2[,sel_idx1]
idx1<-colnames((covariates2))
colnames(covariates2)=c("X1","X2","X3","X4","X5")
heatmap(cor(covariates2))

lmdata1=as.data.frame(cbind(y,covariates2))

### rf_interaction ### X3:X4, X4:X5

set.seed(1)
df <- data.frame(y = y, covariates2)
rf <- randomForest(y ~ ., data = df, ntree = 500)
rf$importance[order(rf$importance, decreasing = T),]

interactions <- min_depth_interactions(rf, colnames(covariates2))
interactions_sorted <- interactions[order(-interactions$occurrences), ]
interactions_sorted_depth <- interactions[order(interactions$mean_min_depth), ]
head(interactions_sorted, 3)
head(interactions_sorted_depth, 3)

### hierNet: NO QUADRATIC ### X3:X4

set.seed(1)
cvfit <- hierNet.path(covariates2, y, diagonal=F)
cvfit <-hierNet.cv(cvfit,covariates2, nfolds=4, y)
bestlam <- cvfit$lamhat
fit <-hierNet(as.matrix(covariates2), y, lam=bestlam, strong = TRUE, diagonal=F)
fit
lmfit<-lm(y ~ X1+X2+X3+X4+X5+X1:X4+X1:X5+X2:X3+X2:X5+X3:X4+X3:X5, data=lmdata1)
summary(lmfit)

### sprinter: QUADRATIC ### X3:X4
set.seed(1)
fit.cv <- cv.sprinter(
  x = covariates2,
  y = y,
  square = T,
  nlam = 100,
  nfold = 4,
  lam_min_ratio = ifelse(nrow(covariates2) < ncol(covariates2), 0.01, 1e-04)
)
coefficients_table <- fit.cv$compact
print(coefficients_table)
lmfit<-lm(y ~ X1+X2+X3+X4+X5+I(X2^2)+I(X3^2)+I(X4^2)+X1:X5+X1:X2+X2:X4+X1:X3+X2:X3+X3:X4, data=lmdata1)
summary(lmfit)

##### 결론: X3:X4 유의함.

###### GROUP 2 ######

data3=data1[group=="mp1" | group=="mp2" | group=="mp3" | group=="mp4",]
covariates3<-data3[,3:ncol(data3)]
y<-data3$pf

hsic_scores <- sapply(1:ncol(covariates3), function(j) {
  dhsic.test(covariates3[, j], y)$statistic
})

ranked_vars <- order(hsic_scores, decreasing = TRUE)
screened<-head(ranked_vars, 25)
covariates3 = as.data.frame(covariates3[,screened])
fit_hsic <- feature.selection(covariates3, y, d = min(5, ncol(covariates3)))
sel_idx  <- fit_hsic$hsic_selected_feature_index

covariates3 = covariates3[,sel_idx1]
idx2<-colnames((covariates3))
colnames(covariates3)=c("X1","X2","X3","X4","X5")
heatmap(cor(covariates3))

lmdata2=as.data.frame(cbind(y,covariates3))

### rf_interaction ### X2:X5
set.seed(1)
df <- data.frame(y = y, covariates3)
rf <- randomForest(y ~ ., data = df, ntree = 500)

interactions <- min_depth_interactions(rf, colnames(covariates3))
interactions_sorted <- interactions[order(-interactions$occurrences), ]
interactions_sorted_depth <- interactions[order(interactions$mean_min_depth), ]
head(interactions_sorted, 3)
head(interactions_sorted_depth, 3)

### hierNet: NO QUADRATIC ### No

set.seed(1)
cvfit <- hierNet.path(covariates3, y, diagonal=F)
cvfit <-hierNet.cv(cvfit,covariates3, nfolds=4, y)
bestlam <- cvfit$lamhat
fit <-hierNet(as.matrix(covariates3), y, lam=bestlam, strong = TRUE, diagonal=F)
fit

lmfit<-lm(y ~ X1+X2+X3+X4+X5+X1:X2+X1:X3+X1:X4+X3:X4+X4:X5, data=lmdata2)
summary(lmfit)

### sprinter: QUADRATIC ### No

set.seed(1)
fit.cv <- cv.sprinter(
  x = covariates3,
  y = y,
  square = T,
  nlam = 100,
  nfold = 4,
  lam_min_ratio = ifelse(nrow(covariates3) < ncol(covariates3), 0.01, 1e-04)
)

coefficients_table <- fit.cv$compact
print(coefficients_table)

lmfit<-lm(y ~ X1+X2+X4+I(X4^2), data=lmdata2)
summary(lmfit)

##### 없음

######  GROUP  3 ######

data4=data1[group=="px1" | group=="px2" | group=="px3" | group=="px4",]
covariates4<-data4[,3:ncol(data4)]
y<-data4$pf

hsic_scores <- sapply(1:ncol(covariates4), function(j) {
  dhsic.test(covariates4[, j], y)$statistic
})

ranked_vars <- order(hsic_scores, decreasing = TRUE)
screened<-head(ranked_vars, 25)
covariates4 = as.data.frame(covariates4[,screened])
fit_hsic <- feature.selection(covariates4, y, d = min(5, ncol(covariates4)))
sel_idx  <- fit_hsic$hsic_selected_feature_index

covariates4 = covariates4[,sel_idx1]
idx3<-colnames((covariates4))
colnames(covariates4)=c("X1","X2","X3","X4","X5")
heatmap(cor(covariates4))

lmdata3=as.data.frame(cbind(y,covariates4))

### rf_interaction ### X1:X5, X1:X2

set.seed(1)
df <- data.frame(y = y, covariates4)
rf <- randomForest(y ~ ., data = df, ntree = 500)
rf$importance

interactions <- min_depth_interactions(rf, colnames(covariates4))
interactions_sorted <- interactions[order(-interactions$occurrences), ]
interactions_sorted_depth <- interactions[order(interactions$mean_min_depth), ]
head(interactions_sorted, 3)
head(interactions_sorted_depth, 3)

### hierNet: NO QUADRATIC ### X2:X4, X2:X5****

set.seed(1)
cvfit <- hierNet.path(covariates4, y, diagonal=F)
cvfit <-hierNet.cv(cvfit,covariates4, nfolds=4, y)
bestlam <- cvfit$lamhat
fit <-hierNet(as.matrix(covariates4), y, lam=bestlam, strong = TRUE, diagonal=F)
fit

lmfit<-lm(y ~ X1+X2+X3+X4+X5+X1:X2+X1:X3+X2:X4+X2:X5+X3:X5+X3:X5+X4:X5, data=lmdata3)
summary(lmfit)

### sprinter: QUADRATIC ### X2:X4, X2:X5****
set.seed(1)
fit.cv <- cv.sprinter(
  x = covariates4,
  y = y,
  square = T,
  nlam = 100,
  nfold = 4,
  lam_min_ratio = ifelse(nrow(covariates4) < ncol(covariates4), 0.01, 1e-04)
)

coefficients_table <- fit.cv$compact
print(coefficients_table)

lmfit<-lm(y ~ X2+X3+X4+I(X1^2)+I(X5^2)+X2:X5+X2:X4, data=lmdata3)
summary(lmfit)

#결론:랜덤포레스트랑 리니어 결과 다르지만. X2:X4, X2:X5가 강력함

######  GROUP  4 ######

data5=data1[group=="sl1" | group=="sl2" | group=="sl3" | group=="sl4",]
covariates5<-data5[,3:ncol(data5)]
y<-data5$pf

hsic_scores <- sapply(1:ncol(covariates5), function(j) {
  dhsic.test(covariates5[, j], y)$statistic
})

ranked_vars <- order(hsic_scores, decreasing = TRUE)
screened<-head(ranked_vars, 25)
covariates5 = as.data.frame(covariates5[,screened])
fit_hsic <- feature.selection(covariates5, y, d = min(5, ncol(covariates5)))
sel_idx  <- fit_hsic$hsic_selected_feature_index

covariates5 = covariates5[,sel_idx1]
idx4<-colnames((covariates5))
colnames(covariates5)=c("X1","X2","X3","X4","X5")
heatmap(cor(covariates5))

lmdata4=as.data.frame(cbind(y,covariates5))

### rf_interaction ### X3:X4, X3:X5

set.seed(1)
df <- data.frame(y = y, covariates5)
rf <- randomForest(y ~ ., data = df, ntree = 500)
rf$importance

interactions <- min_depth_interactions(rf, colnames(covariates5))
interactions_sorted <- interactions[order(-interactions$occurrences), ]
interactions_sorted_depth <- interactions[order(interactions$mean_min_depth), ]
head(interactions_sorted, 3)
head(interactions_sorted_depth, 3)

### hierNet: NO QUADRATIC ### X4:X5

set.seed(1)
cvfit <- hierNet.path(covariates5, y, diagonal=F)
cvfit <-hierNet.cv(cvfit,covariates5, nfolds=4, y)
bestlam <- cvfit$lamhat
fit <-hierNet(as.matrix(covariates5), y, lam=bestlam, strong = TRUE, diagonal=F)
fit

lmfit<-lm(y ~ X1+X2+X3+X4+X5+X1:X4+X1:X5+X2:X4+X3:X5+X4:X5, data=lmdata4)
summary(lmfit)

### sprinter: QUADRATIC ### NULL
set.seed(1)
fit.cv <- cv.sprinter(
  x = covariates5,
  y = y,
  square = T,
  nlam = 100,
  nfold = 4,
  lam_min_ratio = ifelse(nrow(covariates5) < ncol(covariates5), 0.01, 1e-04)
)

coefficients_table <- fit.cv$compact
print(coefficients_table)

lmfit<-lm(y~X1+X2+X3+X5+I(X1^2)+I(X4^2)+X1:X4+X3:X5+X1:X5, data=lmdata4)
summary(lmfit)

#### 유의한 쌍이 없음.
