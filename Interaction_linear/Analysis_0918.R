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
library(RAMP)
library(SIS)
library(iForm)
library(corrplot)

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
#pf<-log(pf)
rawdata2 <- rawdata2[,-(2:5)]
rawdata2 <- as.data.frame(cbind(pf,rawdata2))

rawdata2 <- rbind(rawdata2,rawdata2[85,])
rawdata2 <- rawdata2[order(rawdata2$group),]

ctldata<-(rawdata2[rep(1:8,12),-(1:2)])
trtdata<-(rawdata2[-c(1:8),-(1:2)])
log2FC<-log2(trtdata)-log2(ctldata)

alldata<-rawdata2[-c(1:8),-(1:2)]
alldata<-as.data.frame(cbind(rawdata2[-c(1:8),1:2],log2FC))
group<-alldata$group


###############################################
##########   Screening Meaningful X   #########
###############################################

test1=trtdata[group=="le1" | group=="le2" | group=="le3" | group=="le4",]
pvals1 <- mapply(
  function(col1, col2) wilcox.test(col1, col2)$p.value,
  test1, ctldata
)
p_adj1 <- p.adjust(pvals1, method = "BH")
idx1=(p_adj1<0.01)

test2=trtdata[group=="mp1" | group=="mp2" | group=="mp3" | group=="mp4",]
pvals2 <- mapply(
  function(col1, col2) wilcox.test(col1, col2)$p.value,
  test2, ctldata
)
p_adj2 <- p.adjust(pvals2, method = "BH")
idx2=(p_adj2<0.01)

test3=trtdata[group=="px1" | group=="px2" | group=="px3" | group=="px4",]
pvals3 <- mapply(
  function(col1, col2) wilcox.test(col1, col2)$p.value,
  test3, ctldata
)
p_adj3 <- p.adjust(pvals3, method = "BH")
idx3=(p_adj3<0.01)

test4=trtdata[group=="sl1" | group=="sl2" | group=="sl3" | group=="sl4",]
pvals4 <- mapply(
  function(col1, col2) wilcox.test(col1, col2)$p.value,
  test4, ctldata
)
p_adj4 <- p.adjust(pvals4, method = "BH")
idx4=(p_adj4<0.01)

#############################
########## Fitting ##########
#############################

######  GROUP  1 ######

data1=alldata[group=="le1" | group=="le2" | group=="le3" | group=="le4",]
y1=log(data1$pf)
X1=data1[,-(1:2)]
X1=X1[,idx1]
data1=as.data.frame(cbind(y1,X1))
data1=as.data.frame(scale(data1))

fit1 <- iForm(
  formula      = y1 ~ .,      
  data         = as.data.frame(data1),
  heredity     = "weak",   
  higher_order = FALSE 
)
fit1

fit1_1 <- RAMP(y=y1, X=X1, family="gaussian", tune = 'EBIC', penalty = "SCAD", hier="Weak")
fit1_1

###### GROUP 2 ######
data2=alldata[group=="mp1" | group=="mp2" | group=="mp3" | group=="mp4",]
y2=log(data2$pf)
X2=data2[,-(1:2)]
X2=X2[,idx2]
data2=as.data.frame(cbind(y2,X2))
data2=as.data.frame(scale(data2))

fit2 <- iForm(
  formula      = y2 ~ .,      
  data         = data2,
  heredity     = "weak",   
  higher_order = FALSE
)
fit2

fit2_2 <- RAMP(y=y2, X=X2, family="gaussian", tune = 'EBIC', penalty = "SCAD", hier="Weak")
fit2_2

######  GROUP  3 ######
data3=alldata[group=="px1" | group=="px2" | group=="px3" | group=="px4",]
y3=data3$pf
X3=data3[,-(1:2)]
X3=X3[,idx3]
data3=as.data.frame(cbind(y3,X3))
data3=as.data.frame(scale(data3))

fit3 <- iForm(
  formula      = y3 ~ .,      
  data         = data3,
  heredity     = "weak",   
  higher_order = FALSE 
)
fit3

fit3_3 <- RAMP(y=y3, X=X3, family="gaussian", tune = 'EBIC', penalty = "SCAD", hier="Weak")
fit3_3

######  GROUP  4 ######
data4=alldata[group=="sl1" | group=="sl2" | group=="sl3" | group=="sl4",]
y4=data4$pf
X4=data4[,-(1:2)]
X4=X4[,idx4]
data4=as.data.frame(cbind(y4,X4))
data4=as.data.frame(scale(data4))

fit4 <- iForm(
  formula      = y4 ~ .,      
  data         = data4,
  heredity     = "weak",   
  higher_order = FALSE 
)
fit4

fit4_4 <- RAMP(y=y4, X=X4, family="gaussian", tune = 'EBIC', penalty = "SCAD", hier="Weak")
fit4_4
