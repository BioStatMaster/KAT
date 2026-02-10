#########################
#######  Packages #######
#########################

need <- c("tibble","EstimateGroupNetwork","huge","qgraph","stringr","igraph","ggraph","tidygraph","ggplot2","dplyr","tibble","scales","stringr","purrr")
for (pkg in need) if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
invisible(lapply(need, library, character.only = TRUE))

#########################
####### Load data #######
#########################

##### bak
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rawdata11 <- read.csv("251215_bak15_ms2.csv", header = T, na.strings = "")
rawdata11 <-t(rawdata11)
rawdata11 <- rawdata11[-c(1:56),]
rawdata11 <- rawdata11[1:which(rownames(rawdata11)=="id"),]
colnames(rawdata11) <- rawdata11[NROW(rawdata11),]
rawdata11<-rawdata11[-NROW(rawdata11),]
rawdata11<-as.data.frame(rawdata11)
rawdata11[] <- lapply(rawdata11, function(x) {
  if (is.character(x) | is.factor(x)) {
    suppressWarnings(as.numeric(as.character(x)))
  } else x
})

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rawdata12 <- read.csv("251215_bak15_ms2_cor.csv", header = T, na.strings = "")
rawdata12 <-t(rawdata12)
rawdata12 <- rawdata12[-c(1:56),]
rawdata12 <- rawdata12[1:which(rownames(rawdata12)=="id"),]
colnames(rawdata12) <- rawdata12[NROW(rawdata12),]
rawdata12<-rawdata12[-NROW(rawdata12),]
rawdata12<-as.data.frame(rawdata12)
rawdata12[] <- lapply(rawdata12, function(x) {
  if (is.character(x) | is.factor(x)) {
    suppressWarnings(as.numeric(as.character(x)))
  } else x
})

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rawdata13 <- read.csv("251215_bak15_ms2_dam.csv", header = T, na.strings = "")
rawdata13 <-t(rawdata13)
rawdata13 <- rawdata13[-c(1:56),]
rawdata13 <- rawdata13[1:which(rownames(rawdata13)=="id"),]
colnames(rawdata13) <- rawdata13[NROW(rawdata13),]
rawdata13<-rawdata13[-NROW(rawdata13),]
rawdata13<-as.data.frame(rawdata13)
rawdata13[] <- lapply(rawdata13, function(x) {
  if (is.character(x) | is.factor(x)) {
    suppressWarnings(as.numeric(as.character(x)))
  } else x
})


###### myc

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rawdata21 <- read.csv("251215_myc234_ms2.csv", header = T, na.strings = "")
rawdata21 <-t(rawdata21)
rawdata21 <- rawdata21[-c(1:56),]
rawdata21 <- rawdata21[1:which(rownames(rawdata21)=="id"),]
colnames(rawdata21) <- rawdata21[NROW(rawdata21),]
rawdata21<-rawdata21[-NROW(rawdata21),]
rawdata21<-as.data.frame(rawdata21)
rawdata21[] <- lapply(rawdata21, function(x) {
  if (is.character(x) | is.factor(x)) {
    suppressWarnings(as.numeric(as.character(x)))
  } else x
})

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rawdata22 <- read.csv("251215_myc234_ms2_cor.csv", header = T, na.strings = "")
rawdata22 <-t(rawdata22)
rawdata22 <- rawdata22[-c(1:56),]
rawdata22 <- rawdata22[1:which(rownames(rawdata22)=="id"),]
colnames(rawdata22) <- rawdata22[NROW(rawdata22),]
rawdata22<-rawdata22[-NROW(rawdata22),]
rawdata22<-as.data.frame(rawdata22)
rawdata22[] <- lapply(rawdata22, function(x) {
  if (is.character(x) | is.factor(x)) {
    suppressWarnings(as.numeric(as.character(x)))
  } else x
})

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rawdata23 <- read.csv("251215_myc234_ms2_dam.csv", header = T, na.strings = "")
rawdata23 <-t(rawdata23)
rawdata23 <- rawdata23[-c(1:56),]
rawdata23 <- rawdata23[1:which(rownames(rawdata23)=="id"),]
colnames(rawdata23) <- rawdata23[NROW(rawdata23),]
rawdata23<-rawdata23[-NROW(rawdata23),]
rawdata23<-as.data.frame(rawdata23)
rawdata23[] <- lapply(rawdata23, function(x) {
  if (is.character(x) | is.factor(x)) {
    suppressWarnings(as.numeric(as.character(x)))
  } else x
})

###### npr

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rawdata31 <- read.csv("251215_npr1_ms2.csv", header = T, na.strings = "")
rawdata31 <-t(rawdata31)
rawdata31 <- rawdata31[-c(1:56),]
rawdata31 <- rawdata31[1:which(rownames(rawdata31)=="id"),]
colnames(rawdata31) <- rawdata31[NROW(rawdata31),]
rawdata31<-rawdata31[-NROW(rawdata31),]
rawdata31<-as.data.frame(rawdata31)
rawdata31[] <- lapply(rawdata31, function(x) {
  if (is.character(x) | is.factor(x)) {
    suppressWarnings(as.numeric(as.character(x)))
  } else x
})

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rawdata32 <- read.csv("251215_npr1_ms2_cor.csv", header = T, na.strings = "")
rawdata32 <-t(rawdata32)
rawdata32 <- rawdata32[-c(1:56),]
rawdata32 <- rawdata32[1:which(rownames(rawdata32)=="id"),]
colnames(rawdata32) <- rawdata32[NROW(rawdata32),]
rawdata32<-rawdata32[-NROW(rawdata32),]
rawdata32<-as.data.frame(rawdata32)
rawdata32[] <- lapply(rawdata32, function(x) {
  if (is.character(x) | is.factor(x)) {
    suppressWarnings(as.numeric(as.character(x)))
  } else x
})

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rawdata33 <- read.csv("251215_npr1_ms2_dam.csv", header = T, na.strings = "")
rawdata33 <-t(rawdata33)
rawdata33 <- rawdata33[-c(1:56),]
rawdata33 <- rawdata33[1:which(rownames(rawdata33)=="id"),]
colnames(rawdata33) <- rawdata33[NROW(rawdata33),]
rawdata33<-rawdata33[-NROW(rawdata33),]
rawdata33<-as.data.frame(rawdata33)
rawdata33[] <- lapply(rawdata33, function(x) {
  if (is.character(x) | is.factor(x)) {
    suppressWarnings(as.numeric(as.character(x)))
  } else x
})

###### sobir

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rawdata41 <- read.csv("251215_sobir1_ms2.csv", header = T, na.strings = "")
rawdata41 <-t(rawdata41)
rawdata41 <- rawdata41[-c(1:56),]
rawdata41 <- rawdata41[1:which(rownames(rawdata41)=="id"),]
colnames(rawdata41) <- rawdata41[NROW(rawdata41),]
rawdata41<-rawdata41[-NROW(rawdata41),]
rawdata41<-as.data.frame(rawdata41)
rawdata41[] <- lapply(rawdata41, function(x) {
  if (is.character(x) | is.factor(x)) {
    suppressWarnings(as.numeric(as.character(x)))
  } else x
})

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rawdata42 <- read.csv("251215_sobir1_ms2_cor.csv", header = T, na.strings = "")
rawdata42 <-t(rawdata42)
rawdata42 <- rawdata42[-c(1:56),]
rawdata42 <- rawdata42[1:which(rownames(rawdata42)=="id"),]
colnames(rawdata42) <- rawdata42[NROW(rawdata42),]
rawdata42<-rawdata42[-NROW(rawdata42),]
rawdata42<-as.data.frame(rawdata42)
rawdata42[] <- lapply(rawdata42, function(x) {
  if (is.character(x) | is.factor(x)) {
    suppressWarnings(as.numeric(as.character(x)))
  } else x
})

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rawdata43 <- read.csv("251215_sobir1_ms2_dam.csv", header = T, na.strings = "")
rawdata43 <-t(rawdata43)
rawdata43 <- rawdata43[-c(1:56),]
rawdata43 <- rawdata43[1:which(rownames(rawdata43)=="id"),]
colnames(rawdata43) <- rawdata43[NROW(rawdata43),]
rawdata43<-rawdata43[-NROW(rawdata43),]
rawdata43<-as.data.frame(rawdata43)
rawdata43[] <- lapply(rawdata43, function(x) {
  if (is.character(x) | is.factor(x)) {
    suppressWarnings(as.numeric(as.character(x)))
  } else x
})

###### wt

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rawdata51 <- read.csv("251215_wt_ms2.csv", header = T, na.strings = "")
rawdata51 <-t(rawdata51)
rawdata51 <- rawdata51[-c(1:56),]
rawdata51 <- rawdata51[1:which(rownames(rawdata51)=="id"),]
colnames(rawdata51) <- rawdata51[NROW(rawdata51),]
rawdata51<-rawdata51[-NROW(rawdata51),]
rawdata51<-as.data.frame(rawdata51)
rawdata51[] <- lapply(rawdata51, function(x) {
  if (is.character(x) | is.factor(x)) {
    suppressWarnings(as.numeric(as.character(x)))
  } else x
})

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rawdata52 <- read.csv("251215_wt_ms2_cor.csv", header = T, na.strings = "")
rawdata52 <-t(rawdata52)
rawdata52 <- rawdata52[-c(1:56),]
rawdata52 <- rawdata52[1:which(rownames(rawdata52)=="id"),]
colnames(rawdata52) <- rawdata52[NROW(rawdata52),]
rawdata52<-rawdata52[-NROW(rawdata52),]
rawdata52<-as.data.frame(rawdata52)
rawdata52[] <- lapply(rawdata52, function(x) {
  if (is.character(x) | is.factor(x)) {
    suppressWarnings(as.numeric(as.character(x)))
  } else x
})

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rawdata53 <- read.csv("251215_wt_ms2_dam.csv", header = T, na.strings = "")
rawdata53 <-t(rawdata53)
rawdata53 <- rawdata53[-c(1:56),]
rawdata53 <- rawdata53[1:which(rownames(rawdata53)=="id"),]
colnames(rawdata53) <- rawdata53[NROW(rawdata53),]
rawdata53<-rawdata53[-NROW(rawdata53),]
rawdata53<-as.data.frame(rawdata53)
rawdata53[] <- lapply(rawdata53, function(x) {
  if (is.character(x) | is.factor(x)) {
    suppressWarnings(as.numeric(as.character(x)))
  } else x
})


rawdata6 <- read.csv("metadata_hormone_herbivore.csv", header = T, na.strings = "")
rawdata6$sample <- str_sub(rawdata6$sample, 1, -6)
rownames(rawdata6)<-rawdata6[,1]
rawdata6<-rawdata6[,-c(1:7)] 
rawdata6[] <- lapply(rawdata6, function(x) {
  if (is.character(x) | is.factor(x)) {
    suppressWarnings(as.numeric(as.character(x)))
  } else x
})

rawdata61 <- merge(rawdata6, rawdata51, by = "row.names")
rownames(rawdata61)<-rawdata61[,1]
rawdata61<-rawdata61[,-1]
rowname<-rownames(rawdata61)
rawdata61<-as.data.frame(rawdata61)

rawdata62 <- merge(rawdata6, rawdata52, by = "row.names")
rownames(rawdata62)<-rawdata62[,1]
rawdata62<-rawdata62[,-1]
rowname<-rownames(rawdata62)
rawdata62<-as.data.frame(rawdata62)

rawdata63 <- merge(rawdata6, rawdata53, by = "row.names")
rownames(rawdata63)<-rawdata63[,1]
rawdata63<-rawdata63[,-1]
rowname<-rownames(rawdata63)
rawdata_merged<-as.data.frame(rawdata63)


#############################
########## Fitting ##########
#############################

### Fit GGM with HUGE

### bak15

set.seed(1)
varname11=colnames(rawdata11)
X <- scale(as.matrix(rawdata11))       
hu  <- huge(X, method="glasso", lambda = seq(2,0.05,len=20))
sel <- huge.select(hu, criterion="ebic")
Theta11 <- sel$opt.icov
colnames(Theta11)=rownames(Theta11)=varname11
save(Theta11,file="Theta11.Rdata")


set.seed(1)
varname12=colnames(rawdata12)
X <- scale(as.matrix(rawdata12))       
hu  <- huge(X, method="glasso", lambda = seq(2,0.05,len=20))
sel <- huge.select(hu, criterion="ebic")
Theta12 <- sel$opt.icov
colnames(Theta12)=rownames(Theta12)=varname12
save(Theta12,file="Theta12.Rdata")


set.seed(1)
varname13=colnames(rawdata13)
X <- scale(as.matrix(rawdata13))       
hu  <- huge(X, method="glasso", lambda = seq(2,0.05,len=20))
sel <- huge.select(hu, criterion="ebic")
Theta13 <- sel$opt.icov
colnames(Theta13)=rownames(Theta13)=varname13
save(Theta13,file="Theta13.Rdata")

### myc
set.seed(1)
varname21=colnames(rawdata21)
X <- scale(as.matrix(rawdata21))       
hu  <- huge(X, method="glasso", lambda = seq(2,0.05,len=20))
sel <- huge.select(hu, criterion="ebic")
Theta21<- sel$opt.icov
colnames(Theta21)=rownames(Theta21)=varname21
save(Theta21,file="Theta21.Rdata")


set.seed(1)
varname22=colnames(rawdata22)
X <- scale(as.matrix(rawdata22))       
hu  <- huge(X, method="glasso", lambda = seq(2,0.05,len=20))
sel <- huge.select(hu, criterion="ebic")
Theta22<- sel$opt.icov
colnames(Theta22)=rownames(Theta22)=varname22
save(Theta22,file="Theta22.Rdata")


set.seed(1)
varname23=colnames(rawdata23)
X <- scale(as.matrix(rawdata23))       
hu  <- huge(X, method="glasso", lambda = seq(2,0.05,len=20))
sel <- huge.select(hu, criterion="ebic")
Theta23<- sel$opt.icov
colnames(Theta23)=rownames(Theta23)=varname23
save(Theta23,file="Theta23.Rdata")


### npr
set.seed(1)
varname31=colnames(rawdata31)
X <- scale(as.matrix(rawdata31))       
hu  <- huge(X, method="glasso", lambda = seq(2,0.05,len=20))
sel <- huge.select(hu, criterion="ebic")
Theta31 <- sel$opt.icov
colnames(Theta31)=rownames(Theta31)=varname31
save(Theta31,file="Theta31.Rdata")


set.seed(1)
varname32=colnames(rawdata32)
X <- scale(as.matrix(rawdata32))       
hu  <- huge(X, method="glasso", lambda = seq(2,0.05,len=20))
sel <- huge.select(hu, criterion="ebic")
Theta32 <- sel$opt.icov
colnames(Theta32)=rownames(Theta32)=varname32
save(Theta32,file="Theta32.Rdata")


set.seed(1)
varname33=colnames(rawdata33)
X <- scale(as.matrix(rawdata33))       
hu  <- huge(X, method="glasso", lambda = seq(2,0.05,len=20))
sel <- huge.select(hu, criterion="ebic")
Theta33 <- sel$opt.icov
colnames(Theta33)=rownames(Theta33)=varname33
save(Theta33,file="Theta33.Rdata")

### sobir
set.seed(1)
varname41=colnames(rawdata41)
X <- scale(as.matrix(rawdata41))       
hu  <- huge(X, method="glasso", lambda = seq(2,0.05,len=20))
sel <- huge.select(hu, criterion="ebic")
Theta41 <- sel$opt.icov
colnames(Theta41)=rownames(Theta41)=varname41
save(Theta41,file="Theta41.Rdata")


set.seed(1)
varname42=colnames(rawdata42)
X <- scale(as.matrix(rawdata42))       
hu  <- huge(X, method="glasso", lambda = seq(2,0.05,len=20))
sel <- huge.select(hu, criterion="ebic")
Theta42 <- sel$opt.icov
colnames(Theta42)=rownames(Theta42)=varname42
save(Theta42,file="Theta42.Rdata")


set.seed(1)
varname43=colnames(rawdata43)
X <- scale(as.matrix(rawdata43))       
hu  <- huge(X, method="glasso", lambda = seq(2,0.05,len=20))
sel <- huge.select(hu, criterion="ebic")
Theta43 <- sel$opt.icov
colnames(Theta43)=rownames(Theta43)=varname43
save(Theta43,file="Theta43.Rdata")

### wt
set.seed(1)
varname51=colnames(rawdata51)
X <- scale(as.matrix(rawdata51))       
hu  <- huge(X, method="glasso", lambda = seq(2,0.05,len=20))
sel <- huge.select(hu, criterion="ebic")
Theta51 <- sel$opt.icov
colnames(Theta51)=rownames(Theta51)=varname51
save(Theta51,file="Theta51.Rdata")

set.seed(1)
varname52=colnames(rawdata52)
X <- scale(as.matrix(rawdata52))       
hu  <- huge(X, method="glasso", lambda = seq(2,0.05,len=20))
sel <- huge.select(hu, criterion="ebic")
Theta52 <- sel$opt.icov
colnames(Theta52)=rownames(Theta52)=varname52
save(Theta52,file="Theta52.Rdata")

set.seed(1)
varname53=colnames(rawdata53)
X <- scale(as.matrix(rawdata53))       
hu  <- huge(X, method="glasso", lambda = seq(2,0.05,len=20))
sel <- huge.select(hu, criterion="ebic")
Theta53 <- sel$opt.icov
colnames(Theta53)=rownames(Theta53)=varname53
save(Theta53,file="Theta53.Rdata")

### wt_all_with_hormone
set.seed(1)
varname61=colnames(rawdata61)
X <- scale(as.matrix(rawdata61))       
hu  <- huge(X, method="glasso", lambda = seq(2,0.05,len=20))
sel <- huge.select(hu, criterion="ebic")
Theta61 <- sel$opt.icov
colnames(Theta61)=rownames(Theta61)=varname61
save(Theta61,file="Theta61.Rdata")

set.seed(1)
varname62=colnames(rawdata62)
X <- scale(as.matrix(rawdata62))       
hu  <- huge(X, method="glasso", lambda = seq(2,0.05,len=20))
sel <- huge.select(hu, criterion="ebic")
Theta62 <- sel$opt.icov
colnames(Theta62)=rownames(Theta62)=varname62
save(Theta62,file="Theta62.Rdata")

set.seed(1)
varname63=colnames(rawdata63)
X <- scale(as.matrix(rawdata63))       
hu  <- huge(X, method="glasso", lambda = seq(2,0.05,len=20))
sel <- huge.select(hu, criterion="ebic")
Theta63 <- sel$opt.icov
colnames(Theta63)=rownames(Theta63)=varname63
save(Theta63,file="Theta63.Rdata")
