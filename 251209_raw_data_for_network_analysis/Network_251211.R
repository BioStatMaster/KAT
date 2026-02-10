#########################
#######  Packages #######
#########################

need <- c("tibble","EstimateGroupNetwork","huge","qgraph","stringr","igraph","ggraph","tidygraph","ggplot2","dplyr","tibble","scales","stringr","purrr")
for (pkg in need) if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
invisible(lapply(need, library, character.only = TRUE))

if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_version("QUIC", version = "1.1.1", repos = "https://cran.r-project.org")

#########################
####### Load data #######
#########################

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rawdata1 <- read.csv("251209_bak15_all_features.csv", header = T, na.strings = "")
rawdata1 <-t(rawdata1)
rawdata1 <- rawdata1[-c(1:56),]
colnames(rawdata1) <- rawdata1[NROW(rawdata1),]
rawdata1<-rawdata1[-NROW(rawdata1),]
rawdata1<-as.data.frame(rawdata1)
rawdata1[] <- lapply(rawdata1, function(x) {
  if (is.character(x) | is.factor(x)) {
    suppressWarnings(as.numeric(as.character(x)))
  } else x
})

rawdata2 <- read.csv("251209_myc234_all_features.csv", header = T, na.strings = "")
rawdata2 <-t(rawdata2)
rawdata2 <- rawdata2[-c(1:56),]
colnames(rawdata2) <- rawdata2[NROW(rawdata2),]
rawdata2<-rawdata2[-NROW(rawdata2),]
rawdata2<-as.data.frame(rawdata2)
rawdata2[] <- lapply(rawdata2, function(x) {
  if (is.character(x) | is.factor(x)) {
    suppressWarnings(as.numeric(as.character(x)))
  } else x
})


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rawdata3 <- read.csv("251209_npr1_all_features.csv", header = T, na.strings = "")
rawdata3 <-t(rawdata3)
rawdata3 <- rawdata3[-c(1:56),]
colnames(rawdata3) <- rawdata3[NROW(rawdata3),]
rawdata3<-rawdata3[-NROW(rawdata3),]
rawdata3<-as.data.frame(rawdata3)
rawdata3[] <- lapply(rawdata3, function(x) {
  if (is.character(x) | is.factor(x)) {
    suppressWarnings(as.numeric(as.character(x)))
  } else x
})


rawdata4 <- read.csv("251209_sobir1_all_features.csv", header = T, na.strings = "")
rawdata4 <-t(rawdata4)
rawdata4 <- rawdata4[-c(1:56),]
colnames(rawdata4) <- rawdata4[NROW(rawdata4),]
rawdata4<-rawdata4[-NROW(rawdata4),]
rawdata4<-as.data.frame(rawdata4)
rawdata4[] <- lapply(rawdata4, function(x) {
  if (is.character(x) | is.factor(x)) {
    suppressWarnings(as.numeric(as.character(x)))
  } else x
})


rawdata5 <- read.csv("251209_wt_all_features.csv", header = T, na.strings = "")
rawdata5 <-t(rawdata5)
rawdata5 <- rawdata5[-c(1:56),]
colnames(rawdata5) <- rawdata5[NROW(rawdata5),]
rawdata5<-rawdata5[-NROW(rawdata5),]
rawdata5<-as.data.frame(rawdata5)
rawdata5[] <- lapply(rawdata5, function(x) {
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

rawdata_merged <- merge(rawdata6, rawdata5, by = "row.names")
rownames(rawdata_merged)<-rawdata_merged[,1]
rawdata_merged<-rawdata_merged[,-1]
rowname<-rownames(rawdata_merged)
rawdata_merged<-as.data.frame(rawdata_merged)
View(rawdata_merged)


#############################
########## Fitting ##########
#############################

### Functions ###
compute_hubness_core <- function(W_abs, labels) {
  W_abs <- as.matrix(W_abs); diag(W_abs) <- 0
  g <- igraph::graph_from_adjacency_matrix(W_abs, mode = "undirected", weighted = TRUE, diag = FALSE)
  igraph::V(g)$name <- labels
  w <- igraph::E(g)$weight; eps <- 1e-8; igraph::E(g)$dist <- 1/pmax(w, eps)
  data.frame(
    node        = igraph::V(g)$name,
    degree      = as.numeric(igraph::degree(g)),
    strength    = as.numeric(igraph::strength(g, weights = w)),
    eigencent   = as.numeric(igraph::eigen_centrality(g, directed = FALSE, weights = w)$vector),
    betweenness = as.numeric(igraph::betweenness(g, weights = igraph::E(g)$dist, normalized = TRUE)),
    closeness   = closeness_std_or_harmonic(g),
    stringsAsFactors = FALSE
  )
}


partial_corr_from_precision <- function(Omega, eps = 1e-12) {
  stopifnot(nrow(Omega) == ncol(Omega))
  d <- diag(Omega); d[d < eps] <- eps
  D <- sqrt(d)
  Pc <- -Omega / (D %o% D); diag(Pc) <- 1; Pc
}


closeness_std_or_harmonic <- function(g_w) {
  n <- igraph::gorder(g_w)
  if (igraph::is.connected(g_w)) {
    D <- igraph::distances(g_w, weights = igraph::E(g_w)$dist)
    diag(D) <- 0
    return( (n - 1) / rowSums(D) )
  } else {
    D <- igraph::distances(g_w, weights = igraph::E(g_w)$dist)
    invD <- 1 / D
    diag(invD) <- 0
    invD[!is.finite(invD)] <- 0  # 1/Inf = 0
    return( rowSums(invD) / (n - 1) )
  }
}

### HUGE, 모든 변수

### bak15
set.seed(1)
varname1=colnames(rawdata1)
X <- scale(as.matrix(rawdata1))       
hu  <- huge(X, method="glasso", lambda = seq(2,0.05,len=20))
sel <- huge.select(hu, criterion="ebic")
Theta1 <- sel$opt.icov
colnames(Theta1)=rownames(Theta1)=varname1

### myc
set.seed(1)
varname2=colnames(rawdata2)
X <- scale(as.matrix(rawdata2))       
hu  <- huge(X, method="glasso", lambda = seq(2,0.05,len=20))
sel <- huge.select(hu, criterion="ebic")
Theta2<- sel$opt.icov
colnames(Theta2)=rownames(Theta2)=varname2

### npr
set.seed(1)
varname3=colnames(rawdata3)
X <- scale(as.matrix(rawdata3))       
hu  <- huge(X, method="glasso", lambda = seq(2,0.05,len=20))
sel <- huge.select(hu, criterion="ebic")
Theta3 <- sel$opt.icov
colnames(Theta3)=rownames(Theta3)=varname3

### sobir
set.seed(1)
varname4=colnames(rawdata4)
X <- scale(as.matrix(rawdata4))       
hu  <- huge(X, method="glasso", lambda = seq(2,0.05,len=20))
sel <- huge.select(hu, criterion="ebic")
Theta4 <- sel$opt.icov
colnames(Theta4)=rownames(Theta4)=varname4

### wt_all
set.seed(1)
varname5=colnames(rawdata5)
X <- scale(as.matrix(rawdata5))       
hu  <- huge(X, method="glasso", lambda = seq(2,0.05,len=20))
sel <- huge.select(hu, criterion="ebic")
Theta5 <- sel$opt.icov
colnames(Theta5)=rownames(Theta5)=varname5

### wt_all_with_hormone
set.seed(1)
varname6=colnames(rawdata_merged)
X <- scale(as.matrix(rawdata_merged))       
hu  <- huge(X, method="glasso", lambda = seq(2,0.05,len=20))
sel <- huge.select(hu, criterion="ebic")
Theta6 <- sel$opt.icov
colnames(Theta6)=rownames(Theta6)=varname6

Pcorr<-partial_corr_from_precision(Theta)
W  <- abs(Pcorr); diag(W) <- 0
Scores <- compute_hubness_core(W, varname)

save(Pcorr, W, Scores, file="results_251120.Rdata")
load("results_251120.Rdata")
