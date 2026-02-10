#########################
#######  Packages #######
#########################


need <- c("tibble","EstimateGroupNetwork","huge","qgraph","stringr","igraph","ggraph","tidygraph","ggplot2","dplyr","tibble","scales","stringr","purrr")
for (pkg in need) if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
invisible(lapply(need, library, character.only = TRUE))


#########################
####### Load data #######
#########################


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rawdata1 <- read.csv("annotated_features_DAM.csv", header = T, na.strings = "")
rawdata1 <-t(rawdata1)
colnames(rawdata1) <- rawdata1[3, ]
rawdata1 <- rawdata1[-c(1:16),]
rawdata1<-rawdata1[1:106,]

rawdata1<-as.data.frame(rawdata1)
rownames(rawdata1)
rawname<-str_sub(rownames(rawdata1), 1, -11)
rawname<-str_sub(rawname,10, nchar(rawname))
rownames(rawdata1) <- rawname
View(rawdata1)



rawdata2 <- read.csv("metadata_hormone_herbivore.csv", header = T, na.strings = "")
rawdata2$sample <- str_sub(rawdata2$sample, 1, -6)
rownames(rawdata2)<-rawdata2[,1]
rawdata2<-rawdata2[,-c(1:7)] 
View(rawdata2)


rawdata_merged <- merge(rawdata2, rawdata1, by = "row.names")
rownames(rawdata_merged)<-rawdata_merged[,1]
rawdata_merged<-rawdata_merged[,-1]
# rawdata_merged<-rawdata_merged[-(1:8),] 컨트롤제거
rowname<-rownames(rawdata_merged)
rawdata_merged[] <- lapply(rawdata_merged, function(x) {
  if (is.character(x) | is.factor(x)) {
    suppressWarnings(as.numeric(as.character(x)))
  } else x
})
#rawdata_merged$group <- factor(substr(rownames(rawdata_merged), 1, 1))
rawdata_merged<-as.data.frame(rawdata_merged)
write.csv(rawdata_merged, "data_merge.csv", row.names = FALSE)

View(rawdata_merged)


#############################
########## Fitting ##########
#############################


varname=colnames(rawdata_merged)

### HUGE, 모든 변수
set.seed(1)
X <- scale(as.matrix(rawdata_merged))       
hu  <- huge(X, method="glasso", lambda = seq(2,0.05,len=30))
sel <- huge.select(hu, criterion="ebic")
Theta <- sel$opt.icov

colnames(Theta)=rownames(Theta)=varname


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

Pcorr<-partial_corr_from_precision(Theta)
W  <- abs(Pcorr); diag(W) <- 0
Scores <- compute_hubness_core(W, varname)

save(Pcorr, W, Scores, file="results_251120.Rdata")
load("results_251120.Rdata")
