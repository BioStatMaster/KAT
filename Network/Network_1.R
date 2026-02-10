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
rawdata1 <- read.csv("feature_abundance.csv", header = T, na.strings = "") ###### 여기 데이터 수정
rawdata1 <-t(rawdata1)
colnames(rawdata1) <- rawdata1[3, ]
rawdata1 <- rawdata1[-c(1:2,4:56),]
rawdata1<-rawdata1[-c(1,nrow(rawdata1)),]




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
rawdata_merged$group <- factor(substr(rownames(rawdata_merged), 1, 1))
rawdata_merged<-as.data.frame(rawdata_merged)
write.csv(rawdata_merged, "data_merge.csv", row.names = FALSE)

dim(rawdata_merged)


#############################
########## Fitting ##########
#############################

# Mantel test

### HUGE, 모든 변수
set.seed(1)
X <- scale(as.matrix(rawdata_merged[,-NCOL(rawdata_merged)]))       
hu  <- huge(X, method="glasso", lambda = seq(2,0.05,len=30))
sel <- huge.select(hu, criterion="ebic")
Theta <- sel$opt.icov


# ### HUGE, 모든 변수, 논파라노말
# set.seed(1)
# X2 <- scale(as.matrix(rawdata_merged[,-NCOL(rawdata_merged)]))       
# Z2 <- huge.npn(X2)
# hu2 <- huge(Z2, method = "glasso", lambda = seq(2,0.05,len=30))
# sel2 <- huge.select(hu2, criterion = "ebic")
# Theta2 <- sel2$opt.icov


`%||%` <- function(a,b) if (!is.null(a)) a else b

## ======================= 라벨 보정 ==========================
robust_labels <- function(Theta, labels = NULL, df = NULL) {
  p <- nrow(Theta)
  if (!is.null(labels) && length(labels) == p) return(labels)
  dn <- colnames(Theta)
  if (!is.null(dn) && length(dn) == p) return(dn)
  if (!is.null(df)) {
    cn <- colnames(df)
    if (!is.null(cn)) {
      if (length(cn) - 1 == p) return(cn[-length(cn)])
      if (length(cn) == p)     return(cn)
    }
  }
  warning(sprintf("[labels] 길이 불일치: p=%d, fallback으로 V1..Vp를 사용합니다.", p))
  paste0("V", seq_len(p))
}


################## 정밀도행렬 -> 부분상관행렬 ######################
partial_corr_from_precision <- function(Omega, eps = 1e-12) {
  stopifnot(nrow(Omega) == ncol(Omega))
  d <- diag(Omega); d[d < eps] <- eps
  D <- sqrt(d)
  Pc <- -Omega / (D %o% D)
  diag(Pc) <- 1
  Pc
}

################ 가우시안 그래피컬 모델 인접행렬, 가중인접행렬 ###########
W1 <- partial_corr_from_precision(Theta)
weighted_adj1  <- abs(W1); diag(weighted_adj1) <- 0
adj1 <- ifelse(weighted_adj1>0, 1, 0)


################ 논파라노말 그래피컬 모델 인접행렬, 가중인접행렬 ###########
W2 <- partial_corr_from_precision(Theta2)
weighted_adj2  <- abs(W2); diag(weighted_adj2) <- 0
adj2 <- ifelse(weighted_adj2>0, 1, 0)

save(weighted_adj1, weighted_adj2, adj1, adj2, file = "adjmatrix.RData")
Colnames<-colnames(rawdata_merged[-NCOL(rawdata_merged)])
Rownames<-rownames(rawdata_merged)
save(Colnames,Rownames, file = "meta.RData")

## ====== 허브니스 덴시티 2x3 패널 ======
hubness_density_only_base <- function(Theta,
                                      labels = NULL, df = NULL,
                                      focus_nodes,
                                      metrics = c("degree","strength","eigencent","betweenness","closeness"),
                                      palette = c("red3","blue3","darkgreen"),
                                      label_cex = 1.3,
                                      percent_cex = 1.1,   # ← 레전드 아래 상위% 글자 크기
                                      grid = c(2, 3),      # ← 5개 지표에 적합
                                      legend_pos = "topright") {
  
  ## ---- 헬퍼 ----
  robust_labels <- function(Theta, labels = NULL, df = NULL) {
    p <- nrow(Theta)
    if (!is.null(labels) && length(labels) == p) return(labels)
    dn <- colnames(Theta); if (!is.null(dn) && length(dn) == p) return(dn)
    if (!is.null(df)) {
      cn <- colnames(df)
      if (!is.null(cn)) {
        if (length(cn) - 1 == p) return(cn[-length(cn)])
        if (length(cn) == p)     return(cn)
      }
    }
    warning(sprintf("[labels] 길이 불일치: p=%d → V1..Vp로 대체합니다.", p))
    paste0("V", seq_len(p))
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
  
  compute_hubness_core <- function(W_abs, labels) {
    W_abs <- as.matrix(W_abs); diag(W_abs) <- 0; W_abs <- 0.5*(W_abs + t(W_abs))
    g <- igraph::graph_from_adjacency_matrix(W_abs, mode = "undirected", weighted = TRUE, diag = FALSE)
    igraph::V(g)$name <- labels
    w <- igraph::E(g)$weight; eps <- 1e-8; igraph::E(g)$dist <- 1/pmax(w, eps)
    data.frame(
      node        = igraph::V(g)$name,
      degree      = as.numeric(igraph::degree(g)),
      strength    = as.numeric(igraph::strength(g, weights = w)),
      eigencent   = as.numeric(igraph::eigen_centrality(g, directed = FALSE, weights = w)$vector),
      betweenness = as.numeric(igraph::betweenness(g, weights = igraph::E(g)$dist, normalized = TRUE)),
      closeness   = closeness_std_or_harmonic(g),   # ← 여기! g를 넘김
      stringsAsFactors = FALSE
    )
  }
  add_top_percent <- function(df, ms) {
    n <- nrow(df)
    for (m in ms) { r <- rank(-df[[m]], ties.method="average"); df[[paste0(m,"_top_percent")]] <- 100*(n - r + 1)/n }
    df
  }
  
  ## ---- 준비 ----
  labs <- robust_labels(Theta, labels, df)
  Pc <- partial_corr_from_precision(Theta)
  W  <- abs(Pc); diag(W) <- 0
  stats <- compute_hubness_core(W, labels = labs)
  valid_metrics <- metrics[metrics %in% names(stats)]
  if (!length(valid_metrics)) stop("표시할 지표가 없습니다.")
  stats <- add_top_percent(stats, valid_metrics)
  
  present <- intersect(focus_nodes, stats$node)
  if (!length(present)) warning("관심 노드를 stats에서 찾지 못했습니다.")
  cols <- setNames(rep(palette, length.out = length(present)), present)
  
  ## ---- 패널 레이아웃 ----
  op <- par(no.readonly = TRUE); on.exit(par(op))
  par(mfrow = grid, mar = c(4, 4, 3.5, 1))  # x,y,top,right
  
  ## ---- 패널 루프 ----
  for (m in valid_metrics) {
    vals <- stats[[m]]; vals <- vals[is.finite(vals)]
    can_density <- length(unique(vals)) >= 2
    if (can_density) {
      d <- try(stats::density(vals), silent = TRUE)
      if (inherits(d,"try-error") || any(!is.finite(d$y))) can_density <- FALSE
    }
    
    if (can_density) {
      graphics::plot(d, main = paste("Hubness density —", m), xlab = m, ylab = "", lwd = 2)
      graphics::rug(vals, col = "grey60")
      ymax <- max(d$y, na.rm = TRUE)
    } else {
      graphics::hist(vals, breaks = "FD", col = "grey85", border = "white",
                     main = paste("Hubness histogram —", m), xlab = m, ylab = "")
      usr <- par("usr"); ymax <- usr[4]
    }
    
    ## 점선만(상단 라벨 제거)
    for (fn in present) {
      v <- stats[stats$node == fn, m]
      if (length(v) && is.finite(v)) abline(v = v, lty = 2, col = cols[[fn]])
    }
    
    ## 레전드 (노드명만)
    lg <- legend(legend_pos, legend = present, col = cols[present], lty = 2, bty = "n",
                 cex = 0.95, title = "focus")
    
    ## 레전드 바로 아래에 상위% 표기
    #  - legend() 반환 rect 좌표를 사용해 바로 아래 줄부터 라인별로 출력
    if (length(present)) {
      usr <- par("usr")
      line_h <- strheight("A") * 1.1
      start_x <- lg$rect$left
      start_y <- lg$rect$top - lg$rect$h - line_h * 0.6
      # 지표별 퍼센트 라인
      percent_lines <- vapply(present, function(fn) {
        tp <- stats[stats$node == fn, paste0(m, "_top_percent")]
        sprintf("%s:  %.1f%%", fn, tp)
      }, FUN.VALUE = character(1))
      # 줄바꿈으로 한 번에 출력 (왼쪽 정렬)
      graphics::text(x = start_x, y = start_y,
                     labels = paste(percent_lines, collapse = "\n"),
                     adj = c(0,1), cex = percent_cex, col = "black", xpd = NA)
    }
  }
  
  ## 2x3 그리드에서 남는 칸 채우기
  total_slots <- prod(grid)
  extra <- total_slots - length(valid_metrics)
  if (extra > 0) for (i in seq_len(extra)) plot.new()
  
  invisible(list(
    stats = stats,
    focus_summary = if (length(present))
      stats[stats$node %in% present, c("node", paste0(valid_metrics, "_top_percent")), drop = FALSE] else NULL
  ))
}


####### 결과 ##########


focus_nodes <-colnames(rawdata_merged)[1:3]


res_den_gauss <- hubness_density_only_base(
  Theta = Theta,
  df = if (exists("rawdata_merged")) rawdata_merged else NULL,
  focus_nodes = focus_nodes,
  metrics = c("degree","strength","eigencent","betweenness","closeness"),
  label_cex = 1.3,
  percent_cex = 1.1,
  grid = c(2, 3),         
  legend_pos = "topright" 
)


res_den_gauss2 <- hubness_density_only_base(
  Theta = Theta2,
  df = if (exists("rawdata_merged")) rawdata_merged else NULL,
  focus_nodes = focus_nodes,
  metrics = c("degree","strength","eigencent","betweenness","closeness"),
  label_cex = 1.3,
  percent_cex = 1.1,
  grid = c(2, 3),         
  legend_pos = "topright"  
)
