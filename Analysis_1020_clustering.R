# ============================================================
# Leiden clustering (target ~ 50 clusters) from dissimilarity CSV
# data path rule: ./data/feature_dreams_dissimilarity_matrix.csv
# ============================================================

suppressPackageStartupMessages({
  library(igraph)   # cluster_leiden (>= 1.3 권장)
  library(Matrix)   # nearPD
  library(readr)    # read_csv
})

## ---------- 0) 입력 ----------
infile <- file.path("data", "feature_dreams_dissimilarity_matrix.csv")

# RStudio에서만 WD를 스크립트 폴더로 바꾸고 싶다면(선택):
if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
  if (rstudioapi::isAvailable()) {
    try({
      setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
    }, silent = TRUE)
  }
}

D_raw <- suppressMessages(read_csv(infile, show_col_types = FALSE))
if (!is.null(D_raw[[1]]) && !any(duplicated(D_raw[[1]]))) {
  row_ids <- D_raw[[1]]
  D <- as.matrix(D_raw[,-1, drop=FALSE])
  colnames(D) <- names(D_raw)[-1]
  rownames(D) <- row_ids
} else {
  D <- as.matrix(D_raw)
}
stopifnot(nrow(D) == ncol(D), all(rownames(D) == colnames(D)))
p <- nrow(D)
message(sprintf("Loaded dissimilarity matrix: p = %d features.", p))

## ---------- 1) 유사도 S (0~1 가정: S = 1 - D) ----------
S0 <- 1 - D
S0 <- (S0 + t(S0))/2
diag(S0) <- 1
S0[S0 < 0] <- 0

## ---------- 2) PSD 보정 ----------
# nearPD (느리면 eigen 클리핑으로 대체 가능)
S_psd <- as.matrix(nearPD(S0, corr = FALSE)$mat)
# (S_psd 계산 뒤에 바로 추가)
q <- quantile(S_psd, 0.995)                 # 상위 0.5% 클립 (필요 시 0.99로)
S_psd <- pmin(S_psd, q)
mn <- min(S_psd); mx <- max(S_psd)
if (mx > mn) S_psd <- (S_psd - mn) / (mx - mn)
diag(S_psd) <- 1
## ---------- 3) 그래프 구성: union kNN (+ 선택 SNN) ----------
# 기존: mutual kNN + which(A==1)로 W@x 대입  → 순서 위험
# 교체: union kNN + summary(A) 인덱스로 안전하게 대입

k <- 50
p <- nrow(S_psd)
nn <- t(apply(D + diag(Inf, p), 1, function(x) order(x)[1:k]))

A <- Matrix(0, p, p, sparse = TRUE, dimnames = list(rownames(D), colnames(D)))
for (i in 1:p) A[i, nn[i,]] <- 1
A <- ((A + t(A)) > 0) * 1        # ★ union kNN (합집합)

W  <- A
sA <- summary(A)                 # (i, j, x)
W@x <- S_psd[cbind(sA$i, sA$j)]  # ★ 안전한 순서로 가중치 주입

g <- graph_from_adjacency_matrix(W, mode="undirected", weighted=TRUE, diag=FALSE)

message(sprintf("Graph built: |V|=%d, |E|=%d (union kNN, k=%d).", gorder(g), gsize(g), k))

## ---------- 4) 해상도 자동 튜닝(목표 군집수 ≈ 50) ----------
set.seed(2025)
leiden_tune <- function(g, target=50, low=0.02, high=0.60, maxit=16) {
  best <- NULL
  for (it in 1:maxit) {
    res <- (low + high)/2
    cl  <- cluster_leiden(g, resolution_parameter=res, weights=E(g)$weight)
    kf  <- length(unique(membership(cl)))
    if (is.null(best) || abs(kf - target) < abs(best$k - target)) {
      best <- list(res=res, cl=cl, k=kf)
    }
    if (kf > target) high <- res else if (kf < target) low <- res else break
  }
  best
}

set.seed(2025)
best <- leiden_tune(g, target=50, low=0.02, high=0.60, maxit=16)
final_cl  <- best$cl
final_k   <- best$k
final_res <- best$res
message(sprintf("Leiden tuned: resolution=%.4f -> %d clusters.", final_res, final_k))

## ---------- 5) 결과 정리 ----------
memb <- membership(final_cl)           # 1) 초기 모듈 배정
modules <- split(names(memb), memb)    # 2) 모듈 리스트

## ⬇️⬇️⬇️ 여기 사이에 '작은 모듈 병합 규칙'을 넣는다 ⬇️⬇️⬇️
# (S_psd 필요: feature×feature PSD 유사도 행렬)

# ---- 작은 모듈 병합 규칙 (size < min_size) ----
min_size <- 3   # 원하는 최소 모듈 크기
changed <- TRUE
while (changed) {
  changed <- FALSE
  tab <- table(memb)
  small_ids <- names(tab)[tab < min_size]
  if (length(small_ids) == 0) break
  
  for (sid in small_ids) {
    feats <- modules[[sid]]
    if (length(feats) == 0) next
    
    # 후보: 자신(sid) 제외 모든 모듈
    cand <- setdiff(names(modules), sid)
    if (length(cand) == 0) next
    
    # 각 후보 모듈과의 평균 유사도 계산 (가장 가까운 모듈을 찾음)
    best_to <- NA; best_s <- -Inf
    for (to in cand) {
      # S_psd는 PSD 유사도 행렬이어야 함
      s <- mean(S_psd[feats, modules[[to]], drop = FALSE])
      if (is.finite(s) && s > best_s) { best_s <- s; best_to <- to }
    }
    
    # 병합 수행
    if (!is.na(best_to)) {
      memb[feats] <- as.integer(best_to)
      changed <- TRUE
    }
  }
  
  # memb가 바뀌었으면 modules 갱신
  if (changed) modules <- split(names(memb), memb)
}
# -----------------------------------------------

# 3) 병합 후 결과 테이블 생성 (최종 memb 기준)
res_tbl <- data.frame(
  feature_id = names(memb),
  module_id  = as.integer(memb),
  stringsAsFactors = FALSE
)
module_summary <- do.call(rbind, lapply(names(modules), function(id) {
  feats <- modules[[id]]
  Sg <- S_psd[feats, feats, drop = FALSE]
  n  <- length(feats)
  mean_intra <- if (n > 1) mean(Sg[upper.tri(Sg, diag = FALSE)]) else NA_real_
  data.frame(module_id = as.integer(id), size = n, mean_intra_similarity = mean_intra)
}))
module_summary <- module_summary[order(module_summary$module_id), ]

dir.create("results", showWarnings = FALSE)
write.csv(res_tbl[order(res_tbl$module_id, res_tbl$feature_id), ],
          file.path("results", "leiden_modules_mapping.csv"), row.names = FALSE)
write.csv(module_summary,
          file.path("results", "leiden_modules_summary.csv"), row.names = FALSE)

message("Saved: results/leiden_modules_mapping.csv, results/leiden_modules_summary.csv")

## ---------- 6) 다음 단계 힌트 ----------
# modules 리스트를 사용해 varComp용 K_g를 생성:
# idx <- modules[[g]]; Xg <- X[, idx, drop=FALSE]; Sg <- S_psd[idx, idx, drop=FALSE]
# Kg <- Xg %*% Sg %*% t(Xg) / length(idx)
dir.create("results", showWarnings = FALSE)

pkgs <- c("ggplot2", "Matrix", "igraph", "scales")
opt  <- c("uwot", "ComplexHeatmap", "circlize", "pheatmap", "cluster")
to_install <- setdiff(c(pkgs, opt), rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)

suppressPackageStartupMessages({
  library(ggplot2)
  library(igraph)
  library(scales)
  library(Matrix)
})

if (requireNamespace("uwot", quietly = TRUE)) {
  library(uwot)
  set.seed(2025)
  # 특징-특징 유사도 S_psd를 임베딩 입력으로 사용(대칭/PSD)
  emb <- umap(as.matrix(S_psd), n_neighbors = 30, metric = "cosine", init = "spectral")
  df_umap <- data.frame(UMAP1 = emb[,1], UMAP2 = emb[,2],
                        module = factor(memb[rownames(S_psd)]))
  p2 <- ggplot(df_umap, aes(UMAP1, UMAP2, color = module)) +
    geom_point(size = 0.8, alpha = 0.8) +
    guides(color = "none") +
    labs(title = "UMAP on S_psd (features), colored by module") +
    theme_minimal()
  ggsave("results/vis_umap_spsd.png", p2, width = 7, height = 6, dpi = 250)
}
# 모듈 순서: 크기 큰 순으로
mod_order <- names(sort(table(memb), decreasing = TRUE))
# 재정렬된 피처 순서
feat_order <- unlist(modules[mod_order], use.names = FALSE)

S_re <- S_psd[feat_order, feat_order, drop = FALSE]

# ComplexHeatmap가 있으면 예쁘게
if (requireNamespace("ComplexHeatmap", quietly = TRUE) &&
    requireNamespace("circlize", quietly = TRUE)) {
  library(ComplexHeatmap); library(circlize)
  pal <- colorRamp2(c(0, 0.5, 1), c("#0d0887","#f0f921","#d43d51"))
  png("results/vis_heatmap_Spsd.png", width=1200, height=1000, res=180)
  Heatmap(S_re, name = "S", col = pal, show_row_names = FALSE, show_column_names = FALSE,
          cluster_rows = FALSE, cluster_columns = FALSE,
          column_split = factor(rep(mod_order, lengths(modules[mod_order]))),
          row_split    = factor(rep(mod_order, lengths(modules[mod_order]))))
  dev.off()
} else {
  # fallback: pheatmap
  if (requireNamespace("pheatmap", quietly = TRUE)) {
    png("results/vis_heatmap_Spsd.png", width=1200, height=1000, res=180)
    pheatmap::pheatmap(S_re, show_rownames = FALSE, show_colnames = FALSE,
                       cluster_rows = FALSE, cluster_cols = FALSE)
    dev.off()
  }
}

mod_size <- sort(table(memb), decreasing = TRUE)
df_size <- data.frame(module = factor(names(mod_size), levels = names(mod_size)),
                      size = as.integer(mod_size))
p3 <- ggplot(df_size, aes(x = module, y = size)) +
  geom_col() +
  coord_flip() +
  geom_hline(yintercept = 3, linetype = 2, color = "red") +
  labs(title = "Module size distribution (red dashed = size 3)",
       x = "Module", y = "Size") +
  theme_minimal()
ggsave("results/vis_module_sizes.png", p3, width = 6, height = 8, dpi = 250)

