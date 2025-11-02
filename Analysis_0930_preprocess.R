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
library(MiRKAT)
library(stringr)
#########################
####### Load data #######
#########################
data_dir <- file.path(".", "data")

## 전처리 대상 파일 경로 지정
fa_path  <- file.path(data_dir, "feature_abundance.csv")        # 특성별 풍부도 테이블
meta_path <- file.path(data_dir, "metadata_hormone_herbivore.csv") # 메타데이터 테이블

fa   <- read.csv(fa_path,   header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
meta <- read.csv(meta_path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

## 필요한 열 존재 확인
## stopifnot("메시지" = 조건): 조건이 TRUE가 아니면 에러
stopifnot("'id.x' 컬럼이 feature_abundance.csv에 없습니다." = "id.x" %in% names(fa))

## -------------------------
## (A) feature_abundance 정리
## -------------------------

## 행이름을 feature 식별자(id.x)로 부여
## rownames(x) <- 벡터 : x의 행이름을 지정
rownames(fa) <- fa[["id.x"]]

## 샘플 컬럼 패턴 선택: 'datafile.<무엇이든>.area'
## grepl(pattern, x): x 원소가 pattern(정규식)에 매칭 여부 반환
sample_col_idx <- grepl("^datafile\\..*\\.area$", names(fa))

## 단일 열 'datafile.con_1.mzML.area'가 없으면 경고 (패턴 열들로 계속 진행)
if (!"datafile.con_1.mzML.area" %in% names(fa)) {
  warning("feature_abundance.csv에 'datafile.con_1.mzML.area' 컬럼이 없습니다. 패턴 매칭된 컬럼들로 진행합니다.")
}

## 컬럼명 정리: 앞 'datafile.' 제거, 뒤 '.area' 제거
## sub(pattern, replacement, x): 첫 매칭 부분 치환
clean_names <- names(fa)[sample_col_idx]
clean_names <- sub("^datafile\\.", "", clean_names)  # 접두사 제거
clean_names <- sub("\\.area$", "", clean_names)      # 접미사 제거

## 샘플 열만 추출
## - [, sample_col_idx, drop=FALSE]: 논리 인덱싱으로 해당 열만 보존, 1열이어도 데이터프레임 유지
fa_samples <- fa[, sample_col_idx, drop = FALSE]
## 이름 정리 적용
## names(x) <- 벡터 : 컬럼명 지정
names(fa_samples) <- clean_names

## 전치하여 (행=sample, 열=feature) 구조로 변경
## t(x): 전치, as.matrix로 수치형 행렬 보장
fa_t <- t(as.matrix(fa_samples))

## 전치 결과를 데이터프레임으로 변환
## data.frame(x, check.names=FALSE): 자동명수정 방지
fa_t_df <- data.frame(fa_t, check.names = FALSE)

## 병합키로 쓸 샘플 ID를 열로 승격
## cbind(...): 열 결합, rownames(fa_t_df): 현재 샘플 ID
fa_t_df <- cbind(sample = rownames(fa_t_df), fa_t_df)
rownames(fa_t_df) <- NULL

## -------------------------
## (B) metadata 정리
## -------------------------

## 메타데이터의 샘플 ID 컬럼명 지정 (사용자 제공: "sample")
meta_id_col <- "sample"

## 해당 컬럼 존재 체크
if (!meta_id_col %in% names(meta)) {
  stop(sprintf("메타데이터에 '%s' 컬럼이 없습니다.", meta_id_col))
}

## 필요시 메타의 샘플 ID 형식 정규화(주석 해제해서 사용)
## 예: 메타가 'datafile.xxx.area' 형태라면 동일하게 정리
## meta[[meta_id_col]] <- meta[[meta_id_col]] |>
##   sub("^datafile\\.", "", x=_) |>
##   sub("\\.area$", "", x=_)

## -------------------------
## (C) 병합 (행이 일치하는 경우만)
## -------------------------

## merge(x, y, by, all, sort):
## - x=meta, y=fa_t_df
## - by="sample": 공통 키
## - all=FALSE: 내부조인(겹치는 샘플만)
## - sort=FALSE: 자동 정렬 방지
merged_df <- merge(meta, fa_t_df, by = "sample", all = FALSE, sort = FALSE)

## -------------------------
## (D) 로그/저장
## -------------------------

## message(sprintf(...)): 정보 메시지 출력
message(sprintf("전치 전 feature 수: %d", nrow(fa)))
message(sprintf("전치 후 sample 수 (fa_t_df): %d", nrow(fa_t_df)))
message(sprintf("병합 후 행(샘플) 수: %d", nrow(merged_df)))
message(sprintf("병합 후 열(메타+특성) 수: %d", ncol(merged_df)))

## write.csv(x, file, row.names=FALSE): CSV 저장, 행이름 미포함
out_path <- file.path(data_dir, "merged_metadata_abundance.csv")
write.csv(merged_df, out_path, row.names = FALSE)
message(sprintf("병합 결과를 저장했습니다: %s", out_path))

####
## =========================
## 컬럼 추출 재정의
## =========================

## 제거할 비분석(메타) 컬럼 리스트 정의
## c(...): 벡터 결합
## - 'sl'과 'ls'를 모두 포함해 오타/표기 차이 대응
drop_meta_cols <- c("sl", "ls", "px", "mp", "le", "JA", "SA", "ABA")

## 우선 보존할 핵심 메타 컬럼
keep_meta_cols <- c("sample", "group", "plant")

## names(x): 데이터프레임의 전체 컬럼명 벡터
all_cols <- names(merged_df)

## setdiff(x, y): 집합 차집합 (x에서 y를 뺀 결과)
## - 분석 대상 피처 = 전체 - (핵심 메타 + 제거할 메타)
feature_cols <- setdiff(all_cols, c(keep_meta_cols, drop_meta_cols))

## intersect(x, y): 교집합 (실제 존재하는 컬럼만 사용하도록 안전 장치)
keep_meta_cols <- intersect(keep_meta_cols, all_cols)
feature_cols   <- intersect(feature_cols,   all_cols)

## 데이터 서브셋 만들기
## merged_df[, j, drop=FALSE]:
## - j: 남길 컬럼명 벡터(메타 3개 + 피처 전부)
## - drop=FALSE: 결과가 1열이어도 데이터프레임 유지
df_sel <- merged_df[, c(keep_meta_cols, feature_cols), drop = FALSE]

## 메타/피처 분리
meta_df <- df_sel[, keep_meta_cols, drop = FALSE]
X_df    <- df_sel[, feature_cols,   drop = FALSE]

## =========================
## 피처를 숫자형으로 정리
## =========================

## lapply(X, FUN): 각 열에 함수 적용
## - suppressWarnings(): 경고 억제(숫자 변환 시 NA 생성 경고 숨김)
## - as.numeric(v): 문자→숫자 변환 시도 (변환 불가 값은 NA)
X_df <- data.frame(
  lapply(X_df, function(v) suppressWarnings(as.numeric(v))),
  check.names = FALSE  ## check.names=FALSE: 원래 컬럼명(숫자/특수문자 포함) 보존
)

## =========================
## 점검 출력 (선택)
## =========================

## message(...): 정보 메시지 출력
message(sprintf("최종 보존 메타 컬럼: %s", paste(keep_meta_cols, collapse = ", ")))
message(sprintf("제거된 메타 후보: %s", paste(intersect(drop_meta_cols, all_cols), collapse = ", ")))
message(sprintf("선택된 피처 수: %d", length(feature_cols)))

## head(x): 앞부분 미리보기
## - n=3: 3행만
print(head(df_sel, n = 3))

## 필요 시 현재 컬럼명 확인
## colnames(x): 데이터프레임 컬럼명 반환
print(colnames(df_sel))
write.csv(df_sel, out_path, row.names = FALSE)


# ---------- 0) 패키지 확인/설치 ----------
if (!requireNamespace("varComp", quietly = TRUE)) {
  message("Installing 'varComp' from R-Forge ...")
  install.packages("varComp", repos = "https://R-Forge.R-project.org")
}
library(varComp)

dir.create("results", showWarnings = FALSE)

# ---------- 1) 로드 ----------
X_std   <- readRDS(file.path("results", "preproc_X.rds"))          # n x p, CTRL기준 표준화된 X
df      <- readRDS(file.path("results", "preproc_y_df.rds"))       # data.frame: sample, group, plant
modules <- readRDS(file.path("results", "preproc_modules.rds"))    # list: id -> features
Klist   <- readRDS(file.path("results", "preproc_Klist.rds"))      # list: M{id} -> n x n

# ---------- 2) 기본 점검 ----------
# group 레벨 고정 (참조수준 ctrl)
df$group <- factor(df$group, levels = c("ctrl","le","mp","px","sl"))
stopifnot(all(rownames(X_std) == rownames(df)))
rownames(df) <- df$sample

stopifnot(all(rownames(X_std) == df$sample))

# Klist의 행/열 순서가 df와 같은지 점검
chk <- vapply(Klist, function(K) {
  all(rownames(K) == rownames(df)) && all(colnames(K) == rownames(df))
}, logical(1))
if (!all(chk)) stop("Klist의 행/열 이름이 df의 샘플 순서와 일치해야 합니다.")

# ---------- 3) 적합 ----------
# 고정효과: ctrl을 기준으로 group 조정 (필요 공변량 있으면 +batch +day 등 추가)
fit <- varComp(
  fixed = plant ~ group,
  data  = df,
  varcov = Klist_pruned,
  normalizeTrace = TRUE
)

keep_mod_ids = 1:15
Klist_2   <-  Klist[paste0("M", keep_mod_ids)]

# 요약 저장(고정효과, 분산성분, sigma2 등)
sink(file.path("results", "varcomp_fit_summary.txt"))
cat("=== varComp fit summary ===\n")
print(summary(fit))
cat("\n--- variance components (random only) ---\n")
print(coef(fit, "varComp"))   # tau_hat (각 모듈의 분산성분)
cat("\n--- sigma^2 (residual) ---\n")
print(fit$sigma2)
sink()

# ---------- 4) 전역 검정 (모든 모듈 동시 0?) ----------
# 어떤 깊이의 리스트든 htest를 찾아 p.value를 수집
extract_p <- function(x) {
  if (inherits(x, "htest")) {
    return(x$p.value)
  } else if (is.list(x)) {
    return(unlist(lapply(x, extract_p_any), use.names = FALSE))
  } else {
    return(numeric(0))
  }
}
View(p_each)
glob <- varComp.test(fit, null = integer(0), test = "LinScore")
View(glob)
p_global <- extract_p(glob)[1]
writeLines(sprintf("global_p = %.6g", p_global),
           file.path("results", "varcomp_global_test.txt"))

# ---------- 5) 모듈(클러스터)별 조건부 검정 ----------
G <- length(Klist_2)
p_each <- setNames(numeric(G), names(Klist_2))

for (g in seq_len(G)) {
  null_idx <- setdiff(seq_len(G), g)         # 나머지 모듈은 귀무에 포함
  tg <- varComp.test(fit, null = null_idx,   # 모듈 g만 0인지 조건부 검정
                     test = "LinScore")
  p_each[g] <- extract_p(tg)[1]
}

# FDR (BH) 보정
q_each <- p.adjust(p_each, method = "BH")

# 분산성분 추정치(모듈별 tau_hat)와 분산 기여 비율(PV 비슷한 감각지표)
tau_hat <- coef(fit, "varComp")   # 이름이 Klist 순서와 같아야 함
# tau_hat의 이름이 Klist와 어긋나면 Klist 순서에 맞춰 재정렬
if (!all(names(tau_hat) == names(Klist_2))) {
  tau_hat <- tau_hat[names(Klist_2)]
}

sigma2  <- fit$sigma2
pv_ratio <- tau_hat / (sigma2 + sum(tau_hat))

# 모듈 크기(피처 수) 첨부
mod_size <- vapply(names(Klist_2), function(nm) {
  id <- sub("^M", "", nm)
  length(modules[[id]])
}, integer(1))

# 결과 테이블
res_tbl <- data.frame(
  module_id = sub("^M", "", names(Klist_2)),
  p_value   = as.numeric(p_each),
  q_value   = as.numeric(q_each),
  tau_hat   = as.numeric(tau_hat),
  pv_ratio  = as.numeric(pv_ratio),
  size      = mod_size,
  stringsAsFactors = FALSE
)
res_tbl <- res_tbl[order(res_tbl$q_value, res_tbl$p_value), ]

write.csv(res_tbl, file.path("results", "varcomp_module_tests.csv"), row.names = FALSE)

cat("\n=== Done ===\n",
    sprintf("Global test p = %.6g\n", p_global),
    sprintf("Per-module results saved to: %s\n", file.path("results","varcomp_module_tests.csv")))

################marginal######
# 귀무(랜덤효과 없음) 적합
fit0 <- varComp(fixed = plant ~ group, data = df, varcov = list(), normalizeTrace = TRUE)

# 각 모듈을 단독으로 추가해서 테스트 (added-first)
p_marg <- setNames(numeric(length(Klist)), names(Klist))
for (g in seq_along(Klist)) {
  out <- varComp.test(object = fit0, additional.varcov = list(Klist[[g]]), test = "LinScore")
  p_marg[g] <- extract_p(out)[1]
}
q_marg <- p.adjust(p_marg, "BH")
cbind(p_marg, q_marg)[order(q_marg), ][, , drop=FALSE]

###
# 1) 귀무(랜덤효과 없음) 적합
if (!exists("fit0")) {
  fit0 <- varComp(fixed = plant ~ group, data = df,
                  varcov = list(), normalizeTrace = TRUE)
}

# 2) 각 모듈을 단독으로 추가해서 marginal p 계산
G  <- length(Klist)
nm <- names(Klist)
p_marg <- setNames(numeric(G), nm)

for (g in seq_along(Klist)) {
  out <- varComp.test(object = fit0,
                      additional.varcov = list(Klist[[g]]),
                      test = "LinScore")
  pg <- extract_p(out)
  if (length(pg) == 0)
    stop(sprintf("marginal test: %s p-value 추출 실패", nm[g]))
  p_marg[g] <- pg[1]
}
q_marg <- p.adjust(p_marg, "BH")

# 3) 모듈 크기 붙이기 (modules_pruned가 있으면 우선 사용)
modules_main <- if (exists("modules_pruned")) modules_pruned else modules
mod_size <- vapply(nm, function(nmi){
  id <- sub("^M","", nmi)
  if (!is.null(modules_main[[id]])) length(modules_main[[id]]) else NA_integer_
}, integer(1))

# 4) (선택) 단독 적합으로 tau_hat_marg 추정치도 계산
tau_hat_marg <- rep(NA_real_, G)
for (g in seq_along(Klist)) {
  fit_g <- varComp(fixed = plant ~ group, data = df,
                   varcov = list(Klist[[g]]),
                   normalizeTrace = TRUE)
  # varComp 하나만 있으므로 첫 항목
  th <- as.numeric(coef(fit_g, "varComp")[1])
  tau_hat_marg[g] <- th
}

# 5) 저장 테이블 구성 (varcomp_module_tests.csv와 유사)
marg_tbl <- data.frame(
  module_id = sub("^M","", nm),
  p_value   = as.numeric(p_marg),
  q_value   = as.numeric(q_marg),
  tau_hat   = as.numeric(tau_hat_marg),  # 단독 적합 기반 추정치
  size      = mod_size,
  test_type = "marginal",                # 구분자
  stringsAsFactors = FALSE
)
marg_tbl <- marg_tbl[order(marg_tbl$q_value, marg_tbl$p_value), ]

out_csv <- file.path("results","varcomp_module_tests_marginal.csv")
write.csv(marg_tbl, out_csv, row.names = FALSE)
message("Saved marginal tests to: ", out_csv)

# (옵션) 기존 조건부 결과와 조인하여 하나의 요약 CSV로 저장
cond_csv <- file.path("results","varcomp_module_tests.csv")
if (file.exists(cond_csv)) {
  cond_tbl <- read.csv(cond_csv, stringsAsFactors = FALSE)
  names(cond_tbl)[names(cond_tbl)=="p_value"]  <- "p_value_cond"
  names(cond_tbl)[names(cond_tbl)=="q_value"]  <- "q_value_cond"
  names(cond_tbl)[names(cond_tbl)=="tau_hat"]  <- "tau_hat_cond"
  names(cond_tbl)[names(cond_tbl)=="pv_ratio"] <- "pv_ratio_cond"
  
  # 키: module_id
  merged <- merge(cond_tbl, marg_tbl, by="module_id", all=TRUE, suffixes = c("", "_marg"))
  write.csv(merged, file.path("results","varcomp_module_tests_merged.csv"), row.names = FALSE)
  message("Also wrote merged summary: results/varcomp_module_tests_merged.csv")
}







###############점검##
kern_corr <- function(Klist){
  nm <- names(Klist); G <- length(Klist)
  M <- matrix(NA, G, G, dimnames = list(nm, nm))
  for(i in 1:G) for(j in i:G){
    M[i,j] <- M[j,i] <- cor(as.vector(Klist[[i]]), as.vector(Klist[[j]]))
  }
  M
}
Rk <- kern_corr(Klist)



###################시각화###########################

## ---------- 2) 전역 p-value 출력 ----------
glob_txt <- readLines(file.path("results","varcomp_global_test.txt"))[1]
message(glob_txt)

## ---------- 3) 요약 테이블 가공 ----------
# 모듈 ID는 숫자 형태로 보존
tbl <- tbl %>%
  mutate(module_id = as.integer(module_id),
         neglog10p = -log10(p_value + 1e-300),
         neglog10q = -log10(q_value + 1e-300),
         sig = q_value < alpha)

# 모듈 크기가 없으면 0 처리
if (!"size" %in% names(tbl)) tbl$size <- NA_integer_

## ---------- 4) 요약 그림 1: -log10(p) vs. module size (점 크기=pv_ratio, 색=유의/비유의) ----------
p1 <- ggplot(tbl, aes(x = size, y = neglog10p)) +
  geom_point(aes(size = pmax(pv_ratio, 1e-6), color = sig), alpha = 0.7) +
  scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "firebrick")) +
  scale_size_continuous(name = "pv_ratio", guide = "legend") +
  labs(x = "모듈 크기(# features)", y = expression(-log[10](p)),
       title = "VarComp 모듈 검정 요약: -log10(p) vs. 모듈 크기") +
  theme_minimal(base_size = 12)
ggsave(file.path("results","plot_module_p_vs_size.png"), p1, width = 7, height = 5, dpi = 150)

## ---------- 5) 요약 그림 2: 유의 모듈 상위 20개(또는 sig 전부) 막대(neglog10q) ----------
tbl_sig <- tbl %>% filter(sig) %>% arrange(q_value, p_value)
topN <- min(20, nrow(tbl_sig))
if (topN >= 1) {
  p2 <- ggplot(tbl_sig[1:topN, ],
               aes(x = reorder(paste0("M", module_id), neglog10q), y = neglog10q)) +
    geom_col() +
    coord_flip() +
    labs(x = "모듈 ID", y = expression(-log[10](q)),
         title = sprintf("유의 모듈 Top %d (FDR < %.2f)", topN, alpha),
         subtitle = "막대 길이: -log10(q), 숨은 정보: pv_ratio/size는 CSV 참고") +
    theme_minimal(base_size = 12)
  ggsave(file.path("results","plot_top_sig_modules.png"), p2, width = 6, height = 6, dpi = 150)
} else {
  message(sprintf("FDR < %.2f 유의 모듈이 없습니다.", alpha))
}

## ---------- 6) 유의 모듈의 상위 피처 리스트 산출 ----------
# 아이디어: 검정의 Q_g = r^T K_g r = z^T S_g z (z = X_g^T r)
# feature 수준 기여도: c_j = z_j * (S_g z)_j  (합하면 Q_g)
# 여기서 r은 고정효과(plant ~ group) 잔차를 사용 (P0 y 근사)

# 고정효과만 적합해 잔차 r 취득
lm0 <- lm(plant ~ group, data = df)
r <- resid(lm0)                       # length n

# 안전 장치
stopifnot(all(colnames(X_std) == rownames(S_psd)))

# 모듈별 상위 피처 계산 함수
rank_features_in_module <- function(mod_id, top_k = 20) {
  feats <- modules[[as.character(mod_id)]]
  if (is.null(feats) || length(feats) == 0) return(NULL)
  Xg <- as.matrix(X_std[, feats, drop = FALSE])      # n x p_g
  Sg <- as.matrix(S_psd[feats, feats, drop = FALSE]) # p_g x p_g
  
  # z = X_g^T r
  z <- as.numeric(t(Xg) %*% r)                       # length p_g
  Sz <- as.numeric(Sg %*% z)
  contrib <- z * Sz                                  # length p_g, sum(contrib) ~ Q_g
  
  # (선택) 피처의 within-module 연결도(합유사도)도 함께
  degree <- rowSums(Sg)
  
  out <- data.frame(
    module_id  = mod_id,
    feature_id = feats,
    z_score    = z,
    contrib    = contrib,
    within_deg = degree,
    stringsAsFactors = FALSE
  ) %>% arrange(desc(contrib))
  if (nrow(out) > top_k) out <- out[1:top_k, ]
  out
}

# 유의 모듈만 대상
sig_ids <- tbl %>% filter(sig) %>% pull(module_id)
dir.create(file.path("results","top_features_per_module"), showWarnings = FALSE)

top_list <- lapply(sig_ids, function(mid) {
  res <- rank_features_in_module(mid, top_k = top_k)
  if (!is.null(res) && nrow(res) > 0) {
    # 개별 파일도 저장
    write.csv(res,
              file.path("results","top_features_per_module",
                        sprintf("top_features_M%d.csv", mid)),
              row.names = FALSE)
  }
  res
})
top_all <- do.call(rbind, top_list)
if (!is.null(top_all) && nrow(top_all) > 0) {
  write.csv(top_all, file.path("results","top_features_by_module.csv"), row.names = FALSE)
  message(sprintf("유의 모듈 상위 피처를 저장했습니다: %s",
                  file.path("results","top_features_by_module.csv")))
} else {
  message("유의 모듈이 없거나 상위 피처 산출 결과가 비어 있습니다.")
}

cat("\n=== Summary artifacts written ===\n",
    "- results/plot_module_p_vs_size.png\n",
    if (topN>=1) "- results/plot_top_sig_modules.png\n" else "",
    "- results/top_features_by_module.csv (그리고 results/top_features_per_module/ 개별 파일들)\n", sep="")

















#################################
####################Kernel 제거 ####
# --- 1) 작은 모듈 제거 (size ≤ 10) ---
min_size <- 11L
mod_sizes <- sapply(modules, length)
keep_mod_ids <- names(modules)[mod_sizes >= min_size]      # module_id (문자)
drop_mod_ids <- setdiff(names(modules), keep_mod_ids)

modules_small_removed <- modules[keep_mod_ids]
Klist_small_removed   <- Klist[paste0("M", keep_mod_ids)]

# --- 2) jitter + center (수치 안정) ---
centerK <- function(K){
  n <- nrow(K); H <- diag(n) - matrix(1/n, n, n)
  H %*% K %*% H
}
Klist_fixed <- lapply(Klist_small_removed, function(K){
  centerK(K + diag(1e-6, nrow(K)))  # jitter로 양의정부호 근사 + 센터링
})

# --- 3) 커널 중복/고상관 제거 (r > 0.99) ---
vec <- function(K) as.vector(K)
kn <- names(Klist_fixed)
m  <- length(Klist_fixed)

# 상삼각만 계산하여 상관행렬 구성
pcc <- matrix(NA_real_, m, m, dimnames = list(kn, kn))
for (i in seq_len(m)) {
  vi <- vec(Klist_fixed[[i]])
  pcc[i,i] <- 1
  if (i < m) {
    for (j in (i+1):m) {
      r <- suppressWarnings(cor(vi, vec(Klist_fixed[[j]])))
      pcc[i,j] <- pcc[j,i] <- r
    }
  }
}

# 대표 하나만 남기기(탐욕적): r > 0.99인 커널은 선두 것만 keep
keep <- rep(TRUE, m)
for (i in seq_len(m)) if (keep[i]) {
  dup <- which(pcc[i, ] > 0.99 & keep)
  dup <- setdiff(dup, i)
  if (length(dup)) keep[dup] <- FALSE
}

Klist_pruned   <- Klist_fixed[keep]
kept_mod_ids   <- sub("^M", "", names(Klist_pruned))
modules_pruned <- modules_small_removed[kept_mod_ids]

# --- 리포트 저장 ---
dir.create("results", showWarnings = FALSE)
sink("results/prune_report.txt")
cat("=== PRUNE REPORT ===\n")
cat(sprintf("총 모듈 수(원본): %d\n", length(modules)))
cat(sprintf("작은 모듈 제거(≤10): %d 제거, %d 유지\n",
            length(drop_mod_ids), length(keep_mod_ids)))
cat(sprintf("고상관 커널 제거(r>0.99): %d -> %d\n",
            length(Klist_small_removed), length(Klist_pruned)))
cat("\n유지된 모듈 ID 예시(앞 20개):\n")
print(head(kept_mod_ids, 20))
sink()

saveRDS(modules_pruned, "results/pruned_modules.rds")
saveRDS(Klist_pruned,   "results/pruned_Klist.rds")

cat("Pruning done.\n",
    sprintf("modules: %d -> %d (size>10 & corr pruning)\n",
            length(modules), length(modules_pruned)),
    sprintf("kernels: %d -> %d\n",
            length(Klist), length(Klist_pruned)),
    "Saved: results/pruned_modules.rds, results/pruned_Klist.rds, results/prune_report.txt\n",
    sep = "")




