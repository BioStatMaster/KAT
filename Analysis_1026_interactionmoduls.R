setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

#패키지 수동설치
file.exists("varComp_0.2-0.tar.gz")
install.packages(pkgs = "varComp_0.2-0.tar.gz", dependencies = TRUE)
install.packages(pkgs = "RLRsim", dependencies = TRUE)

#피키지 호출
help(package = "varComp")
library(varComp)

### Oxide/Semiconductor example
install.packages("nlme")
library(nlme)
data(Oxide)
lmef = lme(Thickness~Source, Oxide, ~1|Lot/Wafer)

vcf = varComp(Thickness~Source, Oxide, ~Lot/Wafer)

VarCorr(lmef)
coef(vcf, type = "varComp")    # varComp 추정 분산성분(동일 값 확인)

# (C) 수동으로 K를 구성해 "불균등 분산"을 표현
# - model.matrix(~0+Lot,Oxide): Lot별 원-핫 열벡터 → tcrossprod로 군집 블록 K 만들기
k0 <- tcrossprod(model.matrix(~ 0 + Lot, Oxide))  # Lot 효과 기본 K
k1 <- tcrossprod(Oxide$Source == 1) * k0                   # Source=1일 때 Lot 분산
k2 <- tcrossprod(Oxide$Source == 2) * k0                   # Source=2일 때 Lot 분산
k3 <- tcrossprod(model.matrix(~ 0 + Lot:Wafer, Oxide))     # Lot:Wafer 효과 K

# (C-1) "선호하는" 파라미터화: Source별 Lot 분산을 K1/K2로 나눔
vcf1 <- varComp(
  fixed  = Thickness ~ Source,
  data   = Oxide,
  varcov = list(S1Lot = k1, S2Lot = k2, `Lot:Wafer` = k3)  # 이름을 주면 결과에 라벨 반영
)

# (C-2) "다른" 파라미터화: K0 + K2로 표현(동등하지만 파라미터 해석이 다를 뿐)
vcf2 <- varComp(
  fixed  = Thickness ~ Source,
  data   = Oxide,
  varcov = list(Lot = k0, S2Lot = k2, `Lot:Wafer` = k3)
)

# (C-3) "좋지 않은" 파라미터화(식별성/상관으로 결국 vcf와 동일 결과)
vcf3 <- varComp(
  fixed  = Thickness ~ Source,
  data   = Oxide,
  varcov = list(Lot = k0, S1Lot = k1, `Lot:Wafer` = k3)
)

logLik(vcf);  logLik(vcf1);  logLik(vcf2);  logLik(vcf3)  # 적합값 비교(동등성 확인)

#####Genetics
# 가짜 유전자/처치 데이터
set.seed(2340)
trt <- gl(2, 15)                                 # 2수준 요인(각 15개)
dat <- data.frame(trt = trt)
dat$SNP <- matrix(sample(0:2, 120, replace = TRUE), 30)   # 유전자 행렬(30x4)
dat$Y   <- as.numeric(trt) + rnorm(30) + dat$SNP %*% rnorm(4)

# (A) 기본: ibs(SNP) 한 개 분산성분
(vcf0 = varComp(Y~trt, dat, ~ibs(SNP)))
(vcf00 = varComp(Y~trt, dat, varcov = list(`ibs(SNP)`=IBS(dat$SNP)))) ## same as above
(vcf1 = varComp(Y~trt, dat, ~ibs(SNP):trt)) ## two variance components

######실전###############################
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr)
  library(Matrix)   # nearPD
})

## ---------- 0) 파일 경로 ----------
fp_meta <- file.path("data", "merged_metadata_abundance.csv")
fp_diss <- file.path("data", "feature_dreams_dissimilarity_matrix.csv")
fp_map  <- file.path("data", "leiden_modules_mapping.csv")

## ---------- 1) 메타+풍부도 로드 ----------
meta_raw <- suppressMessages(read_csv(fp_meta, show_col_types = FALSE))

# plant 열명 오타 대비
if (!("plant" %in% names(meta_raw))) {
  if ("plnat" %in% names(meta_raw)) {
    warning("`plant` 열을 찾지 못해 `plnat`를 `plant`로 간주합니다.")
    meta_raw <- meta_raw %>% rename(plant = plnat)
  } else stop("`plant` 열을 찾을 수 없습니다.")
}
stopifnot("sample" %in% names(meta_raw), "`group` 열이 필요합니다." = "group" %in% names(meta_raw))

# sample 문자열화 및 중복 체크
meta_raw <- meta_raw %>% mutate(sample = as.character(sample))
stopifnot(!anyDuplicated(meta_raw$sample))
rownames_df <- meta_raw$sample

## ---------- 2) group 정제/검증 ----------
allowed_levels <- c("ctrl","le","mp","px","sl")
# 공백/대소문자 변동 보정
meta_raw <- meta_raw %>%
  mutate(group = trimws(tolower(as.character(group)))) %>%
  mutate(group = ifelse(group %in% allowed_levels, group, NA_character_))

# 허용 외 레벨 제거(NA)
bad_n <- sum(is.na(meta_raw$group))
if (bad_n > 0) warning(sprintf("허용되지 않은 group 레벨 %d개 샘플 제거(NA로 마킹).", bad_n))
meta_raw <- meta_raw %>% filter(!is.na(group))

# factor 레벨 고정 (순서: ctrl, le, mp, px, sl)
meta_raw$group <- factor(meta_raw$group, levels = allowed_levels)

# 실제 분포 출력 및 저장
grp_counts <- as.data.frame(table(meta_raw$group))
names(grp_counts) <- c("group","n")
message("Observed group counts:")
print(grp_counts, row.names = FALSE)
write.csv(grp_counts, file.path("results","preproc_group_counts.csv"), row.names = FALSE)

## ---------- 3) 어떤 열이 대사체(숫자 ID)인가 ----------
meta_cols <- c("sample","group","plant")
feat_cols <- setdiff(names(meta_raw), meta_cols)

# 숫자형 강제 변환
meta_num <- meta_raw
meta_num[feat_cols] <- lapply(meta_num[feat_cols], as.numeric)

# 반응/공변량 df
df <- meta_num %>% select(all_of(meta_cols))
rownames(df) <- df$sample
y <- df$plant

# 원시 X (n x p)
X_raw <- as.matrix(meta_num[, feat_cols, drop = FALSE])
rownames(X_raw) <- df$sample
stopifnot(!anyDuplicated(colnames(X_raw)))

## ---------- 4) 비유사도 행렬 D 로드 ----------
D_raw <- suppressMessages(read_csv(fp_diss, show_col_types = FALSE))
if (!is.null(D_raw[[1]]) && !anyDuplicated(D_raw[[1]])) {
  row_ids <- as.character(D_raw[[1]])
  D <- as.matrix(D_raw[,-1, drop=FALSE])
  colnames(D) <- names(D_raw)[-1]
  rownames(D) <- row_ids
} else {
  D <- as.matrix(D_raw)
}
stopifnot(nrow(D) == ncol(D))
stopifnot(all(rownames(D) == colnames(D)))

## ---------- 5) 모듈 매핑 로드 ----------
map <- suppressMessages(read_csv(fp_map, show_col_types = FALSE))
stopifnot(all(c("feature_id","module_id") %in% names(map)))
map <- map %>% mutate(feature_id = as.character(feature_id),
                      module_id  = as.integer(module_id))

## ---------- 6) 피처 ID 공통 교집합 ----------
features_common <- Reduce(intersect, list(colnames(X_raw), rownames(D), map$feature_id))
if (length(features_common) < 10)
  stop("공통 feature가 너무 적습니다. 열 이름/매핑/행열 이름 정합을 확인하세요.")

X <- X_raw[, features_common, drop = FALSE]
D <- D[features_common, features_common, drop = FALSE]
map <- map %>% filter(feature_id %in% features_common)

## ---------- 7) 유사도 S = 1 - D (PSD 보정) ----------
S0 <- 1 - D
S0 <- (S0 + t(S0)) / 2
diag(S0) <- 1
S0[S0 < 0] <- 0
S_psd <- as.matrix(nearPD(S0, corr = FALSE)$mat)

## ---------- 8) X 전처리 (CTRL 기반): log1p -> (CTRL 중앙값) 결측 대체 -> (CTRL 기준) z-score ----------
# 8-0) ctrl 인덱스 확인
ctrl_idx <- which(df$group == "ctrl")
if (length(ctrl_idx) < 2L) warning("CTRL 표본이 매우 적습니다. 표준화가 불안정할 수 있습니다.")

# 8-1) log1p
X_log <- log1p(X)  # n x p

# 8-2) 결측 대체값: "CTRL 샘플의 열 중앙값"
ctrl_col_median <- apply(X_log[ctrl_idx, , drop = FALSE], 2, function(v) median(v, na.rm = TRUE))

# 만약 어떤 피처가 CTRL에서 전부 NA라면(희귀 케이스), 전체 표본 중앙값으로 대체
overall_col_median <- apply(X_log, 2, function(v) median(v, na.rm = TRUE))
use_median <- ifelse(is.finite(ctrl_col_median), ctrl_col_median, overall_col_median)

# 모든 샘플에 동일한 대체값 적용(CTRL 기준 일관성)
X_imp <- X_log
for (j in seq_len(ncol(X_imp))) {
  nas <- is.na(X_imp[, j])
  if (any(nas)) X_imp[nas, j] <- use_median[j]
}

# 8-3) z-score 파라미터: "CTRL 샘플에서만" 평균/표준편차 계산
ctrl_means <- colMeans(X_imp[ctrl_idx, , drop = FALSE])
ctrl_sds   <- apply(X_imp[ctrl_idx, , drop = FALSE], 2, sd)

# 분산 0 보호(CTRL에서 변화가 없는 피처)
ctrl_sds[!is.finite(ctrl_sds) | ctrl_sds == 0] <- 1

# 8-4) 모든 샘플에 CTRL 파라미터로 표준화 적용
X_std <- sweep(X_imp, 2, ctrl_means, FUN = "-")
X_std <- sweep(X_std, 2, ctrl_sds,   FUN = "/")

# 형태/이름 복원
X_std <- as.matrix(X_std)
rownames(X_std) <- rownames(df)
colnames(X_std) <- colnames(X)

## ---------- 9) 모듈 리스트 (size >= 3) ----------
modules <- split(map$feature_id, map$module_id)
modules <- lapply(modules, function(v) intersect(v, colnames(X_std)))
min_size <- 3L
modules <- modules[sapply(modules, length) >= min_size]
if (length(modules) < 10)
  warning("유효 모듈이 적습니다. 모듈 크기 임계(min_size)를 낮추거나 클러스터 설정을 조정하세요.")

## ---------- 10) 모듈별 커널 K_g ----------
centerK <- function(K){ n <- nrow(K); H <- diag(n) - matrix(1/n, n, n); H %*% K %*% H }

Klist <- vector("list", length(modules))
names(Klist) <- paste0("M", names(modules))

for (id in names(modules)) {
  feats <- modules[[id]]
  Xg <- X_std[, feats, drop = FALSE]
  Sg <- S_psd[feats, feats, drop = FALSE]
  Kg <- Xg %*% Sg %*% t(Xg) / length(feats)
  Klist[[paste0("M", id)]] <- centerK(Kg)   # trace 정규화는 varComp에서 처리 예정
}

## ---------- 11) 저장 ----------
saveRDS(X_std,   file.path("results","preproc_X.rds"))
saveRDS(df,      file.path("results","preproc_y_df.rds"))    # df$group은 factor(ctrl,le,mp,px,sl)
saveRDS(S_psd,   file.path("results","preproc_S_psd.rds"))
saveRDS(modules, file.path("results","preproc_modules.rds"))
saveRDS(Klist,   file.path("results","preproc_Klist.rds"))

cat("Preprocessing done.\n",
    sprintf("n = %d samples, p = %d features (after intersect)\n", nrow(X_std), ncol(X_std)),
    sprintf("valid modules = %d (min size >= %d)\n", length(modules), min_size),
    "Group counts saved: results/preproc_group_counts.csv\n", sep = "")



#TESTING


# ---------- 2) 첫 모듈 선택 ----------
# Klist는 보통 이름이 "M{module_id}" 형태
if (length(Klist) < 1) stop("Klist 비어있음.")
first_name <- names(Klist)[1]
K_one <- list(Klist[[1]])         # length 1인 리스트로 감싸기
names(K_one) <- first_name

# 모듈 크기/ID 파악(보고용)
mod_id <- sub("^M","", first_name)
mod_size <- if (!is.null(modules[[mod_id]])) length(modules[[mod_id]]) else NA_integer_

# ---------- 3) A안: 단일 모듈 커널로 적합 후 단일 성분 검정 ----------
fit_one <- varComp(
  fixed = plant ~ group,
  data  = df,
  varcov = K_one,          # 오직 첫 모듈 커널만
  normalizeTrace = TRUE
)

# 요약 저장
sink(file.path("results","first_module_fit_summary.txt"))
cat("=== varComp single-kernel fit summary ===\n")
print(summary(fit_one))
cat("\n--- variance component (tau_hat) ---\n")
print(coef(fit_one, "varComp"))
cat("\n--- sigma^2 ---\n")
print(fit_one$sigma2)
sink()

