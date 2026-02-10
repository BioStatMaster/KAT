#########################
#######  Packages #######
#########################

need <- c("glassoFast", "foreach", "doParallel", "stringr", "dplyr")
for (pkg in need) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

#########################
####### Functions #######
#########################

write_log <- function(message) {
  cat(paste0("[", Sys.time(), "] ", message, "\n"), 
      file = "log.txt", append = TRUE)
}

load_and_process <- function(filename) {
  if(!file.exists(filename)) {
    write_log(paste("Error: File not found ->", filename))
    return(NULL)
  }
  
  rawdata <- read.csv(filename, header = TRUE, na.strings = "")
  rawdata <- t(rawdata)
  rawdata <- rawdata[-c(1:56), ]
  colnames(rawdata) <- rawdata[NROW(rawdata), ]
  rawdata <- rawdata[-NROW(rawdata), ]
  rawdata <- as.data.frame(rawdata)
  rawdata[] <- lapply(rawdata, function(x) {
    if (is.character(x) | is.factor(x)) suppressWarnings(as.numeric(as.character(x))) else x
  })
  return(rawdata)
}

# Drton (2014) EBIC 
run_fast_ebic_logged <- function(data, lambda_seq, name_tag, gamma = 0.5) {
  
  X <- scale(as.matrix(data))
  S <- cov(X)
  n <- nrow(X)
  p <- ncol(X) 
  
  write_log(sprintf("[%s] Start Fitting (n=%d, p=%d)", name_tag, n, p))
  
  best_ebic <- Inf
  best_Theta <- NULL
  best_lambda <- NULL
  
  total_steps <- length(lambda_seq)
  
  for (i in 1:total_steps) {
    lam <- lambda_seq[i]
    
    fit <- glassoFast(S, rho = lam)
    Theta <- fit$wi
    
    log_det <- tryCatch(determinant(Theta, logarithm = TRUE)$modulus[1], error = function(e) -Inf)
    tr_S_Theta <- sum(diag(S %*% Theta))
    
    neg_2_log_lik <- -n * (log_det - tr_S_Theta)
    
    upper_tri <- Theta[upper.tri(Theta)]
    n_edges <- sum(abs(upper_tri) > 1e-5)
    
    ebic <- neg_2_log_lik + n_edges * log(n) + 4 * n_edges * gamma * log(p)
    
    
    if (ebic < best_ebic) {
      best_ebic <- ebic
      best_Theta <- Theta
      best_lambda <- lam
      
      write_log(sprintf(" -> [%s] New Best! Lambda=%.4f, EBIC=%.2f, Edges=%d", name_tag, lam, ebic, n_edges))
    }
  }
  
  write_log(sprintf("[%s] Finished. Best Lambda: %.4f, EBIC: %.2f", name_tag, best_lambda, best_ebic))
  
  colnames(best_Theta) <- rownames(best_Theta) <- colnames(data)
  return(best_Theta)
}

#########################
####### Load Data #######
#########################

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

cat("Analysis Log\n", file = "log.txt", append = FALSE)
print("Logs will be saved to 'log.txt' in your working directory.")

files <- list(
  bak15 = "251209_bak15_all_features.csv",
  myc   = "251209_myc234_all_features.csv",
  npr   = "251209_npr1_all_features.csv",
  sobir = "251209_sobir1_all_features.csv",
  wt_all= "251209_wt_all_features.csv"
)

data_list <- list()

for (nm in names(files)) {
  data_list[[nm]] <- load_and_process(files[[nm]])
}


if(file.exists("metadata_hormone_herbivore.csv") && !is.null(data_list[["wt_all"]])) {
  raw6 <- read.csv("metadata_hormone_herbivore.csv", header=T)
  raw6$sample <- str_sub(raw6$sample, 1, -6)
  rownames(raw6) <- raw6[,1]
  raw6 <- raw6[,-c(1:7)]
  raw6[] <- lapply(raw6, function(x) if(is.character(x)) as.numeric(as.character(x)) else x)
  
  rawdata5 <- data_list[["wt_all"]]
  merged <- merge(raw6, rawdata5, by="row.names")
  rownames(merged) <- merged[,1]
  merged <- merged[,-1]
  
  data_list[["wt_all_with_hormone"]] <- as.data.frame(merged)
  write_log("Merged dataset 'wt_all_with_hormone' created.")
}

#############################
###### Parallel Fitting #####
#############################

n_cores <- parallel::detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)

write_log(paste("Starting parallel processing on", n_cores, "cores."))

# Lambda Sequence
lambda_seq <- seq(2, 0.05, length.out = 20)

results <- foreach(nm = names(data_list), .packages = c("glassoFast")) %dopar% {
  
  dat <- data_list[[nm]]
  if(is.null(dat)) return(NULL)
  Theta <- run_fast_ebic_logged(dat, lambda_seq, name_tag = nm)
  return(Theta)
}

stopCluster(cl)
write_log("All parallel jobs finished.")

names(results) <- names(data_list)
Theta1 <- results[["bak15"]]
Theta2 <- results[["myc"]]
Theta3 <- results[["npr"]]
Theta4 <- results[["sobir"]]
Theta5 <- results[["wt_all"]]
Theta6 <- results[["wt_all_with_hormone"]]

print("Done. Check 'log.txt' for details.")