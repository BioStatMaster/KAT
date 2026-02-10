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
load("Theta11.Rdata")
load("Theta12.Rdata")
load("Theta13.Rdata")
load("Theta21.Rdata")
load("Theta22.Rdata")
load("Theta23.Rdata")
load("Theta31.Rdata")
load("Theta32.Rdata")
load("Theta33.Rdata")
load("Theta41.Rdata")
load("Theta42.Rdata")
load("Theta43.Rdata")
load("Theta51.Rdata")
load("Theta52.Rdata")
load("Theta53.Rdata")
load("Theta61.Rdata")
load("Theta62.Rdata")
load("Theta63.Rdata")


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

suffixes <- c("11", "12", "13", 
              "21", "22", "23", 
              "31", "32", "33", 
              "41", "42", "43", 
              "51", "52", "53", 
              "61", "62", "63")

for (suffix in suffixes) {
  
  theta_name <- paste0("Theta", suffix)
  if (exists(theta_name)) {
    
    Theta_current <- get(theta_name)
    Pcorr_current <- partial_corr_from_precision(Theta_current)
    W_current <- abs(Pcorr_current)
    diag(W_current) <- 0
    varname <- colnames(Theta_current)
    Scores_current <- compute_hubness_core(W_current, varname)
    w_name <- paste0("W", suffix)
    scores_name <- paste0("Scores", suffix)
    
    assign(w_name, W_current)
    assign(scores_name, Scores_current)
    
    write.csv(W_current, file = paste0(w_name, ".csv"), row.names = TRUE)
    write.csv(Scores_current, file = paste0(scores_name, ".csv"), row.names = FALSE)
    cat(paste(theta_name, "Done. -> CSV saved\n"))
    
  } else {
    warning(paste(theta_name, "Error. -> Skip"))
  }
}


