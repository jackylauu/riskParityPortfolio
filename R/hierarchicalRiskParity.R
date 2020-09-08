getIVP <- function(Sigma) {
  if(is.null(nrow(Sigma))) return(1) 
  ivp <- 1./diag(Sigma)
  ivp <- ivp/sum(ivp)
  return(ivp)
}

getClusterVar <- function(Sigma, cluster_idx) {
  Sigma_ <- Sigma[cluster_idx, cluster_idx]
  w_ <- getIVP(Sigma_)
  cluster_var <- t(w_) %*% Sigma_ %*% w_
  return(cluster_var[1])
}

bisectClusters <- function(clusters){
  clusters_new <- list()
  for(cl in clusters){
    if(length(cl)<2) next
    n <- length(cl) %/% 2
    len <- length(cl)
    clusters_new <- c(clusters_new, list(cl[1:n]), list(cl[(n+1):len]))
  }
  return(clusters_new)
}

getRecBipart <- function(Sigma, sorted_idx) {
  N <- length(sorted_idx)
  w <- rep(1., N)
  clusters <- bisectClusters(list(sorted_idx))

  while(length(clusters) > 0) {
    for(i in 1:length(clusters)){
      if(i%%2==0 || length(clusters)==1) next

      cl0 <- unlist(clusters[i])
      cl1 <- unlist(clusters[i+1])

      cl_var0 <- getClusterVar(Sigma, cl0)
      cl_var1 <- getClusterVar(Sigma, cl1)
      alpha <- 1 - cl_var0 / (cl_var0 + cl_var1)
      w[cl0] <- w[cl0] * alpha
      w[cl1] <- w[cl1] * (1-alpha)
    }
    clusters <- bisectClusters(clusters)

  }
  names(w) <- rownames(Sigma)
  return(w)
}

#' @export
hierarchicalRiskParity <- function(asset_prices, linkage='single') {
  asset_returns <- diff(log(asset_prices))[-1]
  Sigma <- cov(asset_returns)
  rho <- cor(asset_returns)
  distance <- as.dist(sqrt((1-rho)/2))

  hcluster <- hclust(distance, method=linkage)
  sorted_idx <- hcluster$order

  w <- getRecBipart(Sigma, sorted_idx)

  return(list('w'=w))
}
