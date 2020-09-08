computeClusterInertia <- function(labels, asset_returns) {
  inertia <- 0
  n <- nrow(asset_returns)
  for(label in unique(labels)){
    # removing `as.matrix` increases the speed significantly
    # inertia <- inertia + mean(as.matrix(dist(asset_returns[,labels==label])))
    inertia <- inertia + sum(dist(asset_returns[,labels==label])) * 2 / n^2
  }
  return(log(inertia))
}

computeExpectedInertia <- function(N, D, num_clusters, linkage, num_reference_datasets=5) {
  reference_inertia <- 0
  for(i in 1:num_reference_datasets){
    reference_asset_returns = replicate(N, runif(D))
    reference_corr <- cor(reference_asset_returns)
    reference_dist <- as.dist(sqrt(2 * (1 - reference_corr)))

    reference_cluster_assignments <- dendextend::cutree(hclust(reference_dist, method=linkage), k=num_clusters, order_clusters_as_data=F)

    inertia <- computeClusterInertia(reference_cluster_assignments, reference_asset_returns)
    reference_inertia <- reference_inertia + inertia
  }
  return(reference_inertia/num_reference_datasets)
}

computeNumClusters <- function(rho, linkage, asset_returns, distance) {
  gap_values <- list()
  num_clusters <- 1
  max_clusters <- -Inf
  N <- ncol(asset_returns)
  D <- nrow(asset_returns)
  original_clusters <- hclust(distance, method=linkage)

  for(i in 1:nrow(rho)){
    original_cluster_assignments <- dendextend::cutree(original_clusters, k=num_clusters, order_clusters_as_data=F)
    if(max(original_cluster_assignments)==max_clusters |
       max(original_cluster_assignments)>10){
      break
    }
    max_clusters = max(original_cluster_assignments)
    inertia <- computeClusterInertia(original_cluster_assignments, asset_returns)
    expected_inertia <- computeExpectedInertia(N, D, num_clusters, linkage)

    gap_values <- c(gap_values, expected_inertia - inertia)
    num_clusters <- num_clusters + 1
  }
  return(which.max(gap_values))
}

getIVP <- function(Sigma) {
  if(is.null(nrow(Sigma))) return(1) 
  ivp <- 1./diag(Sigma)
  ivp <- ivp/sum(ivp)
  return(ivp)
}

getClusterVariance <- function(Sigma, cluster_idx) {
  Sigma_ <- Sigma[cluster_idx, cluster_idx]
  w_ <- getIVP(Sigma_)
  cluster_var <- t(w_) %*% Sigma_ %*% w_
  return(cluster_var[1])
}


calculateClusterRiskContribution <- function(risk_measure, Sigma, asset_returns, 
                                             num_clusters, cluster_assignments){
  cluster_contributions <- rep(1, num_clusters)

  for(i in 1:num_clusters){
    cluster_asset_ind <- which(cluster_assignments==i)
    if(risk_measure=='variance'){
      cluster_contributions[i] <- getClusterVariance(Sigma, cluster_asset_ind)
    }
  }
  return(cluster_contributions)
}

getRecursiveBisection <- function(dend, cluster_contributions){
  num_clusters <- nobs(dend)
  weights <- rep(1,num_clusters)

  computeChildClusterWeights <- function(dend, cluster_contributions) {
    left <- dend[[1]]
    right <- dend[[2]]

    left_order <- order.dendrogram(left)
    right_order <- order.dendrogram(right)
    left_contrib <- sum(cluster_contributions[left_order])
    right_contrib <- sum(cluster_contributions[right_order])

    alpha <- right_contrib /(right_contrib+left_contrib)
    weights[left_order] <<- weights[left_order] * alpha
    weights[right_order] <<- weights[right_order] * (1-alpha)

    if(nobs(left)>1)
      computeChildClusterWeights(left, cluster_contributions)
    if(nobs(right)>1)
      computeChildClusterWeights(right, cluster_contributions)
  }
  computeChildClusterWeights(dend, cluster_contributions)
  return(weights)
}


calculateFinalPortfolioWeights <- function(risk_measure, clusters_weights, Sigma,
                                           num_clusters, cluster_assignments) {
  w <- rep(1, nrow(Sigma))
  for(i in 1:num_clusters){
    cluster_inds <- cluster_assignments==i
    Sigma_ <- Sigma[cluster_inds, cluster_inds]
    #cluster_asset_returns <- asset_returns[,cluster_inds]

    parity_weights <- getIVP(Sigma_)

    w[cluster_inds] <- parity_weights * clusters_weights[i]
  }
  return(w)
}

getPortfolioWeights <- function(hcluster, asset_returns, Sigma, num_clusters,
                                risk_measure, cluster_assignments) {
  N <- ncol(Sigma)

  clusters_weights <- rep(1, num_clusters)
  cluster_contributions <- calculateClusterRiskContribution(risk_measure,Sigma,
                                                            asset_returns,num_clusters,
                                                            cluster_assignments)

  dend_depth <- hcluster$height[N - num_clusters]
  clusters_dend <- cut(as.dendrogram(hcluster), h=dend_depth)$upper

  clusters_weights <- getRecursiveBisection(clusters_dend, cluster_contributions)

  w <- calculateFinalPortfolioWeights(risk_measure, clusters_weights, Sigma, 
                                      num_clusters, cluster_assignments)

  names(w) <- names(asset_returns)
  return(w)
}

#' @export
hierarchicalEqualRiskContribution <- function(asset_prices, 
                                              risk_measure=c('variance', 
                                                             'standard-deviation', 
                                                             'equal-weighting',
                                                             'CVaR', 'CDaR'), 
                                              linkage='ward.D2', num_clusters=NULL) {

  risk_measure <- match.arg(risk_measure)
  asset_returns <- diff(log(asset_prices))[-1]
  Sigma <- cov(asset_returns)
  rho <- cor(asset_returns)
  distance <- as.dist(sqrt(2 * (1-rho)))

  hcluster <- hclust(distance, method=linkage)
  if(is.null(num_clusters))
    num_clusters <- computeNumClusters(rho, linkage, asset_returns, distance)

  cluster_assignments <- dendextend::cutree(hcluster, k=num_clusters, order_clusters_as_data=F)
  asset_names <- names(asset_prices)
  cluster_assignments <- cluster_assignments[asset_names]

  w <- getPortfolioWeights(hcluster, asset_returns, Sigma, num_clusters,
                           risk_measure, cluster_assignments)

  return(list('w'=w))
}
