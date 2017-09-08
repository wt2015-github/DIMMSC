#!/usr/bin/env Rscript
# 07/17/2017

# initializing alpha matrix for EM algorithm
#
# @param data a G*C matrix with G genes and C cells
# @param clusters_initial a length C vector of initial clustering member ship
# @param method_alpha_initial method for initializing alpha, "Ronning" (Ronning's method, 1989) or "Weir" (Weir and Hill's method 2002)
# @return a K*G alpha matrix
# @author Zhe Sun <zhs31@pitt.edu>, Ting Wang <tiw33@pitt.edu>, Ark Fang <zhf9@pitt.edu>, Wei Chen <wei.chen@chp.edu>.
# @references Zhe Sun, Ting Wang, Ke Deng, Xiao-Feng Wang, Robert Lafyatis, Ying Ding, Ming Hu, Wei Chen. DIMM-SC: A Dirichlet mixture model for clustering droplet-based single cell transcriptomic data. Bioinformatics 2017.
#' @import dirmult
#' @import stats
EM_initial_alpha <- function(data, clusters_initial, method_alpha_initial)
{
  num_cluster <- length(unique(clusters_initial))
  cluster <- matrix(,num_cluster,1)
  for (i in 1:num_cluster){
    cluster[i,] <-  length(which(clusters_initial==i))/ncol(data)
  }
  cluster <- as.numeric(cluster)
  sort <- rank(cluster,ties.method="random")
  label.new <- clusters_initial
  for(j in 1:num_cluster){
    label.new[which(clusters_initial==j)] <- sort[j]
  }

  location <- list()
  for (m in 1:num_cluster){
    location[[m]] <- which(label.new==m)
  }
  newdata <- list()
  for (m in 1:num_cluster){
    newdata[[m]] <- data[,location[[m]]]
  }

  p <- matrix(,nrow(data),num_cluster)
  for (m in 1:num_cluster){
    sum <- sum(as.vector(newdata[[m]]))
    for ( i in 1:nrow(data)){
      p[i,m] <- sum(newdata[[m]][i,])/sum
    }
  }

  new.data <- list()
  for (m in 1:num_cluster){
    new.data[[m]] <- newdata[[m]][rowSums(newdata[[m]]) != 0, colSums(newdata[[m]]) != 0]
  }
  C <- matrix(,num_cluster,1)
  for (m in 1:num_cluster){
    C[m,1] <- ncol(new.data[[m]])
  }
  G <- matrix(,num_cluster,1)
  for (m in 1:num_cluster){
    G[m,1] <- nrow(new.data[[m]])
  }

  ppp <- list()
  for (m in 1:num_cluster){
    ppp[[m]] <- matrix(,G[m,1],C[m,1])
    for (j in 1:C[m,1]){
      tmp = sum(new.data[[m]][,j])
      ppp[[m]][,j]=new.data[[m]][,j]/tmp
    }
  }

  pp <- list()
  for (m in 1:num_cluster){
    pp[[m]] <- matrix(,G[m,1],1)
    for ( i in 1:G[m,1]){
      pp[[m]][i,1] <- mean(ppp[[m]][i,])
    }
  }

  v <- list()
  for (m in 1:num_cluster){
    v[[m]] <- matrix(,G[m,1],1)
    for ( i in 1:G[m,1]){
      v[[m]][i,1] <- var(ppp[[m]][i,])
    }
    v[[m]][which(v[[m]]==0),1] <- mean(v[[m]],na.rm=T)
  }

  s <- matrix(,num_cluster,1)
  for(m in 1:num_cluster){
    sum <- 0
    for (i in 1:(G[m,1]-1)){
      tmp <- log( ( pp[[m]][i,1]*(1-pp[[m]][i,1])/ v[[m]][i,1] ) -1 )
      sum <- sum+tmp
    }
    s[m,1] <- exp( (1/(G[m,1]-1))*sum )
  }

  new.alpha <- list()
  for(m in 1:num_cluster){
    new.alpha[[m]] <- s[m,1]*p[,m]
    new.alpha[[m]][which(new.alpha[[m]]==0)] <- 0.000001
  }

  new_alpha_R <- matrix(,num_cluster,nrow(data))
  for (m in 1:num_cluster){
    new_alpha_R[m,] <- new.alpha[[m]]
  }

  new.alpha <- list()
  for (m in 1:num_cluster){
    mom <- weirMoM(t(new.data[[m]]), se=FALSE)
    if (mom <= 0) {mom <- 0.005}
    initscalar <- (1 - mom)/mom
    new.alpha[[m]] <- p[,m]*initscalar
    new.alpha[[m]][which(new.alpha[[m]]==0)] <- 0.000001
  }

  new_alpha_W <- matrix(,num_cluster,nrow(data))
  for (m in 1:num_cluster){
    new_alpha_W[m,] <- new.alpha[[m]]
  }

  if(method_alpha_initial == "Ronning"){
    return(new_alpha_R)
  } else{
    return(new_alpha_W)
  }
}

#' Clustering droplet-based single cell transcriptomic data using Dirichlet mixture prior and EM algorithm
#'
#' DIMMSC is for clustering droplet-based single cell transcriptomic data. It uses Dirichlet mixture prior to characterize variations across different clusters. An EM algorithm is used for parameter inference. This package can provide clustering uncertainty.
#' @name DIMMSC
#' @aliases DIMMSC
#' @param data a G*C matrix with G genes and C cells
#' @param K number of desired clusters
#' @param method_cluster_initial method for initializing clusters, "kmeans", "random" or "customized"
#' @param method_alpha_initial method for initializing the alpha matrix for EM algorithm, "Ronning" (Ronning's method, 1989) or "Weir" (Weir and Hill's method 2002)
#' @param maxiter maximum number of iterations
#' @param tol a convergence tolerance for the difference of vector pie between iterations
#' @param lik.tol a convergence tolerance for the difference of log-likelihoods between iterations
#' @param customized_initial a vector of positive integers \{1,2,...,K\} for the customized initial clusters, the length should be equal to the number of cells
#' @return DIMMSC returns a list object containing:
#' @return \itemize{
#'   \item pie: a vector of pie estimates
#'   \item delta: a C*K matrix with probability that each cell belongs to each cluster
#'   \item alpha: a K*G matrix of alpha estimates
#'   \item mem: a vector of clustering label
#'   \item loglik: the final log likelihood after iterations
#'   \item AIC: Akaike information criterion (AIC)
#'   \item BIC: Bayesian information criterion (BIC)
#' }
#' @author Zhe Sun <zhs31@pitt.edu>, Ting Wang <tiw33@pitt.edu>, Ark Fang <zhf9@pitt.edu>, Wei Chen <wei.chen@chp.edu>.
#' @references Zhe Sun, Ting Wang, Ke Deng, Xiao-Feng Wang, Robert Lafyatis, Ying Ding, Ming Hu, Wei Chen. DIMM-SC: A Dirichlet mixture model for clustering droplet-based single cell transcriptomic data. Bioinformatics 2017.
#' @examples
#' # Load the example data data_DIMMSC
#' data("data_DIMMSC")
#'
#' # Run DIMMSC
#' result <- DIMMSC(data=data_DIMMSC, K=3, method_cluster_initial="kmeans",
#'                  method_alpha_initial="Ronning", maxiter=200, tol=1e-4, lik.tol=1e-2)
#'
#' # Plot t-SNE and clusters
#' plot_tsne_clusters(data=data_DIMMSC, cluster=result$mem)
#' @import stats
#' @useDynLib EM_multinomial
#' @export
DIMMSC <- function(data, K=2, method_cluster_initial="kmeans", method_alpha_initial="Ronning", maxiter=200, tol=1e-4, lik.tol=1e-2, customized_initial=NULL)
{
  # Initialize clusters
  cat('Initializing clusters ')
  if(method_cluster_initial == "kmeans"){
    cat('with kmeans clusters.\n')
    kdata <- as.matrix(log2(data+1)) # Normalize data for k-means
    k <- kmeans(t(kdata),K)
    clusters_initial <- k$cluster
  } else if(method_cluster_initial == "random"){
    cat('with random clusters.\n')
    clusters_initial <- sample(seq(K), ncol(data), replace = T)
  } else if(method_cluster_initial == "customized"){
    cat('with customized clusters.\n')
    if(!(length(customized_initial) == ncol(data) & max(as.numeric(names(table(customized_initial)))) == K & min(as.numeric(names(table(customized_initial)))) == 1)){
      stop('Customized initial clusters should be a vector of positive integers {1,2,...,K}, the length should be equal to the number of cells.')
    } else if(!(sum(as.numeric(names(table(customized_initial))) == seq(K)) == K)){
      stop('Customized initial clusters should be a vector of positive integers {1,2,...,K}, the length should be equal to the number of cells.')
    } else{
      clusters_initial <- customized_initial
    }
  } else{
    stop('Method for initializing clusters should be "kmeans", "random" or "customized".')
  }

  # Initialize alpha matrix
  if(!(method_alpha_initial == "Ronning" | method_alpha_initial == "Weir")){
    stop('Method for initializing alpha matrix should be "Ronning" or "Weir".')
  }
  cat('Initializing alpha matrix with',method_alpha_initial,'method.\n')
  alpha <- EM_initial_alpha(data=data, clusters_initial=clusters_initial, method_alpha_initial=method_alpha_initial)

  # DIMMSC clustering
  cat('Performing DIMMSC clustering...\n')
  result <- EM_multinomial(data=data, K=K, alphaInput=alpha, maxiter=maxiter, tol=tol, likTol=lik.tol)

  cat('Analysis is finished.\n')
  return(result)
}

#' Plot t-SNE and clusters
#'
#' plot_tsne_clusters performs log2 normalization, PCA and t-SNE on the data and plot clusters on the t-SNE projection
#' @rdname DIMMSC
#' @param cluster a vector of clustering member ship, e.g. mem of DIMMSC output
#' @return plot_tsne_clusters outputs a figure of t-SNE and clusters
#' @import cellrangerRkit
#' @import stats
#' @export
plot_tsne_clusters <- function(data, cluster){

  if(length(cluster) != ncol(data)){
    stop('the length of clustering member ship is different from the number of cells.\n')
  }

  # log2 normalization, PCA and t-SNE
  cat('log2 normalization, PCA and t-SNE.\n')
  kdata <- as.matrix(log2(data + 1))
  res.pca <- prcomp(t(kdata), center = TRUE, scale. = TRUE)
  res.pca$x <- res.pca$x[,1:10]
  res.tsne <- run_tsne(res.pca)

  # plot t-SNE and clusters
  cat('plot t-SNE projection and clusters.\n')
  tsne_clust <- data.frame(Barcode = as.factor(1:nrow(res.tsne$Y)),TSNE.1 = res.tsne$Y[,1], TSNE.2 = res.tsne$Y[,2], Cluster = cluster)
  visualize_clusters(tsne_clust$Cluster, tsne_clust[c("TSNE.1","TSNE.2")], title="t-SNE and cell clustering", marker_size = 1)
}

