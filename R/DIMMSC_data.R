#!/usr/bin/env Rscript

#' A simulated droplet-based single cell data
#'
#' @docType data
#' @usage data("data_DIMMSC")
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
"data_DIMMSC"
