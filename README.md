# DIMMSC
An R package for clustering droplet-based single cell transcriptomic data

[Homepage @ Github](http://wt2015-github.github.io/DIMMSC/) | [Homepage @ Wei Chen's Lab](http://www.pitt.edu/~wec47/singlecell.html) | [Source Code](https://github.com/wt2015-github/DIMMSC)

## Introduction
**DIMMSC** is an R package for clustering droplet-based single cell transcriptomic data. It uses Dirichlet mixture prior to characterize variations across different clusters. An expectation-maximization algorithm is used for parameter inference. This package can provide clustering uncertainty.

**plot_tsne_clusters** performs log2 normalization, PCA and t-SNE on the data and plot clusters on the t-SNE projection.

## Installation
Install third-party R package [*cellrangerRkit*](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/rkit) from 10X Genomics.

Install DIMMSC with R command
```
install.packages(pkgs = "DIMMSC_*.*.*.tar.gz", repos = NULL, type = "source")

# Or install from GitHub
library(devtools)
install_github("wt2015-github/DIMMSC")
```
Or terminal command
```
R CMD INSTALL DIMMSC_*.*.*.tar.gz
```

## Usage
```
DIMMSC(data, K = 2, method_cluster_initial = "kmeans", method_alpha_initial = "Ronning", maxiter = 200, tol = 1e-04, lik.tol = 0.01, customized_initial = NULL)

plot_tsne_clusters(data, cluster)
```

## Arguments
* *data* : a G*C matrix with G genes and C cells
* *K* : number of desired clusters
* *method_cluster_initial* : method for initializing clusters, "kmeans", "random" or "customized"
* *method_alpha_initial* : method for initializing the alpha matrix for EM algorithm, "Ronning" (Ronning's method, 1989) or "Weir" (Weir and Hill's method 2002)
* *maxiter* : maximum number of iterations
* *tol* : a convergence tolerance for the difference of vector pie between iterations
* *lik.tol* : a convergence tolerance for the difference of log-likelihoods between iterations
* *customized_initial* : a vector of positive integers {1,2,...,K} for the customized initial clusters, the length should be equal to the number of cells
* *cluster* : a vector of clustering member ship, e.g. mem of DIMMSC output

## Values
DIMMSC returns a list object containing:
* *pie* : a vector of pie estimates
* *delta* : a C*K matrix with probability that each cell belongs to each cluster
* *alpha* : a K*G matrix of alpha estimates
* *mem* : a vector of clustering label
* *loglik* : the final log likelihood after iterations
* *AIC* : Akaike information criterion (AIC)
* *BIC* : Bayesian information criterion (BIC)
plot_tsne_clusters outputs a figure of t-SNE and clusters

## Example:
```
# Load the example data data_DIMMSC
data("data_DIMMSC")

# Run DIMMSC
result <- DIMMSC(data=data_DIMMSC, K=3, method_cluster_initial="kmeans", method_alpha_initial="Ronning", maxiter=200, tol=1e-4, lik.tol=1e-2)

# Plot t-SNE and clusters
plot_tsne_clusters(data=data_DIMMSC, cluster=result$mem)
```

## Publications
* Zhe Sun, **Ting Wang**, Ke Deng, Xiao-Feng Wang, Robert Lafyatis, Ying Ding, Ming Hu, Wei Chen. DIMM-SC: A Dirichlet mixture model for clustering droplet-based single cell transcriptomic data. *Bioinformatics* 2017.

## Contact
[Ting Wang](http://wt2015-github.github.io/) ([email](wang9ting@gmail.com)), [Wei Chen](http://www.pitt.edu/~wec47/index.html) ([email](wei.chen@chp.edu)).
