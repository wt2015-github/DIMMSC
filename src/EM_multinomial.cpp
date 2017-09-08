#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h>
#include <vector>
#include <string>
#include <random>
#include <algorithm>

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]


// Estimating clusters for droplet-based single cell data using EM algorithm
//
//  @param data a G*C matrix with G genes and C cells
//  @param K number of desired clusters
//  @param alpha initail alpha matrix coming from the EM_initial_alpha function
//  @param maxiter maximum number of iterations
//  @param tol a convergence tolerance for the difference of vector pie between iterations
//  @param lik.tol a convergence tolerance for the difference of log-likelihoods between iterations
//  @return a list object containing:
//  @return pie: a vector of pie estimates
//  @return delta: a C*K matrix with probability that each cell belongs to each cluster
//  @return alpha: a K*G matrix of alpha estimates
//  @return mem: a vector of clustering label
//  @return loglik: the final log likelihood after iterations
//  @return AIC: Akaike information criterion (AIC)
//  @return BIC: Bayesian information criterion (BIC)
//  @author Zhe Sun <zhs31@pitt.edu>, Ting Wang <tiw33@pitt.edu>, Ark Fang <zhf9@pitt.edu>, Wei Chen <wei.chen@chp.edu>.
//  @references Zhe Sun, Ting Wang, Ke Deng, Xiao-Feng Wang, Robert Lafyatis, Ying Ding, Ming Hu, Wei Chen. DIMM-SC: A Dirichlet mixture model for clustering droplet-based single cell transcriptomic data. Submitted.

using namespace std;

//So...alright no comment is allowed between export line and the function

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec lgammaVec(arma::vec& inVec, int G){
    arma::vec outVec(G);
    for(int g=0; g<G; g++){
        outVec(g) = lgamma(inVec(g));
        if(!isfinite(outVec(g))){ // if inVec(g) is closer than -1e-154 to 0, e.g. -1e-155, -inf will be returned.
            //So we would truncate it to a certain number: here it's lgamma(-1e-154)
            outVec(g) = 354.5891; // lgamma(-1e-154)
        }
    }
    return outVec;
}

arma::vec rdirichlet (int n){
    arma::vec v = arma::randg<arma::vec>(n);
    double sum = arma::sum(v);
    arma::vec out = v/sum;
    return(out);
}

int indexMax( arma::vec& inVec){
    int n = inVec.size();
    int idx = 0;
    double tmpValue = inVec(0);
    if(n>1){
        for(int i=1; i<n; i++){
            if(inVec(i) > tmpValue){
                idx = i;
                tmpValue = inVec(i);
            }
        }
    }
    return(idx);
}

// [[Rcpp::export]]
Rcpp::List EM_multinomial(const arma::mat& data, const int K, const arma::mat& alphaInput, const int maxiter, const double tol, const double likTol){   //data[G,J], alpha[K,G]
    arma::mat alpha=alphaInput;
    arma::vec repTmp(K);
    repTmp.fill(1);
    arma::vec pie(K);

    pie = rdirichlet(K);

    arma::vec pieNew(K);

    int J = data.n_cols;
    int G = data.n_rows;
    double differ = 1.0;
    int iter = 0;
    double logLik = 1.0;
    double dif = 100.0;
    arma::mat delta = arma::zeros<arma::mat>(J,K);

    while( (differ > tol || dif > likTol) && (iter < maxiter) ){
        // E-step compute omega
        arma::mat num1 = arma::zeros<arma::mat>(J,K);
        arma::mat num2 = arma::zeros<arma::mat>(J,K);

        arma::mat lgammaAlpha = arma::zeros<arma::mat>(G,K); // lgamma(alpha[,].t())
        arma::vec lgammaSumAlpha(K); //lgamma(sum(alpha(k,)))
        arma::vec dataAlphaTmp(G); //data[,j]+alpha[,k]
        arma::vec logpie(K); //log(pie(k))
        arma::mat deltaTmp = arma::zeros<arma::mat>(J,K);
        arma::mat alphaT = alpha.t();
        arma::vec alphaSum(K);
        arma::vec dataSum(J);
        delta.fill(0);

        for(int k=0; k<K; k++){
            // arma::vec tmpV = alpha.row(k).t();//always remember to transpose rowvec when you are surely using vec.
            arma::vec tmpV = alphaT.col(k);
            lgammaAlpha.col(k) = lgammaVec(tmpV,G);
            lgammaSumAlpha(k) = lgamma(arma::sum(alpha.row(k)));
            logpie(k) = log(pie(k));
            alphaSum(k) = arma::sum(alpha.row(k));
        }
        for(int j=0; j<J; j++){
            dataSum(j) = arma::sum(data.col(j));
            for(int k=0; k<K; k++){
                dataAlphaTmp = data.col(j) + alphaT.col(k);
                num1(j,k) = arma::sum( lgammaVec(dataAlphaTmp,G) - lgammaAlpha.col(k) );
                num2(j,k) = lgammaSumAlpha(k) - lgamma( alphaSum(k) + dataSum(j) );
                deltaTmp(j,k) = num1(j,k) + num2(j,k) + logpie(k);
            }
        }

            for(int j=0; j<J; j++){
                for(int k=0; k<K; k++){
                    double sumTmp=0.0;
                    for(int m=0; m<K; m++){
                        if(m!=k){
                            sumTmp += exp(arma::as_scalar(deltaTmp(j,m) - deltaTmp(j,k)));
                        }
                    }
                    delta(j,k) = 1.0/(1.0+sumTmp);
                }
            }

            //M-step: update pie and alpha

            for(int k=0; k<K; k++){
                pieNew(k) = arma::sum(delta.col(k))/double(J);
            }

            arma::mat alphaNew = arma::zeros<arma::mat>(K,G);

            for(int k=0; k<K; k++){
                double den = 0.0;
                for(int j=0; j<J; j++){;
                    den += arma::as_scalar(delta(j,k))*dataSum(j)/( dataSum(j)-1+alphaSum(k) );
                }
                for(int g=0; g<G; g++){
                    double tmpNum=0.0;
                    arma::vec dataG(J);
                    dataG = ((data.row(g)+0.000001) / (data.row(g)+0.000001-1+alpha(k,g))).t();
                    // tmpNum += arma::as_scalar(delta.col(k).t()*dataG); //after matrix operation it stays a matrix
                    // below is another option of above, fix, are the results same?
                    for(int j=0; j<J; j++){
                        tmpNum += delta(j,k)*dataG(j);
                    }
                    alphaNew(k,g) = alpha(k,g)*tmpNum/den;
                }
            }

            arma::uvec mem(J);
            for(int j=0; j<J; j++){
                arma::vec tmpVec = delta.row(j).t();
                mem(j) = indexMax(tmpVec); //.index_max(); //elements are unsigned int
            }
            // arma::uvec sort = arma::sort_index(pie);
            arma::vec sort(K); //actually is rank function in R
            arma::uvec sortidx = arma::sort_index(pie);
            int tmpk=0;
            while(tmpk<K){
                for(int k=0; k<K; k++){
                    if(sortidx(k)==tmpk){
                        sort(tmpk)=k;
                        tmpk++;
                        break;
                    }
                }
            }
            arma::uvec res = mem;
            for(int k=0; k<K; k++){
                for(int j=0; j<J; j++){
                    if(mem(j)==k){
                        res(j)=sort(k);
                    }
                }
            }
            mem=res;
            arma::mat num=num1+num2;
            arma::vec lik(J);
            for(int j=0; j<J; j++){
                lik(j) = num(j,mem(j));
            }
            double newLogLik = arma::sum(lik);

            // calculate diff to check convergence
            dif = abs((newLogLik - logLik)/logLik*100.0);

            double sumd=0.0;
            for(int k=0; k<K; k++){
                sumd += pow((arma::as_scalar(pieNew(k) - pie(k))), 2.0);
            }

            differ = sqrt(abs(sumd));
            pie=pieNew;
            alpha=alphaNew;
            logLik = newLogLik;

            for(int k=0; k<K; k++){
                for(int g=0; g<G; g++){
                    if(alpha(k,g)==0){
                        alpha(k,g) = 0.000001;
                    }
                }
            }

            iter++;
    }

        arma::uvec mem(J);
        for(int j=0; j<J; j++){
            for(int k=0; k<K; k++){
                arma::vec tmpVec = delta.row(j).t();
                mem(j) = indexMax(tmpVec); //.index_max() will use the last index when max is tie.
            }
        }
        double AIC, BIC;
        AIC = (-2)*logLik + 2*((K*G)+(K*J));
        BIC = (-2)*logLik + log(G*J)*((K*G)+(K*J));

        Rcpp::List result;

        result["pie"] = pie;
        result["delta"] = delta;
        result["alpha"] = alpha;
        result["mem"] = mem + 1;
        result["loglik"] = logLik;
        result["AIC"] = AIC;
        result["BIC"] = BIC;

        return result;
}
