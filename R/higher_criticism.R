#' This function transfers prior information into the weights
#'
#' The weight is calculated with the equation: w=1/(a x prior_info+b x mean(prior_info)) and then scaled into mean(w)=1. Here we take a=0.95 and b=0.05.
#'
#' @param w a numeric vector which is the original statistic for weight,need to be non-negative
#' @return a numeric vector of the transfered weight with around mean~1, min~0.05.
#'  Missing values in input are transfered into 1
#' @export
#' @examples
#' a=c(NA,exp(rnorm(7,0,1)))
#' trans_w(a)
trans_w=function(w){
  #w transfer to 1/w, with mean(1/0.95*w+0.05*(mean(w)))=1, to avoid zeros, we set lower boundary
  if(max(w,na.rm=TRUE)==min(w,na.rm=TRUE)){
    return(w/w)
  }
  else{
    w0=1/(0.95*w+0.05*(mean(w,na.rm=TRUE)))
    #w0=length(w)/(w*(length(w)-1)+mean(w))
    w0=w0/mean(w0,na.rm=TRUE)
    w0[is.na(w0)]=1
    return(w0)
  }
}

#' This function calculates the weighted-pvalues.
#'
#' @description It calculates the CDF(Cumulative Distribution Function) of weighted pvalues from the original pvalues and their weights and assign new pvalues based on this CDF.
#'
#' @param pm is a statistics matrix of P-values from nrow=n genes(independent tests),ncol=d. Pm are not encouraged to have only 1 rows, if that happend, warning massage will produced and this function returns sqrt(1/pm-1), rather than the double side statistics from linear regression.Thus, we have S(t)~Binomial(n,p).
#' @param w is the weight, if not specify, w=1, if specify w must have the same length as nrow Pm
#'
#' @return numeric vector with each elements is a Higher criticism values calculated from each colum of the Pm
#' @import Rcpp
#' @export
#'
#' @examples
#' pval=matrix(runif(200,0,1),ncol=4,nrow=50)
#' w0=seq(0.5,1.5,length=50)
#' pwval=cal_cdf(pval,w=w0)
#' @references Genovese, C. R., Roeder, K., & Wasserman, L. (2006). False Discovery Control with p-Value Weighting. Biometrika, 93(3), 509–524.
cal_cdf=function(pm,w=1){
  if(w==1 || max(w)==min(w)){
    return(pm)
  }
  pm=as.matrix(pm)
  l=nrow(pm)
  wm=as.matrix(w)%*%as.matrix(t(1/w))
  # Fx=as.matrix(apply(pm,2,function(x){return(rowSums(pmin(x*wm,1))/l)}))
  Fx=colsum_pmin1(pm, wm) #import c code for speeding up.
  return(Fx)
}

#' This function calculates the higher criticism
#  ,the integrated higher criticism and the topsum statistics as well as their variations.
#'
#' For ordered p-values p(1)<p(2)<...<p(n).Define Sn(t)=sum_{i=1}^n 1_{p_{(i)}<=t}.
#' Then define the statistic T as T=(Sn(p(i))-n x p(i))/sqrt(n x p(i) x (1-p(i))).
#' The higher criticism is calculated with HCmax=max T_(i). where 0<i<=(t0ratio x n)
#' @param pm is a statistics matrix of P-values or weighted pvalues, each row represents a gene (independent tests) and each column represents a dataset (e.g. a permutation or an observation). Pm are not encouraged to have only 1 rows, if that happend, warning massage will produced.
#' @param t0ratio is the ratio for the region c(0,t0ratio) of pvalues for statistic calculation.
#' @param filter is the threshold to exclude extremely small pvalues to avoid them driving all signals.default 0.
#'
#' @return a numeric vector with each elements is a Higher criticism values calculated from each colum of the Pm
#' @export
#' @useDynLib wHC
#'
#' @examples
#' pval=matrix(runif(200,0,1),ncol=4,nrow=50)
#' w0=seq(0.5,1.5,length=50)
#' pwval=cal_cdf(pval,w=w0)
#' hc_cal(pwval,t0ratio=0.4)
#' @references Donoho, D., & Jin, J. (2004). Higher Criticism for Detecting Sparse Heterogeneous Mixtures. The Annals of Statistics, 32(3), 962–994.
hc_cal=function(pm,t0ratio=1,filter=0){
  #step1: arrange the matrix
  pm=as.matrix(pm)
  d=ncol(pm)
  n=nrow(pm)
  t0n=max(floor(n*t0ratio),1)
  if(n==1){return(pm)}
  else{
    pm=as.matrix(apply(pm,2,function(x){return(sort(x,na.last = TRUE))}))
    S=apply(pm,2,function(x){return(rank(x,na.last = "keep", ties.method = "max"))})
    preHC=(S-pm*n)/sqrt(pm*(n-pm*n))
    preHCm=preHC
    preHC[pm>=t0ratio]=NA
    preHC[pm<=filter]=NA
    HC=suppressWarnings(apply(preHC,2,function(x){return(max(x,na.rm=TRUE))}))
  }
  return(HC)
}
