#' This function calculates the goodness-of-fit tests
#  ,the integrated the statistics as well as their variations.
#'
#' For ordered p-values p(1)<p(2)<...<p(n).Define \eqn{ Sn(t)=sum_{i=1}^n 1_{p_{(i)}<=t} } .
#' Then define the statistic T as
#' T=(Sn(p(i))-n x p(i))/sqrt(n x p(i) x (1-p(i))) for score statistic. or
#' T=(Sn(p(i))-n x p(i))/sqrt(Sn(p(i))*(n-Sn(p(i)))/n) for wald statistic. or
#' T=nx(p(i)xlog(p(i)xn/Sn(p(i)))+(1-p(i))xlog((1-p(i))xn/(n-Sn(p(i)))))   for log likelihood ratio
#'
#' The max method calculates the Tmax=max T_(i) where 0 < i < t0ration x n
#' The sum method calculates the Tsum=sum T_(i) where 0 < i < t0ration x n
#' The rsum is calculated with Trsum=1/n x sum T_(i). where 0<i<=(i_max),
#' where i_max is the index of the maxium of Ts.
#' The topsum is calculated with Ttopsum=1/n max T_(r).
#' where r belongs to the subset that Ts in the subsets are the top t0ratio proportion among all Ts.
#; The max2, sum2m,rsum2 and topsum2 represents the situations that replacing the statistic T with T^2.They only works on score statistic or wald statistic
#' @param pm is a statistics matrix of P-values or weighted pvalues, each row represents a gene (independent tests) and each column represents a dataset (e.g. a permutation or an observation). Pm are not encouraged to have only 1 rows, if that happend, warning massage will produced.
#' @param t0ratio is the ratio for the region c(0,t0ratio) of pvalues for statistic calculation.
#' @param gof_method is the option for the set based analysis methods, available option includes: "higher_criticism","berk_jones","skat","einmahl_mckeague", see our publication for more details.
#' @param single_statistic is the option for single statistic, it can be c("score","likelihood_ratio","wald"),it would be ignored is gof_method is assigned
#' @param accumulate_option is the option for combining method, it can be c("max","sum","rsum","topsum","max2","sum2","rsum2","topsum2"), it would be ignored is gof_method is assigned, the "max2","sum2","rsum2","topsum2" will be replaced into "max","sum","rsum","topsum" if gof method is assigned.
#' @param filter is the threshold to exclude extremely small pvalues to avoid them driving all signals.default 0.
#' @param weight_option defines the external prior information as weight. It can be "none", "in" and "out". "none" assigns no weight, "in" assigns weight towards each single pvalues and "out" assigns weight towards each statistic T.
#' @param weight is a vector which provides weight towards each genes.
#'
#' @return a numeric vector with each elements is a Higher criticism values calculated from each colum of the Pm
#' @export
#'
#' @examples
#' a=matrix(runif(200,0,1),ncol=4,nrow=50)
#' gof_cal(a)
#' @references Donoho, D., & Jin, J. (2004). Higher Criticism for Detecting Sparse Heterogeneous Mixtures. The Annals of Statistics, 32(3), 962–994.
gof_cal=function(pm,t0ratio=0.1,gof_method=NA,single_statistic="score",accumulate_option="max",filter=0,weight_option="none",weight=1){
  #step0: preparation
  if(!is.na(gof_method)){
    if(gof_method=="higher criticism"){
      single_statistic="score"
      accumulate_option="max"
    }
    if(gof_method=="skat"){
      single_statistic="score"
      accumulate_option="sum2"
    }
    if(gof_method=="berk_jones"){
      single_statistic="likelihood_ratio"
      accumulate_option="max"
    }
    if(gof_method=="einmahl_mckeague"){
      single_statistic="likelihood_ratio"
      accumulate_option="sum"
    }
  }
  if(single_statistic=="likelihood_ratio"){
    if(accumulate_option=="max2"){
      accumulate_option="max"
    }
    if(accumulate_option=="sum2"){
      accumulate_option="sum"
    }
    if(accumulate_option=="rsum2"){
      accumulate_option="rsum"
    }
    if(accumulate_option=="topsum2"){
      accumulate_option="topsum"
    }
  }

  #step1: arrange the matrix
  pm=as.matrix(pm)

  if(weight_option=="in"){
    pm=cal_cdf(pm,weight)
  }

  d=ncol(pm)
  n=nrow(pm)
  t0n=max(floor(n*t0ratio),1,na.rm=TRUE)
  if(n==1){return(pm)}
  else{
    pm=as.matrix(apply(pm,2,function(x){return(sort(x,na.last = TRUE))}))
    S=apply(pm,2,function(x){return(rank(x,na.last = "keep", ties.method = "max"))})
    if(single_statistic=="score"){
      preT=(S-pm*n)/sqrt(pm*(n-pm*n))
    }
    if(single_statistic=="wald"){
      preT=(S-pm*n)/sqrt(S*(n-S)/n)
    }
    if(single_statistic=="likelihood_ratio"){
      preT=n*(pm*log(pm*n/S)+(1-pm)*log((1-pm)*n/(n-S)))
    }


    if(weight_option=="out"){
      preT=preT*weight
    }

    preTm=preT
    preT[pm>=t0ratio]=NA
    preT[pm<=filter]=NA

    if(accumulate_option=="max"){
      T_star=suppressWarnings(apply(preT,2,function(x){return(max(x,na.rm=TRUE))}))
    }
    if(accumulate_option=="max2"){
      T_star=suppressWarnings(apply(preT*preT,2,function(x){return(max(x,na.rm=TRUE))}))
    }
    if(accumulate_option=="sum"){
      T_star=suppressWarnings(apply(preT,2,function(x){return(sum(x,na.rm=TRUE))})/n)
    }
    if(accumulate_option=="sum2"){
      T_star=suppressWarnings(apply(preT*preT,2,function(x){return(sum(x,na.rm=TRUE))})/n)
    }
    if(accumulate_option=="rsum"){
      #calculate rmean
      preTm=apply(preTm,2,function(x){return(cumsum(x)/seq(1:length(x)))})
      preTm[pm>=t0ratio]=NA
      preTm[pm<=filter]=NA
      T_star=suppressWarnings(apply(preTm,2,function(x){return(max(x,na.rm=TRUE))}))
    }
    if(accumulate_option=="rsum2"){
      #calculate rmean
      preTm=apply(preTm*preTm,2,function(x){return(cumsum(x)/seq(1:length(x)))})
      preTm[pm>=t0ratio]=NA
      preTm[pm<=filter]=NA
      T_star=suppressWarnings(apply(preTm,2,function(x){return(max(x,na.rm=TRUE))}))
    }
    if(accumulate_option=="topsum"){ #summing up preT which is over certain t0ratio proportion.
      preT=as.matrix(apply(preT,2,function(x){return(sort(x,decreasing = TRUE,na.last = TRUE))}))
      if((t0n<n)){
        preT[(t0n+1):n,]=NA
      }
      T_star=suppressWarnings(apply(preT,2,function(x){return(sum(x,na.rm=TRUE))}))
    }
    if(accumulate_option=="topsum2"){ #summing up preT which is over certain t0ratio proportion.
      preT=as.matrix(apply(preT*preT,2,function(x){return(sort(x,decreasing = TRUE,na.last = TRUE))}))
      if((t0n<n)){
        preT[(t0n+1):n,]=NA
      }
      T_star=suppressWarnings(apply(preT,2,function(x){return(sum(x,na.rm=TRUE))}))
    }
  }
  return(T_star)
}
