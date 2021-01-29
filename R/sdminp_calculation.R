#' This function calculates the Step Down minP method.
#'
#' It implement the improved step-down minP algorithm by Ge et al, 2003 (box 4)
#' @param res_ob is a numeric vectors of statistics
#' @param res_perm is a matrix with each colum is a permuated statistics with the same length as res_ob
#' @return a numeric vector adjusted pvalues
#' @export
#'
#' @examples
#' ob=rnorm(10,2,2)
#' perm=matrix(rnorm(100,2,2),ncol=10,nrow=10)
#' sdminp(ob,perm)
#' @references Ge, Y., Dudoit, S., & Speed, T. P. (2003). Resampling-based multiple testing for microarray data analysis. Test, 12(1), 1–77. https://doi.org/10.1007/BF02595811
sdminp=function(res_ob,res_perm){
  m=nrow(res_perm)#initialization
  B=ncol(res_perm)#initialization
  #box4. step 0
  print("step 0: calculate the raw-pvalue")
  rawp=apply((res_perm>=res_ob),1,sum)/B
  order_rawp=order(rawp)
  order_rawp_trans=match(seq(1,m,by=1),order_rawp)
  rawp=rawp[order_rawp]
  Tm=res_perm[order_rawp,] #initialization
  Pm=matrix(,ncol=B,nrow=m) #initialization
  Qm=Pm #initialization
  Pm[,B]=1
  Qm=rbind(Qm,1)
  #box4. step 1
  print("step 1-1: calculate the T matrix")
  orderTm=t(apply(Tm,1,function(x){order(x,decreasing = TRUE)}))
  trans_orderTm=t(apply(orderTm,1,function(x){match(seq(1,length(x),by=1),x)}))
  #sort
  for(i in 1:m){
    Tm[i,]=Tm[i,orderTm[i,]]
  }
  print("step 1-2: calculate the P matrix")
  for(j in (B-1):1){
    Pm[,j]=Pm[,(j+1)]
    Pm[(Tm[,j]!=Tm[,(j+1)]),j]=j/B
  }
  #sort back
  for(i in 1:m){
    Pm[i,]=Pm[i,trans_orderTm[i,]]
  }
  #box4. step 2
  print("step 2: calculate the Q matrix")
  for(i in m:1){
    Qm[i,]=pmin(Qm[(i+1),],Pm[i,])
  }
  #box4. step 3
  print("step 3")
  ep=apply((Qm[1:m,]<=rawp),1,sum)/B
  #box4. step 6
  print("step 6")
  for(i in 2:m){
    ep[i]=max(ep[(i-1):i])
  }
  #sort result to original order
  print("sort result to original order")
  ep=ep[order_rawp_trans]
  return(ep)
}



#' This function calculates the Step Down minP and report the results,
#' This is downstream analysis for sdminp
#' It implements the improved step-down minP algorithm based on FDR by Ge et al, 2003 (box 5)
#' It requires another function find_rank
#'
#'
#' @param res_ob is a numeric vectors of statistics
#' @param res_perm is a matrix with each colum is a permuated statistics with the same length as res_ob
#' @param tao0 is a proportion threshold that more than tao0, there will expected to be no significant results.
#' @return a numeric matrix with adjusted pvalues,
#'  the first column is the FDR adjusted pvalues,
#'  the second colum is the corresponding q-values
#'
#' @export
#' @examples
#' ob=rnorm(10,2,2)
#' perm=matrix(rnorm(100,2,2),ncol=10,nrow=10)
#' sdminp_fdr(ob,perm)
#' @references Ge, Y., Dudoit, S., & Speed, T. P. (2003). Resampling-based multiple testing for microarray data analysis. Test, 12(1), 1–77. https://doi.org/10.1007/BF02595811
sdminp_fdr=function(res_ob,res_perm,tao0=0.2){
  m=nrow(res_perm)#initialization
  B=ncol(res_perm)#initialization
  #box4. step 0
  print("step 0: calculate the ob-tao")
  tao=res_ob
  W0=sum(res_ob<=tao0,na.rm = TRUE)
  order_tao=order(tao,decreasing=TRUE)
  order_tao_trans=match(seq(1,m,by=1),order_tao)
  tao=tao[order_tao]

  R=matrix(ncol=1,nrow=m) #initialization
  Rm=matrix(ncol=B,nrow=m) #initialization
  qval=matrix(ncol=1,nrow=m) #initialization
  adj_pval=matrix(ncol=1,nrow=m) #initialization

  #box5. step 1
  print("step 1-1: calculate the R matrix")
  res_obS=sort(res_ob,decreasing = TRUE)
  res_permS=apply(res_perm,2,function(x){sort(x,decreasing = TRUE)})

  R=find_rank(res_obS,tao)
  for(j in 1:B){
    Rm[,j]=find_rank(res_permS[,j],tao)
  }
  print("step 1-2: calculate the W0 vector")
  W0m=as.matrix(apply((res_perm<=tao0),2,sum))

  #box5. step 2
  print("step 2: calculate meanR, meanI meanW ")
  meanR=apply(Rm,1,mean)
  meanI=apply(Rm>0,1,mean)
  meanW0=mean(W0)
  #box5. step 3
  print("step 3: calculate pFDR and FDR ")
  pFDR=W0*meanR/(meanW0*((R+1)+abs(R-1))*meanI/2)
  FDR=W0*meanR/(meanW0*((R+1)+abs(R-1))/2)

  #box5. step 4
  print("step 4 enforcing monotonicity")
  qval[m]=pFDR[m]
  adj_pval[m]=FDR[m]
  for(i in (m-1):1){
    qval[i]=min(qval[(i+1)],pFDR[i])
    adj_pval[i]=min(adj_pval[(i+1)],FDR[i])
  }

  #sort result to original order
  print("sort result to original order")
  qval=qval[order_tao_trans]
  adj_pval=adj_pval[order_tao_trans]
  return(cbind(adj_pval,qval))
}

#' This function calculates the ranks of a vector based on another vector
#'
#' find_rank is used by sdminp_fdr to speed up the calculation
#'
#' @param target a decreacing-sorted vector for measure by the ruler
#' @param ruler has the same length as target, also a decreasing-sorted vector
#' @return a numerica vector, with each element i shows the numbers of values in the target that is bigger than the ruler[i]
#'
#' @export
#' @examples
#' a=c(8,6,4,2)
#' b=c(7,4,3,1)
#' find_rank(a,b)
find_rank=function(target,ruler){
  m=length(ruler)
  res=matrix(ncol=1,nrow=m)
  curRank=1
  for(i in 1:m){
    while((target[curRank]>=ruler[i])&&(curRank<=m)){
      curRank=curRank+1
    }
    res[i]=curRank-1
  }
  return(res)
}


#' This function calculates the impirical raw p-values.
#'
#' It implements the raw-pvalues based on permutation, defined by Ge et al, 2003 (box 1)
#'
#' @param res_ob is a numeric vectors of statistics
#' @param res_perm is a matrix with each colum is a permuated statistics with the same length as res_ob
#' @return a numeric matrix with raw pvalues, defined by Ge et al, 2003 (box 1). NAs will be ignored.
#' @export
#'
#' @examples
#' ob=rnorm(4,2,2)
#' perm=matrix(rnorm(20,2,2),ncol=5,nrow=4)
#' rawp_cal(ob,perm)
#' @references Ge, Y., Dudoit, S., & Speed, T. P. (2003). Resampling-based multiple testing for microarray data analysis. Test, 12(1), 1–77. https://doi.org/10.1007/BF02595811
rawp_cal=function(res_ob,res_perm){
  rawp=apply((res_perm>=res_ob),1,sum)/ncol(res_perm)
  return(rawp)
}

