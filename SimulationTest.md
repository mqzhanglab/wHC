
# The Simulation for power estimation of the Goodness-of-fit Statistics

This is the simulation protocol in support for the paper:

Zhang, Mengqi, et al. "Incorporating external information to improve sparse signal detection in rare‐variant gene‐set‐based analyses." Genetic epidemiology 44.4 (2020): 330-338.
Zhang, Mengqi, et al. "Focused Goodness of Fit Tests for Gene Set Analyses" Submitted in 2021.

The citation is appreciated.

### R package and Functions 

```
if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}

library("wHC")

#Typically, we use "two.sided" or "greater"
norm2pval=function(X0,alter="two.sided"){ #c("two.sided","greater","less")
  if(alter=="two.sided"){
    return(pnorm(-abs(X0)))
  }
  if(alter=="greater"){
    return(pnorm(-(X0)))
  }
  if(alter=="less"){
    return(pnorm((X0)))
  }
} 

RankUp <-function(x,x0){
  return(sum(x>x0)/length(x))
}
```

###  Parameters 


The following parameters provide an example senarios:
```
datasetNum=1000
n=300
pval_method="greater" #"two.sided" #this is for converting norm distributions to pvalues

alpha=1
norm_method=pval_method #this is for generating norm distributions

gof_significant=0.95

r=seq(0.05,1,by=0.05)
Gamma=seq(0.25,1,by=0.05)

mu_seq=sqrt(2*log(n)*r)
pi0_seq=n^-Gamma

single_stat="score"
accumulate_option="max"
```

### Simulation of the situation under H0  
The null situation is used for providing the alpha-level under the null.
Here we take the higher criticism as the example, which is with the parameter 
single_stat="score"
accumulate_option="max"

```
X_u0m=matrix(,ncol=datasetNum,nrow=n)#initialization
p_u0m=matrix(,ncol=datasetNum,nrow=n)#initialization

for(i_d in 1:datasetNum){
  X_u0=rnorm(n, 0, 1)
  p_u0m[,i_d]=norm2pval(X_u0,pval_method)
}

gof_u0=gof_cal(p_u0m,t0ratio=alpha,gof_method=NA,single_statistic=single_stat,accumulate_option=accumulate_option,filter=0,weight_option="none",weight=1)
thres_u0=quantile(gof_u0, probs=gof_significant)
```

### Simulation of the situation under H1
The H1 situation calculate various situation of the mu & pi0 parameter and calculation the method detection power for them. 
Here we take the higher criticism as the example, which is with the parameter 
single_stat="score"
accumulate_option="max"

```
quantile_u1=matrix(,nrow=length(mu_seq),ncol=length(pi0_seq))#initialization
rownames(quantile_u1)=mu_seq
colnames(quantile_u1)=pi0_seq


for(i_mu in 1:length(mu_seq)){
  for(i_pi in 1:length(pi0_seq)){
    p_u1m=matrix(,ncol=datasetNum,nrow=n)#initialization
    #calculate HC
    mu=mu_seq[i_mu]
    pi_0=pi0_seq[i_pi]
    for(i_d in 1:datasetNum){
      cur_mu=mu
      cur_pi_0=pi_0
      flag_u1=c(rep(1,round(cur_pi_0*n)),rep(0,c(n-round(cur_pi_0*n))))
      if(norm_method=="two.sided"){ #H1: pi/2*N(mu,1)+ pi/2*N(-mu,1)+(1-pi)*N(0,1)
        flag_u1[1:round(cur_pi_0*n/2)]=-1
      }
      N_h0=rnorm(n, mean = 0, sd = 1)
      N_h1=rnorm(n, mean = cur_mu, sd = 1)
      #unweighted
      X_u1=flag_u1*N_h1+(1-abs(flag_u1))*N_h0
      p_u1m[,i_d]=norm2pval(X_u1,pval_method)
      
      
    }
    gof_u1=gof_cal(p_u1m,t0ratio=alpha,gof_method=NA,single_statistic=single_stat,accumulate_option=accumulate_option,filter=0,weight_option="none",weight=1)
    
    #single quantile
    quantile_u1[i_mu,i_pi]=RankUp(gof_u1, thres_u0)
    
    print(paste("i_mu=",i_mu," i_pi=",i_pi))
  }
}
```

