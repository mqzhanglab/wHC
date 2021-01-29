#' strat_score_cal_glm calculates the population stratification
#'
#' @param dis vector of disease outcomes (1=case, 0=control) with Ntot subjects.
#' @param pc first 10 principle components(defaut colum=10).
#' @param nstrat_ss number of strata to be separated to, defaut 5
#'
#' @return a numeric matrix with nstrat_ss-1 columns, dis length rows. each elements are 0-1 indicators which shows if this subjects belows to this strata. the last strata are subjects with all other nstrat_ss-1 columns be zero.
#'
#' @seealso \code{\link{glm}} which this function wraps
#' @export
#' @examples
#' pheno=rbinom(100,1,0.5)
#' cur_pc=rbind(matrix(rnorm(500,0,1),ncol=10,nrow=50),matrix(rnorm(500,0.5,1),ncol=10,nrow=50))
#' strata=strat_score_cal_glm(pheno,cur_pc)
#' @references Epstein MP, Allen AS, Satten GA (2007) A simple and improved correction for population stratification in case-control studies. American Journal of Human Genetics 80: 921-930
strat_score_cal_glm=function(dis,pc,nstrat_ss=5){
  Ntot=length(dis)

  #Form stratification score based on 10 principal components
  sscore_pc=glm(dis~pc[,1:10],family=binomial())

  #Assign subjects to 5 strata based on predicted odds
  sslp_pc=sscore_pc$linear.predictors
  qsslp_pc=quantile(sslp_pc,probs=seq(1,(nstrat_ss-1))/nstrat_ss,names=T)
  stratid_pc=matrix(1,ncol=1,nrow=Ntot)
  for(is in 1:(nstrat_ss-1)){
    stratid_pc=(stratid_pc+(sslp_pc>qsslp_pc[is]))
  }
  # Construct indicator functions for each stratum (stratum 5 is baseline category)
  stind_pc=matrix(nrow=Ntot,ncol=(nstrat_ss-1))
  for(k in 1:(nstrat_ss-1)){
    stind_pc[,k]=ifelse(stratid_pc==k,1,0)
  }
  return(stind_pc)
}


#' pval_ss_cal calculates pvalues of each genes between disease and test genotype
#'
#' It featuring for options of adjusting for the stratification-score indicators followed with cmh test, for observation and permutation. If no stratification, this code did simple fisher exact test instead.
#'
#'
#' @param dis_ob 0-1 vector shows the subjects has disease (1) or not (0)
#' @param g matrix with each row demonstrate the genes and each colum demonstrates the subjects. ncol(g)=length(dis_ob)
#' @param strata the output of function strat_score_cal, which are the n rows, k-1 colums matrix for samples stratification information. If ignored, only fisher exact test are applied.
#' @param nperm the number of permutation we should perform. If nperm=0 (default), that means only calculate the observed/input situation. Otherwise, only calculated the permutated situation.
#' @param alter only works if cmh_option is TRUE. The test for p-values calculation. options=c("greater","less","two.sided"): "greater"(default):apply one-side cmh test if mutation are enriched in cases."less":apply one-side cmh test if mutation are enriched in controls. "two.sided": apply two-side cmh test.
#' @param exact only works if cmh_option is TRUE. A logical indicating whether the Mantel-Haenszel test(FALSE) or the exact conditional test(TRUE, default) is applied.
#' @return a numeric vector of p-values
#'
#' @seealso \code{\link{mantelhaen.test}} which this function wraps
#' @seealso \code{\link{cmh_cal}} which this function wraps
#'
#' @export
#' @examples
#' pheno=rbinom(100,1,0.5)
#' cur_pc=rbind(matrix(rnorm(500,0,1),ncol=10,nrow=50),matrix(rnorm(500,0.5,1),ncol=10,nrow=50))
#' strata=strat_score_cal_glm(pheno,cur_pc)
#' geno=matrix(rbinom(2000,1,0.1),nrow=20,ncol=100)
#' pval=pval_ss_cal(pheno,geno,strata,nperm=0, alter="greater",exact_option=TRUE)
#' @references Epstein MP, Allen AS, Satten GA (2007) A simple and improved correction for population stratification in case-control studies. American Journal of Human Genetics 80: 921-930
#'
pval_ss_cal=function(dis_ob,g,strata=NA,nperm=0, alter="greater",exact=TRUE){
  # strata_vector:
  # for example, when we have 5 categories, and when
  # strata = [[1,0,0,0,0];[0,1,0,0,0];[0,0,1,0,0];[0,0,0,0,0]],
  # we have: strata_vector = (2,3,4,0).
  if(is.na(strata)){ #no strata
    k=1
    strata_vector = rep(0L, length(g))
  }
  if(!is.na(strata)){ #strata exists
    k = as.integer(ncol(strata)+1)
    strata_vector = rep(0L, nrow(strata))
    for (ic in 1:ncol(strata)) {
      strata_vector = strata_vector + strata[,ic] * ic
    }
    storage.mode(strata_vector) = "integer"
  }

  if(nperm==0){#that means we calculate observed situation
    pval_ss=matrix(nrow=nrow(g),ncol=1)
    for(ig in 1:nrow(g)){
      pval_ss[ig]=cmh_cal(dis_ob,g[ig,],strata_vector, k, alter_option=alter,exact_option=exact)
      if(ig%%100==0){
        print(c(nperm,ig))
      }
    }
  }
  if(nperm>0){#that means we calculate permutated
    pval_ss=matrix(nrow=nrow(g),ncol=nperm)
    for(ip in 1:nperm){
      # permuatationwith in each strata
      dis_perm=matrix(nrow=length(dis_ob),ncol=1)
      #strata 2-n
      for(is in 1:ncol(strata)){
        dis_ob_s=dis_ob[strata[,is]==1]
        dis_perm[strata[,is]==1]=dis_ob_s[sample.int(length(dis_ob_s))]
      }
      #strata 1
      dis_ob_s=dis_ob[is.na(dis_perm)]
      dis_perm[is.na(dis_perm)]=dis_ob_s[sample.int(length(dis_ob_s))]

      #calculate the pvlues
      for(ig in 1:nrow(g)){
        pval_ss[ig,ip]=cmh_cal(dis_perm,g[ig,],strata_vector, k, alter_option=alter, exact_option=exact)
        #if(ig%%100==0){
        #  print(c(ip,ig))
        #}
      }
    }
  }
  return(pval_ss)
}


#' cmh_cal is the function of the cmh test for single genes with stratas
#'
#' @description cmh_cal belongs to pval_ss_cal
#' when there is only 1 strata, cmh_cal reduces to fisher exact test and ignore the exact_option
#'
#' @param dis is a n-length numeric vector of indicators of phenotypes from n samples.
#' @param cur_g is a n-length numeric vector of indicators of mutations from n samples in a certain gene.
#' @param strata_vector is a n-length numeric vector of categories from the k kinds of strata, it is converted from the result of strat_score_cal_glm, for example, when we have 5 categories, and when strata = [[1,0,0,0,0];[0,1,0,0,0];[0,0,1,0,0];[0,0,0,0,0]], we have: strata_vector = (2,3,4,0).
#' @param k is the kinds of strata, in accordance with content of strata_vector
#' @param alter_option decides the test for p-values calculation. options=c("greater","less","two.sided"):
#' "greater"(default):apply one-side cmh test if mutation are enriched in cases.
#' "less":apply one-side cmh test if mutation are enriched in controls.
#' "two.sided": apply two-side cmh test
#' @param exact_option A logical indicating whether the Mantel-Haenszel test(FALSE) or the exact conditional test (TRUE, default) is applied.If k==1, exact_option will be ignored.
#' @return a p-value from cmh test.
#' @seealso \code{\link{mantelhaen.test}} which this function wraps
#'
#' @export
#' @examples
#' pheno=rbinom(100,1,0.5)
#' cur_pc=rbind(matrix(rnorm(500,0,1),ncol=10,nrow=50),matrix(rnorm(500,0.5,1),ncol=10,nrow=50))
#' strata=strat_score_cal_glm(pheno,cur_pc)
#' cur_geno=rbinom(100,1,0.1)
#' cur_pval=cmh_cal(pheno,cur_geno,strata, alter="greater",exact_option=TRUE)
#' @references William G. Cochran (December 1954). "Some Methods for Strengthening the Common chi-squared Tests". Biometrics. 10 (4): 417–451. doi:10.2307/3001616. JSTOR 3001616.
#' @references Nathan Mantel and William Haenszel (April 1959). "Statistical aspects of the analysis of data from retrospective studies of disease". Journal of the National Cancer Institute. 22 (4): 719–748. doi:10.1093/jnci/22.4.719. PMID 13655060.
#' @references Fisher, R. A. (1922). "On the interpretation of chi-squared from contingency tables, and the calculation of P". Journal of the Royal Statistical Society. 85 (1): 87–94. doi:10.2307/2340521. JSTOR 2340521.
cmh_cal=function(dis,cur_g,strata_vector, k=1, alter_option="greater",exact_option=TRUE){
  # k=ncol(strata)+1
  # #1. make the array for cmh test
  mat = get_2x2xk(dis, cur_g, strata_vector, k)
  #2.do cmh test
  if(k==1){
    res=fisher.test(mat[,,1], alternative = alter_option)$p.value
  }
  if(k>1){
    res=mantelhaen.test(mat, alternative = alter_option, exact = exact_option)$p.value
  }

  return(res)
}

