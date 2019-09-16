

# weighted Higher Criticism Instruction

This package is used for the application of the weighted higher criticism on the given gene sets. It provides codes that is able to calculate the pvalues based on single fisher exact tests or the stratified CMH tests (Epstein, M. P., Allen, A. S., & Satten, G. A. ,2007).

Here's the instructions for a standard workflow.

# Data preparation
The genotype matrix has the format that each row represents a gene and each column represents a subject.The rownames can be the standard gene symbols or the Ensembl ID. Each element of the genotype matrix should be 0 (reference) or 1 (variant).

The phenotype data should be the vector which matches the order of subjects of the genotype matrix, with the element as  0 (control) or 1 (case).

example:
```{r data_preparation}
library("wHC")
pheno=rbinom(100,1,0.5)
geno=matrix(rbinom(2000,1,0.1),nrow=20,ncol=100)
```

If the stratified CMH tests are applied, besides the genotype and phenotype, we still needs the external sample informations, for the population stratification estimation. This information will have each rows represents the subjects and match the order of subjects of the phenotypes. Refer (Epstein, M. P., Allen, A. S., & Satten, G. A. ,2007) for more details.

example:
```{r data_prep}
cur_pc=rbind(matrix(rnorm(500,0,1),ncol=10,nrow=50),matrix(rnorm(500,0.5,1),ncol=10,nrow=50))
```


# pvalue calculation
If the stratified CMH tests are applied,run strat_score_cal_glm to calculate the stratification formation based on generalized linear regression. Refer (Epstein, M. P., Allen, A. S., & Satten, G. A. ,2007) for more details.

example:
```{r pval_strata}
strata=strat_score_cal_glm(pheno,cur_pc)
```

Then pvalues is calculated based on fisher exact tests (if exact_option=TRUE) or cmh test.It can also calculated the permutated pvalues based on mismatching the phenotype labels to genotypes subjects when the parameter nperm, the number of permutation, is greater than 0.

example:
```{r pval_cal}
pval=pval_ss_cal(pheno,geno,strata,nperm=0, alter="greater",exact_option=TRUE)
```


# weight generation
The external information of each individual genes will be grabed from the following functions.Those information can based on the genetic network structures of genes within the given set.

example:
```{r w_centrality}
library("igraph")
interact_m=matrix(rbinom(100000,1,0.3),100,100)
deg=get_centrality(interact_m,w_option="deg",direct_option=FALSE,mode_option="all") #for degree
closn=get_centrality(interact_m,w_option="closn",direct_option=FALSE,mode_option="all") #for closeness
betn=get_centrality(interact_m,w_option="betn",direct_option=FALSE,mode_option="all") # for betweenness
eigen=get_centrality(interact_m,w_option="eigen",direct_option=FALSE,mode_option="all") #for eigenvector centrality
page_rank=get_centrality(interact_m,w_option="pagerank",direct_option=TRUE,mode_option="all") #for pagerank centrality
```

To be more specifically, we extract genes and all of their interactions reported in human beings from biogrid and gets the networks from the MSigDB (gmt format).

example:
```{r w_network}
human_whole=as.matrix(human_whole[,c(2,3,8,9,10,11)])
human_whole=unique(human_whole)
human_whole=as.matrix(human_whole[order(human_whole[,2]),])
human_whole=as.matrix(human_whole[order(human_whole[,1]),])
human_whole[,1]=as.numeric(human_whole[,1])
human_whole[,2]=as.numeric(human_whole[,2])
human_whole=human_whole[,1:2]
res=match_prior_info_centrality(net,human_whole,add_option="none",
report_option=TRUE,w_option="pagerank",direct_option=TRUE,mode_option="all")
```

Otherwise, we can get the gene intolerance information from the database of Petrovski, S. et al, 2013

example:
```{r w_genic_intolerance}
genic_intolerance=get_genic_intolerance()
```

Get gene expression data of specific tissues from GTEx

example:
```{r w_expression}
prior_expression=get_gene_expression(gene_label="symbols",tissue=c("Liver","Lung"),comb="mean")
prior_expression=get_gene_expression()
```

Get the estimated transcripts length

example:
```{r w_gene_length}
prior_length=get_gene_length(gene_label="symbols",comb="median")
```

And we can match those prior information to the given gene set from collection of [MSigDB](http://software.broadinstitute.org/gsea/msigdb)

example:
```{r w_net}
net=net.h.all.v6.1.symbols
#net from MSigDB (software.broadinstitute.org/gsea/msigdb): "h.all.v6.1.symbols.gmt"
prior_gi=get_genic_intolerance()
prior_net=match_prior_info(net,prior_gi)
```



# weighted higher criticism calculation.
The weight can be converted from the external information according to w=1/(a x prior_info+b x mean(prior_info)). It is then scaled into mean(w)=1. Here we take a=0.95 and b=0.05.

example:
```{r w_trans}
w0=trans_w(prior_gi)
```

Then the weighted pvalues is calculated and curved based on the cdf functions.
example:
```{r pressure}
pval=matrix(nrow=w0,ncol=5)
pval=as.matrix(apply(pval,2,function(x){return(runif(length(w0),0,1))}))
pwval=cal_cdf(pval,w=w0)
```

Finally, the higher criticism can be applied to the pvalues. Each colum represents a set of pvalues from the corresponding set of genes, and will return a single higher criticism score through our method.
example:
```{r pressure}
hc_cal(pwval,t0=0.4)
```