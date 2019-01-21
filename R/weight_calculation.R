#' match_prior_info_centrality matches the centrality prior information to the given gene set from collection of MSigDB
#'
#' Different from match_prior_info, the match_prior_info_centrality calculates each interactions based on the true gene sets.
#'
#' @param net dataframe the gmt file from collections of MSigDB, broad institute. each line represents a pathway.please read in with read.csv with header=FALSE and stringAsFactors = FALSE.
#' @param human_whole the 2-column matrix with each line representing the connection from gene in column 1 to gene in column 2
#' @param add_option defines the method of adding up missing values. It can be "none", "mean" or "median". No actions for adding up missing values if "none".
#' @param report_option if TRUE, report current procedures of the path, which is the proportion of sets completed the matching steps.
#' @param w_option the kind of centralities.
#' @param direct_option if it is true, the network will be calculated as directed pathways, parameter especially for pagerank
#' @param mode_option  parameters for centrality calculation, "out" for out-degree, "in" for in-degree or "all" or "total" for the sum of the two.
#'
#' @return a dataframe with the same format as net,which is the gmt files
#' @import igraph
#' @export
#' @seealso \code{\link{igraph}} which this function wraps
#' @examples
#' net=net.h.all.v6.1.entrez
#' #from MSigDB (http://software.broadinstitute.org/gsea/msigdb): "h.all.v6.1.entrez.gmt"
#' human_whole=human_whole_biogird_3.4.147
#' #from BioGRID(https://thebiogrid.org/): "BIOGRID-ORGANISM-Homo_sapiens-3.4.147.tab2.txt"
#' human_whole=as.matrix(human_whole[,c(2,3,8,9,10,11)])
#' human_whole=unique(human_whole)
#' human_whole=as.matrix(human_whole[order(human_whole[,2]),])
#' human_whole=as.matrix(human_whole[order(human_whole[,1]),])
#' human_whole[,1]=as.numeric(human_whole[,1])
#' human_whole[,2]=as.numeric(human_whole[,2])
#' human_whole=human_whole[,1:2]
#' res=match_prior_info_centrality(net,human_whole,add_option="none",
#' report_option=TRUE,w_option="pagerank",direct_option=TRUE,mode_option="all")
#' @references Liberzon, A., Subramanian, A., Pinchback, R., Thorvaldsdóttir, H., Tamayo, P., & Mesirov, J. P. (2011). Molecular signatures database (MSigDB) 3.0. Bioinformatics, 27(12), 1739–1740. https://doi.org/10.1093/bioinformatics/btr260
#' @references Page, L., Brin, S., Motwani, R., & Winograd, T. (1999). The PageRank citation ranking: Bringing order to the web. Stanford InfoLab. Retrieved from http://ilpubs.stanford.edu:8090/422
#' @references White, S., & Smyth, P. (2003). Algorithms for Estimating Relative Importance in Networks. In Proceedings of the Ninth ACM SIGKDD International Conference on Knowledge Discovery and Data Mining (pp. 266–275). New York, NY, USA: ACM. https://doi.org/10.1145/956750.956782
#' @references Chatr-aryamontri, A., Oughtred, R., Boucher, L., Rust, J., Chang, C., Kolas, N. K., … Tyers, M. (2017). The BioGRID interaction database: 2017 update. Nucleic Acids Research, 45(Database issue), D369–D379. https://doi.org/10.1093/nar/gkw1102
match_prior_info_centrality=function(net,human_whole,add_option="none",report_option=TRUE,w_option="deg",direct_option=FALSE,mode_option="all"){
  n=nrow(net)
  set_prior=matrix(nrow=n,ncol=1)
  w=matrix(ncol=1,nrow=0)
  for(i in 1:n){
    #1.match and threshold interactions.
    path_total=unlist(strsplit(net[i,],"\t"))
    gene_ind=path_total[3:length(path_total)]

    gene_id=cbind(match(human_whole[,1],gene_ind),match(human_whole[,2],gene_ind))
    match_gene=1-is.na(gene_id[,1])-is.na(gene_id[,2])
    cur_prior=rep(1,length(match_gene))
    #at least 5 edges
    if(sum(match_gene>0)>4){
      gene_id=gene_id[match_gene>0,] #interact 1 and 2 are both in the network
      #2.#extract Graph Info #direct graph as set
      gene_m=matrix(0,ncol=length(gene_ind),nrow=length(gene_ind))
      gene_m_rev=matrix(0,ncol=length(gene_ind),nrow=length(gene_ind)) #especially for pageRank
      for(i_g in 1:nrow(gene_id)){
        gene_m[gene_id[i_g,1],gene_id[i_g,2]]=1
      }
      cur_prior=get_centrality(gene_m,w_option=w_option,direct_option=direct_option, mode_option=mode_option)
      if(add_option=="median"){
        cur_prior[is.na(cur_prior)]=median(cur_prior,na.rm = TRUE)
      }
      if(add_option=="mean"){
        cur_prior[is.na(cur_prior)]=mean(cur_prior,na.rm = TRUE)
      }
    }
    cur_prior=c(path_total[1:2],cur_prior)
    set_prior[i]=paste(cur_prior,collapse = "\t")
    if(report_option==TRUE && i%%100==0){
      print(paste("complete ",i,"/",n))
    }
  }
  set_prior=as.data.frame(set_prior,col.names=NULL,stringsAsFactors = FALSE)
  return(set_prior)
}


#' get_centrality match the centrality prior information to the given gene set from collection of MSigDB
#'
#' Different from match_prior_info, the match_prior_info_centrality calculates each interactions based on the true gene sets.
#'
#' @param interact_m represents the genetic network. It is the adjacency matrix of a graph with each node represents a gene.
#' @param w_option the kind of centralities.can be "deg" for degree,"closn" for closeness,"betn" for betweenness, "eigen" for eigervector centrality, and "pagerank"
#' @param direct_option if it is true, the network will be calculated as directed pathways, parameter especially for pagerank
#' @param mode_option  parameters for centrality calculation, "out" for out-degree, "in" for in-degree or "all" or "total" for the sum of the two. see \code{\link{igraph}} for more details.
#'
#' @references Page, L., Brin, S., Motwani, R., & Winograd, T. (1999). The PageRank citation ranking: Bringing order to the web. Stanford InfoLab. Retrieved from http://ilpubs.stanford.edu:8090/422
#' @references White, S., & Smyth, P. (2003). Algorithms for Estimating Relative Importance in Networks. In Proceedings of the Ninth ACM SIGKDD International Conference on Knowledge Discovery and Data Mining (pp. 266–275). New York, NY, USA: ACM. https://doi.org/10.1145/956750.956782
get_centrality=function(interact_m,w_option="deg",direct_option=FALSE,mode_option="all"){
  if(direct_option==TRUE){
    net=graph.adjacency(interact_m,mode="directed")
  }
  if(direct_option==FALSE){
    net=graph.adjacency(interact_m,mode="undirected")
  }

  if(w_option=="deg"){
    network_info=degree(net, v = V(net), mode = mode_option, normalized = FALSE)
  }
  if(w_option=="closn"){
    network_info=closeness(net, vids = V(net), mode = mode_option, weights = NULL, normalized = FALSE)
  }
  if(w_option=="betn"){
    network_info=betweenness(net, v = V(net), weights = NULL, normalized = FALSE)
  }
  if(w_option=="eigen"){
    network_info=evcent(net, directed=direct_option, weights=NULL)$vector
  }
  if(w_option=="pagerank"){
    network_info=page_rank(net, algo = "prpack", vids = V(net), directed = direct_option, damping = 0.85, personalized = NULL, weights = NULL, options = NULL)$vector
  }
  return(network_info)
}



#' This function pull down the genic intolerance information
#'
#' The data resources comes from the databased of genic intolerance.
#' (genic-intolerance.org/data/RVIS_Unpublished_ExAC_May2015.txt).
#' The default gene symbol is "CCDS_r15" and "%RVIS_ExAc_0.05%(AnyPopn)"(the percentage of genic intolerance)
#' The genes are labeled with gene symbols.
#' @return the numeric vector representing prior information for each single gene, with genes name corresponding to that in the net, which can be used as input of function match_prior_info
#'
#' @export
#' @examples
#' genic_intolerance=get_genic_intolerance()
#'
#' @references Petrovski, S., Wang, Q., Heinzen, E. L., Allen, A. S., & Goldstein, D. B. (2013). Genic Intolerance to Functional Variation and the Interpretation of Personal Genomes. PLOS Genetics, 9(8), e1003709. https://doi.org/10.1371/journal.pgen.1003709
#'
get_genic_intolerance=function(){
  gene_intolerance_info=read.table("http://genic-intolerance.org/data/RVIS_Unpublished_ExAC_May2015.txt",header = TRUE)
  prior_info=as.numeric(gene_intolerance_info[,7])
  names(prior_info)=gene_intolerance_info[,5]
  return(prior_info)
}


#' This function get gene expression data of specific tissues from GTEx
#'
#' The data resources comes from tpm of genes counts of RNAseq data of GTEx (gtexportal.org/home/datasets).
#' It based on GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz (Retrive Nov (2013))
#' @param gene_label is the gene names of the returning vector or matrix.It can be "ensembl_gene_id" or "symbols"
#' @param tissue a vector of filters the expression level of specific tissues. Use \code{\link{support_gtex_tissue()}} to see supported tissues.
#' @param comb is the operation on combining selected categories of tissues.
#' It can be "median","mean","max","min",and "none", which calculate the median, mean of genes or do nothing on them.
#' @return the numeric vector or matrix(@param mode is "none") representing prior information for each single gene
#'
#' @export
#' @examples
#' prior_expression=get_gene_expression(gene_label="symbols",tissue=c("Brain","Lung"),comb="max")
#' prior_expression=get_gene_expression()
#' @references Lonsdale, John, et al. "The genotype-tissue expression (GTEx) project." Nature genetics 45.6 (2013): 580.
#'
get_gene_expression=function(gene_label="ensembl_gene_id",tissue=support_gtex_tissue(),comb="none"){
  res=as.matrix(whole_tpm_median[,tissue])

  if(gene_label=="ensembl_gene_id"){
    gene_name=whole_tpm_median[,1]
  }
  if(gene_label=="symbols"){
    gene_name=whole_tpm_median[,2]
  }

  if(comb=="median"){
    res=as.numeric(apply(res,1,function(x){return(median(as.numeric(x),na.rm = TRUE))}))
    names(res)=gene_name
  }
  if(comb=="mean"){
    res=as.numeric(apply(res,1,function(x){return(mean(as.numeric(x),na.rm = TRUE))}))
    names(res)=gene_name
  }
  if(comb=="max"){
    res=as.numeric(apply(res,1,function(x){return(max(as.numeric(x),na.rm = TRUE))}))
    names(res)=gene_name
  }
  if(comb=="min"){
    res=as.numeric(apply(res,1,function(x){return(min(as.numeric(x),na.rm = TRUE))}))
    names(res)=gene_name
  }
  if(comb=="none"){
    res=as.matrix(apply(res,2,as.numeric))
    rownames(res)=gene_name
  }
  return(res)
}

#' support_gtex_tissue provides supported option for "tissue" in function get_gene_expression
#'
#' The data resources comes from tpm of genes counts of RNAseq data of "https://gtexportal.org/home/datasets".
#' @seealso \code{\link{get_gene_expression}}
#'
#' @export
#' @examples
#' support_gtex_tissue()
#'
#' @references Lonsdale, John, et al. "The genotype-tissue expression (GTEx) project." Nature genetics 45.6 (2013): 580.
#' @references https://storage.googleapis.com/gtex_analysis_v7/rna_seq_data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz (Retrive Nov (2013)).
support_gtex_tissue=function(){
  return(colnames(whole_tpm_median)[3:ncol(whole_tpm_median)])
}











#' This function estimates the transcripts length.
#'
#' The data resources comes from the ensembl genome browser (useast.ensembl.org/index.html).
#' @param gene_label is the returning gene labels, it can be "ensembl_gene_id" or "symbols"
#' @param comb is the estimation method based on multiple records of transcription lengths.
#' It can be "median","mean","max","min", which calculate the median, mean,max and min of genes or do nothing on them.
#' @return the numeric vector representing prior information for each single gene's length
#' @import biomaRt
#' @seealso \code{\link{biomaRt}}
#' @export
#' @examples
#' prior_length=get_gene_length(gene_label="symbols",comb="median")
#' @references Durinck, Steffen, et al. "Mapping identifiers for the integration of genomic datasets with the R/Bioconductor package biomaRt." Nature protocols 4.8 (2009): 1184.
#'
get_gene_length=function(gene_label="symbols",comb="mean"){
  ensembl = biomaRt::useMart("ensembl")
  human.ensembl = biomaRt::useDataset("hsapiens_gene_ensembl",mart=ensembl)
  ##1. get gene basic information #
  gene_name_list=biomaRt::getBM(attributes = c("ensembl_gene_id","hgnc_symbol"), filters="", values="",mart=human.ensembl)
  ##2. Get gene length #
  gene_name_list_length=biomaRt::getBM(attributes = c("ensembl_gene_id","transcript_length"), filters="", values="",mart=human.ensembl)

  match_index=match(gene_name_list_length$ensembl_gene_id,gene_name_list$ensembl_gene_id)
  match_index=match_index[!is.na(match_index)]
  gene_length_comb=matrix(ncol=1,nrow=nrow(gene_name_list))
  for(ig in 1:nrow(gene_name_list)){
    cur_ig=gene_name_list_length[match_index==ig,]
    if(nrow(cur_ig)>0){ #here the NAs happened since the situation that the ensembl_gene_id shows twice in gene_name_list. Here we take the first appared record since the ID of length list matches the first in the name_list.
      if(comb=="median"){
        gene_length_comb[ig]=median(cur_ig[,2],na.rm = TRUE)
      }
      if(comb=="mean"){
        gene_length_comb[ig]=mean(cur_ig[,2],na.rm = TRUE)
      }
      if(comb=="max"){
        gene_length_comb[ig]=max(cur_ig[,2],na.rm = TRUE)
      }
      if(comb=="min"){
        gene_length_comb[ig]=min(cur_ig[,2],na.rm = TRUE)
      }
    }
  }
  gene_length_comb=as.numeric(gene_length_comb)

  ##3. match gene length#
  if(gene_label=="ensembl_gene_id"){
    nna_flag=(!is.na(gene_name_list[,"ensembl_gene_id"])) & (!is.na(gene_length_comb))
    gene_length_comb=gene_length_comb[nna_flag]
    gene_name_list=gene_name_list[nna_flag,]
    names(gene_length_comb)=gene_name_list[,"ensembl_gene_id"]
  }
  if(gene_label=="symbols"){
    nna_flag=(!is.na(gene_name_list[,"hgnc_symbol"])) & (!is.na(gene_length_comb))
    gene_length_comb=gene_length_comb[nna_flag]
    gene_name_list=gene_name_list[nna_flag,]
    names(gene_length_comb)=gene_name_list[,"hgnc_symbol"]
  }
  return(gene_length_comb)
}











#' match_prior_info match the prior information to the given gene set from collection of MSigDB
#'
#' @param net dataframe the gmt file from collections of Molecular signatures database (MSigDB), broad institute.
#' @param prior_info the numeric vector representing prior information for each single gene, with genes name corresponding to that in the net.
#' @param add_option defines the method of adding up missing values. It can be "none", "mean" or "median". No actions for adding up missing values if "none".
#' @param report_option if TRUE, report current procedures of the path, which is the proportion of sets completed the matching steps.
#'
#' @return a dataframe with the same format as net,which is the gmt files
#'
#' @export
#' @examples
#' net=net.h.all.v6.1.symbols
#' #net from MSigDB (software.broadinstitute.org/gsea/msigdb): "h.all.v6.1.symbols.gmt"
#' prior_gi=get_genic_intolerance()
#' prior_net=match_prior_info(net,prior_gi)
#' @references Liberzon, A., Subramanian, A., Pinchback, R., Thorvaldsdóttir, H., Tamayo, P., & Mesirov, J. P. (2011). Molecular signatures database (MSigDB) 3.0. Bioinformatics, 27(12), 1739–1740. https://doi.org/10.1093/bioinformatics/btr260
match_prior_info=function(net,prior_info,add_option="none",report_option=TRUE){
  n=nrow(net)
  prior_value=as.numeric(prior_info)
  prior_name=names(prior_info)
  set_prior=matrix(,nrow=n,ncol=1)

  w=matrix(ncol=1,nrow=0)
  for(i in 1:n){
    #1.match and threshold interactions.
    path_total=unlist(strsplit(net[i,],"\t"))
    gene_ind=path_total[3:length(path_total)]
    cur_prior=as.numeric(prior_value[match(gene_ind,prior_name)])
    if(add_option=="median"){
      cur_prior[is.na(cur_prior)]=median(cur_prior,na.rm = TRUE)
    }
    if(add_option=="mean"){
      cur_prior[is.na(cur_prior)]=mean(cur_prior,na.rm = TRUE)
    }
    cur_prior=c(path_total[1:2],cur_prior)
    set_prior[i]=paste(cur_prior,collapse = "\t")
    if(report_option==TRUE && i%%100==0){
      print(paste("complete ",i,"/",n))
    }
  }
  set_prior=as.data.frame(set_prior,col.names=NULL,stringsAsFactors = FALSE)
  return(set_prior)
}
