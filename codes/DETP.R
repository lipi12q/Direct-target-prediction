fun_target <- function(querydata){
  
# function used to cacluate association scores between gene expression profiles and gene module pairs
# The caculation is modified from GSEA (Proc Natl Acad Sci U S A. 2005;102(43):15545-50)
  
  nsets<- length(target_dn_symbol[,1])
  genenames<- rownames(querydata)
  n_sample <- ncol(querydata)
  A <-apply(querydata, 2,  function(x){order(x,decreasing = TRUE)})  
  TES <- matrix(ncol = n_sample,nrow = nsets)
  
  for(q in 1:nsets){
    up_gene_order <- match(target_up_symbol[q,], genenames)
    up_gene_order <- up_gene_order[!is.na(up_gene_order)]
    dw_gene_order <- match(target_dn_symbol[q,], genenames)
    dw_gene_order <- dw_gene_order[!is.na(dw_gene_order)]
    for (p in 1:n_sample) {
      gene.list2 <- A[,p]
      upES <- unlist(GSEA.EnrichmentScore2(gene.list = gene.list2, gene.set = up_gene_order, weighted.score.type = 0))
      downES <- unlist(GSEA.EnrichmentScore2(gene.list = gene.list2, gene.set = dw_gene_order, weighted.score.type = 0))
      TES[q,p] <- upES - downES
    }
  }
  colnames(TES) <- colnames(querydata)
  rownames(TES) <- rownames(target_dn_symbol)
  return(TES)
}



#cpp_gene_dat_match2564 <- fread("../CMap_Rproj/cpp_gene_data_match2564_2.csv",header = T,)
#cpp_gene_dat_match2564[1:5,1:4]
#cpp_gene_dat_match2564 <- column_to_rownames(cpp_gene_dat_match2564 ,var = "V1")
#library(dplyr)
#library(tibble)


#cp_geneinfo_beta.txt <- fread("D:/下载/geneinfo_beta (1).txt",sep = "\t",header = T)
#cpp_gene_dat_match2564 <- cpp_gene_dat_match2564[rownames(cpp_gene_dat_match2564) %in% 
                                                   cp_geneinfo_beta.txt$gene_symbol[cp_geneinfo_beta.txt$feature_space != 'inferred'],]


load('D:/data/Cmap_all_coms_geneexp.Rdata') # read gene expression profiles of all CMap compounds. This file is too large and can be obtained on request.
load('D:/data/target_up_symbol.Rdata') # up gene set of target gene module pair
load('D:/data/target_dn_symbol.Rdata') # down gene set of target gene module pair

res <- fun_target(Cmap_all_coms_geneexp)

save(res,file = "Cmap_com_target_DES.Rdata") # DES values between all CMap compounds and targets





fun_target <- function(querydata,target_dn_symbol,target_up_symbol,target_phi){
  
  # cpp_gene_dat_match2564 <- fread("../CMap_Rproj/cpp_gene_data_match2564_2.csv",header = T,)
  # cpp_gene_dat_match2564 <- column_to_rownames(cpp_gene_dat_match2564 ,var = "V1")
  # load('D:/Py/Cmap/CMap_Rproj/target_phi.Rdata')
  # load('D:/Py/Cmap/CMap_Rproj/target_up_symbol.Rdata')
  # load('D:/Py/Cmap/CMap_Rproj/target_dn_symbol.Rdata')
  nsets<- length(target_dn_symbol[,1])
  genenames<- rownames(querydata)
  n_sample <- ncol(querydata)
  A <-apply(querydata, 2,  function(x){order(x,decreasing = TRUE)})  
  TES <- matrix(ncol = n_sample,nrow = nsets)
  
  for(q in 1:nsets){
    up_gene_order <- match(target_up_symbol[q,], genenames)
    up_gene_order <- up_gene_order[!is.na(up_gene_order)]
    dw_gene_order <- match(target_dn_symbol[q,], genenames)
    dw_gene_order <- dw_gene_order[!is.na(dw_gene_order)]
    for (p in 1:n_sample) {
      gene.list2 <- A[,p]
      upES <- unlist(GSEA.EnrichmentScore2(gene.list = gene.list2, gene.set = up_gene_order, weighted.score.type = 0))
      downES <- unlist(GSEA.EnrichmentScore2(gene.list = gene.list2, gene.set = dw_gene_order, weighted.score.type = 0))
      TES[q,p] <- upES - downES
    }
  }
  
  ###2、计算random TESvalue
  rows <- length(genenames)
  nperm <- 1000
  # phi <- matrix(ncol = 1000,nrow = nsets)
  # for(q in 1:nsets){
  #   up_gene_order <- match(target_up_symbol[q,], genenames)
  #   up_gene_order <- up_gene_order[!is.na(up_gene_order)]
  #   dw_gene_order <- match(target_dn_symbol[q,], genenames)
  #   dw_gene_order <- dw_gene_order[!is.na(dw_gene_order)]
  #   for (p in 1:1000) {
  #     reshuffled.gene.labels <- sample(1:rows)#随机的基因列的order
  #     #GSEA.EnrichmentScore2 是GSEA计算函数
  #     upES <- unlist(GSEA.EnrichmentScore2(gene.list = reshuffled.gene.labels, gene.set = up_gene_order, weighted.score.type = 0))
  #     downES <- unlist(GSEA.EnrichmentScore2(gene.list = reshuffled.gene.labels, gene.set = dw_gene_order, weighted.score.type = 0))
  #     phi[q,p] <- upES - downES
  #   }
  # } 
  # target_phi <- phi
  # save(target_phi,file = 'target_phi.Rdata')
  phi <- target_phi[1:8735,]
  ###3、计算 p value  
  pvalue <- matrix(nrow = nsets,ncol = n_sample)
  for(i in 1:n_sample){
    for (j in 1:nsets) {
      pvalue[j,i] <- ifelse(TES[j,i] >= 0, sum(phi[j,] >= TES[j,i])/sum(phi[j,]>=0), 
                            sum(phi[j,] < TES[j,i])/sum(phi[j,]<0))
    } 
  }
  ###3、计算NESvalue and fdr value
  # 
  phi.pos.mean <- apply(phi,1,function(x){
    mean(x[x>=0])
  })
  phi.neg.mean <- apply(phi,1,function(x){
    mean(abs(x[x<0]))
  })
  
  phi.norm.pos <- phi
  phi.norm.pos[phi<0]<-0
  phi.norm.pos <-phi.norm.pos/phi.pos.mean
  phi.norm.neg <- phi
  phi.norm.neg[phi>=0]<-0
  phi.norm.neg <-phi.norm.neg/phi.neg.mean
  phi.norm <- phi
  phi.norm[phi<0] <- phi.norm.neg[phi<0]
  phi.norm[phi>=0] <- phi.norm.pos[phi>=0]
  #r#r#r#r#r#r#r#r#r#r#r#r#r#r
  ES_query <- TES
  ES_query.pos <- ES_query
  ES_query.pos[ES_query<0]<-0
  ES_query.pos <-ES_query.pos/phi.pos.mean
  ES_query.neg <- ES_query
  ES_query.neg[ES_query>=0]<-0
  ES_query.neg <-ES_query.neg/phi.neg.mean
  
  ES_query.norm <- ES_query
  ES_query.norm[ES_query<0] <- ES_query.neg[ES_query<0]
  ES_query.norm[ES_query>=0] <- ES_query.pos[ES_query>=0]
  
  FDRvalue <- matrix(nrow = nsets,ncol = n_sample)
  
  for(i in 1:n_sample){
    for (j in 1:nsets) {
      if(ES_query.norm[j,i] >= 0){
        A <- sum(phi.norm >= ES_query.norm[j,i])/sum(phi.norm>=0)
        B <- sum(ES_query.norm >= ES_query.norm[j,i])/sum(ES_query.norm>=0)
        FDRvalue[j,i] <- A/B 
      }
      else{
        A <- sum(phi.norm <= ES_query.norm[j,i])/sum(phi.norm<0)
        B <- sum(ES_query.norm <= ES_query.norm[j,i])/sum(ES_query.norm<0)
        FDRvalue[j,i] <- A/B
      } 
    }                   
  }
  
  query_target_gsea <- cbind(TES,ES_query.norm, pvalue, FDRvalue)
  rownames(query_target_gsea) <- rownames(target_dn_symbol)
  colnames(query_target_gsea) <- c(paste(rep(colnames(querydata),times=4),rep(c("_CS","_NCS","_Pvalue","_FDR"),each=n_sample),sep = ""))
  gc()
  return(query_target_gsea)
}



#GSEA.EnrichmentScore2 ####

GSEA.EnrichmentScore <- function(gene.list, gene.set, weighted.score.type = 1, correl.vector = NULL) {  
  #
  # Computes the weighted GSEA score of gene.set in gene.list. 
  # The weighted score type is the exponent of the correlation 
  # weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted). When the score type is 1 or 2 it is 
  # necessary to input the correlation vector with the values in the same order as in the gene list.
  #
  # Inputs:
  #   gene.list: The ordered gene list (e.g. integers indicating the original position in the input dataset)  
  #   gene.set: A gene set (e.g. integers indicating the location of those genes in the input dataset) 
  #   weighted.score.type: Type of score: weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted)  
  #  correl.vector: A vector with the coorelations (e.g. signal to noise scores) corresponding to the genes in the gene list 
  #
  # Outputs:
  #   ES: Enrichment score (real number between -1 and +1) 
  #   arg.ES: Location in gene.list where the peak running enrichment occurs (peak of the "mountain") 
  #   RES: Numerical vector containing the running enrichment score for all locations in the gene list 
  #   tag.indicator: Binary vector indicating the location of the gene sets (1's) in the gene list 
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.
  
  tag.indicator <- sign(match(gene.list, gene.set, nomatch=0))    # notice that the sign is 0 (no tag) or 1 (tag) 
  no.tag.indicator <- 1 - tag.indicator 
  N <- length(gene.list) 
  Nh <- length(gene.set) 
  Nm <-  N - Nh 
  if (weighted.score.type == 0) {
    correl.vector <- rep(1, N)
  }
  alpha <- weighted.score.type
  correl.vector <- abs(correl.vector**alpha)
  sum.correl.tag    <- sum(correl.vector[tag.indicator == 1])
  norm.tag    <- 1.0/sum.correl.tag
  norm.no.tag <- 1.0/Nm
  RES <- cumsum(tag.indicator * correl.vector * norm.tag - no.tag.indicator * norm.no.tag)      
  max.ES <- max(RES)
  min.ES <- min(RES)
  if (max.ES > - min.ES) {
    #      ES <- max.ES
    ES <- signif(max.ES, digits = 5)
    arg.ES <- which.max(RES)
  } else {
    #      ES <- min.ES
    ES <- signif(min.ES, digits=5)
    arg.ES <- which.min(RES)
  }
  return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))  
}

GSEA.EnrichmentScore2 <- function(gene.list, gene.set, weighted.score.type = 1, correl.vector = NULL) {  
  #
  # Computes the weighted GSEA score of gene.set in gene.list. It is the same calculation as in 
  # GSEA.EnrichmentScore but faster (x8) without producing the RES, arg.RES and tag.indicator outputs.
  # This call is intended to be used to asses the enrichment of random permutations rather than the 
  # observed one.
  # The weighted score type is the exponent of the correlation 
  # weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted). When the score type is 1 or 2 it is 
  # necessary to input the correlation vector with the values in the same order as in the gene list.
  #
  # Inputs:
  #   gene.list: The ordered gene list (e.g. integers indicating the original position in the input dataset)  
  #   gene.set: A gene set (e.g. integers indicating the location of those genes in the input dataset) 
  #   weighted.score.type: Type of score: weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted)  
  #  correl.vector: A vector with the coorelations (e.g. signal to noise scores) corresponding to the genes in the gene list 
  #
  # Outputs:
  #   ES: Enrichment score (real number between -1 and +1) 
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.
  
  N <- length(gene.list) 
  Nh <- length(gene.set) 
  Nm <-  N - Nh 
  
  loc.vector <- vector(length=N, mode="numeric")
  peak.res.vector <- vector(length=Nh, mode="numeric")
  valley.res.vector <- vector(length=Nh, mode="numeric")
  tag.correl.vector <- vector(length=Nh, mode="numeric")
  tag.diff.vector <- vector(length=Nh, mode="numeric")
  tag.loc.vector <- vector(length=Nh, mode="numeric")
  
  loc.vector[gene.list] <- seq(1, N)
  tag.loc.vector <- loc.vector[gene.set]
  
  tag.loc.vector <- sort(tag.loc.vector, decreasing = F)
  
  if (weighted.score.type == 0) {
    tag.correl.vector <- rep(1, Nh)
  } else if (weighted.score.type == 1) {
    tag.correl.vector <- correl.vector[tag.loc.vector]
    tag.correl.vector <- abs(tag.correl.vector)
  } else if (weighted.score.type == 2) {
    tag.correl.vector <- correl.vector[tag.loc.vector]*correl.vector[tag.loc.vector]
    tag.correl.vector <- abs(tag.correl.vector)
  } else {
    tag.correl.vector <- correl.vector[tag.loc.vector]**weighted.score.type
    tag.correl.vector <- abs(tag.correl.vector)
  }
  
  norm.tag <- 1.0/sum(tag.correl.vector)
  tag.correl.vector <- tag.correl.vector * norm.tag
  norm.no.tag <- 1.0/Nm
  tag.diff.vector[1] <- (tag.loc.vector[1] - 1) 
  tag.diff.vector[2:Nh] <- tag.loc.vector[2:Nh] - tag.loc.vector[1:(Nh - 1)] - 1
  tag.diff.vector <- tag.diff.vector * norm.no.tag
  peak.res.vector <- cumsum(tag.correl.vector - tag.diff.vector)
  valley.res.vector <- peak.res.vector - tag.correl.vector
  max.ES <- max(peak.res.vector)
  min.ES <- min(valley.res.vector)
  ES <- signif(ifelse(max.ES > - min.ES, max.ES, min.ES), digits=5)
  
  return(list(ES = ES))
}



