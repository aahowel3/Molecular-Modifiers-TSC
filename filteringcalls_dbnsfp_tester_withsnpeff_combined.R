#args = commandArgs(trailingOnly=TRUE)

setwd("C:/Users/Owner/OneDrive - Arizona State University/Documents/TSCdata/")

library(VariantAnnotation)
library(biomaRt)
library(dplyr)
library(tidyr)
library(tidyverse)

#combine these two
#vcf <- readVcf("N013_annotated.snpeff.snpsift.dbnsfp4.1_withnarg_1000.vcf","hg19")

vcf <- readVcf("N007.multibreak.tester.vcf","hg19")

######TRYING MULTISPLIT VCFS NOW
df = as.data.frame(info(vcf))
df$rn <- rownames(df)
df2=as.data.frame(rowRanges(vcf),row.names = c(1:nrow(vcf)))
df2$rn <- rownames(df2)
df3=as.data.frame((geno(vcf)$GT))
df3$rn <- rownames(df3)
df4=as.data.frame((geno(vcf)$AD))
df4$rn <- rownames(df4)


df_list <- list(df, df2, df3,df4)

#merge all data frames in list
total = df_list %>% purrr::reduce(full_join, by='rn')

#total = total %>% remove_rownames %>% column_to_rownames(var="Row.names")
total=within(total, GENEINFO<-data.frame(do.call('rbind', strsplit(as.character(GENEINFO), ":", fixed=TRUE))))
mtorgenes=read.csv("mTor_gene_list.csv")
total$mtor=ifelse(total$GENEINFO$X1 %in% mtorgenes$ï..Gene..153.genes., 1, 0)

#total=dplyr::filter(total, total$mtor=="1")

##When the qualitative prediction tools (except for MutationAssessor and Aloft) had D the variant gained one score
total$score_qual = rowSums(sapply(total[c("dbNSFP_SIFT_pred",
                                          "dbNSFP_SIFT4G_pred",
                                          "dbNSFP_Polyphen2_HDIV_pred",
                                          "dbNSFP_Polyphen2_HVAR_pred",
                                          "dbNSFP_LRT_pred",
                                          "dbNSFP_MutationTaster_pred",
                                          "dbNSFP_FATHMM_pred",
                                          "dbNSFP_PROVEAN_pred",
                                          "dbNSFP_MetaSVM_pred",
                                          "dbNSFP_MetaLR_pred",
                                          "dbNSFP_M_CAP_score",
                                          "dbNSFP_PrimateAI_pred",
                                          "dbNSFP_DEOGEN2_pred",
                                          "dbNSFP_BayesDel_addAF_pred",
                                          "dbNSFP_BayesDel_noAF_pred",
                                          "dbNSFP_ClinPred_pred",
                                          "dbNSFP_LIST_S2_pred",
                                          "dbNSFP_fathmm_MKL_coding_pred",
                                          "dbNSFP_fathmm_XF_coding_pred"
                                          
                                          )], grepl, pattern = "D")) 
total=total %>% mutate_at(vars("score_qual"), ~replace_na(.,0))

#An indicator for a deleterious variant of MutationAssessor was "H" and of Aloft was "R" or "D" with high confidence.
#MutationAssessor
#Aloft
#can only use rowsums if greater than one field
#Some ALoFT scores/information are missing in dbNSFP
total$score_mut = ifelse(total$dbNSFP_MutationAssessor_pred == "H", 1, 0)
total=total %>% mutate_at(vars("score_mut"), ~replace_na(.,0))

total$score_mut2 = ifelse(total$dbNSFP_Aloft_pred == "R" | total$dbNSFP_Aloft_pred == "D", 1, 0)
total=total %>% mutate_at(vars("score_mut2"), ~replace_na(.,0))


insert_zero <- function(smpl_vec){
  smpl_vec[is.na(smpl_vec)] <- 0
  baseline=quantile(as.numeric(unlist(smpl_vec)), probs = c(.90),na.rm=TRUE)
  baseline_num=unname(baseline)
  totals=ifelse(grepl("NA",smpl_vec),0,smpl_vec)
  score_quant = ifelse(as.numeric(unlist(totals)) > baseline_num, 1, 0)
  return(score_quant)
}


#fixes your multi zero with two NANA entries but NOT your multi alleleic probelm
insert_zero <- function(smpl_vec){
  totals=ifelse(grepl("NA",smpl_vec),0,smpl_vec)
  baseline=quantile(as.numeric(unlist(smpl_vec)), probs = c(.90),na.rm=TRUE)
  baseline_num=unname(baseline)
  score_quant = ifelse(as.numeric(unlist(totals)) > baseline_num, 1, 0)
  return(score_quant)
}


insert_zero <- function(smpl_vec){
  totals=ifelse(grepl("NA",total$dbNSFP_VEST4_score),0,total$dbNSFP_VEST4_score)
  baseline=quantile(as.numeric(unlist(total$dbNSFP_VEST4_score)), probs = c(.90),na.rm=TRUE)
  baseline_num=unname(baseline)
  score_quant = ifelse(as.numeric(unlist(totals)) > baseline_num, 1, 0)
  return(score_quant)
}



total$score_quant1=apply(total[c('dbNSFP_VEST4_score')], 2, insert_zero)
total$score_quant2=apply(total[c('dbNSFP_REVEL_score')], 2, insert_zero)
total$score_quant3=apply(total[c('dbNSFP_MutPred_score')], 2, insert_zero)
total$score_quant4=apply(total[c('dbNSFP_MVP_score')], 2, insert_zero)
total$score_quant5=apply(total[c('dbNSFP_MPC_score')], 2, insert_zero)
total$score_quant6=apply(total[c('dbNSFP_DANN_score')], 2, insert_zero)
total$score_quant7=apply(total[c('dbNSFP_CADD_phred')], 2, insert_zero)
total$score_quant8=apply(total[c('dbNSFP_Eigen_phred_coding')], 2, insert_zero)
total$score_quant9=apply(total[c('dbNSFP_Eigen_PC_phred_coding')], 2, insert_zero)


###rare alleles table 
unlist_proper <- function(smpl_vec){
  smpl_vec[is.na(smpl_vec)] <- 0
  totals=ifelse(grepl("NA",smpl_vec),0,smpl_vec)
  finals=as.numeric(unlist(totals))
  return(finals)
}

total$gnomAD_exomes=apply(total[c('dbNSFP_gnomAD_exomes_AF')], 2, unlist_proper)

###adding in ref and alts 
r=lapply(total$ALT, function(x) as.list(as.character(x)))
total$combine = as.data.frame(as.matrix(r))


#total=total %>% mutate_at(vars("score_quant"), ~replace_na(.,0))
total$genename = total$GENEINFO$X1 
total$score_final = total$score_qual + total$score_mut + total$score_mut2 + total$score_quant1 + total$score_quant2 +
  total$score_quant3 + total$score_quant4 +total$score_quant5 + total$score_quant6 + total$score_quant7 + total$score_quant8 + total$score_quant9 
#final_table1 = total[c("genename", "seqnames", "start", "Row.names", "ND_N013P01", "ND_N013P02","score_final")]


final_table1 = total[c("genename", "seqnames", "start", "rn", "REF.y", "combine", "gnomAD_exomes",
                       "score_final", "dbNSFP_HGVSc_snpEff","dbNSFP_HGVSp_snpEff")]

final_table2 <- total %>% dplyr:: select(starts_with("ND"))

everything <-cbind(final_table2, final_table1)

ordered_df = everything[with(everything, order(-score_final, gnomAD_exomes)), ]
ordered_df <- apply(ordered_df,2,as.character)
write.csv(ordered_df, file=args[2])


