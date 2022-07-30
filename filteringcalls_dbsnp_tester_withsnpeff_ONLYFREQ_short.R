args = commandArgs(trailingOnly=TRUE)
library(VariantAnnotation)
library(biomaRt)
library(dplyr)
library(tidyr)
library(tidyverse)

#combine these two
vcf <- readVcf(args[1], "hg19")

df = as.data.frame(info(vcf))
df2=as.data.frame(rowRanges(vcf))
df3=as.data.frame((geno(vcf)$GT))

total1 <- merge(df,df2, by="row.names")
total1 = total1 %>% remove_rownames %>% column_to_rownames(var="Row.names")

total=merge(total1,df3,by="row.names")

total=within(total, GENEINFO<-data.frame(do.call('rbind', strsplit(as.character(GENEINFO), ":", fixed=TRUE))))


mtorgenes=read.csv("mTor_gene_list.csv")
total$mtor=ifelse(total$GENEINFO$X1 %in% mtorgenes$gene, 1, 0)
total=dplyr::filter(total, total$mtor=="1")


###rare alleles table 
unlist_proper <- function(smpl_vec){
  smpl_vec[is.na(smpl_vec)] <- 0
  totals=ifelse(grepl("NA",smpl_vec),0,smpl_vec)
  finals=as.numeric(unlist(totals))
  return(finals)
}


total$gnomAD_exomes=apply(total[c('dbNSFP_gnomAD_exomes_AF')], 2, unlist_proper)
total$newfrequency=ifelse(total$gnomAD_exomes < 0.05 & total$gnomAD_exomes > 0, 1, 0)
total=dplyr::filter(total, total$newfrequency=="1")

total$genename = total$GENEINFO$X1 
final_table1 = total[c("genename", "seqnames", "start", "Row.names", "gnomAD_exomes")]
final_table2 <- total %>% dplyr:: select(starts_with("ND"))
everything <-cbind(final_table1, final_table2)

final_table3 <- total %>% dplyr:: select(starts_with("LNTS"))
everything2 <-cbind(everything, final_table3)


ordered_df <- everything2[order(-everything2$gnomAD_exomes),]
ordered_df <- apply(ordered_df,2,as.character)

write.csv(ordered_df, file=args[2])
