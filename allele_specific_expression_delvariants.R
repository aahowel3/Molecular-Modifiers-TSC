args = commandArgs(trailingOnly=TRUE)

library(VariantAnnotation)
library(biomaRt)
library(dplyr)
library(tidyr)
library(tidyverse)

scores=read.csv(args[1])
vcf <- readVcf(args[2], "hg19")


df = as.data.frame(info(vcf))
df2=as.data.frame(rowRanges(vcf))
df3=as.data.frame((geno(vcf)$GT))
#hmm it appears you can get the AD out of them
df4=as.data.frame((geno(vcf)$AD))

total1 <- merge(df,df2, by="row.names")
total1 = total1 %>% remove_rownames %>% column_to_rownames(var="Row.names")

total2=merge(df3, df4, by="row.names")
total2 = total2 %>% remove_rownames %>% column_to_rownames(var="Row.names")

total=merge(total1,total2,by="row.names")

total=within(total, GENEINFO<-data.frame(do.call('rbind', strsplit(as.character(GENEINFO), ":", fixed=TRUE))))



total$selected=ifelse(total$Row.names %in% scores$Row.names, 1, 0)
total=dplyr::filter(total, total$selected=="1")

total$genename = total$GENEINFO$X1 
final_table1 = total[c("genename", "seqnames", "start", "Row.names","REF.y","ALT")]
final_table2 <- total %>% dplyr:: select(starts_with("ND"))

everything <-cbind(final_table1, final_table2)
final_table3 <- total %>% dplyr:: select(starts_with("LNTS"))
everything2 <-cbind(everything, final_table3)

ordered_df=everything[order(match(everything$Row.names, scores$Row.names)), ]
ordered_df <- apply(ordered_df,2,as.character)

write.csv(ordered_df, file=args[3])
 
