#https://bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")
#BiocManager::install("sva")
#https://www.biostars.org/p/452858/

library(EnrichmentBrowser)
library("data.table")
library(DESeq2)
library("tidyverse")
library(statmod)
library(ggrepel)

mergeded=read.csv("htseqs_july2022/expressiondata_july2022.csv",sep = ",")
metadata=read.csv("htseqs_july2022/metadata_LNTS_preDOD_assumptions_minus8_new36_7both.csv",sep=",")

#try putting it into deseq
#will not run without metadata file 
#will need to account for ~family and ~treatment affect 
dds <- DESeqDataSetFromMatrix(countData=mergeded, 
                              colData=metadata, 
                              design=~condition + family, tidy = TRUE)


#even more stringent filterirng 
# at least 3 samples with a count of 10 or higher (not just the row total)
keep <- rowSums(counts(dds) >= 10) >= 3
#keep <- rowSums(counts(dds)) >= 10

dds <- dds[keep,]
nrow(dds)

#actually running the differnetial expression 
ddx = DESeq(dds)

#deleted this section for the ddx results and entrex outputs bc its similar to with the SVA 
#and methodolgically i think we need it since we see such large batch effects and cant back out of that
###########################
#########################################rerun with batch effects s
#check for batch effects 
#https://biodatascience.github.io/compbio/dist/sva.html 
library("sva")
dat  <- counts(ddx, normalized = TRUE)
idx  <- rowMeans(dat) > 1
dat  <- dat[idx, ]
mod  <- model.matrix(~ condition, colData(ddx))
mod0 <- model.matrix(~   1, colData(ddx))
svseq <- svaseq(dat, mod, mod0, n.sv = 2)
norm.cts <- dat[rowSums(dat) > 0,]
fit <- svaseq(norm.cts, mod=mod, mod0=mod0, n.sv=2)


#rerunning deseq with additional surrgate variables to account for batch effects
ddssva <- dds
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]

design(ddssva) <- ~ SV1 + SV2 + condition + family 
#design(ddssva) <- ~ condition + SV1 



ddssva2 = DESeq(ddssva)


#specify what constrasts you want to look at 
res <- results(ddssva2, contrast=c("condition","severe","mild"),alpha=0.05)
summary(res)
resSig <- subset(res, padj < 0.05)


#ANNOTATING AND EXPORTING RESULTS
library("AnnotationDbi")
library("org.Hs.eg.db")
library("EnsDb.Hsapiens.v75")

ens.str <- substr(rownames(resSig), 1, 15)
resSig$symbol <- mapIds(EnsDb.Hsapiens.v75,
                        keys=ens.str,
                        column="SYMBOL",
                        keytype="GENEID",
                        multiVals="first")
resSig$entrez <- mapIds(EnsDb.Hsapiens.v75,
                        keys=ens.str,
                        column="ENTREZID",
                        keytype="GENEID",
                        multiVals="first")

write.csv(resSig, file = "results_10_3_0.05_DGE_SVA_v75.csv")




#does not work 
#enrichment browser only allows for 2 fixed vaariables at a time - condition + group 
#and by adding the SVA variableds in - +condition + group + SV1 + SV2 
#will not take it 
#but i dont know why enrichment browser needs the design of the study since it only looks at genes
#in DE and not in DE set compared to those in gene set and those DE in gene set 
#se <- import(ddssva2, res)

#have to use this one with new res object but old deseq object 
se <- import(ddx, res)


#remaining steps in the enrichmentbrowser are: 
#looks like idMap works just as well as AnnotationDBI 
try <- idMap(se, org = "hsa", from = "ENSEMBL", to = "ENTREZID")
#choose a gene set - could be kegg, GO, etc 
kegg.gs <- getGenesets(org = "hsa", db = "kegg")
go.gs <- getGenesets(org = "hsa", db = "go")
#perform which set-based enrichment analysis (sbea) you choose - ora = overrep analysis
air.res <- sbea(method = "ora", se = try, gs = go.gs, perm = 0, alpha = 0.2)
ggeaout=gsRanking(air.res)
write.csv(ggeaout, file = "results_10_3_0.05_ora_go_v75.csv")

air.res <- sbea(method = "ora", se = try, gs = kegg.gs, perm = 0, alpha = 0.2)
ggeaout=gsRanking(air.res)
write.csv(ggeaout, file = "results_10_3_0.05_ora_kegg_v75.csv")


gsea.all <- sbea(method="gsea", se= try, gs=kegg.gs, perm=100)  
ggeaout=gsRanking(gsea.all)
write.csv(ggeaout, file = "results_10_3_0.05_gsea_go_v75.csv")


gsea.all <- sbea(method="gsea", se= try, gs=go.gs, perm=100)  
ggeaout=gsRanking(gsea.all)
write.csv(ggeaout, file = "results_10_3_0.05_gsea_kegg_v75.csv")


#network analysis ggea
#frequires add gs set - basucally kegg set tho?
data.dir <- system.file("extdata", package="EnrichmentBrowser")
gmt.file <- file.path(data.dir, "hsa_kegg_gs.gmt")
hsa.gs <- getGenesets(gmt.file)
#create grn file
hsa.grn <- compileGRN(org="hsa", db="kegg",kegg.native=TRUE)
#run ggea
ggea.all <- nbea(method="spia", se=try, gs=hsa.gs, grn=hsa.grn, alpha=0.2)
ggeaout=gsRanking(ggea.all)
write.csv(ggeaout, file = "results_10_3_0.05_spia_v75.csv")
