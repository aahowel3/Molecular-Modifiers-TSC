for arg in rna_hcs/ND*.vcf
do
file=$(basename "$arg" .rnaHC_All.snpEff.vcf)
tag=$(basename "$arg")
code=${tag:0:7}
Rscript allele_specific_expression_delvariants.R dbnsfp_scores_4.3/${code}* "$arg" ${file}_allelespecific.txt
done


#for arg in rna_hcs/LNTS*.vcf
#do
#file=$(basename "$arg" .rnaHC_All.snpEff.vcf)
#tag=$(basename "$arg")
#code=${tag:0:10}
#echo "Rscript allele_specific_expression_delvariants.R dbnsfp_scores_4.3/*${code}* "$arg" ${file}_allelespecific.txt"
#done

