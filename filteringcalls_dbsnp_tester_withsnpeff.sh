for arg in dbnsfp_vcfs_4.3/LNTS*dbnsfp4.3a_withnarg.vcf
do
file=$(basename "$arg" _dbnsfp4.3a_withnarg.vcf)
#echo "$arg"
Rscript filteringcalls_dbsnp_tester_withsnpeff.R "$arg" ${file}_deleteriousscores.txt
done
