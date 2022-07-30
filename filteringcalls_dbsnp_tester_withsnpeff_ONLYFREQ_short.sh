for arg in dbnsfp_vcfs_4.3/*dbnsfp4.3a_withnarg.vcf
do
file=$(basename "$arg" _dbnsfp4.3a_withnarg.vcf)
#echo "$arg"
Rscript filteringcalls_dbsnp_tester_withsnpeff_ONLYFREQ_short.R "$arg" ${file}_gnomad_exome_0.05.txt
done
