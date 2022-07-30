while read r; 
do
	base=$(basename "$r" .vcf)
	newfile="${base}".dbnsfp4.3a_withnarg.vcf
	java -jar /packages/easy-build/software/snpEff/4.3t-Java-1.8/snpEff/SnpSift.jar dbnsfp -v -db dbNSFP4.1a.txt.gz -n $r > $newfile

done < list_of_snpeffvcs_locs.txt
