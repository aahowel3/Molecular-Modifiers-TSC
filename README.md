# Molecular-Modifiers-TSC

# Annotating snpEFF GATK vcfs with dbnsfp 
dbnsfp_annotations.sh takes the files listed in list_of_snpeffvcs_locs.txt and annotates them with 30+ functional databases 

filteringcalls_dbsnp_tester_withsnpeff.sh/R takes the output annotated vcfs of dbnsfp_annotations.sh to calculate how many databases annotated them as deleterious 

allele_specific_expression_delvariants.sh/R takes the RNA vcfs listed in list_of_rnas_full.txt and cross references them with the exome variants with high deleterious scores to pull out their assoicated RNA allele count

filteringcalls_dbsnp_tester_withsnpeff_ONLYFREQ_short.sh takes the output annotated vcfs of dbnsfp_annotations.sh to filter for variants with low gnomAD frequency 

# Annotating DeepVariant Calls 
dbnsfp_annotations.sh (in /labs/C4rCD/Naranyan/Sampath/LNTS-DOD/Analysis/dbnsfp_annotations/) takes the files listed in deepvariant_gvcf_list.txt (manually added 35P01 from temp directory to list) and annotates them with 30+ functional databases 
