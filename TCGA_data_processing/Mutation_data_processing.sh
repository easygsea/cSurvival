##################################################
##	       TCGA RNA-seq data processing 	    	##
##  original dataset: Hoadley et al. Cell 2018	##
##	     https://gdc.cancer.gov/node/977		    ##
##################################################

## STEP0 handles poor quality samples and is universal in processing all types of TCGA omics data
## STEP0-i. FLAG BAD QUALITY SAMPLES
## orginal sample quality annotations: https://api.gdc.cancer.gov/data/1a7d7be8-675d-4e60-a105-19d4121bdebf
## downloaded orginal sample quality annotations: https://tau.cmmt.ubc.ca/cSurvival/project_data/977/merged_sample_quality_annotations.tsv
## samples annotated as Do_not_use (n=2614): https://tau.cmmt.ubc.ca/cSurvival/project_data/977/all_bad.tsv
## samples NOT annotated as Do_not_use (n=76707): https://tau.cmmt.ubc.ca/cSurvival/project_data/977/all_good.tsv
awk -F'\t' '$NF=="True"' merged_sample_quality_annotations.tsv > all_bad.tsv
awk -F'\t' '$NF!="True"' merged_sample_quality_annotations.tsv > all_good.tsv

## STEP0-ii. CASES BY STUDY
## available TCGA studies: https://tau.cmmt.ubc.ca/cSurvival/project_data/977/projects
while read -r line; do awk -F'\t' -v x="$line" '$3==x' all_good.tsv | cut -f1|sort -u > $line.testt; done < projects

###################################################################################################################################
## Process solid primary tumor data for non-SKCM & non-LAML, solid primary & metastatic for SKCM, blood-derived primary for LAML ##
###################################################################################################################################

## STEP1. DOWNLOAD ORIGINAL MUTATION CALLS, FILTER LOW QUALITY SAMPLES, EXTRACT PRIMARY TUMOR DATA (01; 01 06 for SKCM; 03 09 for LAML)
## original file: https://api.gdc.cancer.gov/data/1c8cfe5f-e52d-41ba-94da-f15ea1337efc (accessed May 20, 2021)
## downloaded original file on cSurvival server: https://tau.cmmt.ubc.ca/cSurvival/project_data/977/mc3.v0.2.8.PUBLIC.maf.gz
## filtered mutation calls 01 (solid primary tumors): https://tau.cmmt.ubc.ca/cSurvival/project_data/977/mc3.v0.2.8.PUBLIC.good_01.maf
## filtered mutation calls 06 (metastatic tumors): https://tau.cmmt.ubc.ca/cSurvival/project_data/977/mc3.v0.2.8.PUBLIC.good_03.maf
## filtered mutation calls 03 and 09 (blood-derived primary tumors): https://tau.cmmt.ubc.ca/cSurvival/project_data/977/mc3.v0.2.8.PUBLIC.good_0609.maf
## for more info on TCGA Sample Type Code: https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes
awk -F '\t' 'FNR==NR {key[$2]=1}; FNR!=NR {if(!(key[$16])){print}}' all_bad.tsv <(zcat mc3.v0.2.8.PUBLIC.maf.gz) | perl -ne 'if(/^Hugo_Symbol/){print $_;} $line=$_; @ay=split(/\t/,$line); $code=(split(/-/,$ay[15]))[3]; if($code =~ /^01/){print $line;}' > mc3.v0.2.8.PUBLIC.good_01.maf
awk -F '\t' 'FNR==NR {key[$2]=1}; FNR!=NR {if(!(key[$16])){print}}' all_bad.tsv <(zcat mc3.v0.2.8.PUBLIC.maf.gz) | perl -ne 'if(/^Hugo_Symbol/){print $_;} $line=$_; @ay=split(/\t/,$line); $code=(split(/-/,$ay[15]))[3]; if($code =~ /^06/){print $line;}' > mc3.v0.2.8.PUBLIC.good_06.maf
awk -F '\t' 'FNR==NR {key[$2]=1}; FNR!=NR {if(!(key[$16])){print}}' all_bad.tsv <(zcat mc3.v0.2.8.PUBLIC.maf.gz) | perl -ne 'if(/^Hugo_Symbol/){print $_;} $line=$_; @ay=split(/\t/,$line); $code=(split(/-/,$ay[15]))[3]; if($code =~ /^03/ || $code =~ /^09/){print $line;}' > mc3.v0.2.8.PUBLIC.good_0309.maf

## STEP2. ASSIGN MUTATION CALLS TO EACH STUDY, CONVERT INTO CSV FORMAT
## RUN the following R script
# process solid primary tumor data
Rscript Mutation_data_processing.R
# special process on SKCM, >70% are metastases
Rscript Mutation_data_processing.R SKCM
# special process on LAML, blood derived primary
Rscript Mutation_data_processing.R LAML
