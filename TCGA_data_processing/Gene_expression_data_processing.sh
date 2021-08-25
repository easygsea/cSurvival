##################################################
##	   TCGA Gene expression data processing 	  ##
##  original dataset: Hoadley et al. Cell 2018	##
##	     https://gdc.cancer.gov/node/977		    ##
##################################################

## STEP0 handles poor quality samples and is universal in processing all types of TCGA omics data
## STEP0-i. FLAG BAD QUALITY SAMPLES
## orginal sample quality annotations: https://api.gdc.cancer.gov/data/1a7d7be8-675d-4e60-a105-19d4121bdebf (accessed May 20, 2021)
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

## STEP1. DOWNLOAD ORIGINAL GENE EXPRESSION DATA, FILTER LOW QUALITY SAMPLES
## original file: https://api.gdc.cancer.gov/data/9a4679c3-855d-4055-8be9-3577ce10f66e (accessed May 20, 2021)
## downloaded original file on cSurvival server: https://tau.cmmt.ubc.ca/cSurvival/project_data/977/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv
# extract sample IDs
head -n1 EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv|sed 's/\t/\n/g'|sed 's/"//g' > df_gene.all
# filter low quality samples
awk -F '\t' 'FNR==NR {key[$2]=1}; FNR!=NR {if(!(key[$1])){print}}' all_bad.tsv df_gene.all > df_gene.all_good

## STEP2. TRANSPOSE GENE EXPRESSION TABLE, EXTRACT PATIENT IDS
nohup awk -f tst.awk EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv | sed 's/"//g' | awk -F'\t' 'BEGIN {OFS = FS} {if($1 ~ /^TCGA/){split($1,a,"-"); $1=a[1]"-"a[2]"-"a[3]}; print $0}' > EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.t.tsv&

## STEP3. EXTRACT PRIMARY TUMOR DATA (01; 01 06 for SKCM; 03 09 for LAML)
## filtered non-SKCM & non-LAML samples (01 solid primary tumors): https://tau.cmmt.ubc.ca/cSurvival/project_data/977/df_gene.all_good.non_skcmlaml_01
## filtered SKCM (01 primary and 06 metastatic tumors): https://tau.cmmt.ubc.ca/cSurvival/project_data/977/df_gene.all_good.skcm_0106
## filtered LAML (03 09 blood-derived primary tumors): https://tau.cmmt.ubc.ca/cSurvival/project_data/977/df_gene.all_good.laml_0309
## for more info on TCGA Sample Type Code: https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes
# non-SKCM & non-LAML samples (01 solid primary tumors)
awk -F '\t' 'key[$1]; FNR==NR {key[$2]=1}' <(awk -F'\t' '$3!="SKCM" && $3!="LAML"' merged_sample_quality_annotations.tsv) df_gene.all_good > df_gene.all_good.non_skcmlaml
awk -F'-' '$4 ~ /^01/' df_gene.all_good.non_skcmlaml > df_gene.all_good.non_skcmlaml_01
# SKCM samples (01 primary and 06 metastatic tumors)
awk -F '\t' 'key[$1]; FNR==NR {key[$2]=1}' <(awk -F'\t' '$3=="SKCM"' merged_sample_quality_annotations.tsv) df_gene.all_good > df_gene.all_good.skcm
awk -F'-' '$4 ~ /^01|06/' df_gene.all_good.skcm > df_gene.all_good.skcm_0106
# LAML samples (03 09 blood-derived primary tumors)
awk -F '\t' 'key[$1]; FNR==NR {key[$2]=1}' <(awk -F'\t' '$3=="LAML"' merged_sample_quality_annotations.tsv) df_gene.all_good > df_gene.all_good.laml
awk -F'-' '$4 ~ /^03|09/' df_gene.all_good.laml > df_gene.all_good.laml_0309

## STEP4. CHECK DUPLICATE SAMPLES
cut -d'-' -f1-3 df_gene.all_good.non_skcmlaml_01 |sort|uniq -c|sort -k1 -nr|awk '$1>1'|sed 's/ *//'|sed 's/ /\t/' > df_gene.all_good.non_skcmlaml_01dup
cut -d'-' -f1-3 df_gene.all_good.skcm_0106 |sort|uniq -c|sort -k1 -nr|awk '$1>1'|sed 's/ *//'|sed 's/ /\t/' > df_gene.all_good.skcm_0106dup
cut -d'-' -f1-3 df_gene.all_good.laml_0309 |sort|uniq -c|sort -k1 -nr|awk '$1>1'|sed 's/ *//'|sed 's/ /\t/' > df_gene.all_good.laml_0309dup

## STEP5. CALCULATE GEOMETRIC MEANS FOR DUPLICATE SAMPLES (R script)
## RUN the following R script to find
# process solid primary tumor data
Rscript Gene_expression_data_processing.R
# special process on SKCM, >70% are metastases
Rscript Gene_expression_data_processing.R SKCM
# remove duplicated data
awk -F '\t' 'FNR==NR {key["patient_id"]=1; key[$2]=1}; FNR!=NR {if(!key[$1]){print}}' <(cat df_gene.all_good.non_skcmlaml_01dup df_gene.all_good.skcm_0106dup) EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.t.tsv > EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.t_nodup.tsv
# add back geometric means
cat EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.t_nodup.tsv EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.good.non_skcmlaml_01dup.tsv EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.good.skcm_0106dup.tsv | sed 's/^gene/patient/' > EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.good.t_dedup.tsv

## STEP6. ASSIGN DATA BY STUDIES
while read -r line; do awk -F '\t' 'FNR==NR {key["patient_id"]=1; key[$1]=1}; FNR!=NR {if(key[$1]){print}}' $line.testt EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.good.t_dedup.tsv | sed 's/\t/,/g' > ../TCGA-$line/df_gene.csv; done < projects

## STEP7. FILTER BARELY EXPRESSED GENES
nohup sh -c 'for i in ../TCGA-*;do Rscript scale_df.R $i/df_gene;done' &
