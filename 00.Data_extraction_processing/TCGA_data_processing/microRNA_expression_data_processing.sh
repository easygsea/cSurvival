##################################################
##	  TCGA microRNA expression data processing  ##
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

## STEP1. DOWNLOAD ORIGINAL microRNA EXPRESSION DATA, TRANSPOSE, AND FILTER LOW QUALITY SAMPLES
## original file: https://api.gdc.cancer.gov/data/1c6174d9-8ffb-466e-b5ee-07b204c15cf8 (accessed May 20, 2021)
## downloaded original file on cSurvival server: https://tau.cmmt.ubc.ca/cSurvival/project_data/977/pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.csv
# change the delimiter to comma (,) in tst.awk, then transpose table
nohup awk -f tst.awk pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.csv > pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.t.csv&
# extract sample IDs
cut -d',' -f1 pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.t.csv > df_gene_mir.all
# filter low quality samples
awk -F '\t' 'FNR==NR {key[$2]=1}; FNR!=NR {if(!(key[$1])){print}}' all_bad.tsv df_gene_mir.all > df_gene_mir.all_good

## STEP2. EXTRACT PRIMARY TUMOR DATA (01; 01 06 for SKCM; 03 09 for LAML)
## filtered non-SKCM & non-LAML samples (01 solid primary tumors): https://tau.cmmt.ubc.ca/cSurvival/project_data/977/df_gene_mir.all_good.non_skcmlaml_01
## filtered SKCM (01 primary and 06 metastatic tumors): https://tau.cmmt.ubc.ca/cSurvival/project_data/977/df_gene_mir.all_good.skcm_0106
## filtered LAML (03 09 blood-derived primary tumors): https://tau.cmmt.ubc.ca/cSurvival/project_data/977/df_gene_mir.all_good.laml_0309
## for more info on TCGA Sample Type Code: https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes
# non-SKCM & non-LAML samples (01 solid primary tumors)
awk -F '\t' 'FNR==NR {key[$2]=1}; FNR!=NR {if(key[$1]){print $1}}' <(awk -F'\t' '$3!="SKCM" && $3!="LAML"' merged_sample_quality_annotations.tsv) df_gene_mir.all_good > df_gene_mir.all_good.non_skcmlaml
awk -F'-' '$4 ~ /^01/' df_gene_mir.all_good.non_skcmlaml > df_gene_mir.all_good.non_skcmlaml_01
# SKCM samples (01 primary and 06 metastatic tumors)
awk -F '\t' 'FNR==NR {key[$2]=1}; FNR!=NR {if(key[$1]){print $1}}' <(awk -F'\t' '$3=="SKCM"' merged_sample_quality_annotations.tsv) df_gene_mir.all_good > df_gene_mir.all_good.skcm
awk -F'-' '$4 ~ /^01|06/' df_gene_mir.all_good.skcm > df_gene_mir.all_good.skcm_0106
# LAML samples (03 09 blood-derived primary tumors)
awk -F '\t' 'FNR==NR {key[$2]=1}; FNR!=NR {if(key[$1]){print $1}}' <(awk -F'\t' '$3=="LAML"' merged_sample_quality_annotations.tsv) df_gene_mir.all_good > df_gene_mir.all_good.laml
awk -F'-' '$4 ~ /^03|09/' df_gene_mir.all_good.laml > df_gene_mir.all_good.laml_0309

## STEP3. CHECK DUPLICATE SAMPLES
cut -d'-' -f1-3 df_gene_mir.all_good.non_skcmlaml_01 |sort|uniq -c|sort -k1 -nr|awk '$1>1'|sed 's/ *//' > df_gene_mir.all_good.non_skcmlaml_01dup
cut -d'-' -f1-3 df_gene_mir.all_good.skcm_0106 |sort|uniq -c|sort -k1 -nr|awk '$1>1'|sed 's/ *//' > df_gene_mir.all_good.skcm_0106dup
cut -d'-' -f1-3 df_gene_mir.all_good.laml_0309 |sort|uniq -c|sort -k1 -nr|awk '$1>1'|sed 's/ *//' > df_gene_mir.all_good.laml_0309dup

## STEP4. ASSIGN DATA BY STUDIES
while read -r line; do awk -F ',' 'FNR==NR {key["patient_id"]=1; key[$1]=1}; FNR!=NR {if(key[$1]){print}}' $line.testt <(awk -F'-' '$4 ~ /^01|06/ || $1 ~ /^patient_id/' pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.t.csv|perl -ne 'if(/^patient_id/){print $_;next;}if(/^TCGA-23-1023-01R-01R-1564-13/ || /^TCGA-ER-A2NF-06A-11R-A18W-13/){next;}($p,$exp)=split(/,/,$_,2);@ps=split(/-/,$p);print($ps[0]."-".$ps[1]."-".$ps[2].",".$exp);') > ../TCGA-$line/df_mir.csv; done < projects
awk -F ',' 'FNR==NR {key["patient_id"]=1; key[$1]=1}; FNR!=NR {if(key[$1]){print}}' LAML.testt <(awk -F'-' '$4 ~ /^03|09/ || $1 ~ /^patient_id/' pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.t.csv|perl -ne 'if(/^patient_id/){print $_;next;} ($p,$exp)=split(/,/,$_,2);@ps=split(/-/,$p);print($ps[0]."-".$ps[1]."-".$ps[2].",".$exp);') > ../TCGA-LAML/df_mir.csv

## STEP5. FILTER BARELY EXPRESSED microRNAs
nohup sh -c 'for i in ../TCGA-*;do Rscript scale_df.R $i/df_mir;done' &
