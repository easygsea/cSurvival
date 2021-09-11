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

## STEP1. DOWNLOAD ORIGINAL RRPA PROTEIN EXPRESSION DATA, FILTER LOW QUALITY SAMPLES, AND CHECK DUPLICATES
## original file: https://api.gdc.cancer.gov/data/fcbb373e-28d4-4818-92f3-601ede3da5e1 (accessed May 20, 2021)
## downloaded original file on cSurvival server: https://tau.cmmt.ubc.ca/cSurvival/project_data/977/TCGA-RPPA-pancan-clean.txt
# filter low quality samples
awk -F '\t' 'FNR==NR {key[$2]=1}; FNR!=NR {if(!(key[$1])){print}}' all_bad.tsv TCGA-RPPA-pancan-clean.txt > TCGA-RPPA-pancan-clean.good.txt
# check duplicates in NON-SKCM samples -- no duplicate found
awk -F'\t' '$2 != "SKCM"' TCGA-RPPA-pancan-clean.good.txt | awk -F'-' '$4 ~ /^01|03|09/' |cut -d'-' -f1-3|sort|uniq -c|sort -k1,1 -nr|less
# check duplicates in SKCM samples -- TCGA-ER-A2NF found having both primary 01 and metastatic 06 data
awk -F'\t' '$2 == "SKCM"' TCGA-RPPA-pancan-clean.good.txt | awk -F'-' '$4 ~ /^01|06/' |cut -d'-' -f1-3|sort|uniq -c|sort -k1,1 -nr|less

## STEP2. ASSIGN DATA BY STUDIES
while read -r line; do awk -F '\t' 'FNR==NR {key["patient_id"]=1; key[$1]=1}; FNR!=NR {if(key[$1]){print}}' $line.testt <(awk -F'-' '$4 ~ /^01/ || $1 ~ /^patient_id/' TCGA-RPPA-pancan-clean.good.txt|perl -ne ' ($p,$type,$exp)=split(/\t/,$_,3); if($p =~ /^patient_id/){print $p."\t".$exp;next;} @ps=split(/-/,$p);print($ps[0]."-".$ps[1]."-".$ps[2]."\t".$exp);') |sed 's/\t/,/g' > ../TCGA-$line/df_rrpa.csv; done < projects
awk -F '\t' 'FNR==NR {key["patient_id"]=1; key[$1]=1}; FNR!=NR {if(key[$1]){print}}' SKCM.testt <(awk -F'-' '$4 ~ /^01|06/ || $1 ~ /^patient_id/' TCGA-RPPA-pancan-clean.good.txt|perl -ne 'if(/^TCGA-ER-A2NF-06A-21-A241-20/){next;} ($p,$type,$exp)=split(/\t/,$_,3); if($p =~ /^patient_id/){print $p."\t".$exp;next;} @ps=split(/-/,$p);print($ps[0]."-".$ps[1]."-".$ps[2]."\t".$exp);') |sed 's/\t/,/g' > ../TCGA-SKCM/df_rrpa.csv
# no data found for LAML
awk -F '\t' 'FNR==NR {key["patient_id"]=1; key[$1]=1}; FNR!=NR {if(key[$1]){print}}' LAML.testt <(awk -F'-' '$4 ~ /^03|09/ || $1 ~ /^patient_id/' TCGA-RPPA-pancan-clean.good.txt|perl -ne '($p,$type,$exp)=split(/\t/,$_,3); if($p =~ /^patient_id/){print $p."\t".$exp;next;} @ps=split(/-/,$p);print($ps[0]."-".$ps[1]."-".$ps[2]."\t".$exp);') |sed 's/\t/,/g' > ../TCGA-LAML/df_rrpa.csv

## STEP3. FILTER BARELY EXPRESSED PROTEINS
nohup sh -c 'for i in ../TCGA-*;do Rscript scale_df.R $i/df_rrpa;done' &
