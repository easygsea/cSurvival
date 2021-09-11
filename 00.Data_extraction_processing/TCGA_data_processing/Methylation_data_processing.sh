##################################################
##	     TCGA methylation data processing       ##
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

## STEP1. DOWNLOAD ORIGINAL METHYLATION DATA, TRANSPOSE, AND FILTER LOW QUALITY SAMPLES
## we used the "Between-platform normalization for DNA Methylation data - DNA methylation (Merged 27K+450K) Beta Value"
## original file: https://api.gdc.cancer.gov/data/d82e2c44-89eb-43d9-b6d3-712732bf6a53 (accessed May 20, 2021)
## downloaded original file on cSurvival server: https://tau.cmmt.ubc.ca/cSurvival/project_data/977/jhu-usc.edu_PANCAN_HumanMethylation450.betaValue_whitelisted.tsv
# transpose table
nohup awk -f tst.awk jhu-usc.edu_PANCAN_merged_HumanMethylation27_HumanMethylation450.betaValue_whitelisted.tsv > jhu-usc.edu_PANCAN_merged_HumanMethylation27_HumanMethylation450.betaValue_whitelisted.t.tsv &
# extract probe IDs
head -n1 ../jhu-usc.edu_PANCAN_merged_HumanMethylation27_HumanMethylation450.betaValue_whitelisted.t.all_good.tsv |sed 's/\t/\n/g' > met.ids
# filter low quality samples
awk -F '\t' 'FNR==NR {key[$2]=1}; FNR!=NR {if(!(key[$1])){print}}' all_bad.tsv jhu-usc.edu_PANCAN_merged_HumanMethylation27_HumanMethylation450.betaValue_whitelisted.t.tsv > jhu-usc.edu_PANCAN_merged_HumanMethylation27_HumanMethylation450.betaValue_whitelisted.t.all_good.tsv

## STEP2. PROBE ID CONVERSION
## we used annotations by Zhou et al, NAR, 2016: https://zwdzwd.github.io/InfiniumAnnotation (Update Apr-14-2021, accessed May 20, 2021)
## Basic hg38 annotation with suggested overall masking (HM27): http://zhouserver.research.chop.edu/InfiniumAnnotation/20180909/HM27/HM27.hg38.manifest.tsv.gz
## Basic hg38 annotation with suggested overall masking (HM450): http://zhouserver.research.chop.edu/InfiniumAnnotation/20180909/HM450/HM450.hg38.manifest.tsv.gz
## downloaded HM27 file on cSurvival server: https://tau.cmmt.ubc.ca/cSurvival/project_data/977/HM27.hg38.manifest.tsv
## downloaded HM450 file on cSurvival server: https://tau.cmmt.ubc.ca/cSurvival/project_data/977/HM450.hg38.manifest.tsv
# extract essential columns and combine IDs into tab-delimited tables
awk -F'\t' -v OFS="\t" '{ print $5 OFS $1 OFS $2 OFS $3 OFS $21}' HM27.hg38.manifest.tsv |sed 's/\t/|/'|sed 's/\t/:/'|sed 's/\t/-/'|sed 's/\t/|/'|sed 's/|/\t/' > HM27.hg38.manifest.tsv.ids
awk -F'\t' -v OFS="\t" '{ print $5 OFS $1 OFS $2 OFS $3 OFS $21}' HM450.hg38.manifest.tsv |sed 's/\t/|/'|sed 's/\t/:/'|sed 's/\t/-/'|sed 's/\t/|/'|sed 's/|/\t/' > HM450.hg38.manifest.tsv.ids
cat HM27.hg38.manifest.tsv.ids HM450.hg38.manifest.tsv.ids |sort -u -t$'\t' -k1,1 > HM.ids
cat <(echo "patient_id") <(awk -F '\t' 'FNR==NR {key[$1]=1}; FNR!=NR {if(key[$1]){print}}' met.ids HM.ids) > HM.met.ids
# prepare header line with complete ID information
join -t$'\t' met.ids HM.met.ids | sed 's/\t/|/' | sed ':a;N;$!ba;s/\n/\t/g' > HM.met.ids.final
# replace header row by complete information
sed -i -e '1R HM.met.ids.final' -e '1d' jhu-usc.edu_PANCAN_merged_HumanMethylation27_HumanMethylation450.betaValue_whitelisted.t.all_good.tsv

## STEP3. CHECK DUPLICATE SAMPLES
## STEP3a. handle duplicates for non-SKCM samples
# find duplicates
cut -f1 jhu-usc.edu_PANCAN_merged_HumanMethylation27_HumanMethylation450.betaValue_whitelisted.t.all_good.tsv |awk -F'-' '$4 ~ /^01|03|09/'|cut -d'-' -f1-3 |sort|uniq -c|sort -k1 -nr|sed 's/^\s*//g'|sed 's/\s/\t/' > met.dup
# extract data of duplicated samples
mkdir methylation
while read -r line; do grep "^$line" jhu-usc.edu_PANCAN_merged_HumanMethylation27_HumanMethylation450.betaValue_whitelisted.t.all_good.tsv | awk -F'-' '$4 ~ /^01|03|09/' > methylation/$line.tsv; done < <(cut -f2 met.dup)
# calculate geometric means in R
Rscript Methylation_geomean.R
## STEP3b. check if any SKCM tumor has duplicated data -- no duplicate found
awk -F '\t' 'FNR==NR {key[$1]=1}; FNR!=NR {if(key[$2]){print}}' SKCM.testt met.dup

## STEP4. EXTRACT PRIMARY TUMOR DATA (01; 01 06 for SKCM; 03 09 for LAML)
## for more info on TCGA Sample Type Code: https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes
## STEP4a. non-SKCM tumor samples (01 solid and 03 09 blood-derived primary tumors)
# filter primary tumors
awk -F'-' '$1 ~ /^patient_id/ || $4 ~ /^01|03|09/' jhu-usc.edu_PANCAN_merged_HumanMethylation27_HumanMethylation450.betaValue_whitelisted.t.all_good.tsv > jhu-usc.edu_PANCAN_merged_HumanMethylation27_HumanMethylation450.betaValue_whitelisted.t.all_good.010309.tsv
cut -f1 jhu-usc.edu_PANCAN_merged_HumanMethylation27_HumanMethylation450.betaValue_whitelisted.t.all_good.010309.tsv > 1
cut -d "-" -f1-3 1 > 1del
cut -f2- jhu-usc.edu_PANCAN_merged_HumanMethylation27_HumanMethylation450.betaValue_whitelisted.t.all_good.010309.tsv > 2
paste 1del 2 > jhu-usc.edu_PANCAN_merged_HumanMethylation27_HumanMethylation450.betaValue_whitelisted.t.all_good.010309.tsv
rm 1 1del 2
# remove duplicates
awk -F '\t' 'FNR==NR {key[$2]=1}; FNR!=NR {if(!(key[$1])){print}}' met.dup jhu-usc.edu_PANCAN_merged_HumanMethylation27_HumanMethylation450.betaValue_whitelisted.t.all_good.010309.tsv > jhu-usc.edu_PANCAN_merged_HumanMethylation27_HumanMethylation450.betaValue_whitelisted.t.all_good.010309_dedup.tsv
# add back re-calculated geometric means
cat jhu-usc.edu_PANCAN_merged_HumanMethylation27_HumanMethylation450.betaValue_whitelisted.t.all_good.010309_dedup.tsv methylation/*.geomean > jhu-usc.edu_PANCAN_merged_HumanMethylation27_HumanMethylation450.betaValue_whitelisted.t.all_good.010309.tsv

## STEP4b. SKCM samples (01 primary and 06 metastatic tumors)
awk -F'-' '$1 ~ /^patient_id/ || $4 ~ /^01|06/' jhu-usc.edu_PANCAN_merged_HumanMethylation27_HumanMethylation450.betaValue_whitelisted.t.all_good.tsv > jhu-usc.edu_PANCAN_merged_HumanMethylation27_HumanMethylation450.betaValue_whitelisted.t.all_good.0106.tsv

## STEP5. ASSIGN DATA BY STUDIES
while read -r line; do awk -F '\t' 'FNR==NR {key["patient_id"]=1; key[$1]=1}; FNR!=NR {if(key[$1]){print}}' $line.testt jhu-usc.edu_PANCAN_merged_HumanMethylation27_HumanMethylation450.betaValue_whitelisted.t.all_good.010309.tsv |sed 's/\t/,/g' > ../TCGA-$line/df_met.csv; done < projects
awk -F '\t' 'FNR==NR {key["patient_id"]=1; key[$1]=1}; FNR!=NR {if(key[$1]){print}}' SKCM.testt jhu-usc.edu_PANCAN_merged_HumanMethylation27_HumanMethylation450.betaValue_whitelisted.t.all_good.0106.tsv |sed 's/\t/,/g' > ../TCGA-SKCM/df_met.csv

## STEP5. FILTER BARELY DETECTED METHYLATION PROBES
nohup sh -c 'for i in ../TCGA-*;do Rscript scale_df.R $i/df_met;done' &
