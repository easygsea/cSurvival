##################################################
##	         TCGA data quality control       	  ##
##  original dataset: Hoadley et al. Cell 2018	##
##	     https://gdc.cancer.gov/node/977		    ##
##################################################

## FLAG CASES THAT POTENTIALLY JEOPADIZE ANALYSIS ACCURACY
## orginal sample quality annotations: https://api.gdc.cancer.gov/data/1a7d7be8-675d-4e60-a105-19d4121bdebf
## downloaded orginal sample quality annotations: https://tau.cmmt.ubc.ca/cSurvival/project_data/977/merged_sample_quality_annotations.tsv
# Step 1. Use keyword search to extract problematic samples
grep -v 'MDA_RPPA_Core' merged_sample_quality_annotations.tsv | awk -F'\t' '$5 ~ /unacceptable prior treatment/ || $5 ~ /Notification:Prior malignancy/ ||$5 ~ /does not meet/ || $5 ~ /^Redaction/ || $5 ~ /not normal breast/ || $5 ~ /recurrence/ || $5 ~ /origin incorrect/ || ($5 ~ /noncanonical/ && $NF == "True") || ($6 == "FPPE" && $7 ~ /noncanonical/ && $NF == "True") || $5 ~ /normal sample swap/'|grep -v 'RIN value below 7' |awk -F'-' '$6 ~ /^01/'> flagged.tsv

# Step 2. Clean the table by keeping essential columns only and combine samples by cases
sed 's/\t\tFPPE\t/\t/' flagged.tsv | cut -f1,3,5 |sort -u -t$'\t' -k1,1 > flagged_cases.tsv
