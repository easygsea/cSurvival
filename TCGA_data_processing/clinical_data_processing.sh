############################################
## TCGA clinical outcomes data processing ##
## original dataset: Liu et al. Cell 2021 ##
## https://gdc.cancer.gov/about-data/publications/PanCan-Clinical-2018 ##
############################################

## STEP1. DOWNLOAD SOURCE XLSX FILE, EXTRACT ESSENTIAL COLUMNS, TRANSFORM INTO CSV FORMAT
## original data: https://tau.cmmt.ubc.ca/cSurvival/project_data/tcga_survival_o.txt
## transformed data: https://tau.cmmt.ubc.ca/cSurvival/project_data/tcga_survival.csv
cut -f2,3,5,26-33 tcga_survival_o.txt | sed 's/\t/,/g' > tcga_survival.csv

## STEP2. SEPARATE CLINICAL DATA BY STUDY
## list of TCGA studies: https://tau.cmmt.ubc.ca/cSurvival/project_data/list_tcga
while IFS= read -r line; do grep -e $line -e 'patient_id' tcga_survival.csv > TCGA-$line/df_survival_o.csv;done < list_tcga
