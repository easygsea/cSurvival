library("depmap")
library("ExperimentHub")

## create ExperimentHub query object
eh <- ExperimentHub()
query(eh, "depmap")

## retrieve metadata about cancer cell lines
metadata <- depmap::depmap_metadata()

print(metadata)


#START OF CRISPR
crispr <- eh[["EH2261"]]
crispr
result <- data.frame(crispr$depmap_id,crispr$dependency)
#STUDY df_gene,csv
#example <- read.csv('/Applications/Codes/cSurvival/project_data/TCGA-LUAD/TCGA-LUSC/df_gene.csv')