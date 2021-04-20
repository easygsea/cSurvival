library(edgeR)
library(limma)
library(org.Hs.eg.db)

setwd("/Users/jeancheng/Documents/cSurvival/basic_scripts/eVITTA_test")

# ------ gender effect DE analysis ------

# # retrieve expression data
df_gene_gender <- readRDS("df_gene_gender.rds") 
patients <- df_gene_gender$patient_id
genes <- colnames(df_gene_gender)[-1]
df_gene_gender <- transpose(df_gene_gender) %>% .[-1,] %>% as.data.frame()
# reassign patient ids
colnames(df_gene_gender) <- patients
rownames(df_gene_gender) <- genes

# id conversion
# create individual tables using org.Hs
egENS <- toTable(org.Hs.egENSEMBL)
egSYMBOL <- toTable(org.Hs.egSYMBOL)

# bind the tables
id_table <- egENS %>% left_join(egSYMBOL, by = "gene_id")
head(id_table)

# extract the gene ids from df_gene and find their names from id conversion table
gene_ids <- as_tibble_col(rownames(df_gene_gender), column_name = "ensembl_id")
gene_ids_table <- left_join(gene_ids, id_table, by="ensembl_id")

# convert gene ids
df_gene_gender <- df_gene_gender %>% dplyr::mutate(ensembl_id=rownames(df_gene_gender)) %>%
  left_join(gene_ids_table, by = "ensembl_id")

# remove unnecessary columns
df_gene_gender <- df_gene_gender %>%
  dplyr::distinct(symbol,.keep_all = TRUE) %>%
  dplyr::filter(!is.na(symbol))

genes <- df_gene_gender$symbol
rownames(df_gene_gender) <- genes

# convert to numeric matrix
df_gene_gender <- df_gene_gender %>%
  dplyr::select(-c(ensembl_id, gene_id, symbol)) %>%
  dplyr::mutate_all(as.numeric) %>%
  as.matrix(.)

# resave data
saveRDS(df_gene_gender,"df_gene_gender.rds") 

# # retrieve design data
df_design_gender <- readRDS("df_design_gender.rds")
# create the design model
design <- model.matrix(~level.x*level.y, df_design_gender)

# # run DE analysis
# 1) create dgelist
y <- DGEList(counts=df_gene_gender)

# 2) filter low expressing genes
min_n <- min(table(df_design_gender$level))
keep <- rowSums(y$counts>1) >= min_n
y <- y[keep,,keep.lib.sizes=TRUE]

# 3) voom directly on counts, if data are very noisy, as would be used for microarray
v <- voom(y, design, plot=F, normalize.method="quantile")

# 4) DEG analysis
fit <- lmFit(v, design)
fit <- eBayes(fit,trend=TRUE, robust=TRUE)

# results
results <- decideTests(fit)
summary(results)

# available coefficients
coefs <- colnames(design)[-1]

# name the coefficients
lels1 <- levels(df_design_gender$level.x) %>% rev()
lels2 <- levels(df_design_gender$level.y) %>% rev()
names(coefs) <- c(
  paste0(lels1, collapse = " vs. "),
  paste0(lels2, collapse = " vs. "),
  paste0("(",paste0(lels2, collapse = " vs. "),") vs. (",paste0(lels1, collapse = " vs. "),")")
)

# export DEG table
degss = lapply(coefs, function(x){
  topTable(fit, coef=x,sort.by="P",number=Inf)
})
names(degss) <- coefs

saveRDS(degss,"degss.rds")
saveRDS(coefs,"coefs.rds")
