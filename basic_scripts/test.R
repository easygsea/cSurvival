all_file <- "/Users/jeancheng/Documents/cSurvival/project_data/TCGA-BRCA/df_gene.csv"
all_genes <- fread(all_file,sep=",",header=T,nrows=0) %>% names(.)
all_genes <- sapply(all_genes, function(x){
  x <- strsplit(x,"\\|")[[1]]
  if(length(x) == 1){x}else{x[-1]}
}) %>% unname(.)

genes <- "LINC02082|100507661|ENSG00000242268"
genes <- strsplit(genes,"\\|")[[1]]
if(length(genes) == 1){genes <- genes}else{genes <- genes[-1]}
genes

col_to_drop <- (1:length(all_genes))[all_genes != genes] %>% .[-1]
data_o <- fread(all_file,sep=",",header=T,drop = col_to_drop)
exp_scale <- scale(data_o[,2])[,1]
data_o[,2] <- exp_scale

#--------------------

muts <- c(rep("3'UTR",3),rep("Missense_Mutation",7))
stats <- table(muts)

dat <- data.frame(
  Mutation = names(stats),
  Frequency = as.numeric(stats)
) %>% dplyr::arrange(Frequency)

dat$Mutation <- factor(dat$Mutation, levels = unique(dat$Mutation))

fig <- ggplot(data=dat,aes(x=Mutation, y=Frequency)) + geom_bar(stat="identity")
ggplotly(fig)
