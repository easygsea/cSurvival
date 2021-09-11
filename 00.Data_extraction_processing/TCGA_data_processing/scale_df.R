infile <- commandArgs(trailingOnly=TRUE) #args[1]
if(length(infile) > 1){
  print("OOPS! Enter the prefix of a file: e.g. TCGA-ACC/df_gene")
}else{
  test_path <- paste0(infile,".csv")
  
  # function to remove all na per column
  not_all_na <- function(x) any(!is.na(x))
  # function to remove all na per column
  not_any_na <- function(x) all(!is.na(x))
  
  scaled_df <- function(df_path,sep=",",na="all",na.rm=F){
    odir <- paste0(dirname(df_path),"/tmp/")
    if(!dir.exists(odir)){dir.create(odir)}
    genes <- fread(df_path, sep=sep, nrows = 0)
    genes <- colnames(genes)[-1]
    df <- fread(df_path, sep=sep, select = "patient_id")
    if(na == "all"){
      func_df <- not_all_na
    }else if(na == "any"){
      func_df <- not_any_na
    }
    ofile <- paste0(infile,"_scale.csv")
    unlink(ofile)
    ogenes <- "patient_id"
    for(i in seq_along(genes)){
      gene <- genes[i]
      exps <- fread(df_path, sep=sep, select = c("patient_id",gene))
      exps[[gene]] <- as.numeric(exps[[gene]])
      exps[[gene]] <- (exps[[gene]] - mean(exps[[gene]]))/sd(exps[[gene]],na.rm = na.rm)
      if(func_df(exps[[gene]])){
        ogenes <- c(ogenes,gene)
      }
    }
    fwrite(as.list(paste0(ogenes,collapse=",")),ofile, sep = ",", quote = F)
  }
  scaled_df(test_path) 
}
