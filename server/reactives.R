# ------- TCGA data types, e.g expression, snv ---------
data_types <- reactive({
  if(rv$tcga){
    c("Expression"="rna", 
      "Mutation"="snv",
      "CNV"="cnv",
      "miRNA"="mir",
      "Methylation"="met"
      )
  }else if(rv$target){
    c("Expression"="rna", 
      "Mutation"="snv",
      "CNV"="cnv",
      "miRNA"="mir"
    )
  }else if(rv$depmap){
    c("Expression"="rna", 
      "Mutation"="snv",
      "CNV"="cnv",
      "miRNA"="mir",
      "Methylation"="met"
      ,"RRPA"="rrpa")
  }
})
  
