# ------- TCGA data types, e.g expression, snv ---------
data_types <- reactive({
  if(rv$tcga | (!rv$tcga & !rv$target & !rv$depmap)){
    c("Expression"="rna", 
      "Mutation"="snv",
      "CNV"="cnv",
      "miRNA"="mir",
      "Methylation"="met"
      )
  }else if(rv$target){
    
    # get current target projects
    selected_target_projects <- rv$project
    
    # the list that contains all the parameters of target projects
    parameters_target_projects <- (TARGET_existing_data %>%
                                     filter(project_name %in% selected_target_projects))$existing_data
    rv$parameters_target_projects <- parameters_target_projects
    
    # get the overlapped parameter of all the target projects
    overlapped_parameter <- Reduce(intersect, parameters_target_projects)
    
    # name the list to generate the choices for users to select
    overlapped_parameter <- name_project_choices(overlapped_parameter)
    rv$overlapped_parameter <- overlapped_parameter
    overlapped_parameter
  }else if(rv$depmap){
    c("Expression"="rna", 
      "Mutation"="snv",
      "CNV"="cnv",
      # "miRNA"="mir",
      # "Methylation"="met",
      "Proteomics"="pro")
  }
})

# ------- placeholder for gene search field ----------
g_placeholder <- reactive({
  if(rv$depmap){
    paste0("On top, select a gene, cancer (sub)type(s), cell lines, and a ",agene()," to load data ...")
  }else{
    "On top left, select a project(s) and click the confirmation button to load data ..."
  }
})

# ------ agene, help info for DepMap gene/drug selection -------
agene <- reactive({
  if(rv$project == "DepMap-Drug"){"drug"}else{"gene"} #gsub("^DepMap-","",rv$project),
})
depmap_gene_help <- reactive({
  if(rv$project == "DepMap-Drug"){
    "a drug to study if its effects on cell survivals are dependent on certain genomic backgrounds."
  }else if(rv$project == "DepMap-CRISPR"){
    "a gene to study if its CRISPR-Cas9-mediated knockout effects on cell survivals are dependent on certain genomic backgrounds."
  }else if(rv$project == "DepMap-RNAi"){
    "a gene to study if its RNAi-mediated knockdown effects on cell survivals are dependent on certain genomic backgrounds."
  }
})

# -------- methods for survival analysis ----------
surv_methods_r <- reactive({
  if(rv$depmapr){
    c(surv_methods, "Density plot"="dens")
  }else{
    surv_methods
  }
})