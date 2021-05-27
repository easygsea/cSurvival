# ------- TCGA data types, e.g expression, snv ---------
data_types <- reactive({
  dtype <- c("Expression"="rna", 
             "Mutation"="snv",
             "CNV"="cnv",
             "miRNA"="mir",
             "Methylation"="met",
             "RPPA"="rrpa"
  )
  if(!rv$tcga & !rv$target & !rv$depmap){
    dtype
  }else if(rv$tcga){
    if(rv$project != ""){
      # get current target projects
      selected_target_projects <- rv$project
      
      # the list that contains all the parameters of target projects
      parameters_target_projects <- (TCGA_missing_data %>%
                                       filter(project_name %in% selected_target_projects))$missing_data
      overlapped_parameter <- Reduce(intersect, parameters_target_projects)
      dtype[!dtype %in% overlapped_parameter]
    }else{
      dtype
    }
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
    c(surv_methods, "Dependency score distributions"="dens")
  }else{
    surv_methods
  }
})

# -------- denpendency score names & descriptions ----------
dependency_names <- reactive({
  if(rv$project == "DepMap-CRISPR"){
    "Gene effect (CERES)"
  }else if(rv$project == "DepMap-RNAi"){
    "Gene effect (DEMETER2)"
  }else if(rv$project == "DepMap-Drug"){
    "Cell viability"
  }
})

dp_ge_q <- reactive({
  txt <- paste0("According to DepMap (https://depmap.org): The CERES dependency score is based on data from a cell depletion assay."
                ," A <b>lower CERES score</b> indicates a <b>higher likelihood that the gene of interest is essential</b> in a give cell line."
                ," A score of 0 indicates a gene is not essential; correspondingly -1 is comparable to the median of all pan-essential genes."
                ," In our survival curves, values on x-axis are computed as <b>10 ^ CERES</b>; correspondingly 1 means a gene being non-essential in a cell line and -1 is comparable to the median of all pan-essential genes."
  )
  if(rv$project == "DepMap-CRISPR"){
    txt
  }else if(rv$project == "DepMap-RNAi"){
    gsub("CERES","DEMETER2",txt)
  }else if(rv$project == "DepMap-Drug"){
    paste0("Cell viability as measured by logfold change (logFC) relative to DMSO."
           ," In our survival curves, values on x-axis are converted back to fold change by computation as <b>2 ^ logFC</b>.")
  }
})