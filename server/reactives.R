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
    
    # get current target projects
    selected_target_projects <- rv$project
    # the list that contains all the parameters of target projects
    parameters_target_projects <- (TARGET_existing_data %>%
                                     filter(project_name %in% selected_target_projects))$existing_data
    rv$parameters_target_projects <- parameters_target_projects
    print(rv$parameters_target_projects)
    # get the overlapped parameter of all the target projects
    overlapped_parameter <- Reduce(intersect, parameters_target_projects)
    print(overlapped_parameter)
    # name the list to generate the choices for users to select
    overlapped_parameter <- name_project_choices(overlapped_parameter)
    rv$overlapped_parameter <- overlapped_parameter
    overlapped_parameter
    
  }else if(rv$depmap){
    c("Expression"="rna", 
      "Mutation"="snv",
      "CNV"="cnv",
      "miRNA"="mir",
      "Methylation"="met"
      ,"RRPA"="rrpa")
  }
})
  
