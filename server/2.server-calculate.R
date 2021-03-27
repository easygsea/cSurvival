# start analysis when users click the btn
observeEvent(input$confirm,{
  if(input$project == ""){
    shinyalert("Please select a project to begin your analysis")
  }else{
    withProgress(value = 1, message = "Performing analysis. Please wait a minute ...",{
      lapply(1:rv$variable_n, function(x){
        # perform analysis according to input type
        cat_id <- paste0("cat_",x)
        if(input[[cat_id]] == "g"){
          db_id <- paste0("db_",x)
          if(input[[db_id]] != "snv"){
            # perform Surv if expression-like data
            lower_id <- paste0("lower_",x); higher_id <- paste0("upper_",x); step_id <- paste0("step_",x)
            min <- ifelse(is.null(input[[lower_id]]), rv[[lower_id]], input[[lower_id]])
            max <- ifelse(is.null(input[[higher_id]]), rv[[higher_id]], input[[higher_id]])
            step <- ifelse(is.null(input[[step_id]]), rv[[step_id]], input[[step_id]])
            
            # extract gene expression/mutation data
            data <- extract_gene_data(x,input[[db_id]])
            results <- get_info_most_significant_rna(data, min, max, step)
            print(str(results))
          }
        }
      })
    })
  }
})