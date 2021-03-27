#------ start analysis when users click the btn ------
observeEvent(input$confirm,{
  if(input$project == ""){
    shinyalert("Please select a project to begin your analysis")
  }else{
    error_g <- NULL; error_lib <- NULL; error_gs <- NULL
    for(x in 1:rv$variable_n){
      cat_id <- paste0("cat_",x)
      if(input[[cat_id]] == "g"){
        g_ui_id <- paste0("g_",x)
        if(input[[g_ui_id]] == ""){
          error_g <- c(error_g,x)
        }
      }else if(input[[cat_id]] == "gs"){
        gs_db_id <- paste0("gs_db_",x)
        if(input[[gs_db_id]] == ""){
          error_lib <- c(error_lib,x)
        }else{
          gs_lib_id <- paste0("gs_l_",x)
          if(input[[gs_lib_id]] == ""){
            error_gs <- c(error_gs,x)
          }
        }
      }
    }
    
    if(!is.null(error_g)){
      txt <- paste0("Analysis #",error_g," step ",error_g,".3") %>% paste0(., collapse = ", ")
      txt <- paste0("Please select a gene in ",txt)
      shinyalert(txt)
    }else if(!is.null(error_lib)){
      txt <- paste0("Analysis #",error_lib," step ",error_lib,".3a") %>% paste0(., collapse = ", ")
      txt <- paste0("Please select a database in ",txt)
      shinyalert(txt)
    }else if(!is.null(error_gs)){
      txt <- paste0("Analysis #",error_gs," step ",error_gs,".3b") %>% paste0(., collapse = ", ")
      txt <- paste0("Please select a gene set (GS) in ",txt)
      shinyalert(txt)
    }
    
    req(is.null(error_g) & is.null(error_lib) & is.null(error_gs))
    withProgress(value = 1, message = "Performing analysis. Please wait a minute ...",{
      lapply(1:rv$variable_n, function(x){
        # perform analysis according to input type
        cat_id <- paste0("cat_",x)
        if(input[[cat_id]] == "g"){
          db_id <- paste0("db_",x)
          # extract gene expression/mutation data
          data <- extract_gene_data(x,input[[db_id]])
          
          # perform Surv if expression-like data
          if(input[[db_id]] != "snv"){
            iter_id <- paste0("iter_",x)
            yn <- ifelse(is.null(input[[iter_id]]), T, input[[iter_id]] == "iter")
            if(yn){
              lower_id <- paste0("lower_",x); higher_id <- paste0("upper_",x); step_id <- paste0("step_",x)
              min <- ifelse(is.null(input[[lower_id]]), rv[[lower_id]], input[[lower_id]])
              max <- ifelse(is.null(input[[higher_id]]), rv[[higher_id]], input[[higher_id]])
              step <- ifelse(is.null(input[[step_id]]), rv[[step_id]], input[[step_id]])
              
              results <- get_info_most_significant_rna(data, min, max, step)
              
              # extract most significant df
              df <- results[["df"]]
            }else{
              clow_id <- paste0("clow_",x)
              cutoff <- ifelse(is.null(input[[clow_id]]), 49, input[[clow_id]])
              df <- get_df_by_cutoff(data, cutoff)
            }
            
            # perform survival analysis
            km_fit_id <- paste0("km_fit_",x)
            rv[[km_fit_id]] <- cal_surv_rna(df)
            print(rv[[km_fit_id]])
          }
        }
      })
    })
  }
})