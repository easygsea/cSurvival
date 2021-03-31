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
      rv$variable_nr <- rv$variable_n
      df_list <- list()
      rv[["title_all"]] = ""
      rv[["cutoff_all"]] = ""
      for(x in 1:rv$variable_n){
        # perform analysis according to input type
        cat_id <- paste0("cat_",x)
        db_id <- paste0("db_",x)
        gs_mode_id <- paste0("gs_mode_",x)
        g_ui_id <- paste0("g_",x);gs_lib_id <- gs_lib_id <- paste0("gs_l_",x)
        
        if((input[[cat_id]] == "g" & !is.null(input[[db_id]])) | (input[[cat_id]] == "gs" & !is.null(input[[gs_mode_id]]))){
          extract_mode <- ifelse(input[[cat_id]] == "g", input[[db_id]], input[[gs_mode_id]])
          rv[[paste0("title_",x)]] <- ifelse(input[[cat_id]] == "g", input[[g_ui_id]], input[[gs_lib_id]])
          if(rv[["title_all"]] == ""){
            rv[["title_all"]] <- rv[[paste0("title_",x)]]
          }else{
            rv[["title_all"]] <- paste0(c(rv[["title_all"]],rv[[paste0("title_",x)]]),collapse = " vs ")
          }
          # extract gene expression/mutation data
          data <- extract_gene_data(x,extract_mode)

          # perform Surv if expression-like data
          if(input[[db_id]] != "snv" | input[[cat_id]] == "gs"){
            iter_id <- paste0("iter_",x)
            yn <- ifelse(is.null(input[[iter_id]]), T, input[[iter_id]] == "iter")
            if(yn){
              lower_id <- paste0("lower_",x); higher_id <- paste0("upper_",x); step_id <- paste0("step_",x)
              min <- ifelse(is.null(input[[lower_id]]), rv[[lower_id]], input[[lower_id]])
              max <- ifelse(is.null(input[[higher_id]]), rv[[higher_id]], input[[higher_id]])
              step <- ifelse(is.null(input[[step_id]]), rv[[step_id]], input[[step_id]])
              
              results <- get_info_most_significant_rna(data, min, max, step, mode=input[[cat_id]])
              
              # extract most significant df
              df <- results[["df"]]
              rv[[paste0("cutoff_",x)]] <- results[["cutoff"]]
              rv[["cutoff_all"]] <- paste0(rv[["cutoff_all"]],"; #",x,": ",results[["cutoff"]])
            }else{
              clow_id <- paste0("clow_",x)
              cutoff <- ifelse(is.null(input[[clow_id]]), 49, input[[clow_id]])
              df <- get_df_by_cutoff(data, cutoff)
              rv[[paste0("cutoff_",x)]] <- paste0(input[[clow_id]],"%")
            }
            
            # save df
            df_list[[x]] <- df
            rv[[paste0("df_",x)]] <- df

            # no of cases in each group
            lels <- levels(df$level)
            rv[[paste0("lels_",x)]] <- lapply(lels, function(x) table(df$level == x)["TRUE"] %>% unname(.))
            names(rv[[paste0("lels_",x)]]) <- lels

            # perform survival analysis
            cox_id <- paste0("cox_",x)
            rv[[cox_id]] <- cal_surv_rna(df,1)
          }
        }
      }
      
      # generate interaction KM
      if(rv$variable_n > 1){
        # generate interaction df
        df_combined <- Reduce(
          function(x, y) inner_join(x, select(y, patient_id, level), by = "patient_id"),
          df_list
        )
        x_y <- c("x","y")[1:length(df_list)]
        df_combined[["level"]] <- apply(df_combined %>% select(paste0("level.",x_y)),1,paste0,collapse="_")
        rv[["df_all"]] <- df_combined

        # no of cases in each group
        lels <- unique(df_combined$level) %>% sort(.,decreasing = T)
        df_combined$level <- factor(df_combined$level, levels = lels)
        lels <- levels(df_combined$level)
        rv[["lels_all"]] <- lapply(lels, function(x) table(df_combined$level == x)["TRUE"] %>% unname(.))
        names(rv[["lels_all"]]) <- lels
        
        # perform survival analysis
        rv[["cox_all"]] <- cal_surv_rna(df_combined,rv$variable_n)
        # saveRDS(rv[["cox_all"]], "cox_all")
      }
    })
  }
})
