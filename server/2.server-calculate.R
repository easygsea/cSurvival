observeEvent(input$confirm,{
  if(input$project == ""){
    shinyalert("Please select a project to begin your analysis")
  }else{
    #------ 1. check if any errors by user ------
    error_g <- NULL; error_lib <- NULL; error_manual <- NULL; error_gs <- NULL
    for(x in 1:rv$variable_n){
      cat_id <- paste0("cat_",x)
      if(input[[cat_id]] == "g"){
        g_ui_id <- paste0("g_",x)
        if(input[[g_ui_id]] == ""){
          error_g <- c(error_g,x)
        }
      }else if(input[[cat_id]] == "gs"){
        gs_mode_id <- paste0("gs_mode_",x)
        if(input[[gs_mode_id]] == "lib"){
          gs_db_id <- paste0("gs_db_",x)
          if(input[[gs_db_id]] == ""){
            error_lib <- c(error_lib,x)
          }else{
            gs_lib_id <- paste0("gs_l_",x)
            if(input[[gs_lib_id]] == ""){
              error_gs <- c(error_gs,x)
            }
          }
        }else if(input[[gs_mode_id]] == "manual"){
          gs_manual_id <- paste0("gs_m_",x)
          if(rv[[gs_manual_id]] == ""){
            error_manual <- c(error_manual,x)
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
    }else if(!is.null(error_manual)){
      txt <- paste0("Analysis #",error_manual," step ",error_manual,".3") %>% paste0(., collapse = ", ")
      txt <- paste0("Please enter a valid gene list in ",txt)
      txt <- paste0(txt,". Click Submit to confirm your entry")
      shinyalert(txt)
    }else if(!is.null(error_gs)){
      txt <- paste0("Analysis #",error_gs," step ",error_gs,".3b") %>% paste0(., collapse = ", ")
      txt <- paste0("Please select a gene set (GS) in ",txt)
      shinyalert(txt)
    }
    
    req(is.null(error_g) & is.null(error_lib) & is.null(error_manual) & is.null(error_gs))
    
    #------ 2. begin analysis ------
    withProgress(value = 1, message = "Performing analysis. Please wait a minute ...",{
      rv$variable_nr <- rv$variable_n
      if(rv$variable_nr == 1){
        rv$plot_type <- "1"
      }else{
        rv$plot_type <- "all"
      }
      df_list <- list()
      rv[["title_all"]] = ""
      rv[["cutoff_all"]] = ""
      #------ 3. Loop from 1 to rv$variable_n ------
      for(x in 1:rv$variable_n){
        # clear previous RVs
        rv[[paste0("mutations_",x)]] <- NULL
        # perform analysis according to input type
        cat_id <- paste0("cat_",x)
        db_id <- paste0("db_",x)
        gs_mode_id <- paste0("gs_mode_",x)
        g_ui_id <- paste0("g_",x);gs_lib_id <- gs_lib_id <- paste0("gs_l_",x)
        
        #------ 3A. detect type of input data ------
        if((input[[cat_id]] == "g" & !is.null(input[[db_id]])) | (input[[cat_id]] == "gs" & !is.null(input[[gs_mode_id]]))){
          # clear RVs
          rv[[paste0("title_",x)]] <- ""
          # mode of extracting gene expression/mutation data in extract_gene_data
          extract_mode <- ifelse(input[[cat_id]] == "g", input[[db_id]], input[[gs_mode_id]])
          
          # title for the survival plot
          rv[[paste0("title_",x)]] <- ifelse(input[[cat_id]] == "g", input[[g_ui_id]], input[[gs_lib_id]])
          if(rv[["title_all"]] == ""){
            rv[["title_all"]] <- rv[[paste0("title_",x)]]
          }else{
            rv[["title_all"]] <- paste0(c(rv[["title_all"]],rv[[paste0("title_",x)]]),collapse = " vs ")
          }
          # extract gene expression/mutation data
          data <- extract_gene_data(x,extract_mode)

          # check if manually input genes exist in database
          if(extract_mode == "manual"){
            n_col <- ncol(data) - 1
            if(n_col < 1){
              shinyalert(paste0("None of the entered genes in Analysis #",x," 1.3 are found in the expression dataset of the selected cancer project. Please revise your gene list entry in Analysis #",x,". Thank you."))
            }else{
              rv[[paste0("title_",x)]] <- paste0("Manual collection of genes (input n=",length(rv[[paste0("gs_m_",x)]])
                                                 ,", detected n=",n_col,")"
                                                 )
            }
          }
          
          req(rv[[paste0("title_",x)]] != "")

          # ---------- 3B. perform Surv if expression-like data --------
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
              rv[[paste0("cutoff_",x)]] <- paste0("<b>",results[["cutoff"]],"</b>")
              if(rv[["cutoff_all"]] == ""){
                rv[["cutoff_all"]] <- paste0("#",x,": ",rv[[paste0("cutoff_",x)]])
              }else{
                rv[["cutoff_all"]] <- paste0(rv[["cutoff_all"]],", #",x,": ",rv[[paste0("cutoff_",x)]])
              }
            }else{
              clow_id <- paste0("clow_",x)
              cutoff <- ifelse(is.null(input[[clow_id]]), 49, input[[clow_id]])
              df <- get_df_by_cutoff(data, cutoff)
              rv[[paste0("cutoff_",x)]] <- paste0(input[[clow_id]],"%")
            }
            
          # ---------- 3C. survival analysis on SNVs ---------
          }else if(input[[db_id]] == "snv"){
            # non-synonymous
            non_id <- paste0("nonsynonymous_",x)
            nons <- ifelse_rv(non_id)
            # # synonymous
            # syn_id <- paste0("synonymous_",x)
            # syns <- ifelse_rv(syn_id)
            # 
            # # render an error if a mutation class is selected twice
            # error <- 0
            # if(any(nons %in% syns)){
            #   error <- error + 1
            #   ol <- nons[nons %in% syns]
            #   shinyalert(paste0("You have selected ",paste0(ol,collapse = ", ")," in both non- and synonymous mutations in Analysis #",x,". Please refine your selection."))
            # }
            # 
            # req(error == 0)
            
            # create df for survival analysis
            df <- get_df_snv(data, nons)
            rv[[paste0("cutoff_",x)]] <- ""
          }
          
          # save df
          df_list[[x]] <- df
          if(x == 1){rv[[paste0("df_",x)]] <- df}

          # no of cases in each group
          lels <- levels(df$level)
          rv[[paste0("lels_",x)]] <- lapply(lels, function(x) table(df$level == x)["TRUE"] %>% unname(.))
          names(rv[[paste0("lels_",x)]]) <- lels
          
          # perform survival analysis
          cox_id <- paste0("cox_",x)
          rv[[cox_id]] <- cal_surv_rna(df,1)
          
          # save data type to rv
          rv[[paste0("data_type_",x)]] <- extract_mode
        }
      }
      
      # ------------- 3D. perform n=2 interaction Surv ---------
      if(rv$variable_n > 1){
        # generate interaction df
        df_combined <- Reduce(
          function(x, y) inner_join(x, dplyr::select(y, patient_id, level), by = "patient_id"),
          df_list
        )
        
        x_y <- c("x","y")[1:length(df_list)]
        df_combined[["level"]] <- apply(df_combined %>% dplyr::select(paste0("level.",x_y)),1,paste0,collapse="_")
        rv[["df_all"]] <- df_combined

        # no of cases in each group
        lels <- unique(df_combined$level) %>% sort(.,decreasing = T)
        df_combined$level <- factor(df_combined$level, levels = lels)
        rv[["lels_all"]] <- lapply(lels, function(x) table(df_combined$level == x)["TRUE"] %>% unname(.))
        names(rv[["lels_all"]]) <- lels
        
        # perform survival analysis
        rv[["cox_all"]] <- cal_surv_rna(df_combined,rv$variable_n)
        # saveRDS(rv[["cox_all"]], "cox_all")
      }
      
      # ------------- 3E. perform n=1 gender interaction Surv ---------
      if(rv$variable_n == 1){
        # generate interaction df
        df_combined <- df_list[[1]] %>% inner_join(dplyr::select(rv$df_survival, patient_id, gender), by = "patient_id")
        
        df_combined <- dplyr::rename(df_combined, level.x = level, level.y = gender)
        df_combined[["level"]] <- apply(df_combined %>% dplyr::select(paste0("level.",c("x","y"))),1,paste0,collapse="_")
        
        # no of cases in each group
        lels <- unique(df_combined$level) %>% sort(.,decreasing = T)
        if(length(lels) > 2){
          df_combined$level <- factor(df_combined$level, levels = lels)
          lels_y <- unique(df_combined$`level.y`) %>% sort(.,decreasing = T)
          df_combined$`level.y` <- factor(df_combined$`level.y`, levels = lels_y)
          rv[["df_gender"]] <- df_combined
          
          # level tallies
          rv[["lels_gender"]] <- lapply(lels, function(x) table(df_combined$level == x)["TRUE"] %>% unname(.))
          names(rv[["lels_gender"]]) <- lels
          
          # perform survival analysis
          rv[["cox_gender"]] <- cal_surv_rna(df_combined,2)
          rv[["title_gender"]] <- paste0(rv[["title_1"]]," vs Gender")
        }else{
          rv[["df_gender"]] <- ""
        }
      }
    })
  }
})
