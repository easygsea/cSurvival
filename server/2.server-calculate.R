observeEvent(input$confirm,{
  if(rv$projectStatus == "none"){
    shinyalert("Please select a project(s) to begin your analysis")
  }else{
    #------ 0. clear previous data --------
    clear_rds()
    rv$show_ui <- ""
    rv$quantile_graph <- c()
    rv$heatmap_df <- c()

    #------ 1. check if any errors by user ------
    error_g <- NULL; error_lib <- NULL; error_manual <- NULL; error_gs <- NULL
    error_norm_g <- NULL; error_norm_gs_db <- NULL; error_norm_gs_lib <- NULL
    error_norm_g_dup <- NULL; error_norm_gs_lib_dup <- NULL
    for(x in 1:rv$variable_n){
      cat_id <- paste0("cat_",x)
      # save data categories
      rv[[paste0("catr_",x)]] <- input[[cat_id]]
      rv[[paste0("dbr_",x)]] <- ""
      # check if errors in user's input
      if(input[[cat_id]] == "g"){
        g_ui_id <- paste0("g_",x)
        db_id <- paste0("db_",x)
        rv[[paste0("dbr_",x)]] <- input[[db_id]]
        g_ui_norm_id <- paste0("gnorm_",x)
        g_ui_norm_g_id <- paste0("gnorm_g_",x)
        if(input[[g_ui_id]] == ""){
          error_g <- c(error_g,x)
        }else if(input[[db_id]] == "rna"){
          if(input[[g_ui_norm_id]] == "g"){
            if(input[[g_ui_norm_g_id]] == ""){
              error_norm_g <- c(error_norm_g,x)
            }else if(input[[g_ui_norm_g_id]] == input[[g_ui_id]]){
              error_norm_g_dup <- c(error_norm_g_dup,x)
            }
          }else if(input[[g_ui_norm_id]] == "gs"){
            g_ui_norm_gs_db_id <- paste0("gnorm_gs_db_",x)
            g_ui_norm_gs_lib_id <- paste0("gnorm_gs_lib_",x)
            if(input[[g_ui_norm_gs_db_id]] == ""){
              error_norm_gs_db <- c(error_norm_gs_db,x)
            }else if(input[[g_ui_norm_gs_lib_id]] == ""){
              error_norm_gs_lib <- c(error_norm_gs_lib,x)
            }
          }
        }
      }else if(input[[cat_id]] == "gs"){
        gs_mode_id <- paste0("gs_mode_",x)
        if(input[[gs_mode_id]] == "lib"){
          gs_db_id <- paste0("gs_db_",x)
          if(input[[gs_db_id]] == ""){
            error_lib <- c(error_lib,x)
          }else{
            gs_lib_id <- paste0("gs_l_",x)
            g_ui_norm_id <- paste0("gnorm_",x)
            g_ui_norm_g_id <- paste0("gnorm_g_",x)
            g_ui_norm_gs_db_id <- paste0("gnorm_gs_db_",x)
            g_ui_norm_gs_lib_id <- paste0("gnorm_gs_lib_",x)
            if(input[[gs_lib_id]] == ""){
              error_gs <- c(error_gs,x)
            }else if(input[[g_ui_norm_id]] == "g"){
              if(input[[g_ui_norm_g_id]] == ""){
                error_norm_g <- c(error_norm_g,x)
              }
            }else if(input[[g_ui_norm_id]] == "gs"){
              if(input[[g_ui_norm_gs_db_id]] == ""){
                error_norm_gs_db <- c(error_norm_gs_db,x)
              }else{
                if(input[[g_ui_norm_gs_lib_id]] == ""){
                  error_norm_gs_lib <- c(error_norm_gs_lib,x)
                }else if(input[[g_ui_norm_gs_lib_id]] == input[[gs_lib_id]]){
                  error_norm_gs_lib_dup <- c(error_norm_gs_lib_dup,x)
                }
              }
            }
          }
        }else if(input[[gs_mode_id]] == "manual"){
          gs_manual_id <- paste0("gs_m_",x)
          if(length(unique(rv[[gs_manual_id]])) == 1){
            if(rv[[gs_manual_id]] == ""){
              error_manual <- c(error_manual,x)
            }
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
    }else if(!is.null(error_norm_g)){
      txt <- paste0("Analysis #",error_norm_g," step ",error_norm_g,".4b") %>% paste0(., collapse = ", ")
      txt <- paste0("Please select a normalization gene in ",txt,". Or select None in ",paste0(paste0(error_norm_g,".4a"), collapse = ", ")," to analyze gene expression effects without normalization.")
      shinyalert(txt)
    }else if(!is.null(error_norm_gs_db)){
      txt <- paste0("Analysis #",error_norm_gs_db," step ",error_norm_gs_db,".4b") %>% paste0(., collapse = ", ")
      txt <- paste0("Please select a database in ",txt,". Or select None in ",paste0(paste0(error_norm_gs_db,".4a"), collapse = ", ")," to analyze gene expression effects without normalization.")
      shinyalert(txt)
    }else if(!is.null(error_norm_gs_lib)){
      txt <- paste0("Analysis #",error_norm_gs_lib," step ",error_norm_gs_lib,".4c") %>% paste0(., collapse = ", ")
      txt <- paste0("Please select a gene set (GS) in ",txt,". Or select None in ",paste0(paste0(error_norm_gs_lib,".4a"), collapse = ", ")," to analyze gene expression effects without normalization.")
      shinyalert(txt)
    }else if(!is.null(error_norm_g_dup)){
      txt <- paste0("Analysis #",error_norm_g_dup," step ",error_norm_g_dup,".4b") %>% paste0(., collapse = ", ")
      txt <- paste0("Please select a normalization gene different to your gene of interest in ",txt,". Or select None in ",paste0(paste0(error_norm_g_dup,".4a"), collapse = ", ")," to analyze gene expression effects without normalization.")
      shinyalert(txt)
    }else if(!is.null(error_norm_gs_lib_dup)){
      txt <- paste0("Analysis #",error_norm_gs_lib_dup," step ",error_norm_gs_lib_dup,".4b") %>% paste0(., collapse = ", ")
      txt <- paste0("Please select a normalization GS different to your GS of interest in ",txt,". Or select None in ",paste0(paste0(error_norm_gs_lib_dup,".4a"), collapse = ", ")," to analyze gene expression effects without normalization.")
      shinyalert(txt)
    }

    req(is.null(error_g) & is.null(error_lib) & is.null(error_manual) & is.null(error_gs) & is.null(error_norm_g) & is.null(error_norm_gs_db) & is.null(error_norm_gs_lib) & is.null(error_norm_g_dup) & is.null(error_norm_gs_lib_dup))

    if(rv$depmap){
      error_depmap <- 0
      if(input$depmap_gene == ""){
        shinyalert(paste0("Please select a ",agene()," in iv. Select a ",agene()," to load data."))
        error_depmap <- 1
      }
      req(error_depmap == 0)
    }
    rv$depmapr <- rv$depmap
    rv$tcgar <- rv$tcga
    rv$targetr <- rv$target

    #------ 2. survival data processing ------
    # withProgress(value = 1, message = wait_msg("Performing analysis..."),{
      # update survival df
      # filter OS, DFS, or PFS
      rv$df_survival <- rv$df_survival_o
      if(rv$tcga){
        tcga_error <- 0
        if(rv$tcga_stype == "OS"){
          rv$df_survival <- rv$df_survival_o %>% dplyr::mutate(survival_days = OS.time) %>% dplyr::mutate(censoring_status = as.numeric(OS))
        }else if(rv$tcga_stype == "DSS"){
          rv$df_survival <- rv$df_survival_o %>% dplyr::mutate(survival_days = DSS.time) %>% dplyr::mutate(censoring_status = as.numeric(DSS))
        }else if(rv$tcga_stype == "DFI"){
          rv$df_survival <- rv$df_survival_o %>% dplyr::mutate(survival_days = DFI.time) %>% dplyr::mutate(censoring_status = as.numeric(DFI))
        # }else if(rv$tcga_stype == "pss"){
        #   rv$df_survival <- rv$df_survival_o %>% dplyr::filter(new_tumor_event_after_initial_treatment == "YES")
        }else if(rv$tcga_stype == "PFI"){
          rv$df_survival <- rv$df_survival_o %>% dplyr::mutate(survival_days = PFI.time) %>% dplyr::mutate(censoring_status = as.numeric(PFI))
        }
        rv$df_survival <- rv$df_survival %>% dplyr::filter(survival_days != "#N/A") %>%
          dplyr::select(patient_id,survival_days,censoring_status,gender)
        if(nrow(rv$df_survival) == 0){
          shinyalert(paste0("No cases found under category ",vector_names(rv$tcga_stype,tcga_stypes)))
          tcga_error <- 1
        }
        req(tcga_error == 0)
      }
      if(rv$target){
        rv$df_survival <- rv$df_survival %>% dplyr::filter(!is.na(survival_days))
      }
      if(!rv$depmap){
        # re-calculate survival days if censored by time
        if(input$censor_time_ymd != "none"){
          c_time <- ifelse_rv_na("censor_time")
          # calculate censoring time in days
          nnn <- c_time * ymd_unit[[input$censor_time_ymd]]
          # if time record larger than censoring time, convert status to 0
          rv[["df_survival"]][["censoring_status"]] <- ifelse(
            rv[["df_survival"]][["survival_days"]] > nnn,
            0,
            rv[["df_survival"]][["censoring_status"]]
          )
          # if time record larger than censoring time, convert to the time cap
          rv[["df_survival"]][["survival_days"]] <- ifelse(
            rv[["df_survival"]][["survival_days"]] > nnn,
            nnn,
            rv[["df_survival"]][["survival_days"]]
          )
          # mark down as time censored
          ttt <- gsub("s$","",tolower(vector_names(input$censor_time_ymd,ymd_names)))
          rv$censor_time_p <- sprintf(" censored at %s %ss (%s-%s survival)",c_time, ttt,c_time, ttt)
        }else{
          rv$censor_time_p <- ""
        }
        # error on censor time UI
        error_censor <- 0
        if(is.na(input$censor_time)){
          error_censor <- 1
          shinyalert("Please enter a valid censoring time.")
        }else if(input$censor_time < rv$censor_time_min){
          error_censor <- 1
          ymd_name <- vector_names(input$censor_time_ymd,ymd_names)
          shinyalert(paste0("Please enter a longer censoring time. Minimum: ",rv$censor_time_min," ",ymd_name,"."
                            ," Updated to default/previous entry: ",rv$censor_time," ",ymd_name,"."))
        }else if(input$censor_time > rv$censor_time_max){
          error_censor <- 1
          ymd_name <- vector_names(input$censor_time_ymd,ymd_names)
          shinyalert(paste0("Please enter a shorter censoring time. Maximum: ",rv$censor_time_max," ",ymd_name,"."
                            ," Updated to default/previous entry: ",rv$censor_time," ",ymd_name,"."))
        }else{
          rv$censor_time <- input$censor_time
          rv[[paste0("censor_time_",input$censor_time_ymd)]] <- input$censor_time
        }
        updateNumericInput(
          session,"censor_time", NULL,
          value = rv$censor_time
        )
        req(error_censor == 0)

        # automatic adjustment of time units
        d_max <- max(rv[["df_survival"]][["survival_days"]])
        if(d_max < 1800){rv$ymd_int_m<-5;rv$ymd_int_d<-200}
        if(d_max < 1200){rv$ymd_int_m<-3;rv$ymd_int_d<-100}
        if(d_max < 600){rv$ymd_int_m<-2;rv$ymd_int_d<-50}
        if(d_max < 400){rv$ymd_int_m<-1;rv$ymd_int_d<-30}
      }else{
        rv$censor_time_p <- rv$depmap_gene
      }

      # ------- begin analysis after error checking -----
    update_all(m=F)
      show_modal_spinner(
        color = "#BDD5EA", spin = "half-circle",
        text = HTML("<span style='font-size:125%; font-family: sans-serif;'><br>Analysis in progress ...<br><br>This might take a while when permutation is applied.<br><br>Please wait a minute. Thank you.</span>"))
      rv$try_error <- 0; rv$surv_plotted <- ""; rv$gsea_done <- ""; rv[["padj_perm"]] <- NULL
      if(!is.null(input$search_mode)){if(nchar(input$search_mode)>1){rv$search_mode <- input$search_mode}}
      rv$variable_nr <- rv$variable_n
      rv$scatter_gender <- NULL
      if(rv$variable_nr == 1){
        rv$plot_type <- "1"; if(rv$cor_method=="pearson" & !rv$depmap){rv$cor_method <- "kendall"}
        if(rv$depmap){
          rv$cor_method<-"pearson";rv$scatter_log_x<-F
          if(input$cat_1 == "gs"){rv$scatter_log_y <- F}
        }
      }else{
        rv$plot_type <- "all"; if(rv$cor_method!="pearson"){rv$cor_method <- "pearson"}
      }
      df_list <- list()
      rv[["title_all"]] = ""
      rv[["cutoff_all"]] = ""
      #------ 3. Loop from 1 to rv$variable_n ------
      gp_r <- ifelse_rv("risk_gp"); rv$risk_gpr <- gp_r
      updateRV_numeric("min_gp_size",1,90)
      updateRV_numeric("n_perm",10,1000)
      dtypes_tmp <- data_types()
      for(x in 1:rv$variable_n){
        cal_exp <- T; mut_exp_yyy <- F
        iter_mode_value <- paste0("iter_mode_",x)
        assign(iter_mode_value,F)
        rv[[paste0("iter_mode_value",x)]] <- F
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
          title_x <- paste0("title_",x)
          rv[[title_x]] <- ""
          # mode of extracting gene expression/mutation data in extract_gene_data
          extract_mode <- input_mode(x)

          # title for the survival plot
          rv[[title_x]] <- ifelse(input[[cat_id]] == "g", input[[g_ui_id]], ifelse(input[[gs_mode_id]]=="lib",input[[gs_lib_id]],"Manual GS"))
          if(input[[cat_id]] == "g"){rv[[title_x]] <- paste0(rv[[title_x]]," ",tolower(names(dtypes_tmp)[dtypes_tmp == extract_mode]))}
          # check if normalized
          g_ui_norm_id <- paste0("gnorm_",x)
          if(((input[[cat_id]] == "g" & input[[db_id]] == "rna") | (input[[cat_id]] == "gs" & input[[gs_mode_id]] == "lib")) & input[[g_ui_norm_id]] != "none"){
            if(input[[g_ui_norm_id]] == "g"){
              rv[[title_x]] <- paste0(rv[[title_x]]," normalized by ", input[[paste0("gnorm_g_",x)]])
            }else if(input[[g_ui_norm_id]] == "gs"){
              rv[[title_x]] <- paste0(rv[[title_x]]," normalized by ", input[[paste0("gnorm_gs_lib_",x)]])
            }
          }
          # combine titles
          if(rv[["title_all"]] == ""){
            rv[["title_all"]] <- rv[[title_x]]
          }else{
            rv[["title_all"]] <- paste0(c(rv[["title_all"]],rv[[title_x]]),collapse = " vs ")
          }
          # extract gene expression/mutation data
          data <- extract_gene_data(x,extract_mode)
          if(rv$flagged == "y"){
            data <- data %>% dplyr::filter(!patient_id %in% flagged_cases)
          }

          # check if manually input genes exist in database
          if(extract_mode == "manual"){
            n_col <- ncol(data) - 1
            if(n_col < 1){
              shinyalert(paste0("None of the entered genes in Analysis #",x," 1.3 are found in the expression dataset of the selected cancer project. Please revise your gene list entry in Analysis #",x,". Thank you."))
              remove_modal_spinner()
            }else{
              rv[[paste0("title_",x)]] <- paste0("Manual collection of genes (input n=",length(rv[[paste0("gs_m_",x)]])
                                                 ,", detected n=",rv[[paste0("gs_m_len",x)]],")"
                                                 )
            }
          }

          req(rv[[paste0("title_",x)]] != "")

          # ---------- 3B. perform Surv if expression-like data --------
          min_value <- paste0("min_",x); max_value <- paste0("max_",x); step_value <- paste0("step_",x)
          assign(min_value, 0);assign(max_value, 1);assign(step_value, .1)
          if(exp_yyy(extract_mode)){
            iter_id <- paste0("iter_",x)
            yn <- ifelse(is.null(input[[iter_id]]), T, input[[iter_id]] == "iter")
            assign(iter_mode_value,yn)
            rv[[paste0("iter_mode_value",x)]] <- yn
            if(yn){
              lower_id <- paste0("lower_",x); higher_id <- paste0("upper_",x); step_id <- paste0("step_",x)
              min <- rv[[paste0("minF",x)]] <- ifelse(is.null(input[[lower_id]]), rv[[lower_id]], input[[lower_id]])
              max <- rv[[paste0("maxF",x)]] <- ifelse(is.null(input[[higher_id]]), rv[[higher_id]], input[[higher_id]])
              step <- rv[[paste0("stepF",x)]] <- ifelse(is.null(input[[step_id]]), rv[[step_id]], input[[step_id]])
              rv[[paste0("dataF",x)]] <- data

              assign(min_value, min);assign(max_value, max);assign(step_value, step)

              enough_error <- 0
              # check if bivariate
              if(rv$variable_n > 1){
                # ANOVA if bivariate expression-like
                if(x == 1){
                  cal_exp <- F
                }else if(x == 2 & exp_yyy(input_mode(1))){
                  if(exp_iter_yyy(1)){
                    results <- get_info_most_significant_rna(rv[["dataF1"]], rv[["minF1"]], rv[["maxF1"]], rv[["stepF1"]], num=x, data2=data, min2=min, max2=max, step2=step, gp=gp_r)
                  }else{
                    results <- get_info_most_significant_rna(data, min, max, step, data2=rv[["dataF1"]], cat=gp_r)
                  }
                }else if(x == 2 & !exp_yyy(input_mode(1))){
                  results <- get_info_most_significant_rna(data, min, max, step, data2=rv[["dataF1"]], cat=gp_r)
                }
              }else if(rv$variable_n == 1){
                results <- get_info_most_significant_rna(data, min, max, step)
              }
              #extract heatmap_df needed for heatmap
              
              if(cal_exp){
                if(is.null(results)){
                  enough_error <- 1
                  txt <- "No significant result identified"
                  if(rv$variable_n > 0 & if_any_exp()){
                    txt <- paste0(txt," at the defined threshold of at least ",rv$min_gp_size,"% cases in each risk subgroup")
                  }
                  if(input$censor_time_ymd != "none"){
                    txt <- paste0(txt,", at censoring time ",input$censor_time," ",vector_names(input$censor_time_ymd,ymd_names))
                  }
                  rv$try_error <- 1
                  shinyalert(paste0(txt,"."))
                  remove_modal_spinner()
                }

                req(enough_error == 0)
                # print(str(results))

                #CHANGE THE RESULT FUNCTION AND ADD RV$quantile_graph HERE
                if("p_df" %in% names(results)){
                  if(x > 1 & exp_yyy(input_mode(1)) & exp_iter_yyy(1)){
                    rv[["quantile_graph"]][[1]] <- results[["p_df"]][[1]]
                    rv[["quantile_graph"]][[x]] <- results[["p_df"]][[x]]
                  }else{
                    rv[["quantile_graph"]][[x]] <- results[["p_df"]]
                  }
                }
                rv$heatmap_df <- results[["heatmap_df"]]
                # saveRDS(results[["heatmap_df"]], "heatmap_df.rds")
                # View(rv$heatmap_df)
                # extract most significant df
                df <- results[["df"]]
                rv[["padj_perm"]] <- results[["p.adj"]]
                cutoff_n <- paste0("cutoff_",x)
                if(rv$variable_n > 1){
                  cutoff_tmp <- results[["cutoff"]]
                  cutoff_tmp <- ifelse(length(cutoff_tmp)>1, cutoff_tmp[2], cutoff_tmp)
                  rv[[cutoff_n]] <- paste0("<b>",cutoff_tmp,"</b>")
                  if(x==2 & exp_yyy(input_mode(1))){
                    if(exp_iter_yyy(1)){
                      rv[["cutoff_1"]] <- paste0("<b>",results[["cutoff"]][1],"</b>")
                    }
                    rv[["cutoff_all"]] <- paste0("#1: ",rv[["cutoff_1"]],", #",x,": ",rv[[cutoff_n]])
                  }else{
                    rv[["cutoff_all"]] <- paste0("#",x,": ",rv[[cutoff_n]])
                  }
                  if(gp_r == "All"){
                    df_lels <- as.character(df[["level"]])
                    df_lels1 <- sapply(df_lels, function(x) strsplit(x,"_")[[1]][1])
                    df_lels2 <- sapply(df_lels, function(x) strsplit(x,"_")[[1]][2])
                    df[["level.x"]] <- factor(df_lels1); df[["level.x"]] <- relevel(df[["level.x"]], ref = sort(levels(df[["level.x"]]))[2])
                    df[["level.y"]] <- factor(df_lels2); df[["level.y"]] <- relevel(df[["level.y"]], ref = sort(levels(df[["level.y"]]))[2])
                  }
                  # save df
                  if(rv$depmap){
                    df_1 <- df %>% dplyr::select(patient_id,dependency,level.x,value.x)
                    df_2 <- df %>% dplyr::select(patient_id,dependency,level.y,value.y)
                    colnames(df_1) <- c("patient_id","dependency","level","value")
                    colnames(df_2) <- c("patient_id","dependency","level","value")
                  }else{
                    df_1 <- df %>% dplyr::select(patient_id,survival_days,censoring_status,level.x,value.x)
                    df_2 <- df %>% dplyr::select(patient_id,survival_days,censoring_status,level.y,value.y)
                    colnames(df_1) <- c("patient_id","survival_days","censoring_status","level","value")
                    colnames(df_2) <- c("patient_id","survival_days","censoring_status","level","value")
                  }
                  df_list[[1]] <- df_1; df_list[[2]] <- df_2
                  rv[["df_1"]] <- df_1
                  rv[["df_2"]] <- df_2

                  # no of cases in each group
                  lels1 <- levels(df_1$level); lels2 <- levels(df_2$level)
                  rv[["lels_1"]] <- lapply(lels1, function(x) table(df_1$level == x)["TRUE"] %>% unname(.))
                  rv[["lels_2"]] <- lapply(lels2, function(x) table(df_2$level == x)["TRUE"] %>% unname(.))
                  names(rv[["lels_1"]]) <- lels1; names(rv[["lels_2"]]) <- lels2

                  # perform survival analysis
                  if(x==2 & exp_yyy(input_mode(1))){
                    if(exp_iter_yyy(1)){
                      rv[["cox_1"]] <- cal_surv_rna(df_1,1,rv[["minF1"]],rv[["maxF1"]],rv[["stepF1"]])
                    }
                  }
                  rv[["cox_2"]] <- cal_surv_rna(df_2,1,rv[["minF2"]],rv[["maxF2"]],rv[["stepF2"]])

                  # switch to FALSE to prevent errors
                  cal_exp <- F
                }else{
                  rv[[paste0("cutoff_",x)]] <- paste0("<b>",results[["cutoff"]],"</b>")
                  if(rv[["cutoff_all"]] == ""){
                    rv[["cutoff_all"]] <- paste0("#",x,": ",rv[[paste0("cutoff_",x)]])
                  }else{
                    rv[["cutoff_all"]] <- paste0(rv[["cutoff_all"]],", #",x,": ",rv[[paste0("cutoff_",x)]])
                  }
                }
              }
            }else{
              clow_id <- paste0("clow_",x)
              cutoff <- ifelse(is.null(input[[clow_id]]), 49, input[[clow_id]])
              df <- get_df_by_cutoff(data, cutoff)
              rv[[paste0("cutoff_",x)]] <- paste0("<b>",input[[clow_id]],"%</b>")
              if(rv[["cutoff_all"]] == ""){
                rv[["cutoff_all"]] <- paste0("#",x,": ",rv[[paste0("cutoff_",x)]])
              }else{
                rv[["cutoff_all"]] <- paste0(rv[["cutoff_all"]],", #",x,": ",rv[[paste0("cutoff_",x)]])
              }
              rv[[paste0("dataF",x)]] <- df %>% dplyr::select(patient_id,level)

              if(x == 2 & exp_yyy(input_mode(1)) & exp_iter_yyy(1)){
                cal_exp <- T; mut_exp_yyy <- T
                df_list <- cal_conti_cat_interaction(x,gp_r,df_list)
                if(is.null(df_list)){shinyalert(paste0("Not enough data for the selected combination of genomic features w/ the requirement of mininum ",rv$min_gp_size,"% samples in each risk subgroup"))};remove_modal_spinner();req(!is.null(df_list))
                df <- df_list[[1]]; results <- df_list[[3]]; df_list <- df_list[[2]]
                rv[["cutoff_1"]] <- paste0("<b>",results[["cutoff"]],"</b>")
                rv[["cutoff_all"]] <- paste0("#1: ",rv[["cutoff_1"]],", #2: ",rv[["cutoff_2"]])
              }
            }

          # ---------- 3C. survival analysis on SNVs ---------
          }else if(extract_mode == "snv"){
            snv_id <- paste0("snv_method_",x)
            if(is.null(input[[snv_id]])){
              updateSelectizeInput(session,snv_id,selected = "mutect")
            }
            # non-synonymous
            non_id <- paste0("nonsynonymous_",x)
            nons <- ifelse_rv(non_id)
            if(is.null(input[[non_id]])){
              updateSelectizeInput(session,non_id,selected = variant_types_non)
            }
            # synonymous
            syn_id <- paste0("synonymous_",x)
            syns <- ifelse_rv(syn_id)
            if(is.null(input[[syn_id]])){
              updateSelectizeInput(session,syn_id,selected = variant_types_syn)
            }
            # render an error if a mutation class is selected twice
            error <- 0
            if(any(nons %in% syns)){
              error <- error + 1
              ol <- nons[nons %in% syns]
              shinyalert(paste0("You have selected ",paste0(ol,collapse = ", ")," in both Mutated and Other in Analysis #",x,". Please refine your selection."))
              remove_modal_spinner()
            }

            req(error == 0)

            # create df for survival analysis
            df <- get_df_snv(data, nons, syns)
            error_snv <- 0
            if(length(unique(df$level)) < 2){
              shinyalert("Not enough data found for the selected endpoint and parameters")
              error_snv <- 1
              rv$try_error <- 1
              remove_modal_spinner()
            }
            req(error_snv == 0)
            cutoff_x <- paste0("cutoff_",x)
            rv[[cutoff_x]] <- ""
            rv[[paste0("dataF",x)]] <- df %>% dplyr::select(patient_id,level)

            if(x == 2 & exp_yyy(input_mode(1)) & exp_iter_yyy(1)){
              cal_exp <- T; mut_exp_yyy <- T
              df_list <- cal_conti_cat_interaction(x,gp_r,df_list)
              if(is.null(df_list)){shinyalert(paste0("Not enough data for the selected combination of genomic features w/ the requirement of mininum ",rv$min_gp_size,"% samples in each risk subgroup"))};remove_modal_spinner();req(!is.null(df_list))
              df <- df_list[[1]]; df_list <- df_list[[2]]
            }
          # ---------- 3D. survival analysis on CNVs ---------
          }else if(!rv$depmap & extract_mode == "cnv"){
            cnv_mode <- ifelse_rv(paste0("cnv_par_",x))
            rv[[paste0("dataF",x)]] <- data

            results <- get_info_most_significant_cnv(data,cnv_mode)

            # extract most significant df
            df <- results[["df"]]
            tmp <- as.character(df$level)
            names(tmp) <- df$patient_id
            rv[[paste0("mutations_",x)]] <- tmp
            cutoff_x <- paste0("cutoff_",x)
            rv[[cutoff_x]] <- ""
            rv[[paste0("dataF",x)]] <- df %>% dplyr::select(patient_id,level)

            if(x == 2 & exp_yyy(input_mode(1)) & exp_iter_yyy(1)){
              cal_exp <- T; mut_exp_yyy <- T
              df_list <- cal_conti_cat_interaction(x,gp_r,df_list)
              if(is.null(df_list)){shinyalert(paste0("Not enough data for the selected combination of genomic features w/ the requirement of mininum ",rv$min_gp_size,"% samples in each risk subgroup"))};remove_modal_spinner();req(!is.null(df_list))
              df <- df_list[[1]]; df_list <- df_list[[2]]
            }
          }

          if(cal_exp){
            # save df
            df_list[[x]] <- df
            rv[[paste0("df_",x)]] <- df
            # no of cases in each group
            lels <- levels(df$level)
            rv[[paste0("lels_",x)]] <- lapply(lels, function(x) table(df$level == x)["TRUE"] %>% unname(.))
            names(rv[[paste0("lels_",x)]]) <- lels

            # perform survival analysis
            cox_id <- paste0("cox_",x)
            rv[[cox_id]] <- cal_surv_rna(df,1,get(min_value),get(max_value),get(step_value),iter_mode=get(iter_mode_value))
            # saveRDS(rv[[cox_id]], "cox_dm")
          }
          # save data type to rv
          rv[[paste0("data_type_",x)]] <- extract_mode
        }
      }

      # ------------- 3E. perform n=2 interaction Surv ---------
      if(rv$variable_n > 1){
        # generate interaction df
        df_combined <- Reduce(
          function(x, y) inner_join(x, dplyr::select(y, patient_id, level, value), by = "patient_id"),
          df_list
        ) %>% distinct()

        x_y <- c("x","y")[1:length(df_list)]
        df_combined[["level"]] <- apply(df_combined %>% dplyr::select(paste0("level.",x_y)),1,paste0,collapse="_")
        # combine subgroups if indicated
        error_gp <- 0
        if(gp_r != "All"){
          all_lels <- unique(df_combined[["level"]])
          if(!gp_r %in% all_lels){
            if(length(all_lels) < 4){
              gp_rr <- strsplit(gp_r,"_")[[1]]
              if(!(any(grepl(gp_rr[1],all_lels)) & any(grepl(gp_rr[2],all_lels)))){
                error_gp <- 1
              }
            }
          }
          if(error_gp == 0){
            df_combined <- assign_gp(df_combined,gp_r)
            lels <- levels(df_combined$level)
          }
        }else{
          lels <- unique(df_combined$level) %>% sort(.,decreasing = T)
          df_combined$level <- factor(df_combined$level, levels = lels)
        }
        # test if at least 10% (default) of samples in each subgroup
        if(rv$variable_nr > 1){
          lel_min_gps <- table(df_combined$level)
          min_gp_size_n <- rv$min_gp_size / 100 * sum(lel_min_gps)
          lel_min <- min(lel_min_gps)
          if(lel_min < min_gp_size_n){
            lel_min_gps <- lel_min_gps[lel_min_gps < min_gp_size_n]
            lel_min_gps <- sapply(lel_min_gps, function(x){
              paste0(names(lel_min_gps[lel_min_gps == x])," (n=",x,")")
            }) %>% paste0(., collapse = ", ")
            error_gp <- 2
          }
        }
        if(error_gp > 0){
          if(error_gp == 1){
            shinyalert("The selected subgroup does not have enough samples for analysis")
          }else if(error_gp == 2){
            shinyalert(paste0(lel_min_gps," do(es) not meet the minimum number of cases criterion (min n=",min_gp_size_n,", ",rv$min_gp_size,"% of the population)."))
          }
          rv$try_error <- 1
          remove_modal_spinner()
        }
        req(error_gp == 0)
        rv[["df_all"]] <- df_combined

        # no of cases in each group
        rv[["lels_all"]] <- lapply(lels, function(x) table(df_combined$level == x)["TRUE"] %>% unname(.))
        names(rv[["lels_all"]]) <- lels

        # perform survival analysis
        rv[["cox_all"]] <- cal_surv_rna(df_combined,rv$variable_n,0,1,.1,iter_mode=F,gp=gp_r)
        # if(rv$depmap){
        #   p_adj_tmp <- rv[["cox_all"]][["p.adj"]]
        #   for(x in 1:rv$variable_n){
        #     cox_x <- paste0("cox_",x)
        #     if(!is.null(rv[[cox_x]][["p.adj"]])){p_adj_tmp <- ifelse(is.null(p_adj_tmp),0,p_adj_tmp) + rv[[cox_x]][["p.adj"]]}
        #   }
        #   rv[["cox_all"]][["p.adj"]] <- p_adj_tmp
        # }else{
        #   p_adj_tmp_km <- rv[["cox_all"]][["km"]][["p.adj"]]
        #   for(x in 1:rv$variable_n){
        #     cox_x <- paste0("cox_",x)
        #     if(!is.null(rv[[cox_x]][["km"]][["p.adj"]])){p_adj_tmp_km <- ifelse(is.null(p_adj_tmp_km),0,p_adj_tmp_km) + rv[[cox_x]][["km"]][["p.adj"]]}
        #     tmp <- rv[["cox_all"]][["cox"]][["p.adj"]]#[x]
        #     if(is.null(tmp)){tmp <- 0};if(is.na(tmp)){tmp <- 0}
        #     if(!is.null(rv[[cox_x]][["cox"]][["p.adj"]])){rv[["cox_all"]][["cox"]][["p.adj"]] <- tmp + rv[[cox_x]][["cox"]][["p.adj"]]} #[x] correct_p(rv[["cox_all"]][["cox"]][["p"]][x],get(paste0("min_",x)),get(paste0("max_",x)),get(paste0("step_",x)))}
        #   }
        #   rv[["cox_all"]][["km"]][["p.adj"]] <- p_adj_tmp_km
        # }
        # # p adj on pairwise heatmap
        # if(!rv$depmap & rv$km_mul_padj == "padj"){
        #   if(!is.null(rv[["cox_all"]][["km"]][["p.adj"]])){
        #     km2 <- rv[["cox_all"]][["km"]][["stats"]][[1]]
        #     km2$p.value <- km2$p.value + rv[["cox_all"]][["km"]][["p.adj"]]
        #     km2$p.value <- ifelse(km2$p.value > 1, 1, km2$p.value)
        #     rv[["cox_all"]][["km"]][["stats"]][[1]] <- km2
        #   }
        # }else if(rv$depmap & rv$km_mul_padj == "padj"){
        #   if(!is.null(rv[["cox_all"]][["p.adj"]])){
        #     km2 <- rv[["cox_all"]][["stats"]][[1]]
        #     km2$p.value <- km2$p.value + rv[["cox_all"]][["p.adj"]]
        #     km2$p.value <- ifelse(km2$p.value > 1, 1, km2$p.value)
        #     rv[["cox_all"]][["stats"]][[1]] <- km2
        #   }
        # }

        # saveRDS(rv[["cox_all"]], "cox_all")
      }

      # ------------- 3F. perform n=1 gender interaction Surv ---------
      if(rv$variable_n == 1){
        # generate interaction df
        df_combined <- df_list[[1]] %>% inner_join(dplyr::select(rv$df_survival, patient_id, gender), by = "patient_id") %>%
          dplyr::filter(gender != "Unknown")

        # gender types
        gender_cats <- df_combined %>% dplyr::select(gender) %>% unique(.) %>% unlist(.) %>% unname(.)
        rv$scatter_gender <- rv[["genders"]] <- gender_cats[tolower(gender_cats) %in% c("male","female")]

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
          rv[["cox_gender"]] <- cal_surv_rna(df_combined,2,0,1,.1,iter_mode=F)
          # if(rv$depmap){
          #   if(!is.null(rv[["cox_1"]][["p.adj"]])){rv[["cox_gender"]][["p.adj"]] <- ifelse(is.null(rv[["cox_gender"]][["p.adj"]]),0,rv[["cox_gender"]][["p.adj"]]) + rv[["cox_1"]][["p.adj"]]}
          # }else{
          #   if(!is.null(rv[["cox_1"]][["km"]][["p.adj"]])){rv[["cox_gender"]][["km"]][["p.adj"]] <- ifelse(is.null(rv[["cox_gender"]][["km"]][["p.adj"]]),0,rv[["cox_gender"]][["km"]][["p.adj"]]) + rv[["cox_1"]][["km"]][["p.adj"]]}
          #   if(!is.null(rv[["cox_1"]][["cox"]][["p.adj"]])){
          #     rv[["cox_gender"]][["cox"]][["p.adj"]][1] <- rv[["cox_1"]][["cox"]][["p.adj"]]#correct_p(rv[["cox_gender"]][["cox"]][["p"]][1],get("min_1"),get("max_1"),get("step_1"))
          #     # rv[["cox_gender"]][["cox"]][["p.adj"]][2] <- rv[["cox_gender"]][["cox"]][["p.adj"]][2]
          #     # rv[["cox_gender"]][["cox"]][["p.adj"]] <- ifelse(rv[["cox_gender"]][["cox"]][["p.adj"]] > 1, 1, rv[["cox_gender"]][["cox"]][["p.adj"]])
          #   }
          # }
          # # p adj on pairwise heatmap
          # if(!rv$depmap & rv$km_mul_padj == "padj"){
          #   if(!is.null(rv[["cox_1"]][["km"]][["p.adj"]])){
          #     km2 <- rv[["cox_gender"]][["km"]][["stats"]][[1]]
          #     km2$p.value <- km2$p.value + rv[["cox_gender"]][["km"]][["p.adj"]]
          #     km2$p.value <- ifelse(km2$p.value > 1, 1, km2$p.value)
          #     rv[["cox_gender"]][["km"]][["stats"]][[1]] <- km2
          #   }
          # }else if(rv$depmap & rv$km_mul_padj == "padj"){
          #   if(!is.null(rv[["cox_1"]][["p.adj"]])){
          #     km2 <- rv[["cox_gender"]][["stats"]][[1]]
          #     km2$p.value <- km2$p.value + rv[["cox_gender"]][["p.adj"]]
          #     km2$p.value <- ifelse(km2$p.value > 1, 1, km2$p.value)
          #     rv[["cox_gender"]][["stats"]][[1]] <- km2
          #   }
          # }
          rv[["title_gender"]] <- paste0(rv[["title_1"]]," vs Gender")
        }else{
          rv[["df_gender"]] <- unique(df_combined[["level.y"]])
        }
      }
    # })

    # update parameters
    exp_iter_y <- sapply(1:rv$variable_nr,function(x) exp_yyy(input_mode(x)) & exp_iter_yyy(x))
    rv$exp_iter_yyy <- all(exp_iter_y); rv$exp_iter_yyy_any <- any(exp_iter_y)
    if(rv$exp_iter_yyy_any){rv$cox_km <- rv$min_p_kc}else{rv$cox_km <- "km"}
    if(rv$tcga){rv$plot_sstype <- rv$plot_stype}
    else if(rv$target){rv$plot_sstype <- "Overall survival (OS)"}
    else if(rv$depmap){rv$plot_sstype <- dependency_names();rv$depmap_gener <- rv$depmap_gene}
    rv$annot_cells_y <- ""
    rv$annot_data_points_y <- ""
    rv$show_ui <- "yes"
    rv$surv_plotted <- ""
    rv$analysis_no <- rv$analysis_no + 1
    remove_modal_spinner()
  }
})
