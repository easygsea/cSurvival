observeEvent(input$confirm,{
  if(rv$projectStatus == "none"){
    shinyalert("Please select a project(s) to begin your analysis")
  }else{
    #------ 0. clear previous data --------
    clear_rds()
    rv$show_ui <- ""
    rv$quantile_graph <- c()

    #------ 1. check if any errors by user ------
    error_g <- NULL; error_lib <- NULL; error_manual <- NULL; error_gs <- NULL
    error_norm_g <- NULL; error_norm_gs_db <- NULL; error_norm_gs_lib <- NULL
    error_norm_g_dup <- NULL; error_norm_gs_lib_dup <- NULL
    for(x in 1:rv$variable_n){
      cat_id <- paste0("cat_",x)
      if(input[[cat_id]] == "g"){
        g_ui_id <- paste0("g_",x)
        db_id <- paste0("db_",x)
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
    withProgress(value = 1, message = wait_msg("Performing analysis..."),{
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
          rv$censor_time_p <- sprintf(" %s-%s survival",c_time, gsub("s$","",tolower(vector_names(input$censor_time_ymd,ymd_names))))
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
      rv$try_error <- 0; rv$surv_plotted <- ""; rv$gsea_done <- ""
      rv$variable_nr <- rv$variable_n
      rv$scatter_gender <- NULL
      if(rv$variable_nr == 1){
        rv$plot_type <- "1"; if(rv$cor_method=="pearson" & !rv$depmap){rv$cor_method <- "kendall"}
        if(rv$depmap){rv$cor_method<-"pearson";rv$scatter_log_x<-F}
      }else{
        rv$plot_type <- "all"; if(rv$cor_method!="pearson"){rv$cor_method <- "pearson"}
      }
      df_list <- list()
      rv[["title_all"]] = ""
      rv[["cutoff_all"]] = ""
      #------ 3. Loop from 1 to rv$variable_n ------
      for(x in 1:rv$variable_n){
        iter_mode_value <- paste0("iter_mode_",x)
        assign(iter_mode_value,F)
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
          rv[[title_x]] <- ifelse(input[[cat_id]] == "g", input[[g_ui_id]], input[[gs_lib_id]])
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
            }else{
              rv[[paste0("title_",x)]] <- paste0("Manual collection of genes (input n=",length(rv[[paste0("gs_m_",x)]])
                                                 ,", detected n=",n_col,")"
                                                 )
            }
          }

          req(rv[[paste0("title_",x)]] != "")

          # ---------- 3B. perform Surv if expression-like data --------
          min_value <- paste0("min_",x); max_value <- paste0("max_",x); step_value <- paste0("step_",x)
          assign(min_value, 0);assign(max_value, 1);assign(step_value, .1)
          if(extract_mode != "snv" & ((!rv$depmap & extract_mode != "cnv") | rv$depmap)){
            iter_id <- paste0("iter_",x)
            yn <- ifelse(is.null(input[[iter_id]]), T, input[[iter_id]] == "iter")
            assign(iter_mode_value,yn)
            if(yn){
              lower_id <- paste0("lower_",x); higher_id <- paste0("upper_",x); step_id <- paste0("step_",x)
              min <- ifelse(is.null(input[[lower_id]]), rv[[lower_id]], input[[lower_id]])
              max <- ifelse(is.null(input[[higher_id]]), rv[[higher_id]], input[[higher_id]])
              step <- ifelse(is.null(input[[step_id]]), rv[[step_id]], input[[step_id]])
              
              assign(min_value, min);assign(max_value, max);assign(step_value, step)

              enough_error <- 0
              results <- get_info_most_significant_rna(data, min, max, step, mode=input[[cat_id]])
              if(is.null(results)){
                enough_error <- 1
                txt <- "The selected project does not have enough data for the selected gene/locus/gene set"
                if(input$censor_time_ymd != "none"){
                  txt <- paste0(txt,", at censoring time ",input$censor_time," ",vector_names(input$censor_time_ymd,ymd_names))
                }
                shinyalert(paste0(txt,"."))
              }

              req(enough_error == 0)
              
              #CHANGE THE RESULT FUNCTION AND ADD RV$quantile_graph HERE
              if("p_df" %in% names(results)){
                rv[["quantile_graph"]][[x]] <- results[["p_df"]]
              }
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
            }

            req(error == 0)

            # create df for survival analysis
            df <- get_df_snv(data, nons, syns)
            error_snv <- 0
            if(length(unique(df$level)) < 2){
              shinyalert("Not enough data found for the selected endpoint and parameters")
              error_snv <- 1
              rv$try_error <- 1
            }
            req(error_snv == 0)
            rv[[paste0("cutoff_",x)]] <- ""

          # ---------- 3D. survival analysis on CNVs ---------
          }else if(!rv$depmap & extract_mode == "cnv"){
            cnv_mode <- ifelse_rv(paste0("cnv_par_",x))

            results <- get_info_most_significant_cnv(data,cnv_mode)
            
            # extract most significant df
            df <- results[["df"]]
            tmp <- as.character(df$level)
            names(tmp) <- df$patient_id
            rv[[paste0("mutations_",x)]] <- tmp
            rv[[paste0("cutoff_",x)]] <- paste0("<b>",results[["cutoff"]],"</b>")
            if(rv[["cutoff_all"]] == ""){
              rv[["cutoff_all"]] <- paste0("#",x,": ",rv[[paste0("cutoff_",x)]])
            }else{
              rv[["cutoff_all"]] <- paste0(rv[["cutoff_all"]],", #",x,": ",rv[[paste0("cutoff_",x)]])
            }

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
          rv[[cox_id]] <- cal_surv_rna(df,1,get(min_value),get(max_value),get(step_value),iter_mode=get(iter_mode_value))
          # saveRDS(rv[[cox_id]], "cox_dm")

          # save data type to rv
          rv[[paste0("data_type_",x)]] <- extract_mode

        }
      }

      # ------------- 3E. perform n=2 interaction Surv ---------
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
        rv[["cox_all"]] <- cal_surv_rna(df_combined,rv$variable_n,0,1,.1,iter_mode=F)
        if(rv$depmap){
          p_adj_tmp <- rv[["cox_all"]][["p.adj"]]
          for(x in 1:rv$variable_n){
            cox_x <- paste0("cox_",x)
            if(!is.null(rv[[cox_x]][["p.adj"]])){p_adj_tmp <- ifelse(is.null(p_adj_tmp),0,p_adj_tmp) + rv[[cox_x]][["p.adj"]]}
          }
          rv[["cox_all"]][["p.adj"]] <- p_adj_tmp
        }else{
          p_adj_tmp_km <- rv[["cox_all"]][["km"]][["p.adj"]]
          for(x in 1:rv$variable_n){
            cox_x <- paste0("cox_",x)
            if(!is.null(rv[[cox_x]][["km"]][["p.adj"]])){p_adj_tmp_km <- ifelse(is.null(p_adj_tmp_km),0,p_adj_tmp_km) + rv[[cox_x]][["km"]][["p.adj"]]}
            tmp <- rv[["cox_all"]][["cox"]][["p.adj"]][x]
            if(is.null(tmp)){tmp <- 0};if(is.na(tmp)){tmp <- 0}
            if(!is.null(rv[[cox_x]][["cox"]][["p.adj"]])){rv[["cox_all"]][["cox"]][["p.adj"]][x] <- tmp + rv[[cox_x]][["cox"]][["p.adj"]]} # correct_p(rv[["cox_all"]][["cox"]][["p"]][x],get(paste0("min_",x)),get(paste0("max_",x)),get(paste0("step_",x)))}
          }
          rv[["cox_all"]][["km"]][["p.adj"]] <- p_adj_tmp_km
        }
        # p adj on pairwise heatmap
        if(!rv$depmap & rv$km_mul_padj == "padj"){
          if(!is.null(rv[["cox_all"]][["km"]][["p.adj"]])){
            km2 <- rv[["cox_all"]][["km"]][["stats"]][[1]]
            km2$p.value <- km2$p.value + rv[["cox_all"]][["km"]][["p.adj"]]
            km2$p.value <- ifelse(km2$p.value > 1, 1, km2$p.value)
            rv[["cox_all"]][["km"]][["stats"]][[1]] <- km2
          }
        }else if(rv$depmap & rv$km_mul_padj == "padj"){
          if(!is.null(rv[["cox_all"]][["p.adj"]])){
            km2 <- rv[["cox_all"]][["stats"]][[1]]
            km2$p.value <- km2$p.value + rv[["cox_all"]][["p.adj"]]
            km2$p.value <- ifelse(km2$p.value > 1, 1, km2$p.value)
            rv[["cox_all"]][["stats"]][[1]] <- km2
          }
        }

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
          if(rv$depmap){
            if(!is.null(rv[["cox_1"]][["p.adj"]])){rv[["cox_gender"]][["p.adj"]] <- ifelse(is.null(rv[["cox_gender"]][["p.adj"]]),0,rv[["cox_gender"]][["p.adj"]]) + rv[["cox_1"]][["p.adj"]]}
          }else{
            if(!is.null(rv[["cox_1"]][["km"]][["p.adj"]])){rv[["cox_gender"]][["km"]][["p.adj"]] <- ifelse(is.null(rv[["cox_gender"]][["km"]][["p.adj"]]),0,rv[["cox_gender"]][["km"]][["p.adj"]]) + rv[["cox_1"]][["km"]][["p.adj"]]}
            if(!is.null(rv[["cox_1"]][["cox"]][["p.adj"]])){
              rv[["cox_gender"]][["cox"]][["p.adj"]][1] <- rv[["cox_1"]][["cox"]][["p.adj"]]#correct_p(rv[["cox_gender"]][["cox"]][["p"]][1],get("min_1"),get("max_1"),get("step_1"))
              # rv[["cox_gender"]][["cox"]][["p.adj"]][2] <- rv[["cox_gender"]][["cox"]][["p.adj"]][2]
              # rv[["cox_gender"]][["cox"]][["p.adj"]] <- ifelse(rv[["cox_gender"]][["cox"]][["p.adj"]] > 1, 1, rv[["cox_gender"]][["cox"]][["p.adj"]])
            }
          }
          # p adj on pairwise heatmap
          if(!rv$depmap & rv$km_mul_padj == "padj"){
            if(!is.null(rv[["cox_1"]][["km"]][["p.adj"]])){
              km2 <- rv[["cox_gender"]][["km"]][["stats"]][[1]]
              km2$p.value <- km2$p.value + rv[["cox_gender"]][["km"]][["p.adj"]]
              km2$p.value <- ifelse(km2$p.value > 1, 1, km2$p.value)
              rv[["cox_gender"]][["km"]][["stats"]][[1]] <- km2
            }
          }else if(rv$depmap & rv$km_mul_padj == "padj"){
            if(!is.null(rv[["cox_1"]][["p.adj"]])){
              km2 <- rv[["cox_gender"]][["stats"]][[1]]
              km2$p.value <- km2$p.value + rv[["cox_gender"]][["p.adj"]]
              km2$p.value <- ifelse(km2$p.value > 1, 1, km2$p.value)
              rv[["cox_gender"]][["stats"]][[1]] <- km2
            }
          }
          rv[["title_gender"]] <- paste0(rv[["title_1"]]," vs Gender")
        }else{
          rv[["df_gender"]] <- unique(df_combined[["level.y"]])
        }
      }
    })

    # update parameters
    rv$cox_km <- "km"
    if(rv$tcga){rv$plot_sstype <- rv$plot_stype}
    else if(rv$target){rv$plot_sstype <- "Overall survival (OS)"}
    else if(rv$depmap){rv$plot_sstype <- dependency_names();rv$depmap_gener <- rv$depmap_gene}
    rv$annot_cells_y <- ""
    rv$annot_data_points_y <- ""
    rv$show_ui <- "yes"
    rv$surv_plotted <- ""
    rv$analysis_no <- rv$analysis_no + 1
  }
})
