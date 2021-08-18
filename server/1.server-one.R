# ---- UI: censor time cutoff -----
output$ui_censortime <- renderUI({
  req(rv$tcga | rv$target)
  if(rv$censor_time_ymd == "none"){t_w <- "0px";ymd_w <- "90%"}else{t_w <- "45%";ymd_w <- "50%"}
  div(
    HTML(paste0("<h4 style='margin-bottom: 15px;'><b>Censor cases at:</b> ",add_help("censor_time_q"),"</h4>")),
    div(
      style=sprintf("display: inline-block;vertical-align:top;width:%s",t_w),
      uiOutput("ui_censor_time")
    ),
    div(
      style=sprintf("display: inline-block;vertical-align:top;width:%s",ymd_w),
      selectizeInput(
        "censor_time_ymd", NULL,
        choices = c(
          "Years" = "y",
          "Months" = "m",
          "Days" = "d",
          "None" = "none"
        )
        ,selected = rv$censor_time_ymd
      )
    )
    ,bsTooltip("censor_time_q",HTML(paste0(
      "To study a specified time interval. Choose <b>None</b> to study the whole time course"
    )),placement = "right")
  )
})

# the no of time units
output$ui_censor_time <- renderUI({
  req(input$censor_time_ymd != "none")

  numericInput(
    "censor_time", NULL,
    value = rv$censor_time
    ,min = rv$censor_time_min,max = rv$censor_time_max,step = rv$censor_time_step
  )
})

# switch time units
observeEvent(input$censor_time_ymd,{
  ymd <- rv$censor_time_ymd <- input$censor_time_ymd
  req(ymd != "none")
  rv$censor_time <- rv[[paste0("censor_time_",ymd)]]
  rv$censor_time_min <- rv[[paste0("censor_time_min_",ymd)]]
  rv$censor_time_max <- rv[[paste0("censor_time_max_",ymd)]]
  rv$censor_time_step <- rv[[paste0("censor_time_step_",ymd)]]
})

# # correct user's numeric time input
# observeEvent(input$censor_time,{
#   if(!is.na(input$censor_time)){
#     # take intervals
#     n <- floor(input$censor_time)
#     # cap at min and max
#     if(n < rv$censor_time_min){n <- rv$censor_time_min}
#     if(n > rv$censor_time_max){n <- rv$censor_time_max}
#     # update rv
#     rv$censor_time <- rv[[paste0("censor_time_",rv$censor_time_ymd)]] <- n
#   }
# })

#=========================================================================#
####  STEP 0. Freeze project once selected, update gene selection UI   ####
#=========================================================================#
# output option to see if a project(s) is (are) selected or not
output$projectStatus <- reactive({
  rv$projectStatus == "selected"
})
outputOptions(output, "projectStatus", suspendWhenHidden = FALSE)

# extract genes when project selection is confirmed
observeEvent(input$confirm_project,{
  rv$depmap_gene_appear <- "no"
  # check if selected projects are from different studies
  project <- input$project
  study <- sapply(project, function(x){
    strsplit(x, "-")[[1]][1]
  }) %>% unique(.)

  study_length_check <- length(study) > 1
  if(study_length_check){
    shinyalert(paste0("You have selected projects from: ",paste0(study, collapse = ", "),"."
                      ," Please select projects from the same program."))
  }
  req(!study_length_check)

  # # check if exceed maximum no of projects
  if(study == "TARGET" | study == "TCGA"){rv$max_project_n <- 1000}else{rv$max_project_n <- 1}
  project_length_check <- length(input$project) > rv$max_project_n
  if(project_length_check){
    shinyalert(paste0("We support analysis with up to ",rv$max_project_n," project(s) in ",study,"."
                      , " You have selected ",length(input$project),"."
                      , " Please delete unrelated project(s). Thank you."))
  }
  req(!project_length_check)

  # retrieve data
  withProgress(value = 1, message = "Retrieving data from project(s) .... ",{
    rv$project <- input$project
    if(study == "TCGA"){rv$tcga <- T; rv$plot_stype <- vector_names(rv$tcga_stype, tcga_stypes)}else{rv$tcga <- F}
    if(study == "TARGET"){rv$target <- T; rv$plot_stype <- "Overall survival (OS)"}else{rv$target <- F}
    if(study == "DepMap"){rv$depmap <- T; rv$plot_stype <- paste0(gsub("^DepMap-","",project)," dependency")}else{rv$depmap <- F}

    if(study != "DepMap"){
      rv$indir <- paste0(pro_dir,project,"/")
      if(rv$tcga){
        infiles <- paste0(rv$indir,"df_survival_o.csv")
        l <- lapply(infiles, function(x){
          fread(x,sep=",",header=T) #%>%
          # dplyr::select(patient_id,person_neoplasm_cancer_status,new_tumor_event_after_initial_treatment,survival_days,censoring_status,gender)
        })
      }else if(rv$target){
        infiles <- paste0(rv$indir,"df_survival.csv")
        l <- lapply(infiles, function(x){
          fread(x,sep=",",header=T) %>%
            dplyr::select(patient_id,survival_days,censoring_status,gender)
        })
      }
      rv$df_survival_o <- rbindlist(l,use.names = T) %>%
        dplyr::distinct(patient_id, .keep_all = T)
      rv[["ui_parameters"]] <- plot_ui(rv$variable_n)
      if(!is.null(rv$overlapped_parameter)){
        lapply(1:rv$variable_n, function(x){
          db_id <- paste0("db_",x)
          if(!rv[[db_id]] %in% rv$overlapped_parameter){
            rv[[db_id]] <- "rna"
          }
        })
        rv[["ui_parameters"]] <- plot_ui(rv$variable_n)
      }
      update_genes_ui(opt="nil")
    }else{
      rv$indir <- paste0(pro_dir,"DepMap/")
      rv$depmap_path <- paste0(rv$indir,project,".csv")
      # the genes for initial selection
      rv$depmap_genes <- fread(rv$depmap_path, sep = ",", nrows = 0, quote="")
      rv$depmap_genes <- colnames(rv$depmap_genes)[-1]
      # available cell lines
      rv$depmap_ids <- fread(rv$depmap_path, sep = ",", select = "patient_id", quote="") %>% .[["patient_id"]]
      rv$depmap_ccle <- df_ccle %>% dplyr::filter(patient_id %in% rv$depmap_ids)
      rv$cell_lines <- rv$depmap_ccle$CCLE_Name
      names(rv$cell_lines) <- rv$depmap_ccle$patient_id
      # primary cancers
      rv$ccle_cancers <- sort(unique(rv$depmap_ccle$primary_disease))
      rv[["ui_parameters"]] <- plot_ui(rv$variable_n)
    }
  })

  rv$projectStatus <- "selected"
  shinyjs::disable("project")

  # create the modal that display the target project information
  if(rv$target){
    target_project_texts <- ""
    if(length(project) > 1 & length(rv$overlapped_parameter) < 4){
      for(j in seq_along(project)){
        target_project_texts <- paste0(target_project_texts, project[j], ": ", paste0(names(name_project_choices(rv$parameters_target_projects[[j]])), collapse = ", "), ".<br><br>")
      }
      # the modal that displayed the information of target projects' data
      showModal(modalDialog(
        title = h2("Available data"),
        div(
          style="font-size:120%",
          HTML(paste0(target_project_texts
                      , "<b>Common dataset(s) kept for analysis:</b> ", paste(names(rv$overlapped_parameter), collapse = ", "), "."
          ))
        ),
        size = "l", easyClose = T, footer = modalButton("OK")
      ))
    }
  }
})

## reset project
clear_loaded_genes <- function(){
  # update genes placeholder
  lapply(1:rv$variable_n, function(x){
    g_ui_id <- paste0("g_",x)
    updateSelectizeInput(
      session,
      g_ui_id
      ,choices=c()
      ,options = list(
        placeholder = g_placeholder()
      )
      ,server = T
    )
  })
}
observeEvent(input$reset_project,{
  rv$project <- ""; rv$tcga <- T; rv$depmap <- F; rv$target <- F; rv$show_ui <- ""
  rv$depmap_path<-NULL; rv$depmap_genes<-NULL; rv$depmap_ids<-NULL; rv$depmap_ccle<-NULL; rv$cell_lines<-NULL
  rv$ccle_cancer_types<-"";rv$ccle_cancer_subtypes="";rv$depmap_gene=""
  rv$annot_cells_y <- "";rv$annot_data_points_y <- ""
  
  rv[["cox_1"]] <- NULL
  clear_rds()
  shinyjs::enable("project")
  rv$projectStatus <- "none"
  
  rv$surv_plotted <- ""

  updateSelectizeInput(
    session,
    "project",
    selected = ""
  )

  clear_loaded_genes()
})

#======================================================================#
####                          STEP 1. parameters                   ####
#======================================================================#
retrieve_genes_total <- function(){
  lapply(1:rv$variable_n, function(x){
    rv[[paste0("genes",x)]] <- retrieve_genes(x)
    g_ui_id <- paste0("g_",x)

    updateSelectizeInput(
      session,
      g_ui_id,
      choices = rv[[paste0("genes",x)]]
      ,selected = rv[[g_ui_id]]
      ,server = TRUE
      ,options = list(
        placeholder = 'Type to search ...'
      )
    )
  })
}
# ----- 1.1. detect and organize user inputs -------
# capping maximum # of analysis at 2
observeEvent(input$variable_n,{
  req(input$variable_n)
  load <- function(){
    n <- floor(input$variable_n); if(n<1){n <- 1}
    if(n > 2){
      shinyalert("We currently support interaction analysis up to two variables.")
      n <- 2
    }

    rv$variable_n <- n

    if(is.null(rv[[paste0("cat_",rv$variable_n)]])){
      init_rv(rv$variable_n)
    }
    updateNumericInput(
      session,
      "variable_n",
      value = rv$variable_n
    )
    update_all()
    rv[["ui_parameters"]] <- plot_ui(rv$variable_n)
    update_normalization_UI()
  }

  if(rv$variable_n_reached==0){
    load()
  }else{
    withProgress(value = 1, message = "Loading parameters ...",{
      load()
      if(grepl("TCGA|TARGET",rv[["project"]][1]) | (grepl("DepMap",rv[["project"]][1]) & rv$depmap_gene != "")){
        retrieve_genes_total()
      }
    })
  }
  rv$variable_n_reached <- rv$variable_n_reached + 1
})

# main panels for user inputs
output$ui_parameters <- renderUI({
  rv[["ui_parameters"]]
})

# -------- [1A] auto GMT loading ---------
gmt_input_lst <- reactive({
  lapply(1:rv$variable_n, function(x){
    db_id <- paste0("gs_db_",x)
    input[[db_id]]
  })
})

observeEvent(gmt_input_lst(),{
  # lst <- gmt_input_lst()
  array <- 1:rv$variable_n #check_array(lst)
  namespaces <- paste0("gs_db_",array)
  req(req_diff_rv(namespaces))
  rv$show_ui <- ""
  withProgress(value = 1, message = wait_msg("Extracting data from the selected database..."),{
    lapply(array, function(x) {
      gs_db_id <- paste0("gs_db_",x)
      if(rv[[gs_db_id]] != input[[gs_db_id]]){
        # extract GS data
        update_gs_by_db(x)

        # update searchInput
        gs_gene_id <- paste0("gs_lg_",x)
        updateSearchInput(
          session,
          gs_gene_id,
          value = ""
        )

        # update verttext
        gs_gene_genes_id <- paste0("gs_lgg_",x)
        output[[gs_gene_genes_id]] <- NULL
      }
    })
  })
}, ignoreInit = T)

# ------- [1B] verbatimText feedback on GMT genes ---------
lib_input_lst <- reactive({
  lapply(1:rv$variable_n, function(x){
    db_id <- paste0("gs_l_",x)
    input[[db_id]]
  })
})

observeEvent(lib_input_lst(),{
  lst <- lib_input_lst()
  array <- 1:rv$variable_n #check_array(lst)
  namespaces <- paste0("gs_l_",array)
  req(req_diff_rv(namespaces))
  lapply(array, function(x){
    gs_lib_id <- paste0("gs_l_",x)
    gs_lib_genes_id <- paste0("gs_lgs_",x)

    if(!is.na(input[[gs_lib_id]])){
      gs <- rv[[gs_lib_id]] <- input[[gs_lib_id]]

      if(!is.null(gs) & gs != ""){
        genes <- rv[[paste0("gmts",x)]][[gs]]
        rv[[paste0("gs_genes_",x)]] <- genes

        # verbatimTextOutput feedback
        output[[gs_lib_genes_id]] <- renderText({
          paste0("Genes in selected GS (n=",length(genes),"): ", paste0(genes, collapse = " "))
        })
        
        # gene list download handler
        gs_lib_dn_id <- paste0(gs_lib_id,"dn")
        output[[gs_lib_dn_id]] <- downloadHandler(
          filename = function(){paste0(gs,".txt")},
          content = function(file) {
            fwrite(as.list(paste0(sort(genes),collapse = "\n")),file,quote = F)
          }
        )
        
        # official site link
        gs_db_id <- paste0("gs_db_",x)
        gs_lib_link_id <- paste0(gs_lib_id,"link"); gs_lib_link_ui_id <- paste0(gs_lib_link_id,"_ui")
        if(input[[gs_db_id]] %in% names(gmt_links)){
          gs_lib_link_id_q <- paste0(gs_lib_link_id,"_q")
          output[[gs_lib_link_ui_id]] <- renderUI({
            div(id=gs_lib_link_id_q,style="z-index:1000;word-break: break-all;",
                style = "position: absolute; right: -3.5em; top: -0.2em;",
                HTML(paste0("<h4>",link_icon(gs_lib_link_id,paste0(gmt_links[[input[[gs_db_id]]]],tail(strsplit(rv[[gs_lib_id]],"%")[[1]],n=1)),color=d_red),"</h4>"))
                ,bsTooltip(gs_lib_link_id_q,HTML(paste0("Click to visit official ",rv[[gs_db_id]]," website describing ",rv[[gs_lib_id]],".")),placement = "top")
            )
          })
        }else{
          output[[gs_lib_link_ui_id]] <- NULL
        }
      }
    }
  })
}, ignoreInit = T)

# ------- [1C] verbatimText feedback on user-supplied gene(s), filter gss ---------
lg_input_btn_lst <- reactive({
  lapply(1:rv$variable_n, function(x){
    db_id <- paste0("gs_lg_",x,"_search")
    input[[db_id]]
  })
})

lg_input_lst <- reactive({
  lapply(1:rv$variable_n, function(x){
    db_id <- paste0("gs_lg_",x)
    input[[db_id]]
  })
})

observeEvent(lg_input_btn_lst(),{
  lst <- lg_input_lst()
  array <- 1:rv$variable_n #check_array(lst)
  namespaces <- paste0("gs_lg_",array,"_search")
  req(req_diff_rv_btn(namespaces))
  withProgress(value = 1, message = "Filtering gene sets contain the entered gene criteria. Please wait a minute ...",{
    lapply(array, function(x){
      # update the search button value
      lgg_btn_id <- paste0("gs_lg_",x,"_search")
      # saveRDS(input[[lgg_btn_id]], file = paste0(getwd(),"/inc/btn0.rds"))
      if(rv[[lgg_btn_id]] < input[[lgg_btn_id]][1]){
        rv[[lgg_btn_id]] <- input[[lgg_btn_id]][1]

        # the user-supplied GS-filtering genes
        lgg_id <- paste0("gs_lg_",x)
        lgg <- rv[[lgg_id]] <- input[[lgg_id]]

        # check if input is empty
        if(lgg == ""){
          #   shinyalert("Please enter a valid gene.")
          update_gs_by_db(x)
        }else{
          # proceed only rv not equal to input

          # check if valid entry
          genes <- toupper(lgg) %>% gsub(" ","",.) %>% str_split(.,"&") %>% .[[1]] %>% unique()
          no_genes <- sapply(genes, function(x){
            str_split(x,"\\|") %>% .[[1]]
          }) %>% unlist(.) %>% length(.)

          if(genes == ""){
            shinyalert("Please enter valid gene(s).")
          }

          if(no_genes > 10){
            shinyalert("We support evaluation up to 10 genes. Please delete less important genes wherever appropriate.")
          }

          if(genes != "" & no_genes < 11){
            # retrive gmt data
            gmts <- rv[[paste0("gmts",x)]]

            filtered_gmts <- lapply(seq_along(gmts), function(i) {
              gmt <- gmts[i]
              gs <- names(gmt)
              x <- str_split(gs,"_")[[1]]
              words_in_name <- unnest_tokens(tibble(x), word, x)[["word"]] %>% toupper(.)
              txt <- c(words_in_name, gmt[[1]])
              count <- 0
              for(gene in genes){
                gene <- str_split(gene,"\\|")[[1]]
                g_check <- sapply(gene, function(x){
                  return(x %in% txt)
                })
                if(any(g_check)){
                  count <- count + 1
                }
              }

              if(!(count < length(genes))){
                return(gmt)
              }
            })

            filtered_gmts[sapply(filtered_gmts, is.null)] <- NULL
            filtered_gmts <- unlist(filtered_gmts, recursive = F)

            if(length(filtered_gmts) == 0){
              shinyalert(sprintf("Analysis #%s: Unable to detect the entered gene combination in the selected database. Try another database or another gene combination, or double check if your input follows the right format.",x))
            }

            gmt_length <-  length(filtered_gmts)
            if(gmt_length > 0){
              rv[[paste0("gmts_tmp",x)]] <- filtered_gmts
              rv[[paste0("gs_lgg_",x)]] <- paste0("Filtered by: ", paste0(genes, collapse = " "))

              gs_db_id <- paste0("gs_l_",x)

              rv[[paste0("gs_placeholder",x)]] <- sprintf('(Filtered n=%s) Type to search ...',gmt_length)

              # update gene set UI
              updateSelectizeInput(
                session,
                gs_db_id
                ,choices = names(filtered_gmts)
                ,selected="" #rv[[gs_db_id]]
                ,options = list(
                  # `live-search` = TRUE,
                  placeholder = rv[[paste0("gs_placeholder",x)]]
                  ,onInitialize = I(sprintf('function() { this.setValue("%s"); }',rv[[gs_db_id]]))
                )
                ,server = T
              )

              output[[paste0("gs_lgg_",x)]] <- renderText({
                rv[[paste0("gs_lgg_",x)]]
              })

              # output[[paste0("gs_lgg_",x,"_tag")]] <- renderUI({
              #   tags$head(tags$style(
              #     paste0("#gs_lgg_",x,"{",rv$verbTxtStyle1,"}")
              #   ))
              # })
            }
          }
        }
      }
    })
  })
}, ignoreInit = T)

# ------- [1D] reset to all gene sets on clear filtering ---------
lg_input_clearbtn_lst <- reactive({
  lapply(1:rv$variable_n, function(x){
    db_id <- paste0("gs_lg_",x,"_reset")
    input[[db_id]]
  })
})

observeEvent(lg_input_clearbtn_lst(),{
  array <- 1:rv$variable_n
  namespaces <- paste0("gs_lg_",array,"_reset")
  req(req_diff_rv_btn(namespaces))
  withProgress(value = 1, message = "Resetting choices to all gene sets in the selected database ...",{
    lapply(array, function(x) {
      lgg_btn_id <- paste0("gs_lg_",x,"_reset")
      if(rv[[lgg_btn_id]] < input[[lgg_btn_id]][1]){
        # update the search button value
        rv[[lgg_btn_id]] <- input[[lgg_btn_id]][1]
        # update the GS
        update_gs_by_db(x, mode="nil")
      }
    })
  })
}, ignoreInit = T)

# ------- [1E] read and process user input gene list ---------
manual_lst <- reactive({
  lapply(1:rv$variable_n, function(x){
    gs_manual_btn_id <- paste0("add_btn_",x)
    input[[gs_manual_btn_id]]
  })
})

observeEvent(manual_lst(),{
  array <- 1:rv$variable_n
  namespaces <- paste0("add_btn_",array)
  req(req_diff_rv_btn(namespaces))
  withProgress(value = 1, message = "Processing input genes ...",{
    lapply(array, function(x){
      add_btn_id <- paste0("add_btn_",x)

      if(rv[[add_btn_id]] < input[[add_btn_id]][1]){
        # update the Submit button value
        rv[[add_btn_id]] <- input[[add_btn_id]][1]
        # read in the gene list
        gs_manual_id <- paste0("gs_m_",x)
        gs_genes_id <- paste0("gs_mg_",x)


        e_msg <- paste0("Please input a valid gene list in Analysis #",x)

        if(input[[gs_manual_id]] == ""){
          shinyalert(e_msg)
        }else{
          genelist <- as.character(input[[gs_manual_id]])
          genelist = gsub("\"","",genelist)
          if(length(genelist) == 1 & genelist==""){shinyalert(e_msg)}else{
            # genelist = unlist(lapply(genelist, function(x) strsplit(x,'\\s*,\\s*')))
            genelist = strsplit(genelist,"\n") %>% unlist(.)
            genelist = unlist(strsplit(genelist,"\\s*,\\s*"))
            genelist = unlist(strsplit(genelist,"\t"))
            genelist = unlist(strsplit(genelist," "))
            genelist = unique(genelist) %>% toupper(.)

            if(length(genelist) == 1 & genelist==""){
              shinyalert(e_msg)
            }else if(length(genelist)>100){
              shinyalert("We currently support analysis of up to 100 genes. Please revise your input. Thank you.")
            }else{
              genelist <- genelist[genelist!=""]
              rv[[gs_manual_id]] <- genelist
              output[[gs_genes_id]] <- renderText({
                paste0("Input genes (n=",length(genelist),"): ", paste0(genelist, collapse = " "))
              })
              output[[paste0("gs_manual_uploaded",x)]] <- reactive({T})
            }
          }
        }
      }
    })
  })
},ignoreInit = T)

# ------- [1F] update loaded genes when user changed to other mutation callers ---------
snv_lst <- reactive({
  lapply(1:rv$variable_n, function(x){
    snv_id <- paste0("snv_method_",x)
    input[[snv_id]]
  })
})

observeEvent(snv_lst(),{
  req(rv$project != "")
  array <- 1:rv$variable_n
  namespaces <- paste0("snv_method_",array)
  req(req_diff_rv(namespaces))
  withProgress(value = 1, message = "Extracting data from the selected database. Please wait a minute...",{
    lapply(array, function(x) {
      snv_id <- paste0("snv_method_",x)
      snv_diff <- 0
      if(length(rv[[snv_id]]) != length(input[[snv_id]])){
        snv_diff <- 1
      }else if(rv[[snv_id]] != input[[snv_id]]){
        snv_diff <- 1
      }

      if(snv_diff == 1){
        # update SNV calling method stored in RV
        rv[[snv_id]] <- input[[snv_id]]

        # update loaded genes
        update_genes_ui()
      }
    })
  })
},ignoreInit = T)

# ---------- [1G] update loaded genes ---------
## update gene selection UI
genes_lst <- reactive({
  lapply(1:rv$variable_n, function(x){
    input[[paste0("db_",x)]]
  })
})

observeEvent(genes_lst(),{
  array <- 1:rv$variable_n
  # update reset buttons
  lapply(1:array, function(x){rv[[paste0("todefault",x)]] <- 0})
  req(rv$project != "")
  namespaces <- paste0("db_",array)
  req(req_diff_rv(namespaces))
  lapply(array, function(x){
    db_id <- paste0("db_",x)
    if(rv[[db_id]] != input[[db_id]]){
      rv[[db_id]] <- input[[db_id]]
      rv[[paste0("genes",x)]] <- retrieve_genes(x)
      updateSelectizeInput(
        session,
        paste0("g_",x),
        choices = rv[[paste0("genes",x)]]
        ,selected = ""
        ,server = TRUE
        ,options = list(
          placeholder = 'Type to search ...'
        )
      )
    }
  })
},ignoreInit = T)

# ---------- [1H-A] update normalization gene ---------
## update normalization gene selection UI
normg_lst <- reactive({
  lapply(1:rv$variable_n, function(x){
    input[[paste0("gnorm_",x)]]
  })
})

observeEvent(normg_lst(),{
  req(rv$project != "")
  array <- 1:rv$variable_n
  namespaces <- paste0("gnorm_",array)
  req(req_diff_rv(namespaces))
  lapply(array, function(x){
    normg_id <- paste0("gnorm_",x)
    if(rv[[normg_id]] != input[[normg_id]]){
      rv[[normg_id]] <- input[[normg_id]]
      g_ui_id <- paste0("gnorm_g_",x)
      if(rv[[normg_id]] == "g"){
        updateSelectizeInput(
          session,
          g_ui_id,
          choices = rv[[paste0("genes",x)]]
          ,selected = rv[[g_ui_id]]
          ,server = TRUE
          ,options = list(
            placeholder = 'Type to search ...'
          )
        )
      }else{
        updateRV(g_ui_id)
      }
    }
  })
})

update_normalization_UI <- function(){
  lapply(1:rv$variable_n, function(x){
    if(!is.null(rv[[paste0("genes",x)]])){
      g_ui_id <- paste0("gnorm_g_",x)
      updateSelectizeInput(
        session,
        g_ui_id,
        choices = rv[[paste0("genes",x)]]
        ,selected = rv[[g_ui_id]]
        ,server = TRUE
        ,options = list(
          placeholder = 'Type to search ...'
        )
      )
    }
  })
}

# ---------- [1H-B] update normalization GS ---------
normgs_lst <- reactive({
  lapply(1:rv$variable_n, function(x){
    db_id <- paste0("gnorm_gs_db_",x)
    input[[db_id]]
  })
})

observeEvent(normgs_lst(),{
  array <- 1:rv$variable_n #check_array(lst)
  namespaces <- paste0("gnorm_gs_db_",array)
  req(req_diff_rv(namespaces))
  rv$show_ui <- ""
  withProgress(value = 1, message = wait_msg("Extracting data from the selected database..."),{
    lapply(array, function(x) {
      gnorm_gs_db_id <- paste0("gnorm_gs_db_",x)
      if(rv[[gnorm_gs_db_id]] != input[[gnorm_gs_db_id]]){
        # extract GS data
        update_gs_by_db(x,gs_db_id=gnorm_gs_db_id,gs_lib_id = paste0("gnorm_gs_lib_",x), prefix="gnorm_")
      }
    })
  })
}, ignoreInit = T)

# ------- [1H-C] verbatimText feedback on normalization GMT genes ---------
gnorm_lib_input_lst <- reactive({
  lapply(1:rv$variable_n, function(x){
    db_id <- paste0("gnorm_gs_lib_",x)
    input[[db_id]]
  })
})

observeEvent(gnorm_lib_input_lst(),{
  array <- 1:rv$variable_n
  namespaces <- paste0("gnorm_gs_lib_",array)
  req(req_diff_rv(namespaces))
  lapply(array, function(x){
    gs_lib_id <- paste0("gnorm_gs_lib_",x)
    gs_lib_genes_id <- paste0("gnorm_gs_lib_genes_",x)
    
    if(!is.na(input[[gs_lib_id]])){
      gs <- rv[[gs_lib_id]] <- input[[gs_lib_id]]
      
      if(!is.null(gs) & gs != ""){
        genes <- rv[[paste0("gnorm_gmts",x)]][[gs]]
        rv[[paste0("gnorm_gs_genes_",x)]] <- genes
        
        # verbatimTextOutput feedback
        output[[gs_lib_genes_id]] <- renderText({
          paste0("Genes in selected GS (n=",length(genes),"): ", paste0(genes, collapse = " "))
        })

        # gene list download handler
        gs_lib_dn_id <- paste0(gs_lib_id,"dn")
        output[[gs_lib_dn_id]] <- downloadHandler(
          filename = function(){paste0(gs,".txt")},
          content = function(file) {
            fwrite(as.list(paste0(sort(genes),collapse = "\n")),file,quote = F)
          }
        )

        # official site link
        gs_db_id <- paste0("gnorm_gs_db_",x)
        gs_lib_link_id <- paste0(gs_lib_id,"link"); gs_lib_link_ui_id <- paste0(gs_lib_link_id,"_ui")
        if(input[[gs_db_id]] %in% names(gmt_links)){
          gs_lib_link_id_q <- paste0(gs_lib_link_id,"_q")
          output[[gs_lib_link_ui_id]] <- renderUI({
            div(id=gs_lib_link_id_q,style="z-index:1000;word-break: break-all;",
                style = "position: absolute; right: -3.5em; top: -0.2em;",
                HTML(paste0("<h4>",link_icon(gs_lib_link_id,paste0(gmt_links[[input[[gs_db_id]]]],tail(strsplit(rv[[gs_lib_id]],"%")[[1]],n=1)),color=d_red),"</h4>"))
                ,bsTooltip(gs_lib_link_id_q,HTML(paste0("Click to visit official ",rv[[gs_db_id]]," website describing ",rv[[gs_lib_id]],".")),placement = "top")
            )
          })
        }else{
          output[[gs_lib_link_ui_id]] <- NULL
        }
      }
    }
  })
}, ignoreInit = T)

# ----- 1.2. RUN PARAMETERS -------
observe({
  rv[["ui_run_parameters"]] <- plot_run_ui(rv$variable_n)
})

# parameters for KM analysis
output$par_gear <- renderUI({
  if(rv$variable_n > 1){
    db_types <- lapply(1:rv$variable_n, function(x){
      cat_id <- paste0("cat_",x); db_id <- paste0("db_",x)
      cat_name <- ifelse_rv(cat_id)
      x <- ifelse(cat_name=="gs","gs",ifelse_rv(db_id))
      if(x == "cnv"){
        if(rv$depmap){
          risk_hl
        }else{
          c("Gain/Loss","Other")
        }
      }else{
        risk_gp[[x]]
      }
    })
    gps <- apply(expand.grid(db_types[[1]], db_types[[2]]), 1, paste, collapse="_")
    rv$risk_gps <- c("All",gps)
  }
  if_exp <- if_any_exp(); if_iter <- if_any_iter()
  if(if_exp){c_w <- 6;if(if_iter & rv$variable_n > 1){min_per_w <- 4}else{min_per_w <- 12}}else{c_w <- 12}
  if(rv$tcga){c_w_kc <- 6}else{c_w_kc <- 12}
  div(
    if(rv$tcga | if_exp){
      column(
        12,align="center",
        wellPanel(
          style=sprintf("background-color:%s; border: .5px solid #fff;",bcols[1]),
          fluidRow(
            if(rv$tcga){
              column(
                c_w,
                radioGroupButtons(
                  inputId = "flagged",
                  label = HTML(paste0("<h4>If TCGA, exclude flagged cases from analysis? ",add_help("flagged_q"),"</h4>")),
                  choices = c("Yes"="y","No"="n"),
                  selected = rv$flagged,
                  status = "danger", #size = "sm",
                  checkIcon = list(
                    yes = icon("check-square"),
                    no = icon("square-o")
                  )
                )
              )
            }
            ,if(if_exp){
              column(
                c_w_kc,
                radioGroupButtons(
                  "min_p_kc",
                  HTML(paste0("<h4>Method to determine <i>P</i><sub>min</sub>:",add_help("min_p_kc_q"),"</h4>")),
                  choices = c("KM log-rank"="km","Cox PH model"="cox"),
                  selected = rv$min_p_kc,
                  status = "danger", #size = "sm",
                  checkIcon = list(
                    yes = icon("check-square"),
                    no = icon("square-o")
                  )
                )
              )
            },
            if(rv$variable_n > 1){
              div(
                # '(input.cat_1 == "g" & input.db_1 != "snv") | (input.cat_2 == "g" & input.db_2 != "snv")',
                column(
                  12,
                  radioGroupButtons(
                    "risk_gp",
                    HTML(paste0("<h4>Risk subgroup(s) of interest:",add_help("risk_gp_q"),"</h4>")),
                    choices = rv$risk_gps,
                    selected = rv$risk_gp,
                    status = "danger", #size = "sm",
                    checkIcon = list(
                      yes = icon("check-square"),
                      no = icon("square-o")
                    )
                  )
                )
              )
            },
            if(if_exp & if_iter){
              column(
                min_per_w,
                numericInput(
                  "n_perm",
                  HTML(paste0("<h4># of permutations:",add_help("n_perm_q"),"</h4>")),
                  value = rv$n_perm,
                  min = 10, max = 1000, step = 1, width = "180px"
                )
                ,conditionalPanel(
                  'input.n_perm > 100',
                  p(style="color:red;","Increased permutations improves statistical precision. However, longer run time may be expected.")
                )
              )
            },
            if(if_exp & rv$variable_n > 1){
              div(
                column(
                  min_per_w,
                  numericInput(
                    "min_gp_size",
                    HTML(paste0("<h4>Min %:",add_help("min_gp_size_q"),"</h4>")),
                    value = rv$min_gp_size,
                    min = 1, max = 90, step = 1, width = "180px"
                  )
                )
                ,column(
                  min_per_w,
                  selectInput(
                    "search_mode",
                    HTML(paste0("<h4>Search mode:",add_help("search_mode_q"),"</h4>")),
                    choices = c(
                      "Median-anchored greedy" = "heuristic",
                      "Exhaustive" = "exhaustive"
                    ),
                    selected = rv$search_mode
                  )
                  ,conditionalPanel(
                    'input.search_mode == "exhaustive"',
                    p(style="color:red;",HTML("<b>Exhaustive search</b> may need more run time and inflate type II (false negative) errors."))
                  )
                )
              )
            }
          )
          ,bsTooltip("risk_gp_q",HTML("In bivariate outcomes analysis, select the risk group (subgroup) of interest. If <b>All</b> is selected, an ANOVA-like test is done to test if there is significant difference between any of the four subgroups. If a specific subgroup is selected, it is compared against the other subgroups as a whole."),placement = "top")
          ,radioTooltip(id = "risk_gp", choice = "All", title = HTML("ANOVA-like test on all subgroups"))
          ,bsTooltip("min_gp_size_q",HTML("In bivariate outcomes analysis, select the minimum number of cases (default, 10% of the population) to define a subgroup. Only applicable for dynamic iteration on continuous variables (e.g. mRNA gene expression, DNA methylation)."),placement = "top")
          ,bsTooltip("flagged_q",HTML(flagged_exp),placement = "top")
          ,radioTooltip(id = "flagged", choice = "y", title = HTML("Remove flagged cases"))
          ,radioTooltip(id = "flagged", choice = "n", title = HTML("Keep all cases"))
          ,bsTooltip("min_p_kc_q",HTML(paste0(
            "Select the test to determine the optimal cutoff point using the minimum <i>P</i>-value method. "
            ,"Only applicable for continuous variables (e.g. mRNA gene expression, DNA methylation)."))
            ,placement = "top")
          ,radioTooltip(id = "min_p_kc", choice = "km", title = HTML("Kaplan-Meier (KM) log-rank test"))
          ,radioTooltip(id = "min_p_kc", choice = "cox", title = HTML("Cox proportional-hazard (PH) model likelihood ratio test"))
          ,bsTooltip("n_perm_q",HTML("Number of permutations to perform for adjustment on multiple testing arising from assessing a sequence of candidate thresholds with the minimum <i>P</i>-value method in dynamic iteration."),placement = "top")
          ,bsTooltip("search_mode_q",
                     HTML(paste0(
                       "<b>Median-anchored greedy search</b> (heuristic) determines the minimum <i>P</i>-value by finding the percentile in variable 2 that gives the minimum <i>P</i>-value on the median percentile in variable 1, then looking for percentile combinations that give lower <i>P</i>-values via greedy search."
                       ,"<br><br><b>Exhaustive search</b> determines the minimum <i>P</i>-value by testing all percentile combinations."
                       )),placement = "top")
        )
      )
    },
    rv[["ui_run_parameters"]]
  )
})

observeEvent(input$flagged,{rv$flagged <- input$flagged})
observeEvent(input$min_p_kc,{rv$min_p_kc <- input$min_p_kc})

# ------- 1.2a update all run parameters according to analysis #1 --------
observeEvent(input$toall,{
  lapply(2:rv$variable_n, function(x){
    updateRadioGroupButtons(session, paste0("iter_",x), selected = input[["iter_1"]])
    if(input[["iter_1"]] == "iter"){
      updateSliderTextInput(session, paste0("lower_",x), selected = input[["lower_1"]])
      updateSliderTextInput(session, paste0("upper_",x), selected = input[["upper_1"]])
      updateSliderTextInput(session, paste0("step_",x), selected = input[["step_1"]])
    }else{
      updateSliderInput(session, paste0("clow_",x), value = input[["clow_1"]])
    }

  })
})

observeEvent(input$toall_m,{
  lapply(2:rv$variable_n, function(x){
    # updateSelectInput(session, paste0("snv_method_",x), selected = input[["snv_method_1"]])
    # updateRadioGroupButtons(session, paste0("snv_uni_",x), selected = input[["snv_uni_1"]])
    updateSelectizeInput(session, paste0("nonsynonymous_",x), selected = input[["nonsynonymous_1"]])
    updateSelectizeInput(session, paste0("synonymous_",x), selected = input[["synonymous_1"]])
  })
})

# ------- 1.2b reset to default run parameters --------
todefault_lst <- reactive({
  lapply(1:rv$variable_n, function(x){
    db_id <- paste0("todefault",x)
    input[[db_id]]
  })
})

observeEvent(todefault_lst(),{
  array <- 1:rv$variable_n
  namespaces <- paste0("todefault",array)
  req(req_diff_rv_btn(namespaces))
  withProgress(value = 1, message = "Setting to default parameters...",{
    lapply(1:rv$variable_n, function(x){
      todefault_id <- paste0("todefault",x)
      if(rv[[todefault_id]] != input[[todefault_id]][1]){
        rv[[todefault_id]] <- input[[todefault_id]][1]
        cat_id <- paste0("cat_",x); db_id <- paste0("db_",x)
        
        if((input[[cat_id]] == "gs") | (!rv$depmap & input[[cat_id]] == "g" & input[[db_id]] != "snv" & input[[db_id]] != "cnv") | (rv$depmap & input[[cat_id]] == "g" & input[[db_id]] != "snv")){
          updateRadioGroupButtons(session, paste0("iter_",x), selected = "iter")
          updateSliderInput(session, paste0("lower_",x), value = dmin)
          updateSliderInput(session, paste0("upper_",x), value = dmax)
          updateSliderInput(session, paste0("step_",x), value = dstep)
        }else if(input[[cat_id]] == "g" & input[[db_id]] == "snv"){
          updateSelectizeInput(session, paste0("nonsynonymous_",x), selected = variant_types_non)
          updateSelectizeInput(session, paste0("synonymous_",x), selected = variant_types_syn)
        }else if(!rv$depmap & input[[cat_id]] == "g" & input[[db_id]] == "cnv"){
          updateRadioGroupButtons(session, paste0("cnv_par_",x), selected = "auto")
        }
      }
    })
  })
})

# ----- 1.3. TCGA OS, DFI, PFI, DSS -------
output$tcga_pars <- renderUI({
  req(rv$tcga)

  # check if certain outcomes available for selected TCGA project
  if(rv$project != ""){
    code <- tcga_codes %>% dplyr::filter(Study %in% rv$project)
    cords <- code %>% dplyr::select(-Study,-Note)
    if(length(rv$project) > 1){
      cords <- lapply(cords, function(x) if(anyNA(x)){NA}else{"hi"})
      rv$tcga_code <- code
    }
    cords <- codes <- unlist(cords)
    if(length(rv$project) == 1){code <- rv$tcga_code <- code %>% dplyr::select(-Study)}
    # if any NA, render a warning msg on UI
    if(anyNA(cords)){rv$tcga_msg <- code$Note}else{rv$tcga_msg <- ""}
    cords <- c(1:4)[!is.na(cords)]
  }else{
    cords <- c(1:4); rv$tcga_code <- ""; rv$tcga_msg <- ""
  }
  tcga_stypess <- tcga_stypes[cords]
  if(!rv$tcga_stype %in% tcga_stypess){rv$tcga_stype <- tcga_stypess[1]}

  # choice names
  if(rv$project != ""){
    if(length(rv$project) == 1){
      tcga_stypess_names <- sapply(cords, function(x){
        y <- codes[[x]]
        c <- codes_color[[y]]
        i <- codes_icon[[y]]
        return(paste0(names(tcga_stypes)[[x]],span(style=sprintf("color:%s;",c),i)))
      }) %>% unname(.)
    }else{
      tcga_stypess_names <- sapply(cords, function(x){
        return(names(tcga_stypes)[[x]])
      }) %>% unname(.)
    }
  }else{
    tcga_stypess_names <- names(tcga_stypess)
  }

  # render the UI
  column(
    12, id="div_tcga_stype",align="center",
    radioGroupButtons(
      inputId = "tcga_stype",
      label = label_with_help_bttn("If TCGA, select the endpoint to measure survival outcomes: ","tcga_stype_q"),
      choiceNames = tcga_stypess_names,
      choiceValues = unname(tcga_stypess),
      selected = rv$tcga_stype,
      size = "sm",
      checkIcon = list(
        yes = icon("check-square"),
        no = icon("square-o")
      ),
      # status = "primary",
      direction = "horizontal"
    )
    ,uiOutput("tcga_warning")
    ,bsTooltip("tcga_stype_q",HTML(paste0(
      "Curated clinical and survival outcome data by pan-cancer clinical study (Liu et al., <i>Cell</i>, 2018). <b>Click</b> for full assessment table and recommendations."
      # ," To last known disease status, <b>OS</b> assesses all cases"
      # ,"; <b>PFI</b> assesses cases without tumor recurrence"
      # ,"; <b>DFI</b> assesses cases without detectable tumors"
      # ,"; <b>PSS</b> assesses cases with tumor recurrence"
      # ,"; <b>DSS</b> assesses cases with histological evidence of cancer"
      # ,"."
    )), placement = "top")
    ,radioTooltip(id = "tcga_stype", choice = "OS", title = HTML("Assesses all cases: the duration from the time of initial pathological diagnosis till the time of death or loss of followup"))
    ,radioTooltip(id = "tcga_stype", choice = "DSS", title = HTML("Assesses cases with histological evidence of cancer"))
    ,radioTooltip(id = "tcga_stype", choice = "DFI", title = HTML("The length of time after primary treatment for a cancer ends that the patient survives without any signs or symptoms of that cancer"))
    # ,radioTooltip(id = "tcga_stype", choice = "pss", title = HTML("Focus on recurrence cases (with new tumor event) status to last known disease status"))
    ,radioTooltip(id = "tcga_stype", choice = "PFI", title = HTML("The length of time during and after the treatment of cancer, that a patient lives with the disease but it does not get worse"))
  )
})

observeEvent(input$tcga_stype,{rv$tcga_stype <- input$tcga_stype; rv$plot_stype <- vector_names(rv$tcga_stype, tcga_stypes)})

# warning messages about recommendations by Liu 2018
output$tcga_warning <- renderUI({
  req(typeof(rv$tcga_code) == "list")
  if(length(rv$project) == 1){
    # recommendations by Liu, Cell, 2018
    code <- rv[["tcga_code"]][[rv$tcga_stype]]
    
    # beginning wording
    green <- 0
    if(code == "yes"){
      s_msg <- " is recommended"; green <- 1
    }else if(code == "acc"){
      s_msg <- " is accurate and recommended"; green <- 1
    }else if(code == "app"){
      s_msg <- " is approximate and recommended"; green <- 1
    }else if(code == "app|caution"){
      s_msg <- " is approximate and should be <b>used with caution</b>"
    }else if(code == "no"){
      s_msg <- " is <b>not recommended</b>"
    }else if(code == "caution"){
      s_msg <- " should be <b>used with caution</b>"
    }
    
    # assemble msg
    msg <- paste0(rv$plot_stype,s_msg," for ",rv$project)
    if(green == 0){
      msg <- paste0(msg,". Reason: ",rv[["tcga_code"]][["Note"]])
    }else{
      if(rv$tcga_msg != ""){
        msg <- paste0(msg,". Note: ",rv$project," ",rv$tcga_msg)
      }
    }
    
    # render msg
    div(
      p(
        style=sprintf("color:%s;", codes_color[[code]]),
        HTML(paste0(codes_icon[[code]],msg)) #," ",link_icon("liu_link","https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6066282/",color="#363B48",icon="fas fa-book")
      )
      # ,bsTooltip("liu_link",HTML(paste0(
      #   "Assessment by Liu et al, <i>Cell</i>, 2018. Click for full text at PubMed Central."
      # )),placement = "top")
    )
  }else{
    div(
      p(
        style=sprintf("color:%s;", "#939597"),
        HTML(paste0("Reminder: Commonly available clinical endpoints and molecular data are extracted and provided. Click the blue book button above to check recommended use of endpoints for each study."))
      )
    )
    # code_table <- rv$tcga_code
    # code_table[is.na(code_table)] <- "N/A"
    # div(
    #   style = "background:#939597;color:white;",
    #   HTML(paste0("<span style='color:#fff;'><b>Recommended use of clinical endpoints by study:</b></span> ",add_help("liu_rec_q",color = "#F5DF4D"))),
    #   renderTable({code_table},escape = FALSE, spacing = "xs", align = "c", na = "NA")
    # )
  }
})

# ----- 1.b. TCGA OS, DFI, PFI, DSS modal -------
observeEvent(input$tcga_stype_q,{
  tcga_table <- tcga_codes
  tcga_table[["Study"]] <- gsub("^TCGA-","",tcga_table[["Study"]])
  showModal(
    modalDialog(
      title = HTML(paste0("<h3>Recommended use of the endpoints of OS, PFI, DFI, and DSS "
                          ,add_help("liu_rec_q",size="extra-large")
                          ,"</h3><br><div style='font-size:100%'>by <b><a href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6066282/' target='_blank'>TCGA pan-cancer clinical study (Liu et al., <i>Cell</i>, 2018)</a></b></div>")),
      div(
        style="font-size:100%",
        renderTable({
          tcga_table
        },escape = FALSE, spacing = "xs", align = "l", na = "NA", hover = T, striped = F)
        ,bsTooltip("liu_rec_q",HTML(paste0(
          "<b>yes</b> = Recommended;"
          ,"<br><b>no</b> = Not recommended;"
          ,"<br><b>caution</b> = Use with caution;"
          ,"<br><b>acc</b> = Accurate;"
          ,"<br><b>app</b> = Approximate;"
          ,"<br><b>NA</b> = Unavailable."
        )),placement = "bottom")
      ),
      size = "l", easyClose = T, footer = modalButton("OK")
    )
  )
})

# ----- 1.4a. DepMap UI ------
output$depmap_pars <- renderUI({
  req(rv$depmap)
  req(!is.null(rv$depmap_ccle))
  
  div(
    column(
      3,
      pickerInput(
        "ccle_cancer_types",
        HTML(paste0("i. Select cancer type(s) to study:",add_help("ccle_cancer_types_q"))),
        choices = rv$ccle_cancers
        ,selected = rv[["ccle_cancer_types"]]
        ,options = list(
          `actions-box` = TRUE,
          size = 10,
          style = "btn-default",
          `selected-text-format` = "count > 2"
          ,`live-search` = TRUE
        ),
        multiple = TRUE
      )
      ,bsTooltip("ccle_cancer_types_q",HTML(paste0(
        "If DepMap, select a cancer type of interest. Multiple selections are allowed."
      )),placement = "top")
    )
    ,uiOutput("ui_ccle_subtypes")
  )
})

# cancer subtypes
observeEvent(input$ccle_cancer_types,{req(input$ccle_cancer_types!=""); rv$depmap_gene_appear <- "no"; rv$ccle_cancer_types <- input$ccle_cancer_types})
output$ui_ccle_subtypes <- renderUI({
  req(!is.null(input$ccle_cancer_types))
  req(!is.null(rv$depmap_ccle))
  subtypes <- rv$depmap_ccle %>% dplyr::filter(primary_disease %in% rv$ccle_cancer_types) %>%
    .[["Subtype"]] %>% unique() %>% sort()
  div(
    column(
      3,
      pickerInput(
        "ccle_cancer_subtypes",
        HTML(paste0("ii. (Optional) subtypes to study:",add_help("ccle_cancer_subtypes_q"))),
        choices = subtypes
        ,selected = subtypes
        ,options = list(
          `actions-box` = TRUE,
          size = 10,
          style = "btn-default",
          `selected-text-format` = "count > 2"
          ,`live-search` = TRUE
        ),
        multiple = TRUE
      )
      ,bsTooltip("ccle_cancer_subtypes_q",HTML(paste0(
        "Select cancer subtypes. Default: all. Multiple selections are allowed."
      )),placement = "top")
    )
    ,uiOutput("ui_cells")
  )
})

# CCLE cell lines
observeEvent(input$ccle_cancer_subtypes,{req(input$ccle_cancer_subtypes!=""); rv$depmap_gene_appear <- "no"; rv$ccle_cancer_subtypes <- input$ccle_cancer_subtypes})
output$ui_cells <- renderUI({
  req(input$ccle_cancer_subtypes != "")
  req(!is.null(rv$depmap_ccle))
  cells <- rv$depmap_ccle %>% dplyr::filter(primary_disease %in% rv$ccle_cancer_types)
  # execute the following only when newly clicked
  if(length(rv$ccle_cancer_subtypes)==1){
    if(rv$ccle_cancer_subtypes != ""){
      cells <- cells %>%
        dplyr::filter(Subtype %in% rv$ccle_cancer_subtypes)
    }
  }
  cells_names <- cells[["CCLE_Name"]]
  cells <- cells[["patient_id"]]
  names(cells) <- cells_names

  div(
    column(
      3,
      pickerInput(
        "ccle_cells",
        HTML(paste0("iii. (Optional) cell lines:",add_help("ccle_cells_q"))),
        choices = cells
        ,selected = cells
        ,options = list(
          `actions-box` = TRUE,
          size = 10,
          style = "btn-default",
          `selected-text-format` = "count > 1"
          ,`live-search` = TRUE
        ),
        multiple = TRUE
      )
      ,bsTooltip("ccle_cells_q",HTML(paste0(
        "Select cell lines to analyze. Click confirmation button to the right to load data."
      )),placement = "top")
    )
    ,uiOutput("ui_ccle_gene")
  )
})

# DepMap gene/drug selection
output$ui_ccle_gene <- renderUI({
  req(input$ccle_cells != "")
  rv$depmap_gene_appear <- "yes"
  div(
    column(
      3,
      selectizeInput(
        "depmap_gene",
        HTML(paste0("iv. Select a ",agene(),":",add_help("depmap_gene_q"))),
        choices = c()
        ,width = "100%"
        ,options = list(
          `live-search` = TRUE,
          placeholder = "Type to search ..."
          ,onInitialize = I(sprintf('function() { this.setValue(%s); }',""))
        )
      )
      ,bsTooltip("depmap_gene_q",HTML(paste0("To start, select ",depmap_gene_help())),placement = "top")
    )
    # ,conditionalPanel(
    #   'input.depmap_gene != ""',
    #   column(
    #     2,
    #     actionBttn(
    #       "confirm_ccle",strong("Load data!")
    #       ,block = T,style = "simple",color = "warning",size="sm"
    #     )
    #     ,tags$style(type='text/css', "#confirm_ccle { margin-top: 24.2px; height: 33.5px;}")
    #   )
    # )
  )
})

# update available DepMap genes/drugs for selection
observeEvent(list(rv$depmap_gene_appear,input$ccle_cells),{
  req(rv$depmap_gene_appear == "yes")
  updateSelectizeInput(
    session,"depmap_gene",choices = rv$depmap_genes, server = T, selected = rv$depmap_gene
  )
})

# ----- 1.4b. DepMap data processing ------
observeEvent(input$depmap_gene,{
  rv$depmap_gene <- input$depmap_gene
  if(input$depmap_gene == ""){clear_loaded_genes()}
  req(input$depmap_gene != "")

  withProgress(value = 1, message = "Loading data ...",{
    # error if too few cell lines
    error <- 0
    if(length(input$ccle_cells) < 10){
      shinyalert(paste0("Please select at least 10 cell lines to proceed."))
      error <- 1
    }
    req(error == 0)

    # check if input gene has enough data
    error <- 0
    gene <- input$depmap_gene
    df_gene_scale <- fread(rv$depmap_path, sep = ",", select = c("patient_id",gene)) %>%
      dplyr::filter(patient_id %in% input$ccle_cells)
    # df_gene_scale[[gene]] <- (df_gene_scale[[gene]] - mean(df_gene_scale[[gene]])) / sd(df_gene_scale[[gene]])
    # df_gene_scale <- df_gene_scale[!is.na(df_gene_scale[[gene]]),]
    # update depmap RVs
    df_survival_o <- left_join(df_gene_scale, rv$depmap_ccle, by="patient_id") %>%
      dplyr::select(-c(CCLE_Name, primary_or_metastasis, primary_disease, Subtype))
    colnames(df_survival_o) <- c("patient_id","dependency","gender")
    df_survival_o <- df_survival_o %>% dplyr::filter(!is.na(dependency))
    # if(rv$project == "DepMap-Drug"){
    #   df_survival_o[["dependency"]] <- 2^df_survival_o[["dependency"]]
    # }else{
    #   df_survival_o[["dependency"]] <- 10^df_survival_o[["dependency"]]
    # }
    genes_len <- nrow(df_survival_o)
    if(genes_len < 10){
      shinyalert(paste0(gene," only has ",genes_len," data points in the selected cell lines."
                        ," At least 10 are needed. Please select another ",agene()," or adjust cell line choices"))
      error <- 1
    }
    req(error == 0)

    # retrieves molecular data
    retrieve_genes_total()

    rv$df_survival_o <- df_survival_o
    rv$ccle_cells <- input$ccle_cells
    
    # update normalization genes' UI
    update_normalization_UI()
  })
})


# ----- 1.5. confirm to start analysis -------
# the confirm button
output$ui_parameters_confirm <- renderUI({
  column(
    12,align="center",
    add_gear("par_gear"),
    bsButton(
      "confirm",
      strong("Confirm and analyze!")
      # ,block = T
      # ,style = "material-flat"
      ,style = "primary"
      ,size = "large"
    )
  )
})

# ----- NN. DEMO -------
output$ui_demo_header <- renderUI({
  req(rv$demo == "yes")
  column(
    12,align="center",
    wellPanel(
      style = sprintf("background:%s;color:%s;",addalpha(yellow2021,alpha=.5),"black"),
     HTML(paste0("<p><b>You are in demo mode. It will take ~ 10-20 seconds to load example data and complete demo analysis. Then, explore the sample results as the real output.</b></p>"))
     ,HTML(paste0("<u><b><a style='color:",red2021,";font-size:120%;' href='https://tau.cmmt.ubc.ca/cSurvival/'><i class='fas fa-sign-out-alt'></i> Exit Demo Mode</a></b></u>"))
    )
  )
})

