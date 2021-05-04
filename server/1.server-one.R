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
  
  # check if exceed maximum no of projects
  project_length_check <- length(input$project) > rv$max_project_n
  if(study == "TARGET"){rv$max_project_n <- 4}else{rv$max_project_n <- 3}
  if(project_length_check){
    shinyalert(paste0("We support pan-cancer analysis with up to ",rv$max_project_n," projects in ",study,"."
                      , " You have selected ",length(input$project),"."
                      , " Please delete unrelated projects. Thank you."))
  }
  req(!project_length_check)
  
  # retrieve data
  withProgress(value = 1, message = "Retrieving data from project(s) .... ",{
    rv$project <- input$project
    if(study == "TCGA"){rv$tcga <- T}else{rv$tcga <- F}
    if(study == "TARGET"){rv$target <- T}else{rv$target <- F}
    if(study == "DepMap"){rv$depmap <- T}else{rv$depmap <- F}

    rv$indir <- paste0(getwd(),"/project_data/",project,"/")
    infiles <- paste0(rv$indir,"df_survival.csv")
    l <- lapply(infiles, function(x){
      fread(x,sep=",",header=T) %>%
        dplyr::select(patient_id,survival_days,censoring_status,gender)
    })
    rv$df_survival <- rbindlist(l,use.names = T) %>%
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
                      , "<b>Common dataset(s) for pan-cancer analysis:</b> ", paste(names(rv$overlapped_parameter), collapse = ", "), "."
          ))
        ),
        size = "l", easyClose = T, footer = modalButton("OK")
      ))
    }
  }
})

## reset project
observeEvent(input$reset_project,{
  rv$project <- ""
  rv[["cox_1"]] <- NULL
  clear_rds()
  shinyjs::enable("project")
  rv$projectStatus <- "none"
  
  updateSelectizeInput(
    session,
    "project",
    selected = ""
  )
  
  # update genes placeholder
  lapply(1:rv$variable_n, function(x){
    g_ui_id <- paste0("g_",x)
    updateSelectizeInput(
      session,
      g_ui_id
      ,choices=c()
      ,options = list(
        placeholder = g_placeholder
      )
    )
  })
})

#======================================================================#
####                          STEP 1. parameters                   ####
#======================================================================#
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
  }
  
  if(rv$variable_n_reached==0){
    load()
  }else{
    withProgress(value = 1, message = "Loading parameters ...",{
      load()
      if(rv$project != ""){
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
  withProgress(value = 1, message = "Extracting data from the selected database. Please wait a minute...",{
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
        
        output[[gs_lib_genes_id]] <- renderText({
          paste0("Genes in selected GS (n=",length(genes),"): ", paste0(genes, collapse = " "))
        })
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
    db_id <- paste0("add_btn_",x)
    input[[db_id]]
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
            }
          }
        }
      }
    })
  })
},ignoreInit = T)

# ------- [1F] update loaded genes when user changed to other mutation callers ---------
manual_lst <- reactive({
  lapply(1:rv$variable_n, function(x){
    snv_id <- paste0("snv_method_",x)
    input[[snv_id]]
  })
})

observeEvent(manual_lst(),{
  array <- 1:rv$variable_n
  namespaces <- paste0("snv_method_",array)
  req(req_diff_rv(namespaces))
  withProgress(value = 1, message = "Extracting data from the selected database. Please wait a minute...",{
    lapply(array, function(x) {
      snv_id <- paste0("snv_method_",x)
      if(rv[[snv_id]] != input[[snv_id]]){
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
  req(rv$project != "")
  array <- 1:rv$variable_n
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

# ----- 1.2. run parameters -------
observe({
  rv[["ui_run_parameters"]] <- plot_run_ui(rv$variable_n)
})

# parameters for KM analysis
output$par_gear <- renderUI({
  rv[["ui_run_parameters"]]
})

# update all run parameters according to analysis #1
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
    updateSelectInput(session, paste0("snv_method_",x), selected = input[["snv_method_1"]])
    updateSelectizeInput(session, paste0("nonsynonymous_",x), selected = input[["nonsynonymous_1"]])
    updateSelectizeInput(session, paste0("synonymous_",x), selected = input[["synonymous_1"]])
  })
})

# ----- 1.3. confirm to start analysis -------
# the confirm button
output$ui_parameters_confirm <- renderUI({
  column(
    12,align="center",
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


