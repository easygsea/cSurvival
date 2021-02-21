#======================================================================#
####                          STEP 1. parameters                   ####
#======================================================================#
# ----- 1.1. detect and organize user inputs -------
# capping maximum # of analysis at 2
observeEvent(input$variable_n,{
  req(input$variable_n)
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
  lst <- gmt_input_lst()
  array <- check_array(lst)
  namespaces <- paste0("gs_db_",array)
  req(req_diff_rv(namespaces))
  withProgress(value = 1, message = "Extracting data from the selected database. Please wait a minute...",{
    lapply(array, function(x) {
      gs_db_id <- paste0("gs_db_",x)
      gs_lib_id <- paste0("gs_l_",x)
      
      db <- rv[[gs_db_id]] <- isolate(input[[gs_db_id]])
      
      req(db != "")
      gmt_path <- retrieve_gmt_path(db)
      gmt <- gmtPathways(gmt_path)
      
      # update variables
      rv[[paste0("gmts",x)]] <- rv[[paste0("gmts_tmp",x)]] <- gmt
      
      # update gene set UI
      updateSelectInput(
        session,
        gs_lib_id
        ,choices = names(gmt)
        ,selected=rv[[gs_lib_id]]
      )
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
  array <- check_array(lst)
  namespaces <- paste0("gs_l_",array)
  req(req_diff_rv(namespaces))
  lapply(array, function(x){
    gs_lib_id <- paste0("gs_l_",x)
    gs_lib_genes_id <- paste0("gs_lgs_",x)
    
    req(input[[gs_lib_id]])
    gs <- rv[[gs_lib_id]] <- input[[gs_lib_id]]
    
    req(!is.null(gs) & gs != "")
    genes <- rv[[paste0("gmts",x)]][[gs]]
    
    output[[gs_lib_genes_id]] <- renderText({
      paste0("Genes in selected GS (n=",length(genes),"): ", paste0(genes, collapse = " "))
    })
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
  array <- check_array(lst)
  namespaces <- paste0("gs_lg_",array,"_search")
  req(req_diff_rv_btn(namespaces))
  withProgress(value = 1, message = "Filtering gene sets contain the entered gene criteria. Please wait a minute ...",{
    lapply(array, function(x){
      # update the search button value
      lgg_btn_id <- paste0("gs_lg_",x,"_search")
      # saveRDS(input[[lgg_btn_id]], file = paste0(getwd(),"/inc/btn0.rds"))
      rv[[lgg_btn_id]] <- input[[lgg_btn_id]][1]

      # the user-supplied GS-filtering genes
      lgg_id <- paste0("gs_lg_",x)
      lgg <- input[[lgg_id]]
      
      # check if input is empty
      if(lgg == ""){
        shinyalert("Please enter a valid gene.")
      }
      
      req(lgg != "")
      # proceed only rv not equal to input
      # update rv
      rv[[lgg_id]] <- lgg
      
      # check if valid entry
      genes <- toupper(lgg) %>% gsub(" ","",.) %>% str_split(.,"&") %>% .[[1]] %>% unique()
      
      if(genes == ""){
        shinyalert("Please enter valid gene(s).")
      }
      
      req(genes != "")
      
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
        shinyalert("Unable to detect the entered gene(s) in the selected database. Double check if your input follows the right format. Or try another database.")
      }
        
      req(length(filtered_gmts) > 0)
      
      rv[[paste0("gmts_tmp",x)]] <- filtered_gmts
      rv[[paste0("gs_lgg_",x)]] <- paste0("Filter by: ", paste0(genes, collapse = " "))


      # update gene set UI
      updateSelectInput(
        session,
        paste0("gs_l_",x)
        ,choices = names(filtered_gmts)
        ,selected=rv[[paste0("gs_l_",x)]]
      )
      
      output[[paste0("gs_lgg_",x)]] <- renderText({
        rv[[paste0("gs_lgg_",x)]]
      })
      
      # output[[paste0("gs_lgg_",x,"_tag")]] <- renderUI({
      #   tags$head(tags$style(
      #     paste0("#gs_lgg_",x,"{",rv$verbTxtStyle1,"}")
      #   ))
      # })

    })
  })

  
}, ignoreInit = T)

# ------- [1D] reset to all gene sets on clear filtering ---------
# lg_input_clearbtn_lst <- reactive({
#   lapply(1:rv$variable_n, function(x){
#     db_id <- paste0("gs_lg_",x,"_reset")
#     input[[db_id]]
#   })
# })
# 
# observeEvent(lg_input_clearbtn_lst(),{
#   lst <- lg_input_clearbtn_lst()
#   array <- check_array(lst)
#   withProgress(value = 1, message = "Resetting choices to all gene sets in the selected database ...")
#   print("hihihi")
#   print(lst)
#   print("end")
# }, ignoreInit = T)

# ----- 1.2. run parameters -------
# update dynamic rvs
observe({
  rv[["ui_run_parameters"]] <- plot_run_ui(rv$variable_n)
})

# parameters for KM analysis
output$par_gear <- renderUI({
  rv[["ui_run_parameters"]]
})

# update all run parameters according to analysis #1
observeEvent(input$toall,{
  lapply(1:rv$variable_n, function(x){
    updateSliderTextInput(session, paste0("lower_",x), selected = input[["lower_1"]])
    updateSliderTextInput(session, paste0("higher_",x), selected = input[["higher_1"]])
    updateSliderTextInput(session, paste0("step_",x), selected = input[["step_1"]])
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


