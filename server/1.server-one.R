#======================================================================#
####                          STEP 1. parameters                   ####
#======================================================================#
observeEvent(input$variable_n,{
  rv$variable_n <- input$variable_n
  update_all()
  init_rvs()
  rv[["ui_parameters"]] <- plot_ui(rv$variable_n)
})

observe({
  rv[["ui_run_parameters"]] <- plot_run_ui(rv$variable_n)
})

output$ui_parameters <- renderUI({
  rv[["ui_parameters"]]
})

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

# ----- run parameters -------
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

# -------- auto GMT loading ---------
gmt_input_lst <- reactive({
  lapply(1:rv$variable_n, function(x){
    db_id <- paste0("gs_db_",x)
    input[[db_id]]
  })
})

observeEvent(gmt_input_lst(),{
  lst <- gmt_input_lst()
  lst_u <- lst %>% unlist() %>% unique()
  req(!is.null(lst_u) & lst_u != "")
  n <- rv$variable_n
  if(n>1){req(!is.null(lst[[2]]))}
  if(lst_u[1] == ""){array <- 2}else{array <- 1:n}
  withProgress(value = 1, message = "Extracting data from the selected database. Please wait a minute...",{
    for(x in array){
      gs_db_id <- paste0("gs_db_",x)
      gs_lib_id <- paste0("gs_l_",x)
      
      db <- rv[[gs_db_id]] <- isolate(input[[gs_db_id]])
      
      req(db != "")
      gmt_path <- retrieve_gmt_path(db)
      gmt <- gmtPathways(gmt_path)

      # update variables
      rv[[paste0("gmts",x)]] <- gmt

      # update gene set UI
      updatePickerInput(
        session,
        gs_lib_id
        ,choices = names(gmt)
        ,selected=rv[[gs_lib_id]]
      )
    }
  })
})

# ------- verbatimText feedback on GMT genes ---------
lib_input_lst <- reactive({
  lapply(1:rv$variable_n, function(x){
    db_id <- paste0("gs_l_",x)
    input[[db_id]]
  })
})

observe({
  lst <- lib_input_lst()
  lst_u <- lst %>% unlist() %>% unique()
  req(!is.null(lst_u) & lst_u != "")
  n <- rv$variable_n
  if(n>1){req(!is.null(lst[[2]]))}
  if(lst_u[1] == ""){array <- 2}else{array <- 1:n}
  lapply(array, function(x){
    gs_lib_id <- paste0("gs_l_",x)
    gs_lib_genes_id <- paste0("gs_lgs_",x)
    
    req(input[[gs_lib_id]])
    gs <- rv[[gs_lib_id]] <- input[[gs_lib_id]]
    
    req(!is.null(gs) & gs != "")
    genes <- rv[[paste0("gmts",x)]][[gs]]

    output[[gs_lib_genes_id]] <- renderText({
      paste0("(n=",length(genes),") ", paste0(genes, collapse = " "))
    })
  })
})

