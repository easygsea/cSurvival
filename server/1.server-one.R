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
    rv$variable_n <- 2
  }else{
    rv$variable_n <- n
  }
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
  withProgress(value = 1, message = "Extracting data from the selected database. Please wait a minute...",{
    lapply(array, function(x) {
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
    })
  })
})

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


