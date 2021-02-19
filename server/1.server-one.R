# ------- STEP 1. parameters -----------
output$ui_parameters <- renderUI({
  plot_ui(input$variable_n)
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

# -------- auto GMT loading ---------
gmt_input_lst <- reactive({
  lapply(1:input$variable_n, function(x){
    db_id <- paste0("gs_db_",x)
    input[[db_id]]
  })
})

observeEvent(gmt_input_lst(),{
  lapply(1:input$variable_n, function(x){
    db_id <- paste0("gs_db_",x)
    lib_id <- paste0("gs_l_",x)
    lib_gene_id <- paste0("gs_lg_",x)
    db <- input[[db_id]]
    
    req(db != "")
    withProgress(value = 1, message = "Extracting data from the selected database. Please wait a minute...",{
      gmt_path <- retrieve_gmt_path(db)
      gmt <- gmtPathways(gmt_path)
      
      # update variables
      gmts[[x]] <- gmt
      
      # update gene set UI
      updateSelectizeInput(
        session,
        lib_id
        ,choices = names(gmt)
      )
    })
  })
})