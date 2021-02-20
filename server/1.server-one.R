#======================================================================#
####                          STEP 1. parameters                   ####
#======================================================================#
observeEvent(input$variable_n,{
  update_all()
  rv[["ui_parameters"]] <- plot_ui(isolate(input$variable_n))
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

# -------- auto GMT loading ---------
gmt_input_lst <- reactive({
  lapply(1:input$variable_n, function(x){
    db_id <- paste0("gs_db_",x)
    input[[db_id]]
  })
})

observeEvent(gmt_input_lst(),{
  lst <- gmt_input_lst()
  lst_u <- lst %>% unlist() %>% unique()
  req(!is.null(lst_u) & lst_u != "")
  n <- isolate(input$variable_n)
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
  lapply(1:input$variable_n, function(x){
    lib_id <- paste0("gs_l_",x)
    input[[lib_id]]
  })
})
observeEvent(lib_input_lst(),{
  lst <- lib_input_lst()
  lst_u <- lst %>% unlist() %>% unique()
  req(!is.null(lst_u) & lst_u != "")
  n <- isolate(input$variable_n)
  if(n>1){req(!is.null(lst[[2]]))}
  if(lst_u[1] == ""){array <- 2}else{array <- 1:n}
  withProgress(value = 1, message = "Extracting gene information from the selected gene set. Please wait a minute...",{
    for(x in array){
      lib_id <- paste0("gs_l_",x)
      lib_genes_id <- paste0("gs_lgs_",x)
      gs <- rv[[lib_id]] <- isolate(input[[lib_id]])

      req(!is.null(gs) & gs != "")
      genes <- rv[[paste0("gmts",x)]][[gs]]
      
      # removeUI(
      #   selector = paste0("#",lib_genes_id)
      # )
      
      insertUI(
        selector = paste0("#vtxt_anchor",x),
        ui = verbatimTextOutput(lib_genes_id)
      )
      
      output[[lib_genes_id]] <- renderText({
        paste0("(n=",length(genes),") ", paste0(genes, collapse = " "))
      })
    }
  })
})

#======================================================================#
####                    Update rv according to Input               ####
#======================================================================#
# # update these into rv when selections change
update_all <- function(){
  for(x in 1:input$variable_n){
    lst <- list(
      # category to analyze
      cat_id <- paste0("cat_",x)
      # type of db to analyze
      ,db_id <- paste0("db_",x)
      # gene to analyze
      ,g_ui_id <- paste0("g_",x)
      # gene set to analyze
      ,gs_mode_id <- paste0("gs_mode_",x)
      ,gs_db_id <- paste0("gs_db_",x)
      ,gs_lib_id <- paste0("gs_l_",x)
      ,gs_lib_genes_id <- paste0("gs_lgs_",x)
      ,gs_gene_id <- paste0("gs_lg_",x)
      # manual gene input
      ,gs_manual_id <- paste0("gs_m_",x)
      ,gs_genes_id <- paste0("gs_mg_",x)
    )
    
    updateRV(lst)
    
    # req(input[[gs_db_id]] != "")
    # req(rv[[gs_lib_id]] != "")
    # req(rv[[paste0("gmts",x)]])
    # req(length(rv[[paste0("gmts",x)]])>0)
    # 
    # output[[gs_lib_genes_id]] <- renderText({
    #   genes <- rv[[paste0("gmts",x)]][[input[[gs_lib_id]]]]
    #   paste0("(n=",length(genes),") ", paste0(genes, collapse = " "))
    # })
  }
}