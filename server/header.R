# available databases
observeEvent(input$db_download,{
  showModal(
    modalDialog(
      title = h3(HTML("cSurvival source data")),
      fluidRow(
       column(
         12
         
       )
      )
      ,size = "l",
      easyClose = TRUE
      ,footer = modalButton("Dismiss")
    )
  )
})