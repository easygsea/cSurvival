
# ------- parameters for a KM plot -----------
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