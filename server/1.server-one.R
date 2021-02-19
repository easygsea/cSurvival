
# ------- parameters for a KM plot -----------
output$ui_parameters <- renderUI({
  plot_ui(input$variable_n)
})