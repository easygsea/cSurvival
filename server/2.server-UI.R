# Plots area
output$ui_results <- renderUI({
  req(!is.null(rv[["km_fit_1"]]))
  
  box(
    width = 12, status = "danger",
    column(
      # 12,
      
    )
  )
})