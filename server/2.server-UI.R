# Plots area
output$ui_results <- renderUI({
  # req(!is.null(rv[["km_fit_1"]]))
  
  if(rv$variable_n == 1){
    types <- list(
      "KM plot #1" = 1
      ,"Distribution in TCGA"
      ,"Distribution in TARGET"
    )
  }else{
    indi <- as.list(1:rv$variable_n)
    names(indi) <- paste0("KM plot #",1:rv$variable_n)
    types <- list(
      "Interaction KM plot" = "all"
    ) %>% c(.,indi,list(
      "Violin" = "violin"
    ))
  }
  
  box(
    width = 12, status = "danger",
    div(
      id="results_box",
      # src_results,
      # scroll_up_button(),
      column(
        12,align="center",
        radioGroupButtons(
          inputId = "plot_type",
          label = NULL,
          choices = types,
          # size = "sm",
          checkIcon = list(
            yes = icon("check-square"),
            no = icon("square-o")
          ),
          # status = "primary",
          direction = "horizontal"
        )
        # ,tags$hr(style = "border-color: #F5DF4D;")
      )
      ,column(
        8,
        plotOutput("ui_plot")
      )
    )
  )
})