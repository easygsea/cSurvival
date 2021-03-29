# ------------ Plots area ------------
output$ui_results <- renderUI({
  req(!is.null(rv[["cox_1"]]))
  
  if(rv$variable_n == 1){
    types <- list(
      "Surv plot #1" = 1
      ,"Distribution in TCGA"
      ,"Distribution in TARGET"
    )
  }else{
    indi <- as.list(1:rv$variable_n)
    names(indi) <- paste0("Surv plot #",1:rv$variable_n)
    types <- list(
      "Interaction Surv plot" = "all"
    ) %>% c(.,indi,list(
      "Violin" = "violin"
    ))
  }
  
  box(
    width = 12, status = "danger",
    tags$script(HTML(
      "document.getElementById('ui_results').scrollIntoView();"
    )),
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
      7,
      plotOutput("cox_plot",height = "550px")
    )
    ,column(
      5,
      uiOutput("ui_stats")
    )
    ,absolutePanel(
      actionBttn(
        inputId = "up_button", label=NULL, 
        icon = icon("angle-double-up"), style="material-circle", color="default", size="md"
      ),
      tags$script(HTML(
        "document.getElementById('up_button').onclick= function(){
                    document.getElementById('ui_title').scrollIntoView()
                };"
      )),
      right = 20,
      top = 4
    )
  )
})


# --------- 1. display the survival curve ---------
output$cox_plot <- renderPlot({
  withProgress(value = 1, message = "Generating plot ...",{
    if(input$plot_type == "all"){
      fit <- survfit(res.cox)
    }else if(suppressWarnings(!is.na(as.numeric(input$plot_type)))){
      x <- input$plot_type
      df <- rv[[paste0("df_",x)]]
      res <- rv[[paste0("cox_",x)]]
      fig <- res[["fig"]]
    }
    
    # print(fig, risk.table.height = 0.3)
    fig
  })
})