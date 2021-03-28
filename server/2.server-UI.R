# ------------ Plots area ------------
output$ui_results <- renderUI({
  req(!is.null(rv[["cox_fit_1"]]))
  
  if(rv$variable_n == 1){
    types <- list(
      "Surv plot #1" = 1
      ,"Distribution in TCGA"
      ,"Distribution in TARGET"
    )
  }else{
    indi <- as.list(1:rv$variable_n)
    names(indi) <- paste0("KM plot #",1:rv$variable_n)
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
      8,
      box(
        width = 12, status = "warning",
        plotOutput("cox_plot")
      )
    )
    ,column(
      4,
      box(
        width = 12,
        uiOutput("ui_stats")
      )
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
      res <- rv[[paste0("cox_fit_",x)]]
      df <- rv[[paste0("df_",x)]]
      fit <- survfit(res,data=data.frame(level = c("high", "low"), 
                                         age = rep(1, 2)
      ))
    }
    
    ggsurvplot(fit,
               title = "Survival Curves",
               xlab = "Days",
               ylab = "Survival probability",
               pval = TRUE, pval.method = TRUE,    # Add p-value &  method name
               surv.median.line = "hv",            # Add median survival lines
               # legend.title = call_datatype(x),               # Change legend titles
               # legend.labs = c("High", "Low"),  # Change legend labels
               palette = "jco",                    # Use JCO journal color palette
               risk.table = TRUE,                  # Add No at risk table
               cumevents = TRUE,                   # Add cumulative No of events table
               tables.height = 0.15,               # Specify tables height
               tables.theme = theme_cleantable(),  # Clean theme for tables
               tables.y.text = FALSE               # Hide tables y axis text
    )
  })
})