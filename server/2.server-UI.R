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
      h3(rv[["title"]]),
      plotOutput("cox_plot",height = "580px")
      ,div(
        align = "left",
        style = "position: absolute; right: 4.5em; top: 2.5em;",
        # add a id for the gear button in introjs
        div(
          id = "gear_btn",
          style="display: inline-block;vertical-align:top;",
          dropdown(
            uiOutput("plot_gear"),
            circle = TRUE, status = "danger", style = "material-circle",
            size="sm", right = T,
            icon = icon("gear"), width = "300px",
            tooltip = tooltipOptions(title = "Click for advanced plotting parameters", placement = "top")
          )  
        ),
        div(
          id="download_btn",
          style="display: inline-block;vertical-align:top;",
          downloadBttn(
            size = "sm", color = "danger", style = "material-circle",
            outputId = "download_vis", label = NULL
          )
          ,bsTooltip("download_btn","Click to download the plot", placement = "top")
        )
      )
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
      
      # the gene(s)/GS(s) as the title
      rv[["title"]] <- isolate(input[[paste0("g_",x)]])
      
      # no of cases in each group
      rv[["lels"]] <- rv[[paste0("lels_",x)]]
      
      # the cutoff percentile
      rv[["cutoff"]] <- rv[[paste0("cutoff_",x)]]

      # extract statistics
      res <- rv[["res"]] <- rv[[paste0("cox_",x)]]
      stats <- res[[rv[["cox_km"]]]][["stats"]]
      hr <- stats$coefficients[,2]; p <- stats$coefficients[,5]
      rv[["hr"]] <- round(as.numeric(hr), 2)
      rv[["p"]] <- format(as.numeric(p), scientific = T, digits = 3)
      
      # generate figure
      fig <- plot_surv(res)
    }
    
    fig
  })
})

# --------- 1a. plot parameters -------------
output$plot_gear <- renderUI({
  fluidRow(
    column(
      12,
      radioGroupButtons(
        inputId = "cox_km",
        label = HTML(paste0("Select survival analysis method"),add_help(paste0("cox_km_q"))),
        choiceNames = c("Cox proportional-hazards model (Cox)", "Kaplan-Meier logrank (KM)"),
        choiceValues = c("cox","km"),
        selected = rv[["cox_km"]],
        size = "sm",
        checkIcon = list(
          yes = icon("check-square"),
          no = icon("square-o")
        ),
        direction = "horizontal"
      )
      ,bsTooltip("cox_km_q",HTML(paste0("The method for analyzing and summarizing survival data. "
                                        ,"KM describe the survival according to one factor under investigation. "
                                        ,"Cox regression model assesses the effect of several risk factors simultaneously."))
                 ,placement = "top")
    )
    # ,conditionalPanel(
    #   'input.cox_km == "cox',
    #   
    # )
  )
})

observeEvent(input$cox_km,{rv$cox_km <- input$cox_km})

# --------- 2. display the statistics -------------
output$ui_stats <- renderUI({
  req(rv[["hr"]])
  
  col_w <- 12 / length(rv[["lels"]])
  lel1 <- names(rv[["lels"]])[[length(rv[["lels"]])]]
  lel2 <- names(rv[["lels"]])[[1]]
  
  column(12,
    h3("Statistics"),
    boxPad(
      color = "light-blue",
      fluidRow(
        column(
          6,
          descriptionBlock(
            header = rv[["hr"]],
            text = HTML(paste0("HR (hazard ratio)",add_help("hr_q")))
            ,rightBorder = T
          )
        )
        ,column(
          6,
          descriptionBlock(
            header = rv[["p"]],
            text = "P-value"
            ,rightBorder = F
          )
        )
        ,bsTooltip("hr_q",HTML(paste0("Only applicable to Cox regression analysis. A positive HR value indicates that the ",lel1," group have higher risk of death than the ",lel2," group. <i>Vice versa</i>,"
                                      ," a negative HR value indicates a lower risk of death for the ",lel1," as compared to the ",lel2))
                   ,placement = "bottom")
      )
    )
    ,boxPad(
      color = "gray",
      column(
        12, align="center",
        HTML(paste0("Cutoff percentile: <b>",rv[["cutoff"]],"</b>"))
      ),
      tagList(
        lapply(names(rv[["lels"]]), function(x){
          no <- rv[["lels"]][[x]]
          column(
            col_w,
            descriptionBlock(
              header = no,
              text = x
              ,rightBorder = F
            )
          )
        })
      )
      ,renderPrint({print(rv[["res"]][[rv$cox_km]][["stats"]])})
    )
  )
})