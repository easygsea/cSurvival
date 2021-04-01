# ------------ Plots area ------------
output$ui_results <- renderUI({
  req(!is.null(rv[["cox_1"]]))
  
  if(rv$variable_nr == 1){
    types <- list(
      "Surv plot #1" = 1
      ,"Gender effect plot" = "gender"
      ,"Scatter plot" = "scatter"
    )
  }else{
    indi <- as.list(1:rv$variable_nr)
    names(indi) <- paste0("Surv plot #",1:rv$variable_nr)
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
        style = "position: absolute; right: 3.5em; top: 2.5em;",
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
            outputId = "download_plot", label = NULL
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
  ppl <- isolate(input$plot_type)
  withProgress(value = 1, message = "Generating plot ...",{
    if(ppl == "all" | suppressWarnings(!is.na(as.numeric(ppl)))){
      x <- input$plot_type
      
      # # the gene(s)/GS(s) as the title
      # rv[["title"]] <- ifelse(isolate(input[[paste0("cat_",x)]]=="g"),isolate(input[[paste0("g_",x)]]),isolate(input[[paste0("gs_l_",x)]]))
      rv[["title"]] <- rv[[paste0("title_",x)]]
      
      # no of cases in each group
      rv[["lels"]] <- rv[[paste0("lels_",x)]]
      
      # the cutoff percentile
      rv[["cutoff"]] <- rv[[paste0("cutoff_",x)]]
      

      # extract statistics
      res <- rv[["res"]] <- rv[[paste0("cox_",x)]]
    }

    # generate figure
    plot_surv(res,two_rows=ppl)
  })
})

# --------- 1a. plot parameters -------------
output$plot_gear <- renderUI({
  fluidRow(
    column(
      12,
      # survival analysis method
      radioGroupButtons(
        inputId = "cox_km",
        label = HTML(paste0("Survival analysis method:"),add_help(paste0("cox_km_q"))),
        choices = surv_methods,
        selected = rv[["cox_km"]],
        size = "sm",
        checkIcon = list(
          yes = icon("check-square"),
          no = icon("square-o")
        ),
        direction = "horizontal"
      )
      ,bsTooltip("cox_km_q",HTML(paste0("Select the method for analyzing and summarizing survival data. "
                                        ,"Cox regression model assesses the effect of several risk factors simultaneously,"
                                        ," while KM describe the survival according to one factor under investigation. "
                                        ))
                 ,placement = "top")
      
      # median thresholds
      ,checkboxGroupButtons(
        inputId = "median",
        label = HTML(paste0("Draw line(s) at median survival?",add_help("median_q"))),
        choices = c("Horizontal"="h",
                    "Vertical"="v"
                    ),
        selected = rv$median,
        size="s",
        checkIcon = list(
          no = tags$i(class = "fa fa-times",
                      style = "color: crimson"),
          yes = tags$i(class = "fa fa-check",
                       style = "color: green"))
      )
      ,bsTooltip("median_q",HTML(paste0("Select to draw (a) horizontal and/or vertical line(s) at median (50%) survival"
      ))
      ,placement = "top")
      # confidence intervals
      ,materialSwitch(
        inputId = "confi",
        label = HTML(paste0("<b>Plot confidence intervals?</b>",add_help("confi_q"))),
        value = rv$confi, inline = F, width = "100%",
        status = "danger"
      )
      ,bsTooltip("confi_q",HTML(paste0("If TRUE, plots 95% confidence intervals"))
                 ,placement = "top")
    )
    ,conditionalPanel(
      'input.confi',
      column(
        12,
        radioGroupButtons(
          inputId = "confi_opt",
          label = HTML(paste0("Confidence interval style:"),add_help(paste0("confi_opt_q"))),
          choiceNames = c("Ribbon", "Step"),
          choiceValues = c("ribbon","step"),
          selected = rv$confi_opt,
          size = "sm",
          checkIcon = list(
            yes = icon("check-square"),
            no = icon("square-o")
          ),
          direction = "horizontal"
        )
        ,bsTooltip("confi_opt_q",HTML(paste0("Select the confidence interval style. <b>Ribbon</b> plots areas."
                                             ," <b>Step</b> plots boundaries."
        ))
        ,placement = "top")
      )
    )
    # toggle tables for KM plot
    ,conditionalPanel(
      'input.cox_km == "km"',
      column(
        12,
        materialSwitch(
          inputId = "risk_table",
          label = HTML(paste0("<b>Plot survival table?</b>",add_help("risk_table_q"))),
          value = rv$risk_table, inline = F, width = "100%",
          status = "danger"
        )
        ,bsTooltip("risk_table_q",HTML(paste0("If TRUE, plots a table showing the number at risk over time"))
                   ,placement = "top")
        , materialSwitch(
          inputId = "cum_table",
          label = HTML(paste0("<b>Plot cumulative events table?</b>",add_help("cum_table_q"))),
          value = rv$cum_table, inline = F, width = "100%",
          status = "danger"
        )
        ,bsTooltip("cum_table_q",HTML(paste0("If TRUE, plots a table showing the cumulative number of events over time"))
                   ,placement = "top")
      )
    )
  )
})

observeEvent(input$cox_km,{rv$cox_km <- input$cox_km})
observeEvent(input$median,{rv$median <- input$median},ignoreNULL = F)
observeEvent(input$confi,{rv$confi <- input$confi})
observeEvent(input$confi_opt,{rv$confi_opt <- input$confi_opt})
observeEvent(input$risk_table,{rv$risk_table <- input$risk_table})
observeEvent(input$cum_table,{rv$cum_table <- input$cum_table})

# --------- 1b. download plot -------------
output$download_plot <- downloadHandler(
  filename = function(){paste0(rv$cox_km,"_",rv[["title"]],".pdf")},
  content = function(file) {
    pdf(file,onefile = TRUE)
    print(plot_surv(rv[["res"]]),newpage = FALSE)
    dev.off()
    # ggsave(file,print(plot_surv(rv[["res"]]),newpage = FALSE), device = "pdf", width = 10, height = 8, dpi = 300, units = "in")
  }
)

# --------- 2. display the statistics -------------
output$ui_stats <- renderUI({
  req(rv[["res"]])
  
  col_w <- 12 / length(rv[["lels"]])
  lel1 <- names(rv[["lels"]])[[length(rv[["lels"]])]]
  lel2 <- names(rv[["lels"]])[[1]]
  
  res <- rv[["res"]][[rv$cox_km]]
  hr <- res[["hr"]]
  p <- res[["p"]]
  
  
  column(
    12,style="display: inline-block;vertical-align:top; width: 100%;word-break: break-word;",
    h3(paste0("Statistics by ",names(surv_methods)[surv_methods == rv$cox_km])),
    boxPad(
      color = "light-blue",
      fluidRow(
        column(
          6,
          descriptionBlock(
            header = hr,
            text = HTML(paste0("HR (hazard ratio)",add_help("hr_q")))
            ,rightBorder = T
          )
        )
        ,column(
          6,
          descriptionBlock(
            header = p,
            text = "P-value"
            ,rightBorder = F
          )
        )
        ,bsTooltip("hr_q",HTML(paste0("Only applicable to regression analysis by Cox PH model. HR > 1 indicates that the ",lel1," group have higher risk of death than the ",lel2," group. <i>Vice versa</i>,"
                                      ," HR < 1 indicates a lower risk of death for the ",lel1," as compared to the ",lel2))
                   ,placement = "bottom")
      )
    )
    ,boxPad(
      color = "gray",
      if(rv[["cutoff"]] != ""){
        column(
          12, align="center",
          HTML(paste0("Cutoff percentile: <b>",rv[["cutoff"]],"</b>"))
        )
      }
      ,tagList(
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
      ,renderPrint({print(res[["stats"]])})
    )
  )
})