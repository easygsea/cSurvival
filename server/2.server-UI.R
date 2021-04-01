# ------------ Plots area ------------
output$ui_results <- renderUI({
  req(!is.null(rv[["cox_1"]]))
  x <- rv$variable_nr
  if(x == 1){
    types <- list(
      "Survival plot #1" = 1
      ,"Gender effect plot" = "gender"
    )
    
    dtype1 <- rv[["data_type_1"]]
    if(dtype1 == "snv"){
      l_plot <- list("Mutation statistics"="snv_stats")
    }else{
      l_plot <- list("Scatter plot" = "scatter")
    }
    
    types <- c(types, l_plot)
  }else{
    indi <- as.list(1:x)
    names(indi) <- paste0("Survival plot #",1:x)
    types <- list(
      "Interaction Survival plot" = "all"
    ) %>% c(.,indi,list(
      "Violin" = "violin"
    ))
  }
  
  if(rv$cox_km == "km" & (rv$risk_table | rv$cum_table)){
    h_plot <- "725px"
  }else{
    h_plot <- "580px"
  }
  
  # check if to generate survival curves
  surv_yn <- if_surv()

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
        selected = rv$plot_type,
        # size = "sm",
        checkIcon = list(
          yes = icon("check-square"),
          no = icon("square-o")
        ),
        # status = "primary",
        direction = "horizontal"
      )
      # ,tags$hr(style = "border-color: #F5DF4D;")
      ,if(surv_yn){
        div(
          # survival analysis method
          radioGroupButtons(
            inputId = "cox_km",
            label = HTML(paste0("Select survival analysis method:"),add_help(paste0("cox_km_q"))),
            choices = surv_methods,
            selected = rv[["cox_km"]],
            status = "danger",
            # size = "sm",
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
        )
      }
      ,column(
        7, align = "left",
        h3(rv[["title"]]),
        if(surv_yn){
          plotOutput("cox_plot",height = h_plot)
        }else if(rv$plot_type == "scatter"){
          plotlyOutput("scatter_plot", height = "585px")
        }
        ,div(
          align = "left",
          style = "position: absolute; right: 3.5em; top: 1.5em;",
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
        5, align="left",
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

# --------- 1. display the survival curve / scatter plot / mutation statistics ---------
observeEvent(input$plot_type,{
  x <- rv$plot_type <- input$plot_type
  req(if_surv())
  output$cox_plot <- renderPlot({
    withProgress(value = 1, message = "Generating plot ...",{
      # # the gene(s)/GS(s) as the title
      # rv[["title"]] <- ifelse(isolate(input[[paste0("cat_",x)]]=="g"),isolate(input[[paste0("g_",x)]]),isolate(input[[paste0("gs_l_",x)]]))
      rv[["title"]] <- rv[[paste0("title_",x)]]
      
      # no of cases in each group
      rv[["lels"]] <- rv[[paste0("lels_",x)]]
      
      # the cutoff percentile
      rv[["cutoff"]] <- rv[[paste0("cutoff_",x)]]
      
      
      # extract statistics
      res <- rv[["res"]] <- rv[[paste0("cox_",x)]]
      
      # generate survival curve
      plot_surv(res,two_rows=x)
    })
  })
})

# --------- 1a. plot parameters -------------
output$plot_gear <- renderUI({
  if(if_surv()){
    fluidRow(
      column(
        12,
        # median thresholds
        checkboxGroupButtons(
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
  }else if(rv$plot_type == "scatter"){
    fluidRow(
      column(
        12,
        materialSwitch(
          inputId = "scatter_log_x",
          label = HTML(paste0("<b>Log2 transform survival days?</b>",add_help("scatter_log_x_q"))),
          value = rv$scatter_log_x, inline = F, width = "100%",
          status = "danger"
        )
        ,bsTooltip("scatter_log_x_q",HTML(paste0("If TRUE, survival days (x-axis) are log2 transformed"))
                   ,placement = "top")
        ,materialSwitch(
          inputId = "scatter_log_y",
          label = HTML(paste0("<b>Log2 transform gene expression values?</b>",add_help("scatter_log_y_q"))),
          value = rv$scatter_log_y, inline = F, width = "100%",
          status = "danger"
        )
        ,bsTooltip("scatter_log_y_q",HTML(paste0("If TRUE, gene expression values (FPKM, y-axis) are log2 transformed"))
                   ,placement = "top")
      )
    )
  }
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
  
  n_lels <- length(rv[["lels"]])
  col_w <- 12 / n_lels
  lel1 <- names(rv[["lels"]])[[n_lels]]
  lel2 <- names(rv[["lels"]])[[1]]
  
  res <- rv[["res"]][[rv$cox_km]]
  hr <- res[["hr"]]
  p <- res[["p"]]
  
  if(rv$cox_km == "cox"){
    p_w <- 6
  }else{
    p_w <- 12
  }

  if(rv$cox_km == "cox" & rv$plot_type == "all"){
    hr_title <- "HR (hazard ratios)"
    p_title <- "P-values"
    lel1 <- gsub("_"," and/or ",lel1)
    lel2 <- gsub("_"," and/or ",lel2)
  }else{
    hr_title <- "HR (hazard ratio)"
    p_title <- "P-value"
  }
  
  hr_q <- paste0("Only applicable to regression analysis by Cox PH model. HR > 1 indicates that the ",lel1," group have higher risk of death than the ",lel2," group. <i>Vice versa</i>,"
                 ," HR < 1 indicates a lower risk of death for the ",lel1," as compared to the ",lel2)
  
  column(
    12,style="display: inline-block;vertical-align:top; width: 100%;word-break: break-word;",
    h3(paste0("Statistics by ",names(surv_methods)[surv_methods == rv$cox_km])),
    conditionalPanel(
      'input.plot_type == "all" & input.cox_km == "km"',
      selectizeInput(
        "km_mul",
        NULL,
        choices = list(
          "Multiple comparisons test by Holm (1979)" = "holm"
          ,"Multiple comparisons test by Hochberg (1988)" = "hochberg"
          ,"Multiple comparisons test by Hommel (1988)" = "hommel"
          ,"Multiple comparisons test by Bonferroni correction" = "bonferroni"
          ,"Multiple comparisons test by Benjamini & Hochberg (1995)" = "BH"
          ,"Multiple comparisons test by Benjamini & Yekutieli (2001)" = "BY"
          ,"Multiple comparisons test by false discovery rate (FDR)" = "fdr"
        )
        ,selected = rv[["km_mul"]]
      )
    ),
    boxPad(
      color = "light-blue",
      fluidRow(
        if(rv$cox_km == "cox"){
          column(
            6,
            descriptionBlock(
              header = hr,
              text = HTML(paste0(hr_title,add_help("hr_q")))
              ,rightBorder = T
            )
          )
        }
        ,column(
          p_w,
          descriptionBlock(
            header = p,
            text = p_title
            ,rightBorder = F
          )
        )
        ,bsTooltip("hr_q",HTML(hr_q)
                   ,placement = "bottom")
      )
    )
    ,boxPad(
      color = "gray",
      if(rv[["cutoff"]] != ""){
        column(
          12, align="center",
          HTML(paste0("Cutoff percentile: ",rv[["cutoff"]]))
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

# ------------- 3. detailed statistics output -------
observeEvent(input$km_mul,{
  req(input$km_mul != "")
  req(input$km_mul != rv$km_mul)
  rv$km_mul <- input$km_mul
  
  # retrieve df for survival analysis
  df <- rv[[paste0("df_",rv$plot_type)]]
  
  # update statistics
  km2 <- pairwise_survdiff(Surv(survival_days, censoring_status) ~ level, data = df, p.adjust.method = rv$km_mul)
  rv[["res"]][[rv$cox_km]][["stats"]][[2]] <- km2
  
},ignoreInit = T)

# ----------- 4A. expression-survidal days scatter plot ---------------
output$scatter_plot <- renderPlotly({
  df_survival <- rv[["df_1"]] %>% dplyr::select(patient_id,survival_days)
  df <- rv[["exprs_1"]]

  df <- df_survival %>% inner_join(df, by="patient_id")
  colnames(df) <- c("patient_id","survival_days","exp")
   
  if(rv$scatter_log_x){
    df_x <- log2(df$survival_days+1)
    xlab <- "Log2 (survival days + 1)"
  }else{
    df_x <- df$survival_days
    xlab <- "Survival days"
  }
  if(rv$scatter_log_y){
    df_y <- log2(df$exp+1)
    ylab <- "Log2 (FPKM + 1)"
  }else{
    df_y <- df$exp
    ylab <- "Gene expression value (FPKM)"
  }
  fig <- ggplot(df
                ,aes(x=df_x, y=df_y
                     ,text=paste0(
                       "Patient ID: <b>",.data[["patient_id"]],"</b>\n",
                       "Survival days: <b>",.data[["survival_days"]],"</b>\n",
                       "Expression (FPKM): <b>",signif(.data[["exp"]],digits=3),"</b>"
                     )
                )) +
    geom_point(color="#939597") + 
    geom_smooth(method=lm,fill="#F5DF4D",inherit.aes = F,aes(df_x, df_y)) +
    xlab(xlab) +
    ylab(ylab)
  
  suppressWarnings(ggplotly(fig,tooltip = "text"))
})
