# ------------ Plots area ------------
output$ui_results <- renderUI({
  req(rv$try_error < 1)
  req(!is.null(rv[["cox_1"]]))
  x <- rv$variable_nr
  rv[["dtypes"]] <- ""
  
  if(x == 1){
    types <- list(
      "Survival plot #1" = 1
      ,"Gender effect plot" = "gender"
    )
    
    dtype1 <- rv[["data_type_1"]]
    if(dtype1 != "cnv"){
      dtype1_name <- call_datatype_from_rv(dtype1)
      dtype1_scatter <- as.list("scatter")
      names(dtype1_scatter) <- paste0(dtype1_name,"-survival scatter")
      
      if(dtype1 == "snv"){
        l_plot <- list("Mutation statistics"="snv_stats")
      }else{
        l_plot <- dtype1_scatter
      }
      
    }
  }else{
    indi <- as.list(1:x)
    names(indi) <- paste0("Survival plot #",1:x)
    types <- list(
      "Interaction Survival plot" = "all"
    ) %>% c(.,indi)
    
    # check data types
    dtypes <- grep("^data_type_",{names(rv)},value=T)
    dtypes <- sapply(dtypes,function(x) rv[[x]])
    dtypes_names <- call_datatype_from_rv(dtypes)

    # check if both SNV, a single SNV; then decide plot type and add to pre-set plot options
    if(unique(dtypes) == "snv"){
      # l_plot <- list("Mutation statistics"="snv_stats")
    }else if(unique(dtypes) == "cnv"){
      
    }else if(any(c("snv","cnv") %in% dtypes)){
      l_plot <- list("Violin plot"="violin")
    }else{
      l_plot <- as.list("scatter2")
      names(l_plot) <- paste0(paste0(dtypes_names,collapse = "-")," scatter")
    }
    # save to rv
    names(dtypes) <- dtypes_names
    rv[["dtypes"]] <- dtypes
  }
  
  # assemble all types of plots
  types <- c(types, l_plot, "Differential expression & enrichment analysis"="gsea")
  
  if(rv$cox_km == "km" & (rv$risk_table | rv$cum_table)){
    h_plot <- "725px"
  }else{
    h_plot <- "580px"
  }
  
  # check if to generate survival curves
  surv_yn <- if_surv()
  
  # width for the plot area
  if(rv$plot_type != "snv_stats" & rv$plot_type != "gsea"){
    area_w <- 7
  }else{
    area_w <- 12
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
      ,if(surv_yn & typeof(rv[[paste0("df_",input$plot_type)]]) == "list"){
        column(12,align="left",
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
      ,if(typeof(rv[[paste0("df_",input$plot_type)]]) != "list" & rv$plot_type != "scatter" & rv$plot_type != "scatter2" & rv$plot_type != "snv_stats" & rv$plot_type != "gsea"){
        column(
          12, align="center",
          uiOutput("ui_error")
        )
      }else{
        div(
          column(
            area_w, align = "left",
            h3(rv[["title"]]),
            if(surv_yn){
              plotOutput("cox_plot",height = h_plot)
            }else if(rv$plot_type == "scatter" | rv$plot_type == "scatter2"){
              plotlyOutput("scatter_plot", height = "585px")
            }else if(rv$plot_type == "snv_stats"){
              plotlyOutput("snv_stats_plot", height = "585px")
            }else if(rv$plot_type == "gsea"){
              uiOutput("ui_gsea")
            }
            ,div(
              align = "left",
              style = "position: absolute; right: 2.5em; top: 1.5em;",
              div(
                id = "gear_btn",
                style="display: inline-block;vertical-align:top;",
                if(rv$plot_type != "snv_stats" & rv$plot_type != "gsea"){
                  dropdown(
                    uiOutput("plot_gear"),
                    circle = TRUE, status = "danger", style = "material-circle",
                    size="sm", right = T,
                    icon = icon("gear"), width = "300px",
                    tooltip = tooltipOptions(title = "Click for advanced plotting parameters", placement = "top")
                  )  
                }
              )
              ,div(
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
          ,conditionalPanel(
            'input.plot_type != "snv_stats" & input.plot_type != "gsea"',
            column(
              5, align="left",
              uiOutput("ui_stats")
            )
          )
        )
      }
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
observeEvent(list(input$plot_type,rv[["title_1"]]),{
  req(!is.null(input$plot_type))
  req(rv$surv_plotted == "plotted")
  x <- rv$plot_type <- input$plot_type
  # # the gene(s)/GS(s) as the title
  # rv[["title"]] <- ifelse(isolate(input[[paste0("cat_",x)]]=="g"),isolate(input[[paste0("g_",x)]]),isolate(input[[paste0("gs_l_",x)]]))
  if(x == "scatter"){x <- 1}
  rv[["title"]] <- rv[[paste0("title_",x)]]
  
  # --------------- perform differential expression analysis ---------
  if(x == "gsea" & rv$gsea_done==""){
    # a) gene expression matrix
    df_gene <- de_dfgene()
    #    saveRDS(df_gene,"df_gene_all.rds")
    # b) design matrix
    if(rv$variable_nr == 1){
      df_design <- rv[["df_gender"]]
    }else{
      df_design <- rv[["df_all"]]
    }
    #    saveRDS(df_design,"df_design_all.rds")
    # # run DE analysis
    # 1) create dgelist
    y <- DGEList(counts=df_gene)
    
    # 2) filter low expressing genes
    min_n <- min(table(df_design$level))
    keep <- rowSums(y$counts>1) >= min_n
    y <- y[keep,,keep.lib.sizes=TRUE]
    
    # 3) voom directly on counts, if data are very noisy, as would be used for microarray
    v <- voom(y, design, plot=F, normalize.method="quantile")
    
    # 4) DEG analysis
    fit <- lmFit(v, design)
    fit <- eBayes(fit,trend=TRUE, robust=TRUE)
    
    # results
    results <- decideTests(fit)
    summary(results)
    
    # available coefficients
    coefs <- colnames(design)[-1]
    
    # name the coefficients
    lels1 <- levels(df_design_gender$level.x) %>% rev()
    lels2 <- levels(df_design_gender$level.y) %>% rev()
    names(coefs) <- c(
      paste0(lels1, collapse = " vs. "),
      paste0(lels2, collapse = " vs. "),
      paste0("(",paste0(lels2, collapse = " vs. "),") vs. (",paste0(lels1, collapse = " vs. "),")")
    )
    
    # export DEG table
    degss = lapply(coefs, function(x){
      topTable(fit, coef=x,sort.by="P",number=Inf)
    })
    names(degss) <- coefs
    
    # saveRDS(degss,"degss.rds")
    # saveRDS(coefs,"coefs.rds")
    rv$gsea_done <- "yes"
  }
})

output$cox_plot <- renderPlot({
  x <- rv$plot_type
  req(if_surv() & typeof(rv[[paste0("df_",x)]]) == "list")
  withProgress(value = 1, message = "Generating plot ...",{
    # no of cases in each group
    rv[["lels"]] <- rv[[paste0("lels_",x)]]
    
    # the cutoff percentile
    rv[["cutoff"]] <- rv[[paste0("cutoff_",x)]]
    
    
    # extract statistics
    res <- rv[["res"]] <- rv[[paste0("cox_",x)]]
    rv$surv_plotted <- "plotted"
    # generate survival curve
    plot_surv(res,two_rows=x)
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
  }else if(rv$plot_type == "scatter" | rv$plot_type == "scatter2" & !is.null(rv$scatter_gender)){
    fluidRow(
      column(
        12,
        selectizeInput(
          "scatter_gender",
          HTML(paste0("Select gender group:",add_help("scatter_gender_q")))
          ,choices = rv[["genders"]]
          ,selected = rv$scatter_gender
          ,multiple = T
        )
        ,bsTooltip("scatter_gender_q",HTML(paste0("Select gender group(s) to visualize"))
                   ,placement = "top")
        # ,conditionalPanel(
        #   "input.scatter_gender.length > 1",
          ,materialSwitch(
            inputId = "scatter_gender_y",
            label = HTML(paste0("<b>Color data points by gender?</b>",add_help("scatter_gender_y_q"))),
            value = rv$scatter_gender_y, inline = F, width = "100%",
            status = "danger"
          )
          ,bsTooltip("scatter_gender_y_q",HTML(paste0(
            "If TRUE, color scatter points by gender"
          )),placement = "top")
        # )
        ,materialSwitch(
          inputId = "scatter_log_x",
          label = HTML(paste0("<b>Log2 transform x-axis values?</b>",add_help("scatter_log_x_q"))),
          value = rv$scatter_log_x, inline = F, width = "100%",
          status = "danger"
        )
        ,bsTooltip("scatter_log_x_q",HTML(paste0("If TRUE, values along the x-axis are log2 transformed"))
                   ,placement = "top")
        ,materialSwitch(
          inputId = "scatter_log_y",
          label = HTML(paste0("<b>Log2 transform y-axis values?</b>",add_help("scatter_log_y_q"))),
          value = rv$scatter_log_y, inline = F, width = "100%",
          status = "danger"
        )
        ,bsTooltip("scatter_log_y_q",HTML(paste0("If TRUE, values along the y-axis are log2 transformed"))
                   ,placement = "top")
        ,materialSwitch(
          inputId = "scatter_lm",
          label = HTML(paste0("<b>Draw a regression line?</b>",add_help("scatter_lm_q"))),
          value = rv$scatter_lm, inline = F, width = "100%",
          status = "danger"
        )
        ,bsTooltip("scatter_lm_q",HTML(paste0("If TRUE, draw a regerssion line"))
                   ,placement = "top")
        ,conditionalPanel(
          'input.scatter_lm',
          div(
            selectizeInput(
              "lm_method",
              HTML(paste0("Smoothing method:",add_help("lm_method_q"))),
              choices=list(
                "Linear regression model (lm)"="lm"
                ,"Generalized linear model (glm)"="glm"
                ,"Generalized additive model (gam)"="gam"
                ,"Local regression (loess)"="loess"
              )
              ,selected=rv$lm_method
            )
            ,bsTooltip("lm_method_q",HTML("Select the method to draw the regression line")
                       ,placement = "top")
          )
        )
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

observeEvent(input$scatter_gender,{rv$scatter_gender <- input$scatter_gender})
observeEvent(input$scatter_gender_y,{rv$scatter_gender_y <- input$scatter_gender_y})
observeEvent(input$scatter_log_x,{rv$scatter_log_x <- input$scatter_log_x})
observeEvent(input$scatter_log_y,{rv$scatter_log_y <- input$scatter_log_y})
observeEvent(input$scatter_lm,{rv$scatter_lm <- input$scatter_lm})


# --------- 1b. download plot -------------
output$download_plot <- downloadHandler(
  filename = function(){
    if(if_surv()){
      paste0(toupper(rv$cox_km),"_",rv[["title"]],".pdf")
    }else if(rv$plot_type == "scatter"){
      paste0("Scatter_",rv[["title"]],".html")
    }else if(rv$plot_type == "snv_stats"){
      paste0("Mutation_",rv[["title"]],".html")
    }
  },
  content = function(file) {
    if(if_surv()){
      pdf(file,onefile = TRUE)
      print(plot_surv(rv[["res"]]),newpage = FALSE)
      dev.off()
      # ggsave(file,print(plot_surv(rv[["res"]]),newpage = FALSE), device = "pdf", width = 10, height = 8, dpi = 300, units = "in")
    }else if(rv$plot_type == "scatter"){
      saveWidget(as_widget(rv[["scatter_plot"]] ), file, selfcontained = TRUE)
    }else if(rv$plot_type == "snv_stats"){
      saveWidget(as_widget(rv[["snv_stats_fig"]] ), file, selfcontained = TRUE)
    }
  }
)

# --------- 2. display the statistics -------------
output$ui_stats <- renderUI({
  req(rv$plot_type != "snv_stats" & rv$plot_type != "gsea")
  if(if_surv()){
    req(!is.null(rv[["res"]]))
    
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
    stats_title <- paste0("Statistics by ",names(surv_methods)[surv_methods == rv$cox_km])
  }else if(rv$plot_type == "scatter" | rv$plot_type == "scatter2"){
    req(!is.null(rv[["res_scatter"]]))
    stats_title <- "Correlation statistics"
    res <- rv[["res_scatter"]]
    hr <- round(res$estimate, 2)
    hr_title <- "Coefficient"
    p <- format(as.numeric(res$p.value), scientific = T, digits = 3)
    p_title <- "P-value"
    p_w <- 6
  }
  
  column(
    12,style="display: inline-block;vertical-align:top; width: 100%;word-break: break-word;",
    h3(stats_title),
    if(if_surv()){
      conditionalPanel(
        '(input.plot_type == "all" | input.plot_type == "gender") & input.cox_km == "km"',
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
      )
    }else if(rv$plot_type == "scatter" | rv$plot_type == "scatter2"){
      selectizeInput(
        "cor_method",
        NULL,
        choices = list(
          "Pearson's product-moment correlation" = "pearson"
          ,"Kendall's rank correlation tau" = "kendall"
          ,"Spearman's rank correlation rho" = "spearman"
        )
        ,selected = rv[["cor_method"]]
      )
    }
    ,boxPad(
      color = "light-blue",
      fluidRow(
        if(if_surv() & rv$cox_km == "cox"){
          column(
            6,
            descriptionBlock(
              header = hr,
              text = HTML(paste0(hr_title,add_help("hr_q")))
              ,rightBorder = T
            )
            ,bsTooltip("hr_q",HTML(hr_q),placement = "bottom")
          )
        }else if(rv$plot_type == "scatter" | rv$plot_type == "scatter2"){
          column(
            6,
            descriptionBlock(
              header = hr,
              text = hr_title
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
      )
    )
    ,boxPad(
      color = "gray",
      fluidRow(
        if(if_surv()){
          column(12,
            uiOutput("ui_cutoff")
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
          )
        }
        ,if(if_surv()){
          column(
            12,
            renderPrint({print(res[["stats"]])})
          )
        }else if(rv$plot_type == "scatter" | rv$plot_type == "scatter2"){
          column(
            12,
            renderPrint({print(res)})
          )
        }
      )
    )
  )
})

# cutoffs
output$ui_cutoff <- renderUI({
  cutoff <- rv[[paste0("cutoff_",rv$plot_type)]]
  req(cutoff != "")
  if(rv[["data_type_1"]] == "cnv"){
    txt <- "Copy number group: "
  }else{
    txt <- "Cutoff percentile: "
  }
  column(
    12, align="center",
    HTML(paste0(txt,cutoff))
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

# ----------- 4[A]. expression-survival days scatter plot ---------------
output$scatter_plot <- renderPlotly({
  withProgress(value = 1,message = "Updating plot ...",{
    if(rv$plot_type == "scatter"){
      # retrieve survival data
      if(typeof(rv[["df_gender"]]) == "list"){
        req(length(rv$scatter_gender) > 0)
        df_survival <- rv[["df_gender"]] %>% dplyr::select(patient_id,survival_days,level.y) %>%
          dplyr::filter(level.y %in% rv$scatter_gender) 
        genders <- df_survival$`level.y`
        df_survival <- df_survival %>%
          dplyr::select(-level.y)
      }else{
        df_survival <- rv[["df_1"]] %>% dplyr::select(patient_id,survival_days)
      }
      
      df <- rv[["exprs_1"]]
      
      # the unit, e.g. expression (fpkm)
      exp_unit <- input_mode_name("1")

      df_o <- df <- df_survival %>% inner_join(df, by="patient_id")
      if(ncol(df) == 3){
        colnames(df) <- c("patient_id","survival_days","exp")
        exprs <- df$exp
        rv[["gs_no"]] = T; exp_type = exp_unit
        if(rv$scatter_log_y){
          df_y <- log2(df$exp+1)
          ylab <- paste0("Log2 (",exp_unit," + 1)")
        }else{
          df_y <- df$exp
          ylab <- exp_unit
        }
      }else{
        rv[["gs_no"]] = F; exp_type = "mean of Z scores"
        exprs <- df_y <- rowMeans(df[,c(-1,-2)]) %>% unlist(.) %>% unname(.)
        if(rv$scatter_log_y){
          z_min <- min(df_y)
          df_y <- log2(df_y - z_min + 1)
          ylab <- "Log2 (Z score - min(Z scores) + 1)"
        }else{
          ylab <- "Average of gene expression Z scores"
        }
      }
      
      if(rv$scatter_log_x){
        df_x <- log2(df$survival_days+1)
        xlab <- "Log2 (survival days + 1)"
      }else{
        df_x <- df$survival_days
        xlab <- "Survival days"
      }
      
      # calculate correlation
      rv[["res_scatter"]] <- cor.test(df_x, df_y, method = rv$cor_method)
      
      # convert into ranks in necessary
      if(rv$cor_method == "kendall" | rv$cor_method == "spearman"){
        df_x <- rank(df_x,ties.method = "first")
        df_y <- rank(df_y,ties.method = "first")
        xlab <- "Ranks in survival days"; ylab <- paste0("Ranks in ",exp_unit)
      }
      
      # draw the figure
      if(rv$scatter_gender_y & length(rv$scatter_gender)>1){
        fig <- ggplot(df
                      ,aes(x=df_x, y=df_y
                           ,text=paste0(
                             "Patient ID: <b>",.data[["patient_id"]],"</b>\n",
                             "Survival days: <b>",.data[["survival_days"]],"</b>\n",
                             exp_type,": <b>",signif(exprs,digits=3),"</b>"
                           )
                      )) +
          geom_point(aes(color=genders)) + #, shape=genders
          scale_color_manual(values=c("#00BFC4", "#F8766D")) #+ scale_shape_manual(values=c(16, 8))
      }else{
        if(!rv$scatter_gender_y){
          col <- "#939597"
        }else{
          g_val <- as.numeric(genders) %>% unique(.)
          if(g_val == 1){
            col <- "#00BFC4"
          }else{
            col <- "#F8766D"
          }
        }
        fig <- ggplot(df
                      ,aes(x=df_x, y=df_y
                           ,text=paste0(
                             "Patient ID: <b>",.data[["patient_id"]],"</b>\n",
                             "Survival days: <b>",.data[["survival_days"]],"</b>\n",
                             exp_type,": <b>",signif(exprs,digits=3),"</b>"
                           )
                      )) +
          geom_point(color=col)
      }
      fig <- fig + 
        xlab(xlab) +
        ylab(ylab)
      
      # draw a regression line
      if(rv$scatter_lm){
        fig <- fig + geom_smooth(method=rv$lm_method,fill="#F5DF4D",inherit.aes = F,aes(df_x, df_y))
      }
      
      rv[["scatter_plot"]] <- suppressWarnings(ggplotly(fig,tooltip = "text"))
    }else if(rv$plot_type == "scatter2"){
      df <- rv[["exprs_1"]] %>% inner_join(rv[["exprs_2"]], by="patient_id") %>%
        inner_join(dplyr::select(rv$df_survival, patient_id, gender), by="patient_id")
      colnames(df) <- c("patient_id","expa","expb","gender")
      genders <- df$gender
      if(is.null(rv$scatter_gender)){
        rv$scatter_gender <- rv[["genders"]] <- unique(genders)
      }
      df <- df %>% dplyr::filter(gender %in% rv$scatter_gender)
      # the unit, e.g. expression (fpkm)
      exp_unita <- input_mode_name("1")
      exp_unitb <- input_mode_name("2")
      
      # log y, if prompted
      if(rv$scatter_log_y){
        df_y <- log2(df$expb+1)
        ylab <- paste0("Log2 (",exp_unitb," + 1)")
      }else{
        df_y <- df$expb
        ylab <- exp_unitb
      }
      
      # log x, if prompted
      if(rv$scatter_log_x){
        df_x <- log2(df$expa+1)
        xlab <- paste0("Log2 (",exp_unita," + 1)")
      }else{
        df_x <- df$expa
        xlab <- exp_unita
      }
      
      # calculate correlation
      rv[["res_scatter"]] <- cor.test(df_x, df_y, method = rv$cor_method)
      
      # convert into ranks in necessary
      if(rv$cor_method == "kendall" | rv$cor_method == "spearman"){
        df_x <- rank(df_x,ties.method = "first")
        df_y <- rank(df_y,ties.method = "first")
        xlab <- paste0("Ranks in ",exp_unita); ylab <- paste0("Ranks in ",exp_unitb)
      }
      xlab <- paste0(rv[["title_1"]],": ",xlab)
      ylab <- paste0(rv[["title_2"]],": ",ylab)
      
      # figure hovers
      txt <- function(.data){paste0(
        "Patient ID: <b>",.data[["patient_id"]],"</b>\n",
        exp_unita,": <b>",signif(.data[["expa"]],digits=3),"</b>\n",
        exp_unitb,": <b>",signif(.data[["expb"]],digits=3),"</b>"
      )}
      
      # draw the figure
      if(rv$scatter_gender_y & length(rv$scatter_gender)>1){
        fig <- ggplot(df
                      ,aes(x=df_x, y=df_y
                           ,text=txt(.data)
                      )) +
          geom_point(aes(color=genders)) + #, shape=genders
          scale_color_manual(values=c("#00BFC4", "#F8766D")) #+ scale_shape_manual(values=c(16, 8))
      }else{
        if(!rv$scatter_gender_y){
          col <- "#939597"
        }else{
          g_val <- as.numeric(genders) %>% unique(.)
          if(g_val == 1){
            col <- "#00BFC4"
          }else{
            col <- "#F8766D"
          }
        }
        fig <- ggplot(df
                      ,aes(x=df_x, y=df_y
                           ,text=txt(.data)
                      )) +
          geom_point(color=col)
      }
      
      # rename x- and y- axis
      fig <- fig + 
        xlab(xlab) +
        ylab(ylab)
      
      # draw a regression line
      if(rv$scatter_lm){
        fig <- fig + geom_smooth(method=rv$lm_method,fill="#F5DF4D",inherit.aes = F,aes(df_x, df_y))
      }
      
      rv[["scatter_plot"]] <- suppressWarnings(ggplotly(fig,tooltip = "text"))
    }
    
  })
})

# change method of correlation calculation when prompted
observeEvent(input$lm_method,{rv$lm_method <- input$lm_method})
observeEvent(input$cor_method,{rv$cor_method <- input$cor_method})

# ----------- 4[B]. SNV statistics bar plot ---------------
output$snv_stats_plot <- renderPlotly({
  # calculated statistics
  muts <- rv[["mutations_1"]] %>% sapply(., function(x){
    strsplit(x, "\\|")[[1]]
  }) %>% unlist(.)
  stats <- table(muts)

  dat <- data.frame(
    Mutation = names(stats),
    Frequency = as.numeric(stats)
  )
  
  # non-synonymous
  non_id <- "nonsynonymous_1"
  nons <- ifelse_rv(non_id) %>% tolower(.)
  nons_cat <- sapply(dat$Mutation, function(x){
    if(tolower(x) %in% nons){
      "Nonsynonymous"
    }else{
      "Synonymous"
    }
  }) %>% unname(.)

  dat$Category <- nons_cat
  dat <- dat %>% dplyr::arrange(Category,Frequency)
  # patient cases that have each mutation
  Cases <- lapply(dat$Mutation, function(x){
    names(muts)[muts == x] %>% 
      breakvector(.) %>%
      paste0(., collapse = ", ") %>% addlinebreaks(.)
  })
  
  # set Mutation data as factor for ordering in ggplotly
  dat$Mutation <- factor(dat$Mutation, levels = dat$Mutation)
  
  fig <- ggplot(data=dat,aes(x=Mutation, y=Frequency, fill=Category, Cases=Cases)) + 
    geom_bar(stat="identity") +
    xlab("") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  rv[["snv_stats_fig"]] <- ggplotly(fig, tooltip = c("Mutation","Frequency","Category","Cases"))
})

#------------- 5. Error UI -------------
output$ui_error <- renderUI({
  req(!is.null(input$plot_type))
  if(input$plot_type == "gender"){
    div(
      br(),
      p(style="color:gray;font-size:120%;","Unable to assess for gender effect. The selected project(s) only contain(s) ",rv[["df_gender"]]," data.")
      ,br()
    )
  }#else if(rv$plot_type == "1" | rv$plot_type == "2"){
  #   df_name <- paste0("df_",rv$plot_type)
  #   df <- rv[[df_name]]
  #   print(head(df))
  #   lels <- levels(df$level)
  #   i <- rv$plot_type
  #   dtype <- rv[["data_type_"]][[i]]
  #   dtype_name <- call_datatype_from_rv(dtype)
  #   p(style="color:gray;",paste0("Unable to perform survival analysis on the "
  #                                ,dtype_name," data of ", paste0(rv$project, collapse = ", ")
  #                                ,". The selected project(s) only contain(s) ",lels," data."))
  # }
})
