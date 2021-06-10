# ------------ Plots area ------------
output$ui_results <- renderUI({
  req(rv$project != "")
  req(rv$try_error == 0)
  req(!is.null(rv[["cox_1"]]))
  x <- rv$variable_nr
  rv[["dtypes"]] <- ""

  if(is.null(rv[["title"]])){rv[["title"]] <- rv[[paste0("title_",rv$plot_type)]]}

  if(x == 1){
    types <- list(
      "Survival plot #1" = 1
      ,"Percentile tracking" = "track"
      ,"Sex effect plot" = "gender"
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
      "Interaction survival plot" = "all"
      ,"Percentile tracking" = "track"
    ) %>% c(.,indi)

    # check data types
    dtypes <- grep("^data_type_",{names(rv)},value=T)
    dtypes <- sapply(dtypes,function(x) rv[[x]])
    dtypes_names <- call_datatype_from_rv(dtypes)

    # check if both SNV, a single SNV; then decide plot type and add to pre-set plot options
    dtypes_u <- unique(dtypes)
    if(length(dtypes_u) == 1 & (dtypes_u == "snv" | dtypes_u == "cnv")){
      # # l_plot <- list("Mutation statistics"="snv_stats")
      # }else if(dtypes_u == "cnv"){
      #   
      # }
    }else if(any(c("snv","cnv") %in% dtypes)){
      # l_plot <- list("Violin plot"="violin")
    }else{
      l_plot <- as.list("scatter2")
      names(l_plot) <- paste0(paste0(dtypes_names,collapse = "-")," scatter")
    }
    # save to rv
    names(dtypes) <- dtypes_names
    rv[["dtypes"]] <- dtypes
  }


  # assemble all types of plots
  if(exists("l_plot")){types <- c(types, l_plot)}
  types <- c(types, "Transcriptome analysis by eVITTA"="gsea")

  if(rv$cox_km == "km"){
    h_plot <- ifelse(
      rv$risk_table | rv$risk_table,
      ifelse(
        rv$cum_table & rv$risk_table,
        "735px"
        ,"675px"
      )
      ,"550px"
    )
  }else{
    h_plot <- "550px"
  }

  # check if to generate survival curves
  surv_yn <- if_surv()

  # width for the plot area
  if(rv$plot_type != "snv_stats" & rv$plot_type != "gsea"){
    area_w <- 7
  }else{
    area_w <- 12
  }

  # title
  if(rv$depmapr){
    ge <- paste0(dependency_names(),add_help("dp_ge_q")," by ")
  }else{
    ge <- ""
  }
  title_div <- div(
    h3(rv[["title"]]),
    HTML(paste0('<p style="color:gray;font-size:135%;">',paste0(paste0(rv$project, collapse = " & "),": ",ge,rv$censor_time_p),"</p>"))
    ,bsTooltip("dp_ge_q",HTML(dp_ge_q()),placement = "right")
  )

  box(
    width = 12, status = "danger",
    if(rv$show_ui == "yes"){
      tags$script(HTML(
        "document.getElementById('ui_results').scrollIntoView();"
      ))
    },
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
      ,if(surv_yn & typeof(rv[[paste0("df_",input$plot_type)]]) == "list" & !rv$depmapr){
        column(12,align="left",
          # survival analysis method
          radioGroupButtons(
            inputId = "cox_km",
            label = HTML(paste0("Select survival analysis method:"),add_help(paste0("cox_km_q"))),
            choices = surv_methods_r(),
            selected = rv[["cox_km"]],
            status = "danger",
            # size = "sm",
            checkIcon = list(
              yes = icon("check-square"),
              no = icon("square-o")
            ),
            direction = "horizontal"
          )
          ,bsTooltip("cox_km_q",HTML(cox_km_txt)
          ,placement = "right")
        )
      }else if(surv_yn & rv$depmapr){
        div(
          align="left",
          fluidRow(
            column(
              12,
              column(
                12,
                title_div
              ),
              uiOutput("depmap_stats")
            ),
            column(
              12,
              column(
                12,
                tags$hr(style="border-color: light-blue;")
              ),
              column(
                6,
              ),
              column(
                6,align="center",
                selectizeInput(
                  "annot_cells"
                  ,HTML(paste0("(Optional) highlight cell lines in box plot:",add_help("annot_cells_q")))
                  ,choices = c()
                  ,multiple = T
                  ,width = "85%"
                  ,options = list(
                    `live-search` = TRUE,
                    placeholder = "Type to search ..."
                    ,onInitialize = I(sprintf('function() { this.setValue(%s); }',""))
                  )
                )
                ,bsTooltip("annot_cells_q",HTML(paste0(
                  "Select to annotate the cell lines of interest, if any, in the box plot"
                )),placement = "right")
              )
            )
          ),
          fluidRow(
            column(
              12,
              column(
                6,
                plotlyOutput("dens_plot",height = "500px", width = "100%")
                ,div(
                  align = "left",
                  style = "position: absolute; right: 1.5em; top: -3em;",
                  div(
                    id = "gear_btn_dens",
                    style="display: inline-block;vertical-align:top;",
                    if(rv$plot_type != "snv_stats" & rv$plot_type != "gsea"){
                      dropdown(
                        uiOutput("plot_gear_dens"),
                        circle = TRUE, status = "danger", style = "material-circle",
                        size="sm", right = T,
                        icon = icon("gear"), width = "300px",
                        tooltip = tooltipOptions(title = "Click for advanced plotting parameters", placement = "top")
                      )
                    }
                  )
                  ,div(
                    id="download_btn_dens",
                    style="display: inline-block;vertical-align:top;",
                    downloadBttn(
                      size = "sm", color = "danger", style = "material-circle",
                      outputId = "download_plot_dens", label = NULL
                    )
                    ,bsTooltip("download_btn_dens","Click to download plot", placement = "top")
                  )
                )
              ),
              column(
                6,
                plotlyOutput("dens_stats_plot",height = "500px", width = "100%")
                ,div(
                  align = "left",
                  style = "position: absolute; right: 1.5em; top: -3em;",
                  div(
                    id="download_btn_box",
                    style="display: inline-block;vertical-align:top;",
                    downloadBttn(
                      size = "sm", color = "danger", style = "material-circle",
                      outputId = "download_plot_box", label = NULL
                    )
                    ,bsTooltip("download_btn_box","Click to download plot", placement = "top")
                  )
                )
              )
            )
          )
        )
      }
      ,if(typeof(rv[[paste0("df_",input$plot_type)]]) != "list" & rv$plot_type != "scatter" & rv$plot_type != "scatter2" & rv$plot_type != "snv_stats" & rv$plot_type != "gsea" & rv$plot_type != "track"){
        column(
          12, align="center",
          uiOutput("ui_error")
        )
      }else if(rv$plot_type == "track"){
        uiOutput("ui_track")
      }else{
          div(
            column(id="div_surv",
              area_w, align = "left",
              if(!rv$depmapr & rv$plot_type != "gsea"){
                title_div
              },
              if(surv_yn){
                if(!rv$depmapr){
                  conditionalPanel(
                    'input.cox_km != "dens"',
                    plotOutput("cox_plot",height = h_plot)
                  )
                }
              }else if(rv$plot_type == "scatter" | rv$plot_type == "scatter2"){
                plotlyOutput("scatter_plot", height = "585px")
              }else if(rv$plot_type == "snv_stats"){
                plotlyOutput("snv_stats_plot", height = "585px")
              }else if(rv$plot_type == "gsea"){
                uiOutput("ui_gsea")
              }
              ,if((rv$plot_type!="gsea" & !rv$depmapr)|(rv$plot_type=="scatter"|rv$plot_type=="scatter2" & rv$depmapr)){
                div(
                  align = "left",
                  style = "position: absolute; right: 2.5em; top: 1.5em;",
                  div(
                    id = "gear_btn",
                    style="display: inline-block;vertical-align:top;",
                    if(rv$plot_type != "snv_stats" & rv$plot_type != "gsea"){
                      dropdown(
                        inputId = "div_plot_gear",
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
                    ,bsTooltip("download_btn","Click to download plot", placement = "top")
                  )
                )
              }
            )
            ,conditionalPanel(
              'input.plot_type != "snv_stats" & input.plot_type != "gsea"',
              column(
                5, align="left",id="div_ui_stats",
                uiOutput("ui_stats")
              )
            )
          )
      }
      ,conditionalPanel(
        'input.plot_type == "gsea"',
        wellPanel(
          style = paste0("background:#e6f4fc;"),
          HTML(paste0("<h3>Survival-coupled transcriptome analysis by <b>eVITTA</b> ",link_icon("evitta_link","https://tau.cmmt.ubc.ca/eVITTA/"),"</h3>")),
          bsTooltip("evitta_link",HTML("Visit eVITTA webserver"),placement = "top"),
          tags$hr(style="border-color: #939597;margin: 20px;"),
          HTML(paste0("<h4>eVITTA provides modules for differential expression analysis (<b>easyGEO</b>), "
                      ,"functional profiling (<b>easyGSEA</b>), "
                      ,"intersection analysis (<b>easyVizR</b>), and transcriptome pattern visualizations."
                      ,"<br><br>Click button below to start differential expression (DE) analysis by easyGEO.</h4>")),
          br(),
          btn_save_for_geo(id = "btn_jump_to_geo", label = "Start DE analysis and visualization by easyGEO")
          ,htmlOutput("easygeo_iframe")
        )
      )
    )
    ,absolutePanel(
      div(
        style="display: inline-block;vertical-align:top;",
        actionBttn(
          inputId = "up_button", label=NULL,
          icon = icon("angle-double-up"), style="material-circle", color="default", size="md"
        )
      )
      ,div(
        style="display: inline-block;vertical-align:top;",
        actionBttn(
          inputId = "help_button2", label=NULL,
          icon = icon("question"), style="material-circle", color="default", size="md"
        )
      )
      ,tags$script(HTML(
        "document.getElementById('up_button').onclick= function(){
                    document.getElementById('ui_title').scrollIntoView()
                };"
      )),
      right = 20,
      top = 4
    )
  )
})

# ------------ UI for easyGEO -------------
output$easygeo_iframe <- renderUI({
  req(rv$gsea_done == "yes")
  variables_for_geo <- rv$variables_for_geo
  url_easygeo <- paste0('https://tau.cmmt.ubc.ca/eVITTA/easyGEO/',
                "?survival=yes&degss=", variables_for_geo[['degss']], "&coefs=", variables_for_geo[['coefs']])
  # print(url_easygeo)
  div(
    tags$script(HTML(
      "document.getElementById('easygeo_iframe').scrollIntoView();"
    )),
    br(),
    tags$iframe(src = url_easygeo
                     , style="width:100%;",  frameborder="0"
                     ,id="iframe"
                     , height = "800px")
  )
})

# --------- 1. display the survival curve / scatter plot / mutation statistics ---------
observeEvent(input$plot_type,{
  req(!is.null(input$plot_type))
  req(rv$surv_plotted == "plotted")
  if(rv$cox_km == "dens"){if(if_surv(plot_type=input$plot_type)){updateRadioGroupButtons(session,inputId = "cox_km",selected = rv$cox_kmr)}else{rv$cox_km <- rv$cox_kmr}}
  # rv[["res"]] <- NULL
  x <- rv$plot_type <- input$plot_type
  if(x == "scatter"){x <- 1}
  rv[["title"]] <- rv[[paste0("title_",x)]]
  if(!if_surv()){rv$annot_cells_y <- ""}else{
    rv$annot_cells_y <- "yes"
  }
})

observeEvent(list(rv[["title_1"]],rv[["title_all"]]),{
  rv[["title"]] <- rv[[paste0("title_",rv$plot_type)]]
})

extract_plot_data <- function(x){
  # no of cases in each group
  rv[["lels"]] <- rv[[paste0("lels_",x)]]
  
  # the cutoff percentile
  rv[["cutoff"]] <- rv[[paste0("cutoff_",x)]]
  
  
  # extract statistics
  rv[["res"]] <- rv[[paste0("cox_",x)]]
  rv$surv_plotted <- "plotted"
}

output$cox_plot <- renderPlot({
  x <- rv$plot_type
  req(if_surv() & typeof(rv[[paste0("df_",x)]]) == "list")
  withProgress(value = 1, message = "Generating plot ...",{
    extract_plot_data(x)
    # generate survival curve
    plot_surv(rv[["res"]],two_rows=x)
  })
})

# --------- 1a-i. plot parameters -------------
output$plot_gear <- renderUI({
  if(if_surv()){
    if(rv$cox_km == "cox" | rv$cox_km == "km" ){
      fluidRow(
        column(
          12,
          if(rv$tcga | rv$target){
            div(
              # median thresholds
              radioGroupButtons(
                inputId = "ymd",
                label = HTML(paste0("Plot survival in ? ",add_help("ymd_q"))),
                choices = ymd_names,
                selected = rv$ymd,
                size = "sm",
                checkIcon = list(
                  yes = icon("check-square"),
                  no = icon("square-o")
                ),
                direction = "horizontal"
              )
              ,bsTooltip("ymd_q",HTML(paste0("Select the time unit to display on x-axis. Not applicable to DepMap projects"
              ))
              ,placement = "top")
              # fine-tune time intervals
              ,sliderTextInput(
                "ymd_int",
                HTML(paste0("Time inverval on x-axis: ",add_help("ymd_int_q"))),
                choices = rv$ymd_int_range,
                selected = rv$ymd_int
                ,grid=T, force_edges=T
              )
              ,bsTooltip("ymd_int_q",HTML(paste0(
                "Select the # of time units to display on x-axis. Not applicable to DepMap projects"
              )),placement = "top")
            )
          }
          # color scheme
          ,selectInput(
            "palette",
            HTML(paste0("Select color scheme:",add_help("palette_q"))),
            choices = c("Journal of Clinical Oncology palette"="jco","Classic black and red"="br")
            ,selected = rv$palette
          )
          ,bsTooltip("palette_q","Select color palette",placement = "top")
          # confidence intervals
          ,materialSwitch(
            inputId = "confi",
            label = HTML(paste0("<b>Plot confidence intervals?</b> ",add_help("confi_q"))),
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
              label = HTML(paste0("Confidence interval style: "),add_help(paste0("confi_opt_q"))),
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
            ,placement = "top"),
            # median thresholds
            checkboxGroupButtons(
              inputId = "median",
              label = HTML(paste0("Draw line(s) at median survival? ",add_help("median_q"))),
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
            ,materialSwitch(
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
    }
  }else if(rv$plot_type == "scatter" | rv$plot_type == "scatter2" & !is.null(rv$scatter_gender)){
    fluidRow(
      column(
        12,
        selectizeInput(
          "scatter_gender",
          HTML(paste0("Select sex group:",add_help("scatter_gender_q")))
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
            label = HTML(paste0("<b>Color data points by sex?</b>",add_help("scatter_gender_y_q"))),
            value = rv$scatter_gender_y, inline = F, width = "100%",
            status = "danger"
          )
          ,bsTooltip("scatter_gender_y_q",HTML(paste0(
            "If TRUE, color scatter points by sex"
          )),placement = "top")
        # )
        ,materialSwitch(
          inputId = "scatter_log_x",
          label = HTML(paste0("<b>Log2 transform x-axis values?</b>",add_help("scatter_log_x_q"))),
          value = rv$scatter_log_x, inline = F, width = "100%",
          status = "danger"
        )
        ,bsTooltip("scatter_log_x_q",HTML(paste0("If TRUE, values along the x-axis are log2 transformed. Not applicable for DepMap CRISPR-Cas9, RNAi, and drug sensitivity data."))
                   ,placement = "top")
        ,materialSwitch(
          inputId = "scatter_log_y",
          label = HTML(paste0("<b>Log2 transform y-axis values?</b>",add_help("scatter_log_y_q"))),
          value = rv$scatter_log_y, inline = F, width = "100%",
          status = "danger"
        )
        ,bsTooltip("scatter_log_y_q",HTML(paste0("If TRUE, values along the y-axis are log2 transformed. Not applicable for DepMap CRISPR-Cas9, RNAi, and drug sensitivity data."))
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
            materialSwitch(
              inputId = "sm_conf",
              label = HTML(paste0("<b>Draw confidence intervals?</b>",add_help("sm_conf_q"))),
              value = rv$sm_conf, inline = F, width = "100%",
              status = "danger"
            )
            ,bsTooltip("sm_conf_q",HTML(paste0("If TRUE, draw confidence intervals for the regression line"))
                       ,placement = "top"),
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
observeEvent(input$palette,{req(!is.null(input$palette));req(input$palette != "");rv$palette <- input$palette})
observeEvent(input$cox_km,{if(input$cox_km!="dens"){rv$cox_kmr <- input$cox_km};rv$cox_km <- input$cox_km})
observeEvent(input$ymd,{
  req(rv$ymd != input$ymd)
  rv$ymd <- input$ymd
  rv$ymd_int <- rv[[paste0("ymd_int_",input$ymd)]]
  rv$ymd_int_range <- rv[[paste0("ymd_int_range_",input$ymd)]]
},ignoreInit = T)
observeEvent(input$ymd_int,{
  req(rv$ymd_int != input$ymd_int)
  rv$ymd_int <- rv[[paste0("ymd_int_",input$ymd)]] <- input$ymd_int
},ignoreInit = T)
observeEvent(input$median,{rv$median <- input$median},ignoreNULL = F)
observeEvent(input$confi,{rv$confi <- input$confi})
observeEvent(input$sm_conf,{rv$sm_conf <- input$sm_conf})
observeEvent(input$confi_opt,{rv$confi_opt <- input$confi_opt})
observeEvent(input$risk_table,{rv$risk_table <- input$risk_table})
observeEvent(input$cum_table,{rv$cum_table <- input$cum_table})

observeEvent(input$scatter_gender,{rv$scatter_gender <- input$scatter_gender})
observeEvent(input$scatter_gender_y,{rv$scatter_gender_y <- input$scatter_gender_y})
observeEvent(input$scatter_log_x,{rv$scatter_log_x <- input$scatter_log_x})
observeEvent(input$scatter_log_y,{rv$scatter_log_y <- input$scatter_log_y})
observeEvent(input$scatter_lm,{rv$scatter_lm <- input$scatter_lm})

# --------- 1a-ii. plot parameters for density -------------
output$plot_gear_dens <- renderUI({
  fluidRow(
    column(
      12,
      materialSwitch(
        inputId = "dens_fill",
        label = HTML(paste0("<b>Fill density plot?</b>",add_help("dens_fill_q"))),
        value = rv$dens_fill, inline = F, width = "100%",
        status = "danger"
      )
      ,bsTooltip("dens_fill_q",HTML(paste0("If TRUE, fill the density plot with colors"))
                 ,placement = "top")
      ,materialSwitch(
        inputId = "dens_mean",
        label = HTML(paste0("<b>Show mean dependencies</b>",add_help("dens_mean_q"))),
        value = rv$dens_mean, inline = F, width = "100%",
        status = "danger"
      )
      ,bsTooltip("dens_mean_q",HTML(paste0("If TRUE, average dependency scores are plotted for each group"))
                 ,placement = "top")
      # color scheme
      ,selectInput(
        "palette_dens",
        HTML(paste0("Select color scheme:",add_help("palette_dens_q"))),
        choices = c("Journal of Clinical Oncology palette"="jco","Classic black and red"="br")
        ,selected = rv$palette
      )
      ,bsTooltip("palette_dens_q","Select color palette",placement = "top")
    )
  )
})

observeEvent(input$palette_dens,{req(!is.null(input$palette_dens));req(input$palette_dens!="");rv$palette <- input$palette_dens})

# --------- 1b-i. download plot -------------
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

# --------- 1b-ii. download dens plot -------------
output$download_plot_dens <- downloadHandler(
  filename = function(){
    paste0("dpdens_",Sys.time(),".html")
  },
  content = function(file) {
    saveWidget(as_widget(ggplotly(rv[["ggdens"]])), file, selfcontained = TRUE)
  }
)

# --------- 1b-iii. download dens's box plot -------------
output$download_plot_box <- downloadHandler(
  filename = function(){
    paste0("dpbox_",Sys.time(),".html")
  },
  content = function(file) {
    saveWidget(as_widget(ggplotly(rv[["ggbox"]])), file, selfcontained = TRUE)
  }
)

# --------- 2. display the statistics -------------
output$ui_stats <- renderUI({
  # if(rv$depmapr){req(rv$cox_km != "dens")}
  surv_yn <- if_surv(plot_type=input$plot_type)

  if(surv_yn){
    req(!is.null(rv[["res"]]))

    n_lels <- length(rv[["lels"]])
    col_w <- 12 / n_lels
    lel1 <- names(rv[["lels"]])[[n_lels]]
    lel2 <- names(rv[["lels"]])[[1]]

    res <- rv[["res"]][[rv$cox_km]]
    hr <- res[["hr"]]
    p <- res[["p"]]
    p.adj <- res[["p.adj"]]

    if(rv$cox_km == "cox"){
      if(!is.null(p.adj)){p_w <- hr_w <- 4;p_w_r <- T}else{p_w <- hr_w <- 6;p_w_r <- F}
      hf_plot <- ifelse(rv$plot_type == "all","235px","155px")
    }else{
      if(!is.null(p.adj)){p_w <- 6;p_w_r <- T}else{p_w <- 12;p_w_r <- F}
      hf_plot <- "155px"
    }

    if(rv$cox_km == "cox" & rv$plot_type == "all"){
      hr_title <- "HR (hazard ratios)"
      p_title <- "P-values"
      p_adj_title <- "adjusted P-values"
      lel1 <- gsub("_"," and/or ",lel1)
      lel2 <- gsub("_"," and/or ",lel2)
    }else{
      hr_title <- "HR (hazard ratio)"
      p_title <- "P-value"
      p_adj_title <- "adjusted P-value"
    }

    hr_q <- paste0("Only applicable to regression analysis by Cox PH model. HR > 1 indicates that the ",lel1," group have higher risk of death than the ",lel2," group. <i>Vice versa</i>,"
                   ," HR < 1 indicates a lower risk of death for the ",lel1," as compared to the ",lel2)
    stats_title <- paste0("Differential ",firstlower(rv$plot_sstype)," analysis by ",names(surv_methods)[surv_methods == rv$cox_km])
  }else if(rv$plot_type == "scatter" | rv$plot_type == "scatter2"){
    req(!is.null(rv[["res_scatter"]]))
    stats_title <- "Correlation statistics"
    res <- rv[["res_scatter"]]
    hr <- round(res$estimate, 2)
    hr_title <- "Coefficient"
    p <- format(as.numeric(res$p.value), scientific = F, digits = 2)
    p_title <- "P-value"
    p_w <- 6; p_w_r <- F; p.adj <- NULL
  }else{
    stats_title <- ""
  }
  
  req(stats_title != "")
  req(!(rv$depmapr & rv$plot_type != "scatter" & rv$plot_type != "scatter2"))

  column(
    12,style="display: inline-block;vertical-align:top; width: 100%;word-break: break-word;",
    h3(stats_title),
    if(rv$plot_type == "scatter" | rv$plot_type == "scatter2"){
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
        if(surv_yn & rv$cox_km == "cox"){
          column(
            hr_w,
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
            ,rightBorder = p_w_r
          )
        )
        ,if(surv_yn & !is.null(p.adj)){
          column(
            p_w,
            descriptionBlock(
              header = p.adj,
              text = HTML(paste0(p_adj_title,add_help("padj_q")))
              ,rightBorder = F
            )
            ,bsTooltip("padj_q",HTML(padj_q_txt),placement = "bottom")
          )
        }
      )
    )
    ,boxPad(
      color = "gray",
      fluidRow(
        if(surv_yn & !rv$depmapr){
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
        ,if(surv_yn){
          column(
            12,
            conditionalPanel(
              '(input.plot_type == "all" | input.plot_type == "gender") & input.cox_km == "km"',
              selectizeInput(
                "km_mul",
                NULL,
                choices = pairwise_methods
                ,selected = rv[["km_mul"]]
              )
              ,div(
                align="center",
                if(typeof(rv[["res"]][["km"]][["df"]]) == "list"){
                  plotlyOutput("km_hm", height="180px", width = "99.5%")
                }
              )
              ,br()
            ),
            conditionalPanel(
              'input.cox_km == "cox"',
              plotOutput("forest_plot",height = hf_plot)
            ),
            div(
              align="left",
              if(length(res[["stats"]]) == 2){
                renderPrint({print(res[["stats"]][[1]]); cat(paste("\n",sep="\n")); print(res[["stats"]][[2]])})
              }else{
                renderPrint({print(res[["stats"]])})
              }
            )
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

# ------------- 3a. pairwise comparison statistics -------
observeEvent(input$km_mul,{
  req(input$km_mul != "")
  req(input$km_mul != rv$km_mul)
  rv$km_mul <- input$km_mul

  # retrieve df for survival analysis
  df <- rv[[paste0("df_",rv$plot_type)]]

  # # update statistics
  # if(rv$depmapr){
  #   km2 <- pairwise_survdiff(Surv(dependency) ~ level, data = df, p.adjust.method = rv$km_mul)
  # }else{
    km2 <- pairwise_survdiff(Surv(survival_days, censoring_status) ~ level, data = df, p.adjust.method = rv$km_mul)
  # }
  rv[["res"]][[rv$cox_km]][["stats"]][[1]] <- km2
},ignoreInit = T)

# ----------- 3b. pairwise heatmap ----------
plot_heatmap <- function(pvals,mul_methods=rv$km_mul){
  counts <- -log10(pvals)
  counts[is.na(counts)] <- 0
  dat <- expand.grid(y = rownames(counts), x = colnames(counts))
  dat$z <- unlist(as.data.frame(counts),recursive = T)
  pvals <- unlist(as.data.frame(pvals), recursive = T) %>% format(., scientific = F, digits = 2)
  req(length(dat$z)>0)
  
  fig <- plot_ly() %>%
    add_trace(data = dat, x = ~x, y = ~y, z = ~z, type = "heatmap",
              colorscale  = col_scale,zmax = 3,zmin=0,
              colorbar = list(
                title = list(text="-log10(P)", side = "right")
                ,len = 1.2),
              text = pvals,
              hovertemplate = paste('<b>%{x}</b> vs <b>%{y}</b><br>',
                                    'P-value: <b>%{text}</b>'
              )
    ) %>% layout(
      # title = "Pariwise comparisons",
      xaxis = list(title = paste0("Pariwise comparisons adjusted by ",mul_methods), showticklabels = T),
      yaxis = list(title = "", showticklabels = T)
      # ,margin = list(l=200)
    ) %>%
    add_annotations(x = dat$x, y = dat$y,
                    text = as.character(pvals),
                    showarrow = FALSE, xref = 'x', yref = 'y', font=list(color='black')
                    ,ax = 20, ay = -20
    )
  
  fig
}

output$km_hm <- renderPlotly({
  req(length(rv[["res"]][[rv$cox_km]][["stats"]]) == 2)
  pvals <- rv[["res"]][[rv$cox_km]][["stats"]][[1]]$p.value
  req(is.numeric(pvals))
  plot_heatmap(pvals)
})

# ----------- 4[A]. expression-survival days scatter plot ---------------
output$scatter_plot <- renderPlotly({
  withProgress(value = 1,message = "Updating plot ...",{
    if(rv$plot_type == "scatter"){
      # retrieve survival data
      if(typeof(rv[["df_gender"]]) == "list"){
        req(length(rv$scatter_gender) > 0)
        if(rv$depmapr){
          df_survival <- rv[["df_gender"]] %>% dplyr::select(patient_id,dependency,level.y) %>%
            dplyr::filter(level.y %in% rv$scatter_gender)
        }else{
          df_survival <- rv[["df_gender"]] %>% dplyr::select(patient_id,survival_days,level.y) %>%
            dplyr::filter(level.y %in% rv$scatter_gender)
        }
        genders <- df_survival$`level.y`
        df_survival <- df_survival %>%
          dplyr::select(-level.y)
      }else{
        if(rv$depmapr){
          df_survival <- rv[["df_1"]] %>% dplyr::select(patient_id,dependency)
        }else{
          df_survival <- rv[["df_1"]] %>% dplyr::select(patient_id,survival_days)
        }
      }

      df <- rv[["exprs_1"]]

      # the unit, e.g. expression (fpkm)
      exp_unit <- input_mode_name("1")

      df_o <- df <- df_survival %>% inner_join(df, by="patient_id")
      if(ncol(df) == 3){
        gene_name <- colnames(df)[3]
        colnames(df) <- c("patient_id","survival_days","exp")
        exprs <- df$exp
        rv[["gs_no"]] = T; exp_type = exp_unit
        if(rv$scatter_log_y){
          if(input$db_1 == "crispr" | input$db_1 == "rnai" | input$db_1 == "drug"){
            df_y <- df$exp
            ylab <- exp_unit
          }else{
            df_y <- log2(df$exp+1)
            ylab <- paste0("Log2 (",exp_unit," + 1)")
          }
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

      if(rv$depmapr){
        df[["patient_id"]] <- paste0(translate_cells(df[["patient_id"]]),"|",df[["patient_id"]])
        pid <- "Cell line"
        ltitle <- dep_name <- dependency_names()
        # if(rv$project == "DepMap-Drug"){
        #   df[["survival_days"]] <- log2(df[["survival_days"]])
        # }else{
        #   df[["survival_days"]] <- log10(df[["survival_days"]])
        # }
      }else{
        pid <- "Patient ID"
        ltitle <- dep_name <- "Survival days"
      }

      if(rv$scatter_log_x){
        if(rv$depmapr){
          df_x <- df$survival_days
          xlab <- ltitle
        }else{
          df_x <- log2(df$survival_days+1)
          xlab <- "Log2 (survival days + 1)"
        }
      }else{
        df_x <- df$survival_days
        if(rv$depmapr){
          xlab <- dep_name
        }else{
          xlab <- ltitle
        }
      }
      
      # if DepMap, add gene names on to axises
      if(rv$depmapr){
        ylab <- paste0(gene_name,": ",ylab)
        xlab <- paste0(rv$depmap_gener,": ",xlab)
      }

      # calculate correlation
      rv[["res_scatter"]] <- cor.test(df_x, df_y, method = rv$cor_method)

      # convert into ranks in necessary
      if(rv$cor_method == "kendall" | rv$cor_method == "spearman"){
        df_x <- rank(df_x,ties.method = "first")
        df_y <- rank(df_y,ties.method = "first")
        xlab <- paste0("Ranks in ",firstlower(dep_name)); ylab <- paste0("Ranks in ",firstlower(exp_unit))
      }

      # draw the figure
      if(rv$scatter_gender_y & length(rv$scatter_gender)>1){
        fig <- ggplot(df
                      ,aes(x=df_x, y=df_y
                           ,text=paste0(
                             pid,": <b>",.data[["patient_id"]],"</b>\n",
                             dep_name,": <b>",.data[["survival_days"]],"</b>\n",
                             exp_type,": <b>",signif(exprs,digits=3),"</b>"
                           )
                      )) +
          geom_point(aes(color=genders)) + #, shape=genders
          theme_classic() +
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
                             dep_name,": <b>",.data[["survival_days"]],"</b>\n",
                             exp_type,": <b>",signif(exprs,digits=3),"</b>"
                           )
                      )) +
          theme_classic() +
          geom_point(color=col)
      }
      fig <- fig +
        xlab(xlab) +
        ylab(ylab)

      # draw a regression line
      if(rv$scatter_lm){
        fig <- fig + geom_smooth(method=rv$lm_method,fill="#F5DF4D",inherit.aes = F,aes(df_x, df_y),se=rv$sm_conf,size=0.5)
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
        if(input$db_2 == "crispr" | input$db_2 == "rnai" | input$db_2 == "drug"){
          df_y <- df$expb
          ylab <- exp_unitb
        }else{
          df_y <- log2(df$expb+1)
          ylab <- paste0("Log2 (",exp_unitb," + 1)")
        }
      }else{
        df_y <- df$expb
        ylab <- exp_unitb
      }

      # log x, if prompted
      if(rv$scatter_log_x){
        if(input$db_1 == "crispr" | input$db_1 == "rnai" | input$db_1 == "drug"){
          df_x <- df$expa
          xlab <- exp_unita
        }else{
          df_x <- log2(df$expa+1)
          xlab <- paste0("Log2 (",exp_unita," + 1)")
        }
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
          theme_classic() +
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
          theme_classic() +
          geom_point(color=col)
      }

      # rename x- and y- axis
      fig <- fig +
        xlab(xlab) +
        ylab(ylab)

      # draw a regression line
      if(rv$scatter_lm){
        fig <- fig + geom_smooth(method=rv$lm_method,fill="#F5DF4D",inherit.aes = F,aes(df_x, df_y),se=rv$sm_conf,size=0.5)
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
      "Mutated"
    }else{
      "Other"
    }
  }) %>% unname(.)

  dat$Category <- nons_cat
  dat <- dat %>% dplyr::arrange(Category,Frequency)
  # patient cases that have each mutation
  Cases <- lapply(dat$Mutation, function(x){
    xx <- names(muts)[muts == x]
    # convert to CCLE cell ids if depmap
    if(rv$depmapr){xx <- paste0(translate_cells(xx),"|",xx)}
    xx <- xx %>%
      breakvector(.,max=20) %>%
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
    if(rv$depmapr){pname <- "cell lines";pcont <- "contain"}else{pname <- "project(s)";pcont <- "contain(s)"}
    div(
      br(),
      p(style="color:gray;font-size:120%;",paste0("Unable to assess for sex effect. The selected ",pname," only ",pcont," ",rv[["df_gender"]]," data."))
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

# ------ 6. forest plot -------
output$forest_plot <- renderPlot({
  if(rv$plot_type == "all" | rv$plot_type == "gender"){
    maint <- "Hazard ratios"
  }else{
    maint <- "Hazard ratio"
  }
  res <- rv[["res"]]
  req(!is.null(res[["cox"]][["cox_fit"]]))
  ggforest(res[["cox"]][["cox_fit"]], main = maint, data = res[["cox"]][["cox_df"]], fontsize = 0.7)
})

# ------ 7a. stats & heatmap, if dependency ---------
output$depmap_stats <- renderUI({
  req(rv$project != "")
  x <- rv$plot_type
  req(if_surv() & typeof(rv[[paste0("df_",x)]]) == "list")
  extract_plot_data(x)
  rv$annot_cells_y <- "yes"
  # stats names
  x_numeric <- suppressWarnings(!is.na(as.numeric(x)))
  if(x_numeric){
    stats_name <- "Mann-Whitney U test,"
  }else{
    stats_name <- "Kruskal-Wallis rank sum test, overall"
  }
  
  p_title <- paste0(stats_name," P-value = ",rv[["res"]][["p"]])
  if(!is.null(rv[["res"]][["p.adj"]])){p_title <- paste0(p_title,", adjusted P-value",add_help("padj_dp_q")," = ",rv[["res"]][["p.adj"]])}
  
  column(
    12,style="display: inline-block;vertical-align:top; width: 100%;word-break: break-word;",
    boxPad(
      color = "light-blue",
      fluidRow(
        column(
          12,align="center",
          HTML(paste0("<h4>",p_title,"</h4>"))
          ,bsTooltip("padj_dp_q",HTML(padj_q_txt),placement = "bottom")
        )
        ,if(!x_numeric){
          div(
            column(
              12,align="center",
              selectizeInput(
                "km_mul_dp",
                NULL,
                choices = pairwise_methods
                ,selected = rv[["km_mul_dp"]]
              )
            ),
            uiOutput("dm_stats_render")
          )
        }
      )
    )
  )
})

# depmap heatmap
observeEvent(input$km_mul_dp,{
  req(input$km_mul_dp != "")
  rv$km_mul_dp <- input$km_mul_dp
  
  # retrieve df for survival analysis
  df <- rv[[paste0("df_",rv$plot_type)]]
  levl <- length(unique(df$level))
  
  # the height of the heatmap
  if(levl > 3){
    dp_h <- "218px"
  }else if(levl == 3){
    dp_h <- "198px"
  }else if(levl == 2){
    dp_h <- "178px"
  }else if(levl < 2){
    dp_h <- "158px"
  }
  
  # update statistics
  surv_diff <- pairwise.wilcox.test(df$dependency, df$level, p.adjust.method = rv$km_mul_dp)
  rv[["res"]][["fit"]] <- surv_diff
  
  output$dm_stats_render <- renderUI({
    div(
      column(
        6,
        renderPrint({print(surv_diff)})
      ),
      column(
        6,
        plotlyOutput("dp_hm", height=dp_h, width = "100%")
      )
    )
  })
},ignoreInit = F)

output$dp_hm <- renderPlotly({
  pvals <- rv[["res"]][["fit"]]$p.value
  req(is.numeric(pvals))
  dims <- dim(pvals)
  req(dims[1] * dims[2] >= 1)
  plot_heatmap(pvals,mul_methods=rv$km_mul_dp)
})

# ------ 7b. density plot, if dependency ---------
lel_colors <- function(n_lel){
  # adjust colors, if applicable
  if(rv$palette == "br"){
    col_alt <- c("#939597","#F0A1BF") # grey pink
    if(n_lel > 2){
      c("black",col_alt[1:(n_lel-2)],"red")
    }else{
      c("black","red")
    }
  }else{
    pal_jco("default")(n_lel)
  }
}
observeEvent(input$dens_fill,{rv$dens_fill <- input$dens_fill})
observeEvent(input$dens_mean,{rv$dens_mean <- input$dens_mean})
output$dens_plot <- renderPlotly({
  req(rv$depmapr)
  df <- retrieve_dens_df()
  dep_name <- dependency_names()
  n_lel <- length(levels(df$Level))
  c_values <- lel_colors(n_lel)
  
  if(rv$dens_fill){
    p <- ggplot(df, aes(x=.data[[dep_name]], fill=Level)) + geom_density(alpha=0.4) + #, ..scaled..
      scale_fill_manual(values=c_values)
  }else{
    p <- ggplot(df, aes(x=.data[[dep_name]], color=Level)) + geom_density() +
      scale_color_manual(values=c_values)
  }
  
  if(rv$dens_mean){
    mu <- df %>% dplyr::group_by(Level) %>% dplyr::summarise(grp.mean = mean(.data[[dep_name]], na.rm=T))
    p <- p +
      geom_vline(data=mu, aes(xintercept=grp.mean, color=Level), linetype="dashed") +
      scale_color_manual(values=c_values)
  }
  p <- p +
    # geom_histogram(aes(y=..density..), alpha=0.5, position="identity", binwidth=0.02) +
    labs(title=paste0(dep_name," distribution"),x=dep_name, y = "Density") +
    theme_classic()
  rv[["ggdens"]] <- p
  ggplotly(p)
})

# -------- 7c. stats of density plot ---------
output$dens_stats_plot <- renderPlotly({
  req(rv$project != "")
  req(!is.null(rv[["dens_df"]]))
  withProgress(value=1,message = "Generating plots ...",{
    df <- rv[["dens_df"]]
    dep_name <- dependency_names()
    
    pos <- position_jitter(0.1, seed = 2)
    lels_len <- length(levels(df$Level))
    cols <- lel_colors(lels_len)
    p <- ggplot(df, aes(x=Level, y=.data[[dep_name]], color=Level, Line=Cell)) + geom_boxplot() + coord_flip() +
      scale_color_manual(values=cols) +
      geom_jitter(shape=16, position=pos) +
      theme_classic() +
      labs(title="Cell line distribution",x="", y = dep_name)
    
    if(length(input$annot_cells) > 0){
      if(input[["annot_cells"]][1] != ""){
        p <- p + 
          # geom_jitter(position=pos, aes(color = paste0(df$Level,", ",ifelse(Cell %in% input$annot_cells, "Highlighted", "Not highlighted")))) +
          # scale_color_manual(values = c(cols,c(rbind(cols,rep("red",lels_len)))))
          geom_jitter(position=pos, aes(shape = ifelse(Cell %in% input$annot_cells, "Highlighted", "Not highlighted"), size = ifelse(Cell %in% input$annot_cells, "Highlighted", "Not highlighted"))) + #data = subset(df, Cell %in% input$annot_cells)
          scale_shape_manual(values=c(8,16)) + scale_size_manual(values=c(3,1)) +
          labs(color = "", shape = "", size = "")
          #   geom_text(
        #   data = subset(df, Cell %in% input$annot_cells),
        #   aes(label = Cell),
        #   size = 3
        # )
      }
    }
    rv[["ggbox"]] <- p
    ggplotly(p,tooltip = c("Line","x","y"))
  })
})

# highlight cells
observeEvent(list(rv$annot_cells_y,rv[["dens_df"]]),{
  req(rv$annot_cells_y == "yes")
  updateSelectizeInput(session,"annot_cells",choices = rv[["dens_df"]][["Cell"]],server = T)
})

# --------- 8. P-value tracking -------------
output$ui_track <- renderUI({
  req(nrow(rv$quantile_graph >= 1))
  
  
  column(
    width = (12 / rv$variable_nr),
    plotlyOutput("quantile_graph")
  )
})

#Quantile Plot Output
output$quantile_graph <- renderPlotly({
  #Check at lease some rows are in the quantile graph
  req(nrow(rv$quantile_graph >= 1))
   
  #NEW ONE
  fig <- plot_ly(rv$quantile_graph, x = rv$quantile_graph$quantile)
  fig <- fig %>% add_trace(y = ~rv$quantile_graph$p_value,type = 'scatter',#color =~p_value,
                           line = list(color = 'rgb(173,173,173)', width = 2),
                           
                           name = 'P Value',
                           marker=list(
                             color=~p_value,
                             # colorbar=list(
                             #   title='Colorbar'
                             # ),
                             colorscale=col_scale,#'YlOrRd',#custom_colorscale,
                             cmid = 0.5,
                             reversescale =TRUE
                           ),
                           text = rv$quantile_graph$expression,
                           
                           #colors = brewer.pal("YlOrRd "),
                           #rev('YlOrRd'),#brewer.pal(length(rv$quantile_graph$p_value),"YlOrRd "),
                           name = '',mode = 'lines+markers', hovertemplate = paste(
                             #"P value is : %{y:.3f}<br>",
                             "%{y:.3f}<br>",
                             "Quantile(in %) : %{x:.0f}<br>",
                             "ExpressionS : %{text:.3f}<br>"
                           )) %>%
    add_trace(y = rv$quantile_graph$hr, name = 'Hazard Ratio',mode = 'lines+markers',type = 'scatter',
              line = list(color = 'rgb(0,88,155)', width = 2),
              marker=list(
                symbol = 'diamond',
                color=rv$quantile_graph$hr,
                colorscale='RdBu',#'RdBu',#col_scale_hr,
                #'YlOrRd',#custom_colorscale,
                #cmin = min(rv$quantile_graph$hr),
                #cmax = max(rv$quantile_graph$hr)-0.2,
                cmid = 1,
                reversescale =FALSE
              ),#hovertemplate = '',
              yaxis = "y2"
    )%>%
    layout(title = 'P Values of Different Quantiles',
           xaxis = list(title = 'Quantile(%)'),
           yaxis = list (title = 'P Value'),
           yaxis2 = list(overlaying = "y",
                         side = "right",
                         title = "Harzard Ratio"),
           hovermode = "x unified"
    )
  fig

})



