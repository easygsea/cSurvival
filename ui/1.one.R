bodyOne <- tabItem(tabName = "one",
    fluidRow(
        box(
            width = 12, status = "danger",
            column(
              12,
              fluidRow(
                column(
                  10,
                  selectizeInput(
                    "project",
                    HTML(paste0("<h4><b>To start, select project(s) to analyze:</b>",add_help("project_q"),"</h4>"))
                    ,choices = projects[grepl("TCGA|TARGET|DepMap",names(projects))]
                    ,width = "100%"
                    ,multiple = T
                    ,options = list(
                      `live-search` = TRUE,
                      placeholder = "Type to search ..."
                      ,onInitialize = I(sprintf('function() { this.setValue(%s); }',"['TCGA-LUAD','TCGA-LUSC']"))
                    )
                  )
                  ,bsTooltip("project_q",HTML(paste0(
                    "Select a project(s) from DepMap, TARGET or TCGA. Multiple selections (maximum 5) are allowed for pan-cancer analysis."
                    ," Click button below to confirm your selection and proceed."
                  )),placement = "right")
                )
                ,column(
                  2,#align="right",
                  numericInput(
                    "variable_n",
                    HTML(paste0("<h4><b>No. of analysis:</b>",add_help("variable_n_q"),"</h4>")),
                    value = 1,
                    min = 1, max = 2, step = 1
                  )
                  # ,tags$style(type='text/css', "#variable_n { margin-top: 10px;}"),
                  ,bsTooltip("variable_n_q",HTML(paste0(
                    "1 to analyze a single gene, locus, or gene set."
                    ," 2 to analyze interactions between genes, loci and/or gene sets."
                  )),placement = "right")
                )
                ,column(
                  12,
                  conditionalPanel(
                    'input.project != "" & !output.projectStatus',
                    actionBttn(
                      "confirm_project",
                      "Click to confirm project(s) selection to proceed!"
                      ,block = T,style = "simple",color = "warning",size="md"
                    )
                    # ,tags$style(type='text/css', "#confirm_project { margin-top: 43.5px;}"),
                  )
                  ,conditionalPanel(
                    'output.projectStatus',
                    actionButton(
                      "reset_project",
                      strong("Click to reset project(s) selection")
                      ,width = "100%"
                    )
                    # ,tags$style(type='text/css', "#reset_project { margin-top: 43.5px;}"),
                  )
                )
              )
            ,br()

              ,fluidRow(
                uiOutput("ui_parameters")
              )
              ,fluidRow(
                column(
                  12, align="center"
                  ,uiOutput("ui_parameters_confirm")
                  ,add_gear("par_gear")
                )
              )
            )

        )
    )

    ,fluidRow(
      uiOutput("ui_results")
    )
)
