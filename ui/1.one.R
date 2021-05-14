# the panel to confirm project(s) selection
confirm_panel <- conditionalPanel(
  'input.project != "" & !output.projectStatus',
  actionBttn(
    "confirm_project",
    tags$b("Confirm selection!")
    ,block = T,style = "simple",color = "warning",size="sm"
  )
  ,tags$style(type='text/css', "#confirm_project { margin-top: 43.5px; height: 33.5px;}"),
)

# the panel to reset project(s) selection
reset_panel <- conditionalPanel(
  'output.projectStatus',
  actionButton(
    "reset_project",
    strong("Reset selection")
    ,width = "100%"
  )
  ,tags$style(type='text/css', "#reset_project { margin-top: 43.5px;}"),
)

# assemble the UI
bodyOne <- tabItem(tabName = "one",
    fluidRow(
        box(
            width = 12, status = "danger",
            column(
              12,
              fluidRow(
                column(
                  6,
                  selectizeInput(
                    "project",
                    HTML(paste0("<h4><b>To start, select project(s) to analyze:</b>",add_help("project_q"),"</h4>"))
                    ,choices = projects[grepl("TCGA|TARGET|DepMap",names(projects))]
                    ,width = "100%"
                    ,multiple = T
                    ,options = list(
                      `live-search` = TRUE,
                      placeholder = "Type to search ..."
                      ,onInitialize = I(sprintf('function() { this.setValue(%s); }',"")) #['TCGA-LUAD','TCGA-LUSC']
                    )
                  )
                  ,bsTooltip("project_q",HTML(paste0(
                    "Select a project(s) from DepMap, TARGET or TCGA. Multiple selections (maximum 3) are allowed for TARGET."
                    ," Click button to the right to confirm your selection and proceed."
                  )),placement = "right")
                )
                ,column(
                  2,
                  confirm_panel
                  ,reset_panel
                  # ,tags$style(type='text/css', "#variable_n { margin-top: 10px;}"),
                )
                ,column(
                  2,
                  uiOutput("ui_censortime")
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
                    ," 2 to analyze interactions and relationships between genes, loci and/or gene sets."
                  )),placement = "right")
                )
              )
            ,br()
            
            ,fluidRow(
              # TCGA only disease-free survival and progression-free survival
              uiOutput("tcga_pars")
              # DepMap only cell line selection
              ,uiOutput("depmap_pars")
            )
            
            # control widgets for individual analysis
            ,fluidRow(
              uiOutput("ui_parameters")
            )
            
            # the confirm button
            ,fluidRow(
              uiOutput("ui_parameters_confirm")
            )
            ,br()
          )
        )
    )

    ,fluidRow(
      uiOutput("ui_results")
    )
)
