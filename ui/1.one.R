bodyOne <- tabItem(tabName = "one",
    fluidRow(
        box(
            width = 12, status = "danger", 
            column(
              12,
              fluidRow(
                column(
                  8,
                  selectizeInput(
                    "project",
                    h4(strong("Cancer type by project"))
                    ,choices = projects[grepl("TCGA|TARGET",names(projects))]
                    ,width = "100%"
                    ,options = list(
                      `live-search` = TRUE,
                      placeholder = "Type to search ..."
                      ,onInitialize = I(sprintf('function() { this.setValue("%s"); }',"TCGA-BRCA"))
                    )
                  )
                )
                ,column(
                  2,
                  conditionalPanel(
                    'input.project != ""',
                    actionButton(
                      "reset_project",
                      "Reset project"
                      ,width = "100%"
                    )
                    ,tags$style(type='text/css', "#reset_project { margin-top: 45px;}"),
                  )
                )
                ,column(
                  2,#align="right",
                  numericInput(
                    "variable_n",
                    h4(strong("No. of analysis")),
                    value = 1,
                    min = 1, max = 2, step = 1
                  )
                  # ,tags$style(type='text/css', "#variable_n { margin-top: 10px;}"),
                )
              )
              
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