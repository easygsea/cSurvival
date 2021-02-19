bodyOne <- tabItem(tabName = "one",
    fluidRow(
        box(
            width = 12, status = "danger", 
            column(
                8,
                div(id="project_div",selectizeInput(
                    "project",
                    HTML("Cancer type by project &nbsp")
                    ,choices = ""
                    ,width = "100%"
                ))
                ,br()
            )
            ,column(
                4,align="right",
                div(id="variable_n_div",numericInput(
                    "variable_n",
                    HTML("Number of analysis &nbsp"),
                    value = 1,
                    min = 1, max = 2, step = 1
                ))
            )
            
            ,uiOutput("ui_parameters")
            ,uiOutput("ui_parameters_confirm")
        )
    )
)