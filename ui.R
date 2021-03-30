source("ui/1.one.R")
source("ui/css_addons.R")

sidebar <- dashboardSidebar(
    disable = T
)

loadMsg = "cSurival"

shinyUI(
    dashboardPage(
        
        title="cSurival",

        dashboardHeader(title = div(id="ui_title",HTML("<b>nSurvival</b>: multivariate cancer survival analysis"))
                        ,titleWidth = "100%"
        )
        # skin = "black",
        ,sidebar
        ,dashboardBody(
            tags$head(
                tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
            ),
            # theme = shinytheme("flatly"),
            use_waiter(), # dependencies
            waiter_show_on_load(tagList(spin_three_bounce(),h4(loadMsg)), color = "#00589b"), # shows before anything else
            disconnectMessage(text = "Your session has timed out. Please refresh page and start again. For bug report, email us at jcheng@cmmt.ubc.ca. Thank you for your support."),

            useShinyalert(),  # Set up shinyalert
            rintrojs::introjsUI(), # introjs
            useShinyjs(), # Set up shinyjs
            
            # apply specific css adjustments additionally
            css_addons,
            
            bodyOne
            
            
            # ,tags$footer(HTML("<b>Taubert Lab</b> | BC Children's Hospital Research Institute | Centre for Molecular Medicine and Therapeutics | University of British Columbia. 2019-2020. All Rights Reserved."),
            #             align = "left", style = "
            #   position:absolute;
            #   bottom:0;
            #   width:100%;
            #   height:30px;
            #   color: white;
            #   padding: 5px;
            #   background-color: #3179ae;
            #   z-index: 1000;")


        )

    )
)
