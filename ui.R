source("ui/1.one.R")
source("ui/css_addons.R")

sidebar <- dashboardSidebar(
    disable = T
)

# loading message
loadMsg = "cSurival"

# home buttons style
db_style <- "color: #fff; background-color: transparent; border-color: #fff; margin-top:8px; margin-right:8px; border-radius:2rem; border:0.125rem solid #fff;"

# assemble the UI
shinyUI(
    dashboardPage(
        
        title="cSurival",

        dashboardHeader(title = div(id="ui_title",align="left",HTML("<span style='color:#F5FBEF;'><b>cSurvival</b><sup><img src='android-chrome-512x512.png' height='13px' title='Logo by Jiaming (Caitlyn) Xu' style='margin-top: -3px;'></img>v1.0.0</sup></span>: a mechanistic cancer survival database"))
                        ,titleWidth = "80%"
                        ,tags$li(class = "dropdown", actionButton("db_download", NULL,icon("database"),style=db_style))
                        ,tags$li(class = "dropdown", actionButton("db_help", NULL,icon("question"),style=db_style))
                        ,tags$li(class = "dropdown", actionButton("db_demo", NULL,icon=icon("play"),style=db_style))
        )
        # skin = "black",
        ,sidebar
        ,dashboardBody(
            tags$head(
                includeHTML(("google-analytics.html")),
                tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
                ,tags$link(rel = "shortcut icon", href = "favicon.ico")
                ,tags$link(rel = "apple-touch-icon", sizes = "180x180", href = "apple-touch-icon.png")
                ,tags$link(rel = "android-chrome-icon", sizes = "512x512", href = "android-chrome-512x512.png")
                ,tags$link(rel = "android-chrome-icon", sizes = "192x192", href = "android-chrome-192x192.png")
                ,tags$link(rel = "icon", type = "image/png", sizes = "32x32", href = "favicon-32x32.png")
                ,tags$link(rel = "icon", type = "image/png", sizes = "16x16", href = "favicon-16x16.png")
            ),
            # theme = shinytheme("flatly"),
            use_waiter(), # dependencies
            waiter_show_on_load(tagList(spin_orbiter(),h4(loadMsg)), color = "#2D4059"), # shows before anything else
            disconnectMessage(text = "Your session has timed out. Please refresh page and start again. For bug report, email us at jcheng@cmmt.ubc.ca. Thank you for your support."),

            useShinyalert(),  # Set up shinyalert
            rintrojs::introjsUI(), # introjs
            useShinyjs(), # Set up shinyjs
            
            # apply specific css adjustments additionally
            css_addons,
            
            bodyOne
            
            ,bsTooltip("db_download","Download source data")
            ,bsTooltip("db_help","Help")
            ,bsTooltip("db_demo","Example runs")
            
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
