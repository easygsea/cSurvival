plot_ui <- function(n){
  # automatically adjust column width according to # of analysis selected
  if(n == 1){col_w <- 12}else{col_w <- 6}
  
  # create the UI list
  ui <- lapply(1:n, function(x){
    # color number for the wellpanel
    bcol_n <- x %% 2; if(bcol_n == 0){bcol_n <- 2}
    # category to analyze
    cat_id <- paste0("cat_",x)
    # gene to analyze
    g_ui_id <- paste0("g_",x)
    gs_ui_id <- paste0("gs_",x)
    # the UI
    column(
      col_w,
      # tags$hr(style="border: .5px solid lightgrey; margin-top: 0.5em; margin-bottom: 0.5em;"),
      wellPanel(
        style = paste0("background-color: ", bcols[bcol_n], "; border: .5px solid #fff;"),
        prettyRadioButtons(
          cat_id,
          paste0("Analysis #",x),
          choices = c("Gene" = "g", "Gene set" = "gs")
          ,selected = "g"
          ,status = "warning"
          ,icon = icon("check")
          ,shape = "curve"
          # ,plain = T
          # ,outline = T
          ,animation = "jelly"
          ,inline = T
        )
        ,conditionalPanel(
          condition = sprintf("input.%s=='g'", cat_id),
          selectizeInput(
            g_ui_id,
            paste0(x,".1. Select your gene of interest:")
            ,choices=c()
          )
        )
        ,conditionalPanel(
          condition = sprintf("input.%s=='gs'", cat_id),
          selectizeInput(
            gs_ui_id,
            paste0(x,".1. Select your gene set of interest:")
            ,choices=c()
          )
        )
      )
      
    )
    
  })
  
  # create the UI
  do.call(tagList, ui)
}