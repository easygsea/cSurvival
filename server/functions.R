plot_ui <- function(n){
  if(n == 1){col_w <- 12}else{col_w <- 6}
  ui <- lapply(1:n, function(x){
    bcol_n <- x %% 2
    if(bcol_n == 0){bcol_n <- 2}
    column(
      col_w,
      # tags$hr(style="border: .5px solid lightgrey; margin-top: 0.5em; margin-bottom: 0.5em;"),
      wellPanel(
        style = paste0("background-color: ", bcols[bcol_n], "; border: .5px solid #fff;"),
        prettyRadioButtons(
          paste0("cat_",x),
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
      )
      
    )
    
  })
  do.call(tagList, ui)
}