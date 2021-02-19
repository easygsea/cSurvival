css_addons <- 
  tags$head(
    tags$style(HTML(paste0(
      # file upload in GMT area
      "
      #variable_n_div label {display: table-cell; text-align: center; vertical-align: middle;}
      #variable_n_div .form-group { display: table-row;}
      
      #project_div label {display: table-cell; text-align: center; vertical-align: middle;}
      #project_div .form-group { display: table-row;}
      "
    )))
  )