css_addons <- 
  tags$head(
    tags$style(HTML(paste0(
      "
      #variable_n_div label {display: table-cell; text-align: center; vertical-align: middle;}
      #variable_n_div .form-group { display: table-row;}
      
      #project_div label {display: table-cell; text-align: center; vertical-align: middle;}
      #project_div .form-group { display: table-row;}
      
      #evitta {font-size:120%;}
      
      #confirm_project {display: inline-block;}
      #reset_project {display: inline-block;}
      
      #censor_time {height: 50%;}
      "
    )))
  )