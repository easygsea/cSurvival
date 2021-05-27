# sets max upload size to 50MB, modal appears when it exceeds 10MB,
# see 1.1.server-run.R line 325 for more details
options(shiny.maxRequestSize=50*1024^2) 

server <- function(input, output, session) {
  waiter_hide() # will hide *on_load waiter
    
  source("server/rv.R", local = TRUE)
  source("server/reactives.R", local = TRUE)
  source("server/variables.R", local = TRUE)
  source("server/functions.R", local = TRUE)
  source("server/functions-ui.R", local = TRUE)
  source("server/functions-vis.R", local = TRUE)
  source("server/functions-calculate.R", local = TRUE)
  source("server/header.R", local = TRUE)
  source("server/help.R", local = TRUE)
  source("server/1.server-one.R", local = TRUE)
  source("server/2.server-calculate.R", local = TRUE)
  source("server/2.server-UI.R", local = TRUE)
  source("server/3.server-eVITTA.R", local = TRUE)
  
  # delete the cSurvival variables we have save for easygeo
  onStop(fun = function(){
    clear_rds()
  })
}
