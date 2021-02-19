# sets max upload size to 50MB, modal appears when it exceeds 10MB,
# see 1.1.server-run.R line 325 for more details
options(shiny.maxRequestSize=50*1024^2) 

server <- function(input, output, session) {
    waiter_hide() # will hide *on_load waiter
    
  source("server/functions.R", local = TRUE)
  source("server/1.server-one.R", local = TRUE)
    
}
