# the function that call rintrojs
call_introjs <- function(file_name) {
  rintrojs::introjs(session, options = list(showStepNumbers=FALSE,
                                            steps = file_name)
  )
}
# the events trigger by pressing the help button on header
observeEvent(input$db_help, {
  call_introjs(rbind(intros$R_pre))
})

observeEvent(input$help_button2, {
  call_introjs(rbind(intros$surv))
})