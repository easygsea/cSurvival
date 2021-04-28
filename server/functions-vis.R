#======================================================================#
####                     clear RDS when prompted                      ####
#======================================================================#
clear_rds <- function(){
  for(i in seq_along(isolate(rv$variables_for_geo))){
    unlink(
      paste0(surv_dir, isolate(rv$variables_for_geo[[i]]), ".rds"))
  }
}