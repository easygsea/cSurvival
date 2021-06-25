#======================================================================#
####                     clear RDS when prompted                    ####
#======================================================================#
clear_rds <- function(){
  for(i in seq_along(isolate(rv$variables_for_geo))){
    unlink(
      paste0(surv_dir, isolate(rv$variables_for_geo[[i]]), ".rds"))
  }
}

# change heximal colors to 90% transparency
addalpha <- function(colors, alpha=0.35) {
  r <- col2rgb(colors, alpha=T)
  # Apply alpha
  r[4,] <- alpha*255
  r <- r/255.0
  return(rgb(r[1,], r[2,], r[3,], r[4,]))
}