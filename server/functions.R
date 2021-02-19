# update RV according to changes in input
updateRV <- function(id_list){
  for (id in id_list){
    if(is.null(input[[id]])==F){ 
      rv[[id]] <- input[[id]]
    }
  }
}