library(stringr)

#* only get data on list
#* @filter data_check
function(req){
  cat(as.character(Sys.time()), "-",
      req$QUERY_STRING, "\n")
  plumber::forward()
}
#* Get the data inside of project_data
#* @param file:[str] subproject name  + file name
#* @serializer contentType list(type="text/csv")
#* @get /project_data/<file>
function(file = "DepMap_ccle.csv", res){
  #Set work directory to current file
  setwd("/Applications/Codes/cSurvival_dev")
  #!!! subject to change
  subproject = unlist(str_split(file, "_", n = 2))[1]
  filename = unlist(str_split(file, "_", n = 2))[2]
  filename <- paste0("./www/project_data/",subproject,"/",filename)
  plumber::include_file(filename, res, "text/csv")
}

# unlist(strsplit("DepMap+ccle_df.csv","+", fixed = TRUE))[2]




# #* Get the data inside of project_data
# #* @param subproject:[str] subproject folder name
# #* @param file:[str] file name
# #* @serializer contentType list(type="text/csv")
# #* @get /project_data/<subproject>/<file>
# function(subproject = "DepMap",file = "ccle.csv", res){
#   #Set work directory to current file
#   setwd("/Applications/Codes/cSurvival_dev")
#   #!!! subject to change
#   filename <- paste0("./www/project_data/",subproject,"/",file)
#   # data <- read.csv(filename)
#   res$setHeader("foo", "bar")
#   plumber::include_file(filename, res, "text/csv")
#   #plumber::as_attachment(data, paste0(subproject,"_",file))
# }
