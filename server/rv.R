rv <- reactiveValues(
  variable_n_reached=0
  ,variable_n = 1
  ,verbTxtStyle1 = "box-shadow: 0 0 .2em #F5DF4D;color: red;"
  ,verbTxtStyle2 = ""
  
  ,tcga=T # check if TCGA project
  
  ,cox_km = "cox"
  ,median = NULL
  ,confi = T,confi_opt = "ribbon"
  ,risk_table = T,cum_table=T
)
