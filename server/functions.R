# assign values to dynamic RVs when initialized
init_rv <- function(x){
  # category to analyze
  rv[[paste0("cat_",x)]] <- "g"
  # type of db to analyze
  rv[[paste0("db_",x)]] <- "rna"
  # gene to analyze
  rv[[paste0("g_",x)]] <- ""
  # gs mode
  rv[[paste0("gs_mode_",x)]] <- "lib"
  # gs db to analyze
  rv[[paste0("gs_db_",x)]] <- ""
  # gs lib to analyze
  rv[[paste0("gs_l_",x)]] <- ""
  # gs genes in feedback
  rv[[paste0("gs_lgs_",x)]] <- ""
  # gs gene to search
  rv[[paste0("gs_lg_",x)]] <- ""
  # manual gene input
  rv[[paste0("gs_m_",x)]] <- ""
  # feedback on manual gene input
  rv[[paste0("gs_mg_",x)]] <- ""
  # lower bound for quantile loop
  rv[[paste0("lower_",x)]] <- .15
  # upper bound for quantile loop
  rv[[paste0("upper_",x)]] <- .85
  # step size
  rv[[paste0("step_",x)]] <- .01
}

# update RV according to changes in input
updateRV <- function(id_list){
  for (id in id_list){
    x <- isolate(input[[id]])
    if(is.null(x)==F){ 
      rv[[id]] <- x
    }
  }
}

