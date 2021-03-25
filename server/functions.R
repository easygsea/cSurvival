#======================================================================#
####                       general functions                        ####
#======================================================================#
# assign values to dynamic RVs when initialized
# btn0_rds <- readRDS(paste0(getwd(),"/inc/btn0.rds"))
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
  rv[[paste0("gs_lg_",x,"_search")]] <- 0 # btn0_rds
  rv[[paste0("gs_lg_",x,"_reset")]] <- 0
  rv[[paste0("gs_lgg_",x)]] <- ""
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
  rv[[paste0("snv_method_",x)]] <- "mutect"
}

# update these into rv when selections change
update_all <- function(){
  for(x in 1:rv$variable_n){
    lst <- dyn_list(x)
    updateRV(lst)
  }
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

# req for every element in a list to be both non-null & non-""
req_lst <- function(lst){
  lapply(lst, function(x){
    req(x)
    req(x!="")
  })
}

# req rv not equal to input value
req_diff_rv <- function(namespaces){
  !all(
    sapply(namespaces, function(x){
      rv[[x]] == input[[x]]
    })
  )
}

# req rv larger than or equal to input value, applies to btn
req_diff_rv_btn <- function(namespaces){
  !all(
    sapply(namespaces, function(x){
      rv[[x]] >= input[[x]][1]
    })
  )
}

# req GS filter input content exist
req_filter_on <- function(namespaces, filter="", target="rv", mode="equal"){ # namespace = paste0("gs_lg_",x),
  !all(
    sapply(namespaces, function(x){
      # user's input
      if(target=="rv"){
        if(mode == "equal"){
          rv[[x]] == filter
        }else{
          rv[[x]] != filter
        }
      }else if(target == "input"){
        if(!is.null(input[[x]])){
          if(mode == "equal"){
            input[[x]] == filter
          }else{
            input[[x]] != filter
          }
        }else{
          FALSE
        }
      }
    })
  )
}

# # specific function to handle the bug when second panel is initiated but not responding to UI update
# check_array <- function(lst){
#   # lst_u <- lst %>% unlist() %>% unique()
#   # req(lst_u != "")
#   n <- rv$variable_n
#   if(n>1){
#     # req(!is.null(lst[[2]]))
#     req(lst[[2]] != "")
#   }
#   if(lst[1] == ""){array <- 2}else{array <- 1:n}
#   return(array)
# }

# ------- return all GSs when a db is selected -------
update_gs_by_db <- function(x){
  gs_db_id <- paste0("gs_db_",x)
  gs_lib_id <- paste0("gs_l_",x)
  
  db <- rv[[gs_db_id]] <- isolate(input[[gs_db_id]])
  
  if(db != ""){
    gmt_path <- retrieve_gmt_path(db)
    gmt <- gmtPathways(gmt_path)
    
    # update variables
    rv[[paste0("gmts",x)]] <- rv[[paste0("gmts_tmp",x)]] <- gmt
    
    # update placeholder
    rv[[paste0("gs_placeholder",x)]] <- sprintf('(Total n=%s) Type to search ...',length(gmt))
    
    # update gene set UI
    updateSelectizeInput(
      session,
      gs_lib_id
      ,choices = names(gmt)
      ,selected=rv[[gs_lib_id]]
      ,options = list(
        # `live-search` = TRUE,
        placeholder = rv[[paste0("gs_placeholder",x)]]
        ,onInitialize = I(sprintf('function() { this.setValue("%s"); }',rv[[gs_db_id]]))
      )
    )
  }
  
}