#======================================================================#
####                       General functions                        ####
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
  rv[[paste0("gs_mg_",x)]] <- ""
  rv[[paste0("add_btn_",x)]] <- 0
  # feedback on manual gene input
  rv[[paste0("gs_mg_",x)]] <- ""
  # lower bound for quantile loop
  rv[[paste0("lower_",x)]] <- .15
  # upper bound for quantile loop
  rv[[paste0("upper_",x)]] <- .85
  # step size
  rv[[paste0("step_",x)]] <- .01
  # parameters for SNV mutation analysis
  rv[[paste0("snv_method_",x)]] <- "mutect"
  rv[[paste0("nonsynonymous_",x)]] <- variant_types_non
  # rv[[paste0("synonymous_",x)]] <- variant_types_syn
  rv[[paste0("iter_",x)]] <- "iter"
  rv[[paste0("clow_",x)]] <- 50
  rv[[paste0("cnv_par_",x)]] <- "auto"
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

# pass value rv if input is.null
ifelse_rv <- function(id){
  if(is.null(input[[id]])){
    rv[[id]]
  }else{
    input[[id]]
  }
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

#======================================================================#
####                       Data handling                        ####
#======================================================================#
# rbind a list of dfs by common columns only
rbind_common <- function(df_list){
  Reduce(function(df_1,df_2){
    common_genes <- intersect(colnames(df_1), colnames(df_2))
    df_1 <- select(df_1, all_of(common_genes))
    df_2 <- select(df_2, all_of(common_genes))
    rbindlist(list(df_1, df_2))
  }, df_list)
}

# break vectors if too long
breakvector <- function(x, max=60){
  if(length(x)>max){
    c(x[1:max],"...")
  }else{
    x
  }
}

# add line breaks into a string
addlinebreaks <- function(x, max=50, lbtype="<br>"){
  x = gsub(paste0('(.{1,',max,'})(\\s|$)'), paste0('\\1',lbtype), x)
  return(x)
}

# call the data type
call_datatype <- function(x){
  ddd <- c(data_types(),data_types_gs)
  names(ddd)[match(input[[paste0("db_",x)]], ddd)]
}

call_datatype_from_rv <- function(x){
  ddd <- c(data_types(),data_types_gs)
  names(ddd)[match(x, ddd)]
}

# return all GSs when a db is selected
update_gs_by_db <- function(x, mode="nil"){
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
    
    if(mode == "nil"){gg=""}else{gg=rv[[gs_lib_id]]}
    # update gene set UI
    updateSelectizeInput(
      session,
      gs_lib_id
      ,choices = names(gmt)
      ,selected=gg
      ,options = list(
        # `live-search` = TRUE,
        placeholder = rv[[paste0("gs_placeholder",x)]]
        ,onInitialize = I(sprintf('function() { this.setValue("%s"); }',rv[[gs_db_id]]))
      )
    )
  }
}

# retrieve genes from a project
retrieve_genes <- function(x){
  db_id <- paste0("db_",x)
  method <- ifelse(is.null(input[[paste0("snv_method_",x)]]),"mutect",input[[paste0("snv_method_",x)]])

  if(is.null(input[[db_id]])){
    infiles <- paste0(rv$indir,"df_gene_scale.csv")
  }else if(input[[db_id]] == "rna"){
    infiles <- paste0(rv$indir,"df_gene_scale.csv")
  }else if(input[[db_id]] == "snv"){
    infiles <- paste0(rv$indir,"df_snv_class_",method,".csv")
  }else if(input[[db_id]] == "cnv"){
    infiles <- paste0(rv$indir,"df_cnv.csv")
  }else if(input[[db_id]] == "mir"){
    infiles <- paste0(rv$indir,"df_mir_scale.csv")
  }
  
  l <- lapply(infiles, function(x){
    info <- fread(x,sep=",",header=T,nrows = 0)
    return(info)
  })
  
  # df_gene <- rbindlist(l, fill = T, use.names = T)
  df_gene <- rbind_common(l)
  
  # the genes
  genes <- names(df_gene) %>% .[-1]
  
  # save into rv$snv_genes
  if(input[[db_id]] == "snv"){
    rv[[paste0("snv_genes_",x)]] <- lapply(l, function(x){
      names(x) %>% .[-1]
    })
  }
  
  # return the genes
  return(genes)
}

# update genes in the UI accordingly
update_genes_ui <- function(opt="hi"){
  lapply(1:rv$variable_n, function(x){
    rv[[paste0("genes",x)]] <- retrieve_genes(x)
    
    g_ui_id <- paste0("g_",x)
    if(opt == "nil"){gg <- ""}else{gg <- rv[[g_ui_id]]}
    
    updateSelectizeInput(
      session,
      g_ui_id,
      choices = rv[[paste0("genes",x)]]
      ,selected = gg
      ,server = TRUE
      ,options = list(
        placeholder = 'Type to search ...'
      )
    )
  })
}
