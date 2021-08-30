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
  rv[[paste0("lower_",x)]] <- dmin
  # upper bound for quantile loop
  rv[[paste0("upper_",x)]] <- dmax
  # step size
  rv[[paste0("step_",x)]] <- dstep
  # parameters for SNV mutation analysis
  rv[[paste0("snv_method_",x)]] <- "mutect"
  rv[[paste0("nonsynonymous_",x)]] <- variant_types_non
  rv[[paste0("synonymous_",x)]] <- variant_types_syn
  rv[[paste0("iter_",x)]] <- "iter"
  rv[[paste0("clow_",x)]] <- 50
  rv[[paste0("cnv_par_",x)]] <- "auto"
  rv[[paste0("snv_uni_",x)]] <- "int"
  rv[[paste0("todefault",x)]] <- 0
  rv[[paste0("gnorm_",x)]] <- "none"
  rv[[paste0("gnorm_g_",x)]] <- ""
  rv[[paste0("gnorm_gs_db_",x)]] <- ""
  rv[[paste0("gnorm_gs_lib_",x)]] <- ""
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

# update RV only when acceptable numeric range
updateRV_numeric <- function(id, min, max){
  id_value <- ifelse_rv(id)
  update_id_value <- function(){updateNumericInput(session,id,value = rv[[id]])}
  if(!is.na(id_value)){if(id_value != rv[[id]] & id_value >= min & id_value <= max){rv[[id]] <- id_value}else{update_id_value()}}else{update_id_value()}
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
  suppressWarnings(
    !all(
      sapply(namespaces, function(x){
        sapply(seq_along(rv[[x]]), function(i){
          rv[[x]][i] == input[[x]][i]
        })
      }) %>% unlist(.)
    )
  )
}

# req rv larger than or equal to input value, applies to btn
req_diff_rv_btn <- function(namespaces){
  suppressWarnings(
    !all(
      sapply(namespaces, function(x){
        rv[[x]] >= input[[x]][1]
      })
    )
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

# pass value rv if input is.na
ifelse_rv_na <- function(id){
  if(is.na(input[[id]])){
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

# waiting message for withProgress if too long
wait_msg <- function(msg){
  paste0(
    msg,
    " This might take a while. Please wait a minute. Thank you."
  )
}

vector_names <- function(x, vector){
  names(vector)[vector == x]
}

#======================================================================#
####                       Data handling                        ####
#======================================================================#
# lower case first letter
firstlower <- function(x) {
  substr(x, 1, 1) <- tolower(substr(x, 1, 1))
  x
}

# function to remove all na per column
not_all_na <- function(x) any(!is.na(x))
# function to remove all na per column
not_any_na <- function(x) all(!is.na(x))
# paste without NA
paste_na <- function(...,sep="\\|") {
  L <- list(...)
  L <- lapply(L,function(x) {x[is.na(x)] <- ""; x})
  ret <-gsub(paste0("(^",sep,"|",sep,"$)"),"",
             gsub(paste0(sep,sep),sep,
                  do.call(paste,c(L,list(sep=sep)))))
  is.na(ret) <- ret==""
  ret
}

# find common mutations
common_mut <- function(row, mode="int"){
  muts <- sapply(row, function(x){
    strsplit(x, "\\|")
  })
  if(mode == "int"){
    muts <- Reduce(intersect, muts) 
  }else if(mode == "uni"){
    muts <- Reduce(union, muts)
  }
  muts <- muts[!is.na(muts)]
  if(identical(muts,character(0))){
    return(NA)
  }else{
    return(paste0(muts,collapse = "|"))
  }
}

# rbind a list of dfs by common columns only
rbind_common <- function(df_list){
  Reduce(function(df_1,df_2){
    common_genes <- intersect(colnames(df_1), colnames(df_2))
    df_1 <- dplyr::select(df_1, all_of(common_genes))
    df_2 <- dplyr::select(df_2, all_of(common_genes))
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

# check input data type/mode
input_mode <- function(x){
  cat_id <- paste0("cat_",x)
  db_id <- paste0("db_",x)
  gs_mode_id <- paste0("gs_mode_",x)
  ifelse(input[[cat_id]] == "g", input[[db_id]], input[[gs_mode_id]])
}

input_mode_name <- function(x, if_crispr=F){
  inmode <- input_mode(x)
  # special mode for normalized counts
  cat_id <- paste0("cat_",x)
  db_id <- paste0("db_",x)
  gs_mode_id <- paste0("gs_mode_",x)
  g_ui_norm_id <- paste0("gnorm_",x)
  if(input[[g_ui_norm_id]] != "none" & ((input[[cat_id]] == "g" & input[[db_id]] == "rna")| (input[[cat_id]] == "gs" & input[[gs_mode_id]] == "lib"))){
    inmode <- input[[g_ui_norm_id]]
  }
  y <- names(input_mode_names)[input_mode_names == inmode]
  if(inmode == "rna"){
    if(rv$tcga){
      y <- gsub("UQ-FPKM","upper quartile normalized RSEM",y)
    }else if(rv$depmap){
      y <- gsub("UQ-FPKM","TPM",y)
    }
  }
  if(if_crispr){
    if(rv$depmapr){
      y <- names(input_mode_names)[input_mode_names == rv$project]
    }
  }
  return(y)
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
update_gs_by_db <- function(x, mode="nil", gs_db_id = paste0("gs_db_",x), gs_lib_id = paste0("gs_l_",x), prefix=""){
  db <- rv[[gs_db_id]] <- isolate(input[[gs_db_id]])
  
  if(db != ""){
    gmt_path <- retrieve_gmt_path(db)
    gmt <- gmtPathways(gmt_path)
    
    # update variables
    rv[[paste0(prefix,"gmts",x)]] <- gmt
    if(prefix == ""){
      rv[[paste0(prefix,"gmts_tmp",x)]] <- gmt
    }
    
    # update placeholder
    placeholder_id <- paste0(prefix,"gs_placeholder",x)
    rv[[placeholder_id]] <- sprintf('(Total n=%s) Type to search ...',length(gmt))
    
    if(mode == "nil"){gg=""}else{gg=rv[[gs_lib_id]]}
    # update gene set UI
    updateSelectizeInput(
      session,
      gs_lib_id
      ,choices = names(gmt)
      ,selected=gg
      ,server = T
      ,options = list(
        # `live-search` = TRUE,
        placeholder = rv[[placeholder_id]]
        ,onInitialize = I(sprintf('function() { this.setValue("%s"); }',rv[[gs_db_id]]))
      )
    )
  }
}

# retrieve genes from a project
retrieve_genes <- function(x){
  db_id <- paste0("db_",x)
  snv_id <- paste0("snv_method_",x)
  # if(is.null(input[[snv_id]])){
  #   method <- "mutect"
  # }else{method <- input[[snv_id]]}
  dbt <- rv[[db_id]]
  if(is.null(dbt)){
    infiles <- paste0(rv$indir,"df_gene_scale.csv")
  }else if(dbt == "rna"){
    infiles <- paste0(rv$indir,"df_gene_scale.csv")
  }else if(dbt == "snv"){
    if(rv$tcga){
      # infiles <- paste0(rv$indir,"df_snv_class_",method,".csv")
      infiles <- paste0(rv$indir,"df_snv_class_977.csv")
    }else{
      infiles <- paste0(rv$indir,"df_snv_class",".csv")
    }
  }else if(dbt == "cnv"){
    infiles <- paste0(rv$indir,"df_cnv.csv")
  }else if(dbt == "mir"){
    infiles <- paste0(rv$indir,"df_mir_scale.csv")
  }else if(dbt == "met"){
    infiles <- paste0(rv$indir,"df_met_scale.csv")
  }else if(dbt == "pro"){
    infiles <- paste0(rv$indir,"df_proteomic_scale.csv")
  }else if(dbt == "rrpa"){
    infiles <- paste0(rv$indir,"df_rrpa_scale.csv")
  }else if(dbt == "crispr"){
    infiles <- paste0(rv$indir,"DepMap-CRISPR.csv")
  }else if(dbt == "rnai"){
    infiles <- paste0(rv$indir,"DepMap-RNAi.csv")
  }else if(dbt == "drug"){
    infiles <- paste0(rv$indir,"DepMap-Drug.csv")
  }
  
  l <- lapply(infiles, function(x){
    info <- fread(x,sep=",",header=T,nrows = 0)
    return(info)
  })
  
  # df_gene <- rbindlist(l, fill = T, use.names = T)
  df_gene <- rbind_common(l)
  
  # the genes
  genes <- names(df_gene) %>% .[-1] %>% sort()
  
  # save into rv$snv_genes
  if(!is.null(input[[db_id]])){
    if(input[[db_id]] == "snv"){
      rv[[paste0("snv_genes_",x)]] <- lapply(l, function(x){
        names(x) %>% .[-1]
      })
    }
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


# the button ui that will save cSurvival's variables for easygeo
btn_save_for_geo <- function(id, label){
  bsButton(
    inputId = id,
    label = tags$b(label),
    style = "primary"
    #,onclick ="location.href='http://tau.cmmt.ubc.ca/eVITTA/';target='_blank'"
    
  )
}

# save cSurvival's variables to easyGEO 's variables folder for future read
# input: the rv you would like to save
# output: the random string we have generated
save_csurvival_variable <- function(rv){
  random_string <- ids::random_id(bytes = 8)
  ofile <- paste0(surv_dir,random_string, ".rds")
  saveRDS(object = rv, 
          file = ofile)
  system(paste0("chmod -R a+rwX ",ofile))
  print(random_string)
  # url <- paste0('https://tau.cmmt.ubc.ca/eVITTA/easyGSEA/',"?data_head_o=", random_string)
  # runjs(paste0("window.open('", url, "','_blank');"))
  return(random_string)
}


# name a vector based on the value to form the choices of "Select type of molecular data:"
name_project_choices <- function(overlapped_parameter){
  if(!is.null(overlapped_parameter)){
    names(overlapped_parameter)[which(overlapped_parameter=="rna")] <- "Expression"
    names(overlapped_parameter)[which(overlapped_parameter=="snv")] <- "Mutation"
    names(overlapped_parameter)[which(overlapped_parameter=="cnv")] <- "CNV"
    names(overlapped_parameter)[which(overlapped_parameter=="mir")] <- "miRNA"
    names(overlapped_parameter)[which(overlapped_parameter=="met")] <- "Methylation"
    names(overlapped_parameter)[which(overlapped_parameter=="rrpa")] <- "RRPA"
  }
  overlapped_parameter
}

#======================================================================#
####                           Calculations                         ####
#======================================================================#
# format p values
format_p <- function(p, max = 0.0001){
  if(!is.numeric(p)){
    p
  }else{
    if(p < max){
      format(p, scientific = T, digits = 2)
    }else{
      format(p, scientific = F, digits = 2)
    }
  }
}
