#======================================================================#
####                     general UI functions                       ####
#======================================================================#
# radiotooltip
radioTooltip <- function(id, choice, title, placement = "top", trigger = "hover", options = NULL){
  
  options = shinyBS:::buildTooltipOrPopoverOptionsList(title, placement, trigger, options)
  options = paste0("{'", paste(names(options), options, sep = "': '", collapse = "', '"), "'}")
  bsTag <- shiny::tags$script(shiny::HTML(paste0("
    $(document).ready(function() {
      setTimeout(function() {
        $('input', $('#", id, "')).each(function(){
          if(this.getAttribute('value') == '", choice, "') {
            opts = $.extend(", options, ", {html: true});
            $(this.parentElement).tooltip('destroy');
            $(this.parentElement).tooltip(opts);
          }
        })
      }, 500)
    });
  ")))
  htmltools::attachDependencies(bsTag, shinyBS:::shinyBSDep)
}

# a link icon, need to wrap in HTML
# example: HTML(paste0("Visit eVITTA at ",link_icon("evitta_link","https://tau.cmmt.ubc.ca/eVITTA/")))
link_icon <- function(id, link, title="Click to visit", icon="fas fa-external-link-alt", color="#00c0ef", style=""){
  sprintf(
    '<a href="%s" target="_blank"><i class="%s" id="%s" style = "color:%s"></i></a>'
  ,link,icon,id,color)
}

# # add a gear button
add_gear <- function(
  id, left="6.9em", top="1em", title="Click for advanced run parameters", up = F, width = "80%"
  , placement="top"
){
  div(
    style=sprintf("position: relative; align: center; left: %s; top: %s;",left, top),
    dropdownButton(
      inputId = paste0("div_",id),
      circle = TRUE, status = "info",
      size = "xs",
      icon = icon("gear"),# class = "opt"),
      up = up, width = width,
      tooltip = tooltipOptions(title = title, placement = placement),
      
      fluidRow(
        uiOutput(id)
      )
    )
  )
}

# # download link with download btn
dlink <- function(link, basename){
  paste0("<a href='",link,"' download='",basename,"' target='_blank'> <i class='fa fa-download'> </i>",basename,"</a><br/>")
}

# # red title
red_title <- function(txt,color = "#CF5C78", size="120%"){
  strong(style = sprintf("color:%s;font-size:%s",color,size),txt)
}

# # default hr
default_hr <- function(border="0.5px",line="solid", color="#F0EEE9", margin="0em"){
  tags$hr(style=sprintf("border: %s %s %s; margin: %s;",border,line,color,margin))
}

# # strong h4
strong_h4 <- function(txt, h="4"){
  HTML(paste0("<h",h,"><b>",txt,"</b></h",h,">"))
}

# # grey h4
grey_h4 <- function(txt,color="#939597",h="4"){
  HTML(paste0('<h',h,' style="color:',color,';">',txt,'</h',h,'>'))
}

#======================================================================#
####          function to generate dynamic 1 2, or more G/GS UIs     ####
#======================================================================#
extract_color <- function(x, cols=bcols){
  bcol_n <- x %% 2; if(bcol_n == 0){bcol_n <- 2}
  return(cols[bcol_n])
}

plot_ui <- function(n){
  # automatically adjust column width according to # of analysis selected
  if(n == 1){col_w <- 12}else{col_w <- 6}
  
  # create the UI list
  ui <- lapply(1:n, function(x){
    # color number for the wellpanel
    col <- extract_color(x)
    # category to analyze
    cat_id <- paste0("cat_",x); cat_id_q <- paste0(cat_id,"_q")
    # type of db to analyze
    db_id <- paste0("db_",x); db_id_q <- paste0(db_id,"_q")
    # gene to analyze
    g_ui_id <- paste0("g_",x); g_ui_id_q <- paste0(g_ui_id,"_q")
    # gene normalization options
    g_ui_norm_id <- paste0("gnorm_",x); g_ui_norm_id_q <- paste0(g_ui_norm_id,"_q")
    g_ui_norm_g_id <- paste0("gnorm_g_",x); g_ui_norm_g_id_q <- paste0(g_ui_norm_g_id,"_q")
    g_ui_norm_gs_db_id <- paste0("gnorm_gs_db_",x); g_ui_norm_gs_db_id_q <- paste0(g_ui_norm_gs_db_id,"_q")
    g_ui_norm_gs_lib_id <- paste0("gnorm_gs_lib_",x); g_ui_norm_gs_lib_id_q <- paste0(g_ui_norm_gs_db_id,"_q")
    g_ui_norm_gs_lib_genes_id <- paste0("gnorm_gs_lib_genes_",x)
    g_ui_norm_gs_lib_dn_id <- paste0(g_ui_norm_gs_lib_id,"dn"); g_ui_norm_gs_lib_dn_id_q <- paste0(g_ui_norm_gs_lib_dn_id,"_q")
    g_ui_norm_gs_lib_link_id <- paste0(g_ui_norm_gs_lib_id,"link"); g_ui_norm_gs_lib_link_ui_id <- paste0(g_ui_norm_gs_lib_link_id,"_ui")
    # gene set to analyze
    gs_mode_id <- paste0("gs_mode_",x); gs_mode_id_q <- paste0(gs_mode_id,"_q")
    gs_db_id <- paste0("gs_db_",x); gs_db_id_q <- paste0(gs_db_id,"_q")
    gs_lib_id <- paste0("gs_l_",x); gs_lib_id_q <- paste0(gs_lib_id,"_q")
    gs_lib_dn_id <- paste0(gs_lib_id,"dn"); gs_lib_dn_id_q <- paste0(gs_lib_dn_id,"_q")
    gs_lib_link_id <- paste0(gs_lib_id,"link"); gs_lib_link_ui_id <- paste0(gs_lib_link_id,"_ui")
    gs_lib_genes_id <- paste0("gs_lgs_",x) # verbatimTextOutput on gs genes
    gs_gene_id <- paste0("gs_lg_",x); gs_gene_id_q <- paste0(gs_gene_id,"_q") # gene to search
    gs_gene_genes_id <- paste0("gs_lgg_",x) # verbatimTextOutput on input genes to filter GS
    # manual gene input
    gs_manual_id <- paste0("gs_m_",x); gs_manual_id_q <- paste0(gs_manual_id,"_q")
    gs_manual_btn_id <- paste0("add_btn_",x)
    gs_genes_id <- paste0("gs_mg_",x)

    if(!rv[[db_id]] %in% data_types()){rv[[db_id]] <- data_types()[1]}
    # the UI
    column(
      col_w,
      # tags$hr(style="border: .5px solid lightgrey; margin-top: 0.5em; margin-bottom: 0.5em;"),
      wellPanel(
        style = paste0("background-color: ", col, "; border: .5px solid #fff;"),
        h4(paste0("Analysis #",x), align = "center"),
        prettyRadioButtons(
          cat_id,
          label = HTML(paste0(x,".1. Select data category:",add_help(cat_id_q))),
          choices = c("Gene or locus" = "g", "Gene set (GS)" = "gs")
          ,selected = rv[[cat_id]]
          ,status = "default"
          ,icon = icon("check", class="icon_red")
          ,shape = "curve"
          # ,plain = T
          # ,outline = T
          ,animation = "jelly"
          ,inline = T
        )
        ,conditionalPanel(
          condition = sprintf("input.%s=='g'", cat_id),
          radioGroupButtons(
            inputId = db_id,
            label = HTML(paste0(x,".2. Select type of molecular data:",add_help(db_id_q))),
            choices = data_types(),
            status = "danger",
            selected = rv[[db_id]],
            checkIcon = list(
              yes = tags$i(class = "fa fa-check-square", 
                           style = "color: white"),
              no = tags$i(class = "fa fa-square-o", 
                          style = "color: white"))
          )
          ,div(id = "div_g",
            selectizeInput(
              g_ui_id,
              HTML(paste0(x,".3. Select your gene of interest:",add_help(g_ui_id_q)))
              ,choices=c()
              ,selected=rv[[g_ui_id]]
              ,width = "100%"
              ,options = list(
                placeholder = g_placeholder()
                ,onInitialize = I(sprintf('function() { this.setValue("%s"); }',rv[[g_ui_id]]))
              )
            )
          )
        )
        ,conditionalPanel(
          condition = sprintf("input.%s=='gs'", cat_id),
          radioButtons(
            gs_mode_id,
            HTML(paste0(x,".2. Mode of analysis:",add_help(gs_mode_id_q)))
            ,choices = c("Library"="lib","Manual"="manual")
            ,selected = rv[[gs_mode_id]]
            ,inline = T
          )
          ,conditionalPanel(
            condition = sprintf("input.%s=='lib'", gs_mode_id),
            fluidRow(id=paste0("div_gs_db_",x),
              column(
                4,
                selectizeInput(
                  gs_db_id,
                  HTML(paste0(x,".3a. Select database:"),add_help(gs_db_id_q))
                  ,choices=gmt_dbs
                  ,selected=rv[[gs_db_id]]
                  ,options = list(
                    # `live-search` = TRUE,
                    placeholder = 'Type to search ...'
                    ,onInitialize = I(sprintf('function() { this.setValue("%s"); }',rv[[gs_db_id]]))
                  )
                )
              )
              ,conditionalPanel(
                condition = sprintf("input.%s != ''", gs_db_id),
                column(
                  4,
                  selectizeInput(
                    gs_lib_id,
                    HTML(paste0(x,".3b. Select gene set (GS):"),add_help(gs_lib_id_q))
                    ,choices=names(rv[[paste0("gmts_tmp",x)]])
                    ,selected=rv[[gs_lib_id]]
                    ,options = list(
                      # `live-search` = TRUE,
                      placeholder = rv[[paste0("gs_placeholder",x)]]
                      ,onInitialize = I(sprintf('function() { this.setValue("%s"); }',rv[[gs_lib_id]]))
                    )
                  )
                )
                
              )
              ,conditionalPanel(
                condition = sprintf("input.%s != ''", gs_db_id),
                column(
                  4,
                  searchInput(
                    gs_gene_id,
                    HTML(paste0("(Optional) gene set filtering by gene combination:", add_help(gs_gene_id_q))),
                    value = rv[[gs_gene_id]],
                    placeholder = "Enter genes in HUGO symbol format",
                    btnSearch = icon("search"),
                    btnReset = icon("remove"),
                    width = "100%"
                  )
                  
                )
              )
            )
            ,fluidRow(
              column(
                    8,
                    conditionalPanel(
                      condition = sprintf("input.%s != ''", gs_lib_id),
                      div(id=gs_lib_dn_id_q,style="z-index:1000;",
                        style = "position: absolute; right: -1.5em; top: 0.2em;",
                        downloadBttn(
                          gs_lib_dn_id,NULL
                          ,size = "xs", color = "danger", style = "material-circle",
                        )
                      )
                      ,uiOutput(gs_lib_link_ui_id)
                      ,bsTooltip(gs_lib_dn_id_q,HTML("Click to download the gene list. Adjust and re-upload using the <b>Manual</b> mode if necessary."),placement = "top")
                      ,span(verbatimTextOutput(gs_lib_genes_id), style = rv$verbTxtStyle1)
                    )
              )
              ,column(
                  4,
                  conditionalPanel(
                    condition = sprintf("input.%s != ''", gs_gene_id),
                    span(verbatimTextOutput(gs_gene_genes_id), style = rv$verbTxtStyle2)
                )
                
              )
            )
          )
          ,conditionalPanel(
            condition = sprintf("input.%s=='manual'", gs_mode_id),
            textAreaInput(
              gs_manual_id,
              HTML(paste0(x,".3. Enter your genes:"),add_help(gs_manual_id_q))
              ,value = rv[[gs_manual_id]]
              ,placeholder = "Type to enter..."
            )
            ,span(verbatimTextOutput(gs_genes_id), style = paste0(rv$verbTxtStyle1))
            ,bsButton(gs_manual_btn_id,tags$strong("Submit"),icon=icon("upload"),style = "primary")
          )
        )
        ,conditionalPanel(
          sprintf("(input.%s== 'g' & input.%s== 'rna' & input.%s!='') | (input.%s== 'gs' & input.%s=='lib' & input.%s!='') | (input.%s== 'gs' & input.%s=='manual' & output.%s)", cat_id,db_id,g_ui_id, cat_id,gs_mode_id,gs_lib_id, cat_id,gs_mode_id,paste0("gs_manual_uploaded",x)),
          radioGroupButtons(
            g_ui_norm_id,
            HTML(paste0(x,".4a. (Optional) normalize gene expression by:",add_help(g_ui_norm_id_q))),
            choices = c("None"="none","Gene"="g","Gene set (GS)"="gs"),
            selected = rv[[g_ui_norm_id]],
            size = "sm",
            checkIcon = list(
              yes = icon("check-square"),
              no = icon("square-o")
            ),
            direction = "horizontal"
          )
          ,conditionalPanel(
            sprintf("input.%s == 'g'", g_ui_norm_id),
            selectizeInput(
              g_ui_norm_g_id,
              HTML(paste0(x,".4b. Select a normalization gene:",add_help(g_ui_norm_g_id_q)))
              ,choices=c()
              ,selected=rv[[g_ui_norm_g_id]]
              ,width = "100%"
              ,options = list(
                placeholder = "Type to search ..."
                ,onInitialize = I(sprintf('function() { this.setValue("%s"); }',rv[[g_ui_norm_g_id]]))
              )
            )
          )
          ,conditionalPanel(
            sprintf("input.%s == 'gs'", g_ui_norm_id),
            fluidRow(
              column(
                4,
                selectizeInput(
                  g_ui_norm_gs_db_id,
                  HTML(paste0(x,".4b. Select database:"),add_help(g_ui_norm_gs_db_id_q))
                  ,choices=gmt_dbs
                  ,selected=rv[[g_ui_norm_gs_db_id]]
                  ,options = list(
                    # `live-search` = TRUE,
                    placeholder = 'Type to search ...'
                    ,onInitialize = I(sprintf('function() { this.setValue("%s"); }',rv[[g_ui_norm_gs_db_id]]))
                  )
                )
              )
              ,conditionalPanel(
                condition = sprintf("input.%s != ''", g_ui_norm_gs_db_id),
                column(
                  8,
                  selectizeInput(
                    g_ui_norm_gs_lib_id,
                    HTML(paste0(x,".4c. Select gene set (GS):"),add_help(g_ui_norm_gs_lib_id_q))
                    ,choices=names(rv[[paste0("gnorm_gmts",x)]])
                    ,selected=rv[[g_ui_norm_gs_lib_id]]
                    ,options = list(
                      # `live-search` = TRUE,
                      placeholder = rv[[paste0("gnorm_gs_placeholder",x)]]
                      ,onInitialize = I(sprintf('function() { this.setValue("%s"); }',rv[[g_ui_norm_gs_lib_id]]))
                    )
                  )
                )
              )
              ,conditionalPanel(
                condition = sprintf("input.%s != ''", g_ui_norm_gs_lib_id),
                column(
                  11,
                  div(id=g_ui_norm_gs_lib_dn_id_q,style="z-index:1000;",
                      style = "position: absolute; right: -1.5em; top: 0.2em;",
                      downloadBttn(
                        g_ui_norm_gs_lib_dn_id,NULL
                        ,size = "xs", color = "danger", style = "material-circle",
                      )
                  )
                  ,uiOutput(g_ui_norm_gs_lib_link_ui_id)
                  ,bsTooltip(g_ui_norm_gs_lib_dn_id_q,HTML("Click to download the gene list"),placement = "top")
                  ,span(verbatimTextOutput(g_ui_norm_gs_lib_genes_id), style = rv$verbTxtStyle1)
                )
              )
            )
          )
        )

        # tooltip for data category
        ,bsTooltip(cat_id_q, HTML(cat_id_q_txt)
                   ,placement = "right")
        ,bsTooltip(db_id_q, HTML(db_id_q_txt)
                   ,placement = "right")
        ,bsTooltip(g_ui_id_q, HTML(g_id_txt)
                   ,placement = "right")
        ,bsTooltip(gs_mode_id_q, HTML("Select <b>Library</b> to analyze a pathway, a biological process, a cellular location, a transcriptional factor, a drug, or a gene\\'s interacting partners.<br>Alternatively, select <b>Manual</b> to enter your own list of genes.")
                   ,placement = "right")
        ,bsTooltip(gs_db_id_q, HTML("For full list of available gene set databases, visit <b>easyGSEA User Guide</b>.")
                   ,placement = "right")
        ,bsTooltip(gs_lib_id_q, HTML("Search for keywords (e.g. glycolysis, chemokine, tor signaling) and select the one of interest.")
                   ,placement = "right")
        ,bsTooltip(gs_gene_id_q, HTML(paste0("To filter out gene sets that contains your gene(s) of interest, in HUGO symbol format, delimited by \"&\" (and) or \"|\" (or). | is evaluated before &. Example: MYC&TP53|BCL2&BRCA1|BRCA2 is evaluated as MYC&(TP53|BCL2)&(BRCA1|BRCA2), which means MYC and (TP53 or BCL2) and (BRCA1 or BRCA2)."
                                             ," A maximum of 10 genes are supported."))
                   ,placement = "top")
        ,bsTooltip(gs_manual_id_q,HTML("Newline-, space- or comma-delimited")
                   ,placement = "right")
        ,radioTooltip(id = db_id, choice = "rna", title = HTML("Gene expression level quantified by RNA-seq"))
        ,radioTooltip(id = db_id, choice = "snv", title = HTML("Simple Nucleotide Variation (SNV)"))
        ,radioTooltip(id = db_id, choice = "cnv", title = HTML("Copy Number Variation (CNV)"))
        ,radioTooltip(id = db_id, choice = "mir", title = HTML("microRNA expression level"))
        ,radioTooltip(id = db_id, choice = "met", title = HTML("DNA methylation level"))
        ,radioTooltip(id = db_id, choice = "rrpa", title = HTML("Reverse-phase protein array (RPPA)"))
        ,radioTooltip(id = db_id, choice = "pro", title = HTML("Normalized protein expression data by mass spectrometry"))
        ,radioTooltip(id = db_id, choice = "crispr", title = HTML("Gene effect measured by CRISPR-Cas9"))
        ,radioTooltip(id = db_id, choice = "rnai", title = HTML("Gene effect measured by RNA interference"))
        ,radioTooltip(id = db_id, choice = "drug", title = HTML("Drug sensitivity assays by PRISM compound repurposing screening"))
        ,bsTooltip(g_ui_norm_id_q,HTML(
          paste0("To study gene ratios.<br>"
            ,"<b>None</b>: No normalization"
            ,"<br><b>Gene</b>: Normalize expression by expression of another gene"
            ,"<br><b>Gene set (GS)</b>: Normalize expression by expression of another GS"))
          ,placement = "right")
        ,radioTooltip(id = g_ui_norm_id, choice = "none", title = HTML("No normalization"))
        ,radioTooltip(id = g_ui_norm_id, choice = "g", title = HTML("Normalize against a gene"))
        ,radioTooltip(id = g_ui_norm_id, choice = "gs", title = HTML("Normalize against a GS"))
        ,bsTooltip(g_ui_norm_g_id_q,HTML("Search and select a gene for normalization purpose")
                   ,placement = "right")
        ,bsTooltip(g_ui_norm_gs_db_id_q,HTML("For full list of available gene set databases, visit <b>easyGSEA User Guide</b>.")
                   ,placement = "right")
        ,bsTooltip(g_ui_norm_gs_lib_id_q,HTML("Search for keywords (e.g. glycolysis, chemokine, tor signaling) and select the one of interest.")
                   ,placement = "right")
      )
      
    )
    
  })
  
  # create the UI
  do.call(tagList, ui)
}

#======================================================================#
####           function to generate dynamic run parameter UI     ####
#======================================================================#
plot_run_ui <- function(n){
  if(n == 1){col_w <- 12}else{col_w <- 6}
  
  ui <- lapply(1:n, function(x){
    
    iter_id <- paste0("iter_",x); iter_id_q <- paste0(iter_id,"_q")
    
    lower_id <- paste0("lower_",x); lower_id_q <- paste0(lower_id,"_q")
    higher_id <- paste0("upper_",x); higher_id_q <- paste0(higher_id,"_q")
    step_id <- paste0("step_",x); step_id_q <- paste0(step_id,"_q")
    
    clow_id <- paste0("clow_",x); clow_id_q <- paste0(clow_id,"_q")

    col <- extract_color(x)
    
    # snv_id <- paste0("snv_method_",x); snv_id_q <- paste0(snv_id,"_q")
    # snv_uni_id <- paste0("snv_uni_",x); snv_uni_id_q <- paste0(snv_uni_id,"_q")
    non_id <- paste0("nonsynonymous_",x); non_id_q <- paste0(non_id,"_q")
    syn_id <- paste0("synonymous_",x); syn_id_q <- paste0(syn_id,"_q")
    
    cnv_id <- paste0("cnv_par_",x); cnv_id_q <- paste0(cnv_id,"_q")
    
    cat_id <- paste0("cat_",x); db_id <- paste0("db_",x)
    check_inputs <- function(){
      if(is.null(input[[cat_id]]) & is.null(input[[db_id]])){
        return(rv[[cat_id]] == "g" & rv[[db_id]] != "snv")
      }else if(is.null(input[[cat_id]])){
        return(rv[[cat_id]] == "g" & input[[db_id]] != "snv")
      }else if(is.null(input[[db_id]])){
        return(input[[cat_id]] == "g" & rv[[db_id]] != "snv")
      }else{
        return((input[[cat_id]] == "g" & input[[db_id]] != "snv") | input[[cat_id]] == "gs")
      }
    }
    datatype <- call_datatype(x)
    
    # # preselected SNV callers
    # snv_pre_a <- ifelse(
    #   length(rv[[snv_id]]) == 1,
    #   paste0('"',rv[[snv_id]],'"'),
    #   paste0("['",paste0(rv[[snv_id]],collapse = "','"),"']")
    # )
    
    # render the UI
    column(
      col_w,align="center",
      # tags$hr(style="border: .5px solid lightgrey; margin-top: 0.5em; margin-bottom: 0.5em;"),
      wellPanel(
        style = paste0("background-color: ", col, "; border: .5px solid #fff;"),
        h4(paste0("Advanced run parameters for Analysis #",x), align = "center"),
        h4(paste0("(",datatype,")")),
        tags$hr(style="border-color: #c2bfb5;"),
        if((!rv$depmap & check_inputs() & ifelse_rv(db_id) != "cnv")|(rv$depmap & check_inputs())){
          div(
            radioGroupButtons(
              inputId = iter_id,
              label = HTML(paste0("Select method to stratify high- and low- ",datatype," groups"),add_help(iter_id_q)),
              choiceNames = c("Dynamic iteration", "Manual cutoffs"),
              choiceValues = c("iter","manual"),
              selected = rv[[iter_id]],
              size = "sm",
              checkIcon = list(
                yes = icon("check-square"),
                no = icon("square-o")
              ),
              # status = "primary",
              direction = "horizontal"
            ),
            conditionalPanel(
              sprintf('input.%s == "iter"',iter_id),
              sliderInput(
                lower_id,
                label = HTML(paste0("Start percentile:",add_help(lower_id_q)))
                ,value = rv[[lower_id]]
                # ,choices = c(.05, .1, .15, .2, .25, .3, .35, .4, .45, .5)
                # ,grid = TRUE
                ,min=.05,max=.5,step=.01
              )
              ,sliderInput(
                higher_id,
                label = HTML(paste0("End percentile:",add_help(higher_id_q)))
                ,value = rv[[higher_id]]
                # ,choices = c(.5, .55, .6, .65, .7, .75, .8, .85, .9, .95)
                # ,grid = TRUE
                ,min=.5,max=.95,step=.01
              )
              ,sliderInput(
                step_id,
                label = HTML(paste0("Step size:",add_help(step_id_q)))
                ,value = rv[[step_id]]
                # ,choices = c(.01, .02, .03, .05, .1, .15, .2, .25)
                # ,grid = TRUE
                ,min=.01,max=.25,step=.005
              )
            )
            ,conditionalPanel(
              sprintf('input.%s == "manual"',iter_id),
              sliderInput(
                clow_id,
                HTML(paste0("Cutoff percentile (%):",add_help(clow_id_q))),
                value = rv[[clow_id]],
                min = 10,max = 90,step=.5
              )
            )
            ,bsTooltip(iter_id_q,HTML(paste0("<b>Dynamic iteration</b>: Determine the optimal cutoff by searching for the percentile yielding the lowest P-value"
                                             ,"<br><b>Manual cutoffs</b>: Manually enter the cutoffs for high- and low-",datatype," groups"))
                       ,placement = "top")
            ,bsTooltip(lower_id_q,HTML("The percentile to start iteration")
                       ,placement = "right")
            ,bsTooltip(higher_id_q,HTML("The percentile to end iteration")
                       ,placement = "right")
            ,bsTooltip(step_id_q,HTML(paste0("Step size to iterate to find the optimum cutoff"))
                       ,placement = "right")
            ,bsTooltip(clow_id_q,HTML("Cases &le; the cutoff will be classified as low, while those &gt; the cutoff will be classified as high")
                       ,placement = "right")
          )
        }else if(input[[cat_id]] == "g" & input[[db_id]] == "cnv"){
          div(
            radioGroupButtons(
              inputId = cnv_id,
              label = HTML(paste0("Select group to analyze:"),add_help(cnv_id_q)),
              choiceNames = c("Automatic", "Copy number gain", "Copy number loss"), #,"Copy number gain and loss"
              choiceValues = c("auto","gain","loss"), #,"both"
              selected = rv[[cnv_id]],
              size = "sm",
              checkIcon = list(
                yes = icon("check-square"),
                no = icon("square-o")
              ),
              # status = "primary",
              direction = "horizontal"
            )
            ,bsTooltip(cnv_id_q,HTML(paste0("<b>Automatic</b>: Automatically determines whether copy number gain and/or loss results in more significant survival difference"
                                            ,"<br><b>Copy number gain</b>: To compare cases with copy number gain with the rest of the population"
                                            ,"<br><b>Copy number loss</b>: To compare cases with copy number loss with the rest of the population"
                                            ,"<br><b>Copy number gain and loss</b>: To compare cases with copy number gain and loss, respectively, with the rest of the population"
            ))
            ,placement = "top")
          )
        }else{
          
          div(
            # if(rv$tcga){
            #   div(
            #     # mutation caller options
            #     selectizeInput(
            #       snv_id
            #       ,label = HTML(paste0("Select somatic mutation caller(s):",add_help(snv_id_q)))
            #       ,choices = snv_algorithms
            #       ,selected = rv[[snv_id]]
            #       ,multiple = T
            #       ,options = list(
            #         # `live-search` = TRUE,
            #         placeholder = 'Type to search ...'
            #         ,onInitialize = I(sprintf('function() { this.setValue(%s); }',snv_pre_a)))
            #     )
            #     ,bsTooltip(snv_id_q
            #                ,HTML("Algorithm(s) used for calling somatic nucleotide variations. Multiple selections are allowed")
            #                      ,placement = "top")
            #     # intersect or union
            #     ,conditionalPanel(
            #       sprintf('input.%s.length > 1',snv_id),
            #       radioGroupButtons(
            #         inputId = snv_uni_id,
            #         label = HTML(paste0("Select method to handle results by different callers: ",add_help(snv_uni_id_q))),
            #         choices = c(
            #           "Intersect" = "int"
            #           ,"Union" = "uni"
            #         ),
            #         selected = rv[[snv_uni_id]],
            #         size = "sm",
            #         checkIcon = list(
            #           yes = icon("check-square"),
            #           no = icon("square-o")
            #         ),
            #         # status = "primary",
            #         direction = "horizontal"
            #       )
            #     )
            #     ,bsTooltip(snv_uni_id_q,HTML(paste0(
            #       "When multiple callers are selected, select <b>Intersect</b> to extract consensus results for analysis, or <b>Union</b> to find a positive hit in any algorithm."
            #     )),placement = "top")
            #   )
            # }
            
            # non-silent variants classifications
            selectizeInput(
              non_id
              ,label = HTML(paste0("Select variants of interest (Mutated): ",add_help(non_id_q)))
              ,choices = variant_types
              ,selected = rv[[non_id]]
              ,multiple = T
            )
            ,bsTooltip(non_id_q,HTML("By default, variants to be classified as High/Moderate variant consequences are selected. Adjust to suit the purpose of your study. For more information about variant classification: http://uswest.ensembl.org/Help/Glossary?id=535")
                       ,placement = "top")
            
            # silent variants classifications
            ,selectizeInput(
              syn_id
              ,label = HTML(paste0("Select control group (Other): ",add_help(syn_id_q)))
              ,choices = variant_types
              ,selected = rv[[syn_id]]
              ,multiple = T
            )
            ,bsTooltip(syn_id_q,HTML("By default, wild-type (WT) and variants to be classified as Low/No variant consequences are selected. Adjust to suit the purpose of your study. For more information, visit http://uswest.ensembl.org/Help/Glossary?id=535")
                       ,placement = "top")
          )
        }
        ,div(
          style="display: inline-block;vertical-align:top;",
          bsButton(paste0("todefault",x), strong("Reset to default"), icon = icon("redo"), style = "primary")
        )
        ,if(x == 1 & rv$variable_n > 1){
          if(req_filter_on(paste0("db_",2:rv$variable_n),filter="snv",target="input")){
            div(
              style="display: inline-block;vertical-align:top;",
              bsButton("toall", strong("Apply to all"), icon = icon("globe"), style = "primary")
            )
          }else if(req_filter_on(paste0("db_",2:rv$variable_n),filter="snv",target="input",mode="unequal")){
            div(
              style="display: inline-block;vertical-align:top;",
              bsButton("toall_m", strong("Apply to all"), icon = icon("globe"), style = "primary")
            )
          }
        }
      )
    )
  })

  do.call(tagList, ui)
  
}




#======================================================================#
####           function to deynamically plot the percentile graphs     ####
#======================================================================#

#helper function for vertical line
vline <- function(x = 0, color = "pink") {
  list(
    type = "line", 
    y0 = 0, 
    y1 = 1,
    yref = "paper",
    x0 = x, 
    x1 = x, 
    line = list(color = color, dash="dashdot")
  )
}


single_plot <- function(quantile_df, index){
  #Due to a bug in Plotly, I am using the solution mentioned in this website:
  #https://stackoverflow.com/questions/55251470/missing-data-when-supplying-a-dual-axis-multiple-traces-to-subplot
  #If index == 1, this is first graph and we need to pass in normal y axis
  cutoff <- quantile_df$quantile[which.min(quantile_df$p_value)]
  if(index == 1){
    yaxes = c("y","y2")
    HR_Axis <- list(overlaying = yaxes[1],
                    side = "right", title = "Hazard ratio")
    P_Axis <- list(side = "left", title = "P-value")
  }
  else{
    yaxes = c("y2","y3")
    HR_Axis <- list(overlaying = yaxes[2],
                    side = "right", title = "Hazard ratio")
    P_Axis <- list(side = "left", title = "P-value")
  }
  subtitle <- paste0(rv[[paste0("title_",index)]]) #"Subplot of ",
  
  fig <- plot_ly(quantile_df, x = quantile_df$quantile)
  fig <- fig %>% add_trace(y = ~quantile_df$p_value,type = 'scatter',#color =~p_value,
                           line = list(color = 'rgb(173,173,173)', width = 2),
                           yaxis = yaxes[1],
                           name = 'P-value',
                           marker=list(
                             color=~p_value,
                             colorscale=col_scale,
                             cmid = 0.5,
                             reversescale =TRUE
                           ),
                           text = quantile_df$expression,
                           name = '',mode = 'lines+markers', hovertemplate = paste0(
                             "P-value: %{y:.8f}<br>",
                             "Quantile (in %): %{x:.0f}<br>",
                             "Expression: %{text:.3f}<br>"
                           )) %>%
    add_annotations(
      text = subtitle,
      x = 0.5,
      y = 1,
      yref = "paper",
      xref = "paper",
      xanchor = "middle",
      yanchor = "top",
      showarrow = FALSE,
      font = list(size = 20))
    # )%>%
    # add_annotations(
    #   text = "minimum p value cut off",
    #   x = quantile_df$quantile[which.min(quantile_df$p_value)],
    #   y = 0.1,
    #   yref = "paper",
    #   xanchor = "middle",
    #   yanchor = "top",
    #   showarrow = FALSE,
    #   font = list(size = 15)
    # )
  
  #%>%
    #add_segments(x = quantile_df$quantile[which.min(quantile_df$p_value)], xend = quantile_df$quantile[which.min(quantile_df$p_value)], y = 0, yend = 1)
  #Add a DepMap Check
  if(!all(is.na(quantile_df$hr))){
    fig <- fig %>%
      add_trace(y = quantile_df$hr, name = 'Hazard ratio',mode = 'lines+markers',type = 'scatter',
                line = list(color = 'rgb(212, 235, 242)', width = 2),
                marker=list(
                  symbol = 'diamond',
                  color=quantile_df$hr,
                  colorscale='RdBu',
                  cmid = 1,
                  reversescale =FALSE
                ),
                yaxis = yaxes[2]
      )
  }
    
  #different layout setting for first graph and second graph
  if(index == 1){
    #omit second axis if all hr values = null
    if(all(is.na(quantile_df$hr))){
      fig <- fig %>%
        add_trace(y = c(0,1), x = c(cutoff, cutoff),type = 'scatter', mode = 'lines',#color =~p_value,
                  line = list(color = 'pink', dash="dashdot"),
                  yaxis = yaxes[1],
                  name = 'Minimum P-value cutoff') %>%
        layout(title = 'P-Values of Different Percentiles',
               xaxis = list(title = 'Quantile (%)'),
               yaxis = P_Axis,
               hovermode = "x unified"#,
               #shapes = list(vline(quantile_df$quantile[which.min(quantile_df$p_value)]))
        )
    }else{
      fig <- fig %>%
        add_trace(y = c(0,1), x = c(cutoff, cutoff),type = 'scatter', mode = 'lines',#color =~p_value,
                  line = list(color = 'pink', dash="dashdot"),
                  yaxis = yaxes[1],
                  name = 'Minimum P-value cutoff') %>%
        layout(title = 'P-Values and Harzard Ratios of Different Percentiles',
               xaxis = list(title = 'Quantile (%)'),
               yaxis = P_Axis,
               yaxis2 = HR_Axis,
               hovermode = "x unified"#,
               #shapes = list(vline(quantile_df$quantile[which.min(quantile_df$p_value)]))
        )
    }
  }else{
    fig <- fig %>%
      add_trace(y = c(0,1), x = c(cutoff, cutoff),type = 'scatter', mode = 'lines',#color =~p_value,
                line = list(color = 'pink', dash="dashdot"),
                yaxis = yaxes[1],
                name = 'Minimum P-value cutoff',
                showlegend = FALSE) %>%
      layout(title = 'P-Values and Harzard Ratios of Different Percentiles',
             xaxis = list(title = 'Quantile (%)'),
             yaxis2 = P_Axis,
             yaxis3 = HR_Axis,
             hovermode = "x unified"#,
             #shapes = list(vline(quantile_df$quantile[which.min(quantile_df$p_value)]))
      )
  }
  
  return(fig)
}


assemble_percentile_plot <- function(quantile_df_list){
  fig_list <- c()
  i <- 0
  for(index in 1:length(quantile_df_list)){
    if(!is.null(rv[["quantile_graph"]][[index]])){
      i <- i + 1
      fig_list[[i]] <- single_plot(quantile_df_list[[index]], index = index)
    }
  }
  
  return(fig_list)
}

#heatmap----
pvalue_heatmap <- function(heatmap_df){
  #prepare dat needed for plotly
  dat = heatmap_df
  #annotation and hr columns are for texts on top of heatmap graph on hovertext
  dat$annotation = as.character(lapply(dat$annotation,function(x){format_heatmap_p(x)}))
  dat$hr = as.character(round(dat$hr,2))
  dat[is.na(dat)] <- "NA"
  #p value and log p value columns are for plotly ploting and color scaling
  dat$p_value = round(dat$p_value,4)
  #p value column is for texts on top of heatmap graph
  dat$log_p_value <- -log10(dat$p_value)
  dat$log_p_value[is.na(dat$log_p_value)] <- 0
  #sometimes this value could become Inf, and thus make the heatmap's distribution unclear,
  #so add this line to replace all exceeding-threshold values as threshold values
  dat$log_p_value[dat$log_p_value > heatmap_maximum_thershold] <- heatmap_maximum_thershold
  #axis names:
  x_axis_title = str_remove(rv[["title_1"]], fixed(" expression"))
  y_axis_title = str_remove(rv[["title_2"]], fixed(" expression"))
  
  req(length(dat$log_p_value)>0)
  selected_color = col_scale_o
  #default color scale
  selected_color_no <- c(0, 0.16666666, 0.33333333333, 0.5, 0.66666666, 1)
  if(rv$tracking_heatmap_color != "white"){
    if(rv$tracking_heatmap_color != "default"){
      selected_color <- sequential_hcl(5, palette = rv$tracking_heatmap_color) %>% rev(.)
    }
    selected_color <- col2rgb(selected_color) #adjust_transparency(selected_color,alpha=rv$hm_alpha)
    selected_color <- sapply(1:ncol(selected_color), function(x) {x <- c(selected_color[,x],rv$hm_alpha); paste0("rgba(",paste0(x, collapse = ", "),")")})
    if(rv$hm_text_color == "black"){
      selected_color <- c(paste0("rgba(255, 255, 255, ",rv$hm_alpha,")"), selected_color)
    }else if(rv$hm_text_color == "white"){
      selected_color <- c(paste0("rgba(190, 190, 190, ",rv$hm_alpha,")"), selected_color)
    }
    selected_color <- lapply(1:length(selected_color), function(i) list(selected_color_no[i], selected_color[i]))
  }else if(rv$tracking_heatmap_color == "white"){
    selected_color <- rep("rgb(255, 255, 255)",6)
    selected_color <- lapply(1:length(selected_color), function(i) list(selected_color_no[i], selected_color[i]))
  }
  fig <- plot_ly(hoverinfo = "text",hovertext = paste0("P value: ", '<b>',dat$annotation,'</b>',
                                  "<br> Hazard Ratio: ", '<b>', dat$hr,'</b>',
                                  "<br>", x_axis_title, " quantile: ", '<b>',dat$Q1,'%','</b>',
                                  "<br>", y_axis_title, " quantile: ", '<b>',dat$Q2,'%','</b>'
                                  )) %>%
    add_trace(data = dat, x = ~Q1, y = ~Q2, z = ~log_p_value, type = "heatmap",
              colorscale  = selected_color,zmax = heatmap_maximum_thershold,zmin=0, #max(dat$log_p_value)
              colorbar = list(
                title = list(text="-log10(P)", side = "right")
                ,len = 1.2
              )
              #text = ~annotation,
              #hovertemplate = paste('%{yaxis.title.text} quantile: <b>%{x}%</b> and %{yaxis.title.text} quantile: <b>%{y}%</b><br>',
              #                      'Corresponding P-value: <b>%{text}</b>','<extra></extra>'
              #)
    )%>%
    layout(
      title = "Quantile P value Heatmap",
      xaxis = list(title = x_axis_title, showticklabels = T, showgrid = F),
      yaxis = list(title = y_axis_title, showticklabels = T, showgrid = F)
      # ,margin = list(l=200)
    )
  # cols <- ifelse(dat$annotation == "NA","#000000","#ffffff")
  if(rv$tracking_heatmap_annotation){
    if(rv$hm_text_color == "black"){
      tcol <- "#000000"
    }else if(rv$hm_text_color == "white"){
      tcol <- "#ffffff"
    }
    fig <- fig %>%
      add_annotations(x = dat$Q1, y = dat$Q2,
                      text = dat$annotation,
                      showarrow = FALSE, xref = 'x', yref = 'y', font=list(color=tcol, size = rv$tracking_heatmap_text_size)
                      ,ax = 20, ay = -20
      )
  }
  fig
}