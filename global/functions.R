# add help buttons to labels (need to wrap again in HTML)
# example of use: label=HTML("Label here", add_help("id1", style="padding:1px 1px 1px 1px;") )
add_help <- function(id, color="#00c0ef", size="medium", style=""){
  out <- paste0("<i class='fa fa-question-circle'
                style = 'color:",color,";
                font-size:",size,";",style,"'
                id='",id,"'></i>")
  
  HTML(out)
}

# LABELS WITH CLICKABLE BS BUTTON 
# construct a label with a clickable help bs button
label_with_help_bttn <- function(label_text, bttn_id, bttn_status="info", bttn_style="", icon_i=icon("book-open")){
  p(style="margin-block-end: 2px;",
    label_text,
    tags$style(type = "text/css", paste0("#",bttn_id,"{display: inline-block;width: 17px;height: 17px;padding: 0;border-radius: 50%;vertical-align: text-top;margin-left: 3px;font-size: 10px;padding-top: 1px;",bttn_style,"}")),
    bsButton(bttn_id, label = "", icon = icon_i, style = bttn_status, size = "extra-small"))
}

# sort a named list by its first-layer elements
sort_list <- function(alist, na.last = TRUE, decreasing=FALSE){
  alist[order(names(alist), decreasing=decreasing, na.last = na.last)]
}
