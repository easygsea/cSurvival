# ------- call rintrojs --------
call_introjs <- function(file_name) {
  rintrojs::introjs(session, options = list(showStepNumbers=FALSE,
                                            steps = file_name)
  )
}
# the events trigger by pressing the help button on header
observeEvent(input$db_help, {
  showModal(
    modalDialog(
      size = "l",easyClose = TRUE
      ,footer = modalButton("Dismiss"),
      div(
        div(align="center",
          actionBttn(
            "db_help_nav",
            "Click for intro tours"
          )
        ),
        hr(),
        tags$iframe(
          src = "https://tau.cmmt.ubc.ca/cSurvival/help.html",
          style="width:100%;",  frameborder="0",id="iframe", height = "800px"
        )
      )
    )
  )
})

observeEvent(input$db_help_nav,{
  removeModal()
  
  df <- intros$R_pre
  
  if(rv$tcga){
    df <- rbind(df,intros$P_TCGA)
  }
  
  df <- rbind(df,intros$P_main)
  
  if(length(input[["cat_1"]]) > 0){
    if(input[["cat_1"]] == "g"){
      df <- rbind(df,intros$P_g)
    }else if(input[["cat_1"]] == "gs"){
      df <- rbind(df,intros$P_gs)
      if(input[["gs_mode_1"]] == "lib"){
        df <- rbind(df,intros$P_gs_l)
      }else if(input[["gs_mode_1"]] == "manual"){
        df <- rbind(df,intros$P_gs_m)
      }
    }
  }
  
  df <- rbind(df,intros$P_confirm)
  
  call_introjs(df)
})

observeEvent(input$help_button2, {
  call_intro_post()
})

call_intro_post <- function(){
  df <- intros$R_post
  if(rv$plot_type == "all" | rv$plot_type == "1" | rv$plot_type == "2" | rv$plot_type == "gender"){
    if(rv$tcgar | rv$targetr){
      df <- rbind(df,intros$surv1,intros$V_stats)
    }else if(rv$depmapr){
      df <- rbind(df,intros$depmap1)
    }
  }else if(rv$plot_type == "track"){
    df <- rbind(df,intros$V_track)
  }else if(rv$plot_type == "scatter" | rv$plot_type == "scatter2"){
    df <- rbind(df,intros$V_data_points,intros$V_plot,intros$V_plot_gear,intros$V_scatter,intros$V_stats)
  }else if(rv$plot_type == "violin"){
    df <- rbind(df,intros$V_data_points,intros$V_plot,intros$V_plot_gear,intros$V_stats)
  }
  
  call_introjs(df)
}

# ------- run demos --------
observeEvent(input$db_demo,{
  showModal(
    modalDialog(
      title = h3("Explore the sample output that performs interactively in the same way as real output"),
      div(
        align="center",
        h3(strong("Select an example run:")),
        actionBttn("tcga_demo","TCGA|TARGET",color = "danger",style = "simple",size = "lg")
        ,actionBttn("depmap_demo","DepMap",color = "warning",style = "simple",size = "lg")
      )
      ,br()
      ,size = "l",
      easyClose = TRUE
      ,footer = modalButton("Dismiss")
    )
  )
})

observeEvent(input$tcga_demo,{
  init_tcga()
  removeModal()
  load_demo_modal(time=10000)
})

observeEvent(input$depmap_demo,{
  init_depmap()
  removeModal()
  load_demo_modal()
})

observeEvent(input$welcome_modal,{
  removeModal()
  call_intro_post()
})

# --------- demo functions ----------
init_tcga <- function(){
  rv$demo <- "yes"
  project <- "TCGA-LUAD"
  
  # update TCGA UI
  updateSelectInput(session,"project",selected = project)
  updateNumericInput(session,"variable_n",value = 2)
  shinyjs::delay(1000, shinyjs::click("confirm_project"))

  # update GS selection
  updatePrettyRadioButtons(session, "cat_1", selected = "gs")
  shinyjs::delay(3500, updatePrettyRadioButtons(session, "cat_2", selected = "gs"))
  shinyjs::delay(2000, updateSelectizeInput(session,"gs_db_1", selected = "WikiPathways"))
  shinyjs::delay(2000, updateSelectizeInput(session,"gs_db_2", selected = "WikiPathways"))
  shinyjs::delay(4500, updateSelectizeInput(session,"gs_l_1", selected = "WP_NRF2-ARE_regulation%WP4357"))
  shinyjs::delay(5000, updateSelectizeInput(session,"gs_l_2", selected = "WP_HIF1A_and_PPARG_regulation_of_glycolysis%WP2456"))
  
  # mimic clicking run button
  shinyjs::delay(8000, shinyjs::click("confirm"))
}


init_depmap <- function(){
  rv$demo <- "yes"
  project <- "DepMap-CRISPR"
  
  # update TCGA UI
  updateSelectInput(session,"project",selected = project)
  updateNumericInput(session,"variable_n",value = 2)
  shinyjs::delay(1000, shinyjs::click("confirm_project"))
  
  # update cell lines selection
  shinyjs::delay(1500, updatePickerInput(session,"ccle_cancer_types",selected = "Lung Cancer"))
  
  # update CRISPR gene selection
  rv$depmap_gene <- "NFE2L2|4780"
  # shinyjs::delay(6000, updateSelectizeInput(session,"depmap_gene",selected = "NFE2L2|4780"))
  
  # update GS selection
  updatePrettyRadioButtons(session, "cat_1", selected = "gs")
  shinyjs::delay(2000, updateSelectizeInput(session,"gs_db_1", selected = "WikiPathways"))
  shinyjs::delay(6500, updateSelectizeInput(session,"gs_l_1", selected = "WP_NRF2-ARE_regulation%WP4357"))
  shinyjs::delay(2000, updateRadioGroupButtons(session,"db_2", selected = "snv"))
  shinyjs::delay(1800, rv[["g_2"]] <- "KEAP1|9817")

  # mimic clicking run button
  shinyjs::delay(10000, shinyjs::click("confirm"))
}

load_demo_modal <- function(time=11000){
  shinyjs::delay(
    time,
    showModal(modalDialog(
      title = tags$h3("Welcome to cSurvival demo session"),
      tags$h4("Explore the sample output that performs interactively in the same way as real output.")
      ,br()
      ,tags$h4("Click OK to follow the intro tour."),
      size = "m",easyClose = FALSE,footer = actionButton("welcome_modal",label = "OK")))
  )
}