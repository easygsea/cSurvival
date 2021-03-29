library(survival) # to do the survival analysis
library(survminer) # to plot the survival analysis nicer

setwd("/Users/jeancheng/Documents/cSurvival/")

cox_fit <- readRDS("basic_scripts/cox_fit")
data <- readRDS("basic_scripts/data")
lels <- unique(data$level) %>% sort(.,decreasing = T)
data$level <- factor(data$level, levels = lels)
new_data <- with(data,data.frame(level = c("low", "high")))

km.fit <- survfit(Surv(survival_days, censoring_status) ~ level, data = data)
km.surv <- ggsurvplot(km.fit, data=data, risk.table = TRUE, palette = "jco")

cox.fit <- survfit(cox_fit,newdata=new_data)

cox.surv <- ggsurvplot(cox.fit,data=new_data,
           title = "Survival Curves",
           xlab = "Days",
           ylab = "Survival probability",
           conf.int=TRUE,
           # pval = TRUE,    # Add p-value
           # surv.median.line = "hv",            # Add median survival lines
           # legend.title = call_datatype(x),               # Change legend titles
           legend.labs = c("Low", "High"),  # Change legend labels
           palette = "jco"                    # Use JCO journal color palette
           # risk.table = T,                  # Add No at risk table
           # cumevents = TRUE,                   # Add cumulative No of events table
           # tables.height = 0.15,               # Specify tables height
           # tables.theme = theme_cleantable(),  # Clean theme for tables
           # tables.y.text = FALSE               # Hide tables y axis text
)
cox.surv$table <- km.surv$table
