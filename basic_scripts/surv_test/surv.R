library(survival) # to do the survival analysis
library(survminer) # to plot the survival analysis nicer

setwd("/Users/jeancheng/Documents/cSurvival/")

# ------- 1. basic cox and km fit ---------
cox_fit <- readRDS("basic_scripts/cox_fit")
data <- readRDS("basic_scripts/data")
lels <- unique(data$level) %>% sort(.,decreasing = T)
data$level <- factor(data$level, levels = lels)
new_data <- with(data,data.frame(level = c("low", "high")))

km.stats <- survdiff(Surv(survival_days, censoring_status) ~ level, data = data)
p <- 1 - pchisq(km.stats$chisq, length(km.stats$n) - 1)
km.fit <- survfit(Surv(survival_days, censoring_status) ~ level, data = data)
km.surv <- ggsurvplot(km.fit, data=data, 
                      xlab = "Days",
                      ylab = "Survival probability",
                      legend.labs = c("Low", "High"),  # Change legend labels
                      risk.table = TRUE, 
                      cumevents = F,                   # Add cumulative No of events table
                      ggtheme = theme_survminer(
                        base_size = 18,
                        font.main = c(20, "plain", "black"),
                        font.submain = c(18, "plain", "black"),
                        font.x = c(18, "plain", "black"),
                        font.y = c(18, "plain", "black"),
                        font.caption = c(18, "plain", "black"),
                        font.tickslab = c(16, "plain", "black"),
                        # legend = c("top", "bottom", "left", "right", "none"),
                        font.legend = c(18, "plain", "black")
                      ),
                      palette = "jco")

km.surv

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
           ggtheme = theme_survminer(
             base_size = 18,
             font.main = c(20, "plain", "black"),
             font.submain = c(18, "plain", "black"),
             font.x = c(18, "plain", "black"),
             font.y = c(18, "plain", "black"),
             font.caption = c(18, "plain", "black"),
             font.tickslab = c(16, "plain", "black"),
             # legend = c("top", "bottom", "left", "right", "none"),
             font.legend = c(18, "plain", "black")
           ),
           palette = "jco"                    # Use JCO journal color palette
           # risk.table = T,                  # Add No at risk table
           # cumevents = TRUE,                   # Add cumulative No of events table
           # tables.height = 0.15,               # Specify tables height
           # tables.theme = theme_cleantable(),  # Clean theme for tables
           # tables.y.text = FALSE               # Hide tables y axis text
)

# ggsave("test.pdf",print(cox.surv))
pdf("test.pdf",onefile = TRUE)
print(cox.surv,newpage = FALSE)
dev.off()

# cox.surv$table <- km.surv$table + theme_cleantable(
#   base_size = 18,
#   font.main = c(18, "plain", "black"),
#   font.submain = c(18, "plain", "black"),
#   font.caption = c(18, "plain", "black"),
#   font.tickslab = c(16, "plain", "black"),
#   # legend = c("top", "bottom", "left", "right", "none"),
#   font.legend = c(18, "plain", "black")
# )
# print(cox.surv,risk.table.height = 0.3)

# ----------- 2. categorical interaction analysis --------
df1 <- readRDS("basic_scripts/df_1")
df2 <- readRDS("basic_scripts/df_2")
df_list <- list(df1,df2)
df_combined <- Reduce(
  function(x, y, ...) inner_join(x, select(y, patient_id, level), by = "patient_id"), 
  df_list
)

x_y <- c("x","y")[1:length(df_list)]

df_combined[["level"]] <- apply(df_combined %>% select(paste0("level.",x_y)),1,paste0,collapse="_")
lels <- unique(df_combined$level) %>% sort(.,decreasing = T)
df_combined$level <- factor(df_combined$level, levels = lels)
lels <- levels(df_combined$level)

res_all <- readRDS("basic_scripts/cox_all")
res_cox <- res_all[["cox"]]
res_km <- res_all[["km"]]

res.cox <- coxph(Surv(survival_days, censoring_status) ~ level.x + level.y + level.x*level.y, data =  df_combined)
summary(res.cox)

# ------------- 3. continuous variable analysis ---------------
df_exp1 <- readRDS("basic_scripts/df_exp1")

