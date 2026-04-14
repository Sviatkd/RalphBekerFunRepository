library(readxl)

traffic_flow_borough <- read_excel("traffic-flow-borough.xlsx", 
                                   sheet = "Traffic Flows - Cars")
View(traffic_flow_borough)
library(tidyverse)
library(mgcv)
library(stargazer)
traffic_panel <- traffic_flow_borough %>%
  # 1. Remove the empty first row if it's all NA
  filter(!is.na(`LA Code`)) %>%
  
  # 2. Pivot all year columns (from "1993" to the last column)
  pivot_longer(
    cols = `1993`:`2024`, 
    names_to = "Year", 
    values_to = "Traffic_Flow"
  ) %>%
  
  # 3. Clean the Year column (remove "[note X]" and convert to numeric)
  mutate(
    Year = str_extract(Year, "\\d{4}"), # Extracts just the 4 digits
    Year = as.numeric(Year)
  )
traffic_panel <- traffic_panel %>%
  group_by(`Local Authority`) %>%
  arrange(Year) %>%
  mutate(
    log_diff = log(Traffic_Flow) - log(lag(Traffic_Flow))
  ) %>%
  ungroup()

# View the result
print(traffic_panel)
ggplot(traffic_panel,mapping=aes(x=Year,y=log_diff,colour=`Local Authority`))+stat_summary(fun=mean,geom="point")+stat_summary(fun=mean,geom="line")+geom_vline(xintercept = 2019,linetype="dashed")+ ggtitle("Means with a policy breakpoint(Traffic Volume)")    
unique(traffic_panel$`Local Authority`)    
intervention <- c(2019,2020,2021,2022,2023,2024)
treatment <- c("London")
control <- c("North West","West Midlands")
panel_dfindlr <- traffic_panel %>%
  mutate(
    intervention = ifelse(Year %in% intervention, 1,0),
    control = ifelse(`Local Authority` %in% control, 1, 0),
    treatment= ifelse(`Local Authority` %in% treatment,1,0),
    did= intervention*treatment)
treatment_control_df <- panel_dfindlr %>%
  arrange(Year) %>%
  mutate(d_Traffic = Traffic_Flow - lag(Traffic_Flow))%>%
  filter(treatment== 1,Year <= 2024& Year > 2000)
  
treatment_control_df$Covid_Shock <- ifelse(treatment_control_df$Year %in% c(2020, 2021), 1, 0)
ggplot(treatment_control_df,mapping=aes(x=Year,y=log_diff,colour=factor(`Local Authority`)))+stat_summary(fun=mean,geom="point")+stat_summary(fun=mean,geom="line")+geom_vline(xintercept = 2019,linetype="dashed")+ ggtitle("Means with a policy breakpoint(Traffic Flows)")
ggplot(treatment_control_df,mapping=aes(x=Year,y=log_diff,colour=factor(interaction(treatment,intervention))))+geom_smooth(method="loess",se=FALSE)+geom_jitter()+geom_vline(xintercept = 2019,linetype="dashed")+ggtitle("Treatment and Control Group Local Regression (Traffic Flows)")
ggplot(treatment_control_df,mapping=aes(x=Year,y=log_diff,colour=factor(treatment)))+stat_summary(fun=mean,geom="point")+stat_summary(fun=mean,geom="line")+geom_vline(xintercept = 2019,linetype="dashed")+ ggtitle("Means with a policy breakpoint(Traffic Flows)")
breakmodel <- lm(Traffic_Flow~factor(intervention)*Year+Covid_Shock,data=treatment_control_df)
breakmodelldf <- lm(log_diff~factor(intervention)*Year+Covid_Shock,data=treatment_control_df)
nobreakmodel <- lm(Traffic_Flow~Year+Covid_Shock,data=treatment_control_df)
eventmodel <- lm(Traffic_Flow~factor(treatment)*factor(Year)+factor(`Local Authority`)+factor(Year),data=treatment_control_df)
IQR(treatment_control_df$Traffic_Flow,na.rm = TRUE)
hc3_vcov <- vcovHC(breakmodel, type = "HC3")
cluster_vcov <- vcovCR(eventmodel,type = "CR2",cluster = treatment_control_df$`Local Authority`)
#Event Plot Code#
all_coefs <- names(coef(eventmodel))
qtr_coefs <- all_coefs[grep("factor\\(treatment\\)1:factor\\(Year\\)", all_coefs)]
qtr_coefs <- qtr_coefs[order(gsub(".*\\)", "", qtr_coefs))]
coefplot(eventmodel, 
         coefficients = qtr_coefs, vcov = cluster_vcov,
         innerCI=0,outerCI=0.95, horizontal = TRUE) + 
  
  # The "Thinning" Logic
  scale_y_discrete(labels = function(x) {
    clean <- gsub("factor\\(Year\\)|factor\\(treatment\\)1:factor\\(Year\\)", "", x)
  })
theme_gray()
coeftest(breakmodel,vcov. = hc3_vcov)
BIC(breakmodel,nobreakmodel)

