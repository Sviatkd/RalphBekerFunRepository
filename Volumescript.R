library(readxl)
library(readxl)
tra8902 <- read_excel("C:/Users/alexa/Downloads/tra8902-miles-by-local-authority-and-selected-vehicle-type.xlsx", 
                      sheet = "TRA8902", skip = 4)
View(tra8902_miles_by_local_authority_and_selected_vehicle_type)
library(tidyverse)
library(mgcv)
library(stargazer)
traffic_panel_full <- tra8902 %>%
  # 1. Pivot the year columns (from '1993' to the last column)
  # We use matches("\\d{4}") to catch any column starting with 4 digits
  pivot_longer(
    cols = matches("^\\d{4}"), 
    names_to = "Year", 
    values_to = "Traffic_Volume"
  ) %>%
  
  # 2. Clean the data
  mutate(
    # Clean the Year column (remove [note XX])
    Year = as.numeric(str_extract(Year, "\\d{4}")),
    
    # Convert Traffic_Volume to numeric (important since it's currently <chr>)
    # Note: This will turn any non-numeric placeholders into NA
    Traffic_Volume = as.numeric(Traffic_Volume)
  ) %>%
  
  # 3. Optional: Remove rows where Traffic_Volume is NA
  filter(!is.na(Traffic_Volume))
traffic_panel_full <- traffic_panel_full %>%
  group_by(`Local Authority`, Vehicle) %>%
  arrange(Year) %>%
  mutate(
    log_diff = log(Traffic_Volume) - log(lag(Traffic_Volume))
  ) %>%
  ungroup()
# Look at the first few rows
treatment <- c("London All")
control <- c("Greater Manchester ITA All","Tyne and Wear ITA All","Glasgow City")
panel_dfindlr <- traffic_panel_full %>%
  mutate(
    intervention = ifelse(Year %in% intervention, 1,0),
    control = ifelse(`Local Authority` %in% control, 1, 0),
    treatment= ifelse(`Local Authority` %in% treatment,1,0),
    did= intervention*treatment)
treatment_control_df <- panel_dfindlr %>%
  filter(treatment + control == 1,Year < 2024& Year > 2000)
ggplot(treatment_control_df, aes(x = log(Traffic_Volume), fill = factor(intervention), color = factor(intervention))) +
  geom_density(alpha = 0.3) +
  # Use facet_wrap to create the side-by-side view
  facet_wrap(~treatment, labeller = as_labeller(c("0" = "Control", "1" = "Treatment"))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_minimal() +
  labs(
    title = "Comparison of Densities(KDE)",
    subtitle = "Side-by-side density plots of Treatment vs Control",
    x = "log(Volume)",
    y = "Density",
    fill = "Intervention",
    color = "Intervention"
  ) +
  theme(legend.position = "bottom")
ggplot(treatment_control_df, aes(x = log(Traffic_Volume), fill = factor(intervention), color = factor(intervention))) +
  geom_density(alpha = 0.3) +
  # Use facet_wrap to create the side-by-side view
  facet_wrap(~`Local Authority`) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_minimal() +
  labs(
    title = "Comparison of Densities(KDE)Traffic Volume",
    subtitle = "Side-by-side density plots before and after ULEZ",
    x = "log(Volume)",
    y = "Density",
    fill = "Intervention",
    color = "Intervention"
  ) +
  theme(legend.position = "bottom")
ggplot(treatment_control_df,mapping=aes(x=Year,y=log_diff,colour=factor(`Local Authority`)))+stat_summary(fun=mean,geom="point")+stat_summary(fun=mean,geom="line")+geom_vline(xintercept = 2019,linetype="dashed")+ ggtitle("Means with a policy breakpoint(Traffic Volume)")
ggplot(treatment_control_df,mapping=aes(x=Year,y=log_diff,colour=factor(interaction(treatment,intervention))))+geom_smooth(method="loess",se=FALSE)+geom_jitter()+geom_vline(xintercept = 2019,linetype="dashed")+ggtitle("Treatment and Control Group Local Regression (Traffic Volume)")
ggplot(treatment_control_df,mapping=aes(x=Year,y=log_diff,colour=factor(treatment)))+stat_summary(fun=mean,geom="point")+stat_summary(fun=mean,geom="line")+geom_vline(xintercept = 2019,linetype="dashed")+ ggtitle("Means with a policy breakpoint(Traffic Volume)")
breakmodel <- lm(Traffic_Volume~factor(intervention)*Year+Covid_Shock,data=treatment_control_df)
nobreakmodel <- lm(Traffic_Volume~Year+Covid_Shock,data=treatment_control_df)
eventmodel <- lm(Traffic_Volume~factor(treatment)*factor(Year)+factor(`Local Authority`)+factor(Year),data=treatment_control_df)
IQR(treatment_control_df$Traffic_Volume,na.rm = TRUE)
hc3_vcov <- vcovHC(eventmodel, type = "HC1")
cluster_vcov <- vcovCR(eventmodel,type = "CR2",cluster = treatment_control_df$`Local Authority`)
#Event Plot Code#
all_coefs <- names(coef(eventmodel))
qtr_coefs <- all_coefs[grep("factor\\(treatment\\)1:factor\\(Year\\)", all_coefs)]
qtr_coefs <- qtr_coefs[order(gsub(".*\\)", "", qtr_coefs))]
coefplot::coefplot(eventmodel, 
         coefficients = qtr_coefs, vcov = hc3_vcov,
         innerCI=0,outerCI=0.95, horizontal = TRUE) + 
  
  # The "Thinning" Logic
  scale_y_discrete(labels = function(x) {
    clean <- gsub("factor\\(Year\\)|factor\\(treatment\\)1:factor\\(Year\\)", "", x)
  })
theme_gray()
coeftest(breakmodel,vcov. = hc3_vcov)
BIC(breakmodel,nobreakmodel)
coeftest(eventmodel,vcov. = hc1_vcov)
stargazer(eventmodel,se=list(sqrt(diag(cluster_vcov))),type = "text")
