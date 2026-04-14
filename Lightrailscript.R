library(readxl)
library(readxl)
library(tidyverse)
library(mgcv)
library(stargazer)
library(lmtest)
library(sandwich)
library(coefplot)
library(clubSandwich)
lrt0101 <- read_excel("lrt0101.xlsx", sheet = "LRT0101", 
                      skip = 7)
View(lrt0101)
lrt_panel <- lrt0101 %>%
  pivot_longer(
    cols = -`Year ending March`, 
    names_to = "transport_system", 
    values_to = "passenger_journeys",
    # This force-converts everything to character so they can be combined
    values_transform = list(passenger_journeys = as.character)
  ) %>%
  mutate(
    # Clean the system names
    transport_system = str_remove_all(transport_system, " \\[note \\d+\\]"),
    # Convert [w] to NA, then turn the whole column into numbers
    passenger_journeys = as.numeric(na_if(passenger_journeys, "[w]"))
  )

# Check the first few rows
head(lrt_panel)
intervention <- c(2019,2020,2021,2022,2023)
treatment <- c("London Trams","Docklands Light Railway","London Underground")
control <- c("Manchester Metrolink","Glasgow Subway","Sheffield Supertram","Tyne and Wear Metro","Edinburgh Trams")
ggplot(lrt_panel,mapping = aes(x=`Year ending March`,y=passenger_journeys,colour = factor(transport_system)))+stat_summary(fun=mean,geom="line")+geom_vline(xintercept = 2019,linetype="dashed")+ ggtitle("Means with a policy breakpoint(journeys)")
ggplot(lrt_panel,mapping = aes(x=`Year ending March`,y=passenger_journeys,colour = factor(transport_system)))+geom_point()+ geom_point()+geom_smooth(method = "gam",se=FALSE)+geom_vline(xintercept = 2019,linetype="dashed")
panel_dfindlr <- lrt_panel %>%
  mutate(
    intervention = ifelse(`Year ending March` %in% intervention, 1,0),
    control = ifelse(transport_system %in% control, 1, 0),
    treatment= ifelse(transport_system %in% treatment,1,0),
    did= intervention*treatment)
treatment_control_df <- panel_dfindlr %>%
  filter(treatment + control == 1,`Year ending March` < 2024&`Year ending March` > 2000)
ggplot(treatment_control_df,aes(x=`Year ending March`,y=log(passenger_journeys),colour = factor(treatment)))+geom_smooth(method = "loess",se=TRUE)+geom_vline(xintercept = 2019,linetype="dashed")+geom_jitter()+ggtitle("Semiparametric Conditional Means with a policy breakpoint(logjourneys)")
ggplot(treatment_control_df,aes(x=`Year ending March`,y=log(passenger_journeys),colour = factor(treatment)))+stat_summary(fun=mean,geom="line")+geom_vline(xintercept = 2019,linetype="dashed")+  ggtitle("Means with a policy breakpoint(logjourneys)")
ggplot(treatment_control_df,aes(x=`Year ending March`,y=log(passenger_journeys),colour = factor(interaction(treatment,intervention))))+geom_smooth(method = "loess",se=FALSE)+geom_vline(xintercept = 2019,linetype="dashed")+geom_jitter()+ggtitle("Semiparametric Conditional Means with a policy breakpoint(logjourneys)")
ggplot(treatment_control_df,aes(x=`Year ending March`,y=passenger_journeys,colour = factor(interaction(treatment,intervention))))+geom_smooth(method = "loess",se=FALSE)+geom_vline(xintercept = 2019,linetype="dashed")+geom_jitter()+ggtitle("Semiparametric Conditional Means with a policy breakpoint(logjourneys)")
eventmodel <- lm(log(passenger_journeys)~factor(treatment)*factor(`Year ending March`)+factor(transport_system),data=treatment_control_df)
eventmodelnl <- lm(passenger_journeys~factor(treatment)*factor(`Year ending March`)+factor(transport_system)+factor(`Year ending March`),data=treatment_control_df)
IQR(log(treatment_control_df$passenger_journeys),na.rm = TRUE)
hc3_vcov <- vcovHC(eventmodel, type = "HC3")
cluster_vcov <- vcovCR(eventmodel,type = "CR2",cluster = treatment_control_df$transport_system)
#Event Plot Code#
all_coefs <- names(coef(eventmodel))
qtr_coefs <- all_coefs[grep("factor\\(treatment\\)1:factor\\(`Year ending March`\\)", all_coefs)]
qtr_coefs <- qtr_coefs[order(gsub(".*\\)", "", qtr_coefs))]
coefplot::coefplot(eventmodel, 
         coefficients = qtr_coefs, vcov = cluster_vcov,innerCI=0,outerCI=1.96,horizontal = TRUE) + 
  scale_y_discrete(labels = function(x) {
    clean <- gsub("factor\\(`Year ending March`\\)|factor\\(treatment\\)1:factor\\(`Year ending March`\\)", "", x)
  })
  theme_gray()
  stargazer(eventmodel,se=list(sqrt(diag(cluster_vcov))),omit = c("^factor\\(`Year ending March`\\)","factor\\(transport_system\\)",type = "text"))
  
            