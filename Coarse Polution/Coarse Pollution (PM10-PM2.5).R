rm(list = ls())

# =============================================================================
# PM10-PM2.5 (Non-Combustion Coarse)
# =============================================================================

library(readxl)
library(tidyverse)
library(stargazer)
library(lmtest)
library(sandwich)
library(clubSandwich)
library(aod)
library(modelsummary)
library(broom)


### DATA

FullDataset <- read_csv("Full_Air_Quality_Dataset_Merged.csv")

pm25quarterly <- FullDataset %>%
  mutate(yr_qtr = quarter(date, with_year = TRUE)) %>%
  group_by(yr_qtr, Station) %>%
  summarise(mean_pm25 = mean(pm25, na.rm = TRUE), .groups = "drop") %>%
  arrange(Station, yr_qtr)

pm10quarterly <- FullDataset %>%
  mutate(yr_qtr = quarter(date, with_year = TRUE)) %>%
  group_by(yr_qtr, Station) %>%
  summarise(mean_pm10 = mean(pm10, na.rm = TRUE), .groups = "drop") %>%
  arrange(Station, yr_qtr)

no2quarterly <- FullDataset %>%
  mutate(yr_qtr = quarter(date, with_year = TRUE)) %>% 
  group_by(yr_qtr, Station) %>%
  summarise(mean_value = mean(no2, na.rm = TRUE), .groups = "drop") %>%
  arrange(Station, yr_qtr)


# Merge and compute the difference index
coarse_quarterly <- left_join(pm10quarterly, pm25quarterly,
                              by = c("yr_qtr", "Station")) %>%
  mutate(
    pmcoarse_idx = mean_pm10 - mean_pm25   # non-combustion proxy
  )

# Add City labels
coarse_quarterly <- coarse_quarterly %>%
  mutate(City = case_when(
    str_detect(Station, regex("Birmingham",  ignore_case = TRUE)) ~ "Birmingham",
    str_detect(Station, regex("Glasgow",     ignore_case = TRUE)) ~ "Glasgow",
    str_detect(Station, regex("Leeds",       ignore_case = TRUE)) ~ "Leeds",
    str_detect(Station, regex("London",      ignore_case = TRUE)) ~ "London",
    str_detect(Station, regex("Manchester",  ignore_case = TRUE)) ~ "Manchester",
    str_detect(Station, regex("Newcastle",   ignore_case = TRUE)) ~ "Newcastle",
    str_detect(Station, regex("Sheffield",   ignore_case = TRUE)) ~ "Sheffield",
    str_detect(Station, regex("Camden",      ignore_case = TRUE)) ~ "London",
    str_detect(Station, regex("Salford",     ignore_case = TRUE)) ~ "Salford",
    str_detect(Station, regex("Sunderland",  ignore_case = TRUE)) ~ "Sunderland",
    str_detect(Station, regex("Westminster", ignore_case = TRUE)) ~ "London",
    TRUE ~ "Other"
  ))

# Intervention / treatment / control flags
intervention_periods <- c(2019.2,2019.3,2019.4,
                          2020.1,2020.2,2020.3,2020.4,
                          2021.1,2021.2,2021.3,2021.4,
                          2022.1,2022.2,2022.3,2022.4)

treatment_stations <- c(
  "camden-- euston road",
  "london-bloomsbury",
  "LondonC",
  "LondonMR"
  # LondonW excluded for PM10 / PM2.5 – no observations
)

control_stations <- c(
  "glasgow-anderston", "glasgow-byres road", "glasgow-high street",
  "glasgow-townhead",  "GlasgowK",
  "manchester-oxford road", "ManchesterPic",
  "NewcastleC",
  "salford-eccles", "salford-m60"
)

treatment_control_df <- no2quarterly %>%
  mutate(
    intervention = ifelse(yr_qtr %in% intervention_periods, 1, 0),
    control      = ifelse(Station %in% control_stations,    1, 0),
    treatment    = ifelse(Station %in% treatment_stations,  1, 0),
    did          = intervention * treatment
  ) %>%
  filter(treatment + control == 1,
         yr_qtr < 2023.1) %>%
  mutate(
    Station = as.factor(Station),
    yr_qtr  = as.factor(yr_qtr),
    did     = as.factor(did)
  )

panel_coarse <- coarse_quarterly %>%
  mutate(
    intervention = ifelse(yr_qtr %in% intervention_periods, 1, 0),
    control      = ifelse(Station %in% control_stations,    1, 0),
    treatment    = ifelse(Station %in% treatment_stations,  1, 0),
    did          = intervention * treatment
  )

# Balanced panel: treatment + control only, up to 2022 Q4
# Pre-period starts 2015 Q2
tc_coarse <- panel_coarse %>%
  filter(treatment + control == 1,
         yr_qtr >= 2015.2,
         yr_qtr <  2023.1)

# Do the same for PM2.5 on its own (long-format, mean_value = mean_pm25)
tc_pm25 <- tc_coarse %>%
  select(yr_qtr, Station, City,
         mean_value = mean_pm25,
         intervention, control, treatment, did)

# And for the coarse index
tc_pmcoarse <- tc_coarse %>%
  select(yr_qtr, Station, City,
         mean_value = pmcoarse_idx,
         intervention, control, treatment, did)

# Factor-ise for modelling
for (df_name in c("tc_pm25", "tc_pmcoarse")) {
  df <- get(df_name)
  df$Station      <- as.factor(df$Station)
  df$yr_qtr       <- as.factor(df$yr_qtr)
  df$City         <- as.factor(df$City)
  df$did          <- as.factor(df$did)
  assign(df_name, df)
}

### Event Study + Heterogeneous Station Models for PM2.5 and PMcoarse

# ── NO2 ──────────────────────────────────────────────────
reducedmodel <- lm(mean_value ~ did + yr_qtr + Station,
                   data = treatment_control_df)

hetmodel <- lm(mean_value ~ did * Station + yr_qtr,
               data = treatment_control_df)

eventmodel <- lm(mean_value ~ treatment * yr_qtr + Station,
                 data = treatment_control_df)

hcr_vcov_het_no2 <- vcovHC(hetmodel, type = "HC3")

cluster_vcovr <- vcovCR(
  eventmodel,
  type = "CR2",
  cluster = treatment_control_df$Station
)

# ── PM2.5 ────────────────────────────────────────────────────────────────────
eventmodel_pm25  <- lm(mean_value ~ treatment * yr_qtr + Station, data = tc_pm25)
hetmodel_pm25    <- lm(mean_value ~ did * Station + yr_qtr,        data = tc_pm25)

hcr_vcov_het_pm25 <- vcovHC(hetmodel_pm25, type = "HC3")

cluster_vcov_event_pm25 <- vcovCR(eventmodel_pm25, type = "CR2",
                                  cluster = tc_pm25$City)
cluster_vcov_het_pm25   <- vcovCR(hetmodel_pm25,   type = "CR2",
                                  cluster = tc_pm25$City)

# ── PMcoarse (PM10 − PM2.5) ──────────────────────────────────────────────────
eventmodel_pmc  <- lm(mean_value ~ treatment * yr_qtr + Station, data = tc_pmcoarse)
hetmodel_pmc    <- lm(mean_value ~ did * Station + yr_qtr,        data = tc_pmcoarse)

hcr_vcov_het_pmc <- vcovHC(hetmodel_pmc, type = "HC3")

cluster_vcov_event_pmc <- vcovCR(eventmodel_pmc, type = "CR2",
                                 cluster = tc_pmcoarse$City)
cluster_vcov_het_pmc   <- vcovCR(hetmodel_pmc,   type = "CR2",
                                 cluster = tc_pmcoarse$City)


### Visuals for PMcoarse (PM10 − PM2.5)

# 3a. Overall time trend with policy break
ggplot(tc_pmcoarse,
       aes(x = as.numeric(as.character(yr_qtr)),
           y = mean_value, colour = factor(Station))) +
  stat_summary(fun = mean, geom = "line") +
  geom_vline(xintercept = 2019.2, linetype = "dashed") +
  labs(title = "PM10 − PM2.5 Index (Non-Combustion Proxy)",
       x = "Year.Quarter", y = "PMcoarse Index",
       colour = "Station") +
  theme_minimal()

# 3b. Density plots: treatment vs control, pre vs post
ggplot(tc_pmcoarse,
       aes(x = mean_value,
           fill  = factor(intervention),
           colour = factor(intervention))) +
  geom_density(alpha = 0.3) +
  facet_wrap(~ treatment,
             labeller = as_labeller(c("0" = "Control", "1" = "Treatment"))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_minimal() +
  labs(title    = "Comparison of Densities (KDE) – PM10 − PM2.5 Index",
       subtitle = "Treatment vs Control, Pre (pink) / Post (teal) intervention",
       x = "Mean Value", y = "Density",
       fill = "Intervention", colour = "Intervention") +
  theme(legend.position = "bottom")

# 3c. Station-level density for treated stations
ggplot(filter(tc_pmcoarse, treatment == 1),
       aes(x = mean_value,
           fill  = factor(intervention),
           colour = factor(intervention))) +
  geom_density(alpha = 0.3) +
  facet_wrap(~ Station) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_minimal() +
  labs(title    = "Densities – PM10 − PM2.5, Treated London Stations",
       x = "Mean Value", y = "Density",
       fill = "Intervention", colour = "Intervention") +
  theme(legend.position = "bottom")

tc_pmcoarse <- tc_pmcoarse %>%
  mutate(yr_qtr_num = as.numeric(as.character(yr_qtr)))

# 3d. Smooth counterfactual plot (LOESS) – PMcoarse
control_loess_pmc <- loess(mean_value ~ yr_qtr_num,
                           data = filter(tc_pmcoarse, treatment == 0),
                           span = 0.75)

qtr_grid_pmc <- seq(min(tc_pmcoarse$yr_qtr_num),
                    max(tc_pmcoarse$yr_qtr_num),
                    by = 0.1)

ctrl_trend_pmc <- data.frame(yr_qtr_num = qtr_grid_pmc) %>%
  mutate(pred_val = predict(control_loess_pmc,
                            newdata = data.frame(yr_qtr_num = yr_qtr_num)),
         drift    = pred_val - first(pred_val))

treated_starts_pmc <- tc_pmcoarse %>%
  filter(treatment == 1,
         as.numeric(as.character(yr_qtr)) < 2019.2) %>%
  group_by(Station) %>%
  summarise(start_mean = mean(mean_value, na.rm = TRUE))

loess_cf_pmc <- treated_starts_pmc %>%
  crossing(ctrl_trend_pmc) %>%
  mutate(counterfactual_loess = start_mean + drift)

ggplot() +
  geom_smooth(data = filter(tc_pmcoarse, treatment == 1),
              aes(x = as.numeric(as.character(yr_qtr)),
                  y = mean_value, colour = Station),
              method = "loess", se = FALSE, linewidth = 1) +
  geom_line(data = loess_cf_pmc,
            aes(x = yr_qtr_num,
                y = counterfactual_loess, colour = Station),
            linetype = "dashed", alpha = 0.7, linewidth = 0.8) +
  geom_vline(xintercept = 2019.2, linetype = "dotted", linewidth = 1) +
  labs(title    = "PM10 − PM2.5: Treated Stations vs Smooth Counterfactuals",
       subtitle = "Solid = Actual LOESS | Dashed = Control trend shifted to station baseline",
       x = "Year.Quarter", y = "PMcoarse Index",
       caption  = "Divergence after dotted line = estimated policy impact") +
  theme_minimal()

# 3e. Event study coefficient plot – PMcoarse, City-level CR2
all_coefs_pmc  <- names(coef(eventmodel_pmc))
qtr_coefs_pmc  <- all_coefs_pmc[grep("treatment:yr_qtr", all_coefs_pmc)]

modelplot(eventmodel_pmc,
          vcov       = cluster_vcov_event_pmc,
          coef_omit  = "^(?!.*treatment:)",
          color      = "darkblue") +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = "treatment:yr_qtr2019.2",
             linetype = "dashed", colour = "red") +
  labs(title = "PM10 − PM2.5 Event Study (2019.2), CR2 City-Clustered SE") +
  theme_classic() +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 90))

# Same plot for PM2.5 itself
modelplot(eventmodel_pm25,
          vcov       = cluster_vcov_event_pm25,
          coef_omit  = "^(?!.*treatment:)",
          color      = "darkgreen") +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = "treatment:yr_qtr2019.2",
             linetype = "dashed", colour = "red") +
  labs(title = "PM2.5 Event Study (2019.2), CR2 City-Clustered SE") +
  theme_classic() +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 90))


# Mega modelsummary Table

# Model list
models_mega <- list(
  "NO2 Het"        = hetmodel,
  "NO2 Event"      = eventmodel,
  "PM2.5 Het"      = hetmodel_pm25,
  "PM2.5 Event"    = eventmodel_pm25,
  "PMcoarse Het"   = hetmodel_pmc,
  "PMcoarse Event" = eventmodel_pmc
)

# Vcov list in same order
vcovs_mega <- list(
  hcr_vcov_het_no2,
  cluster_vcovr,
  hcr_vcov_het_pm25,
  cluster_vcov_event_pm25,
  hcr_vcov_het_pmc,
  cluster_vcov_event_pmc
)

# Final table
modelsummary(
  models_mega,
  vcov = vcovs_mega,
  coef_omit = "Station|yr_qtr",
  gof_map = c("nobs", "r.squared", "adj.r.squared", "aic", "bic"),
  stars = TRUE,
  title = "ULEZ Effect on Air Pollution – All Pollutants"
)

