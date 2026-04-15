# ============================================================
# Cycling Demand - Difference-in-Differences (DiD) Analysis
# Data: TRA0403 - Pedal cycle traffic (billion vehicle miles)
# Treatment: London (ULEZ introduced 2019)
# Control: Scotland, North West
# Policy cutoff: 2019
# SE Method: HC3 (Jackknife) - preferred for small samples
#            with high leverage points and noisy data
# ============================================================

# --- 1. Install & Load Libraries ---
if (!require("readODS"))      install.packages("readODS")
if (!require("tidyverse"))    install.packages("tidyverse")
if (!require("lmtest"))       install.packages("lmtest")
if (!require("sandwich"))     install.packages("sandwich")
if (!require("clubSandwich")) install.packages("clubSandwich")
if (!require("stargazer"))    install.packages("stargazer")
if (!require("ggplot2"))      install.packages("ggplot2")
if (!require("skimr"))        install.packages("skimr")

library(readODS)
library(tidyverse)
library(lmtest)
library(sandwich)
library(clubSandwich)
library(stargazer)
library(ggplot2)
library(skimr)

# ============================================================
# 2. Set Working Directory & Load Data
# ============================================================

setwd(".")   # Change to your folder path if needed

raw <- read_ods("tra0403-miles-pedal-cycle-traffic-by-region.ods",
                sheet = "TRA0403",
                skip  = 3)

# ============================================================
# 3. Clean & Reshape Data
# ============================================================

colnames(raw)[1] <- "Region"
colnames(raw)[2] <- "Subarea"
colnames(raw)[3] <- "Notes"
colnames(raw)[4] <- "Units"
colnames(raw)[5:ncol(raw)] <- as.character(1993:2024)

# Keep only "All" summary rows for London, North West, Scotland
regions_keep <- c("London", "North West", "Scotland")

df_wide <- raw %>%
  filter(Region %in% regions_keep,
         grepl("All", Subarea, ignore.case = TRUE)) %>%
  select(Region, all_of(as.character(1993:2024)))

# Force all year columns to numeric BEFORE pivoting (fixes mixed type error)
df_wide <- df_wide %>%
  mutate(across(as.character(1993:2024), as.numeric))

# Pivot to long format
df_long <- df_wide %>%
  pivot_longer(cols      = as.character(1993:2024),
               names_to  = "Year",
               values_to = "CycleTraffic_bn") %>%
  mutate(Year            = as.integer(Year),
         CycleTraffic_bn = as.numeric(CycleTraffic_bn))

# ============================================================
# 4. DiD Variable Construction
# ============================================================

df_did <- df_long %>%
  mutate(
    treated   = if_else(Region == "London", 1L, 0L),
    post      = if_else(Year >= 2019, 1L, 0L),
    did       = treated * post,
    log_cycle = log(CycleTraffic_bn)
  )

# ============================================================
# 5. Summary Statistics
# ============================================================

cat("\n========== FULL SUMMARY (skimr) ==========\n")
skim(df_did)

cat("\n========== SUMMARY STATS BY REGION ==========\n")
df_did %>%
  group_by(Region) %>%
  summarise(
    Mean = round(mean(CycleTraffic_bn, na.rm = TRUE), 3),
    SD   = round(sd(CycleTraffic_bn,   na.rm = TRUE), 3),
    Min  = round(min(CycleTraffic_bn,  na.rm = TRUE), 3),
    Max  = round(max(CycleTraffic_bn,  na.rm = TRUE), 3),
    N    = n()
  ) %>%
  print()

cat("\n========== SUMMARY STATS PRE vs POST 2019 ==========\n")
df_did %>%
  group_by(Region, post) %>%
  summarise(
    Mean = round(mean(CycleTraffic_bn, na.rm = TRUE), 3),
    SD   = round(sd(CycleTraffic_bn,   na.rm = TRUE), 3),
    N    = n(),
    .groups = "drop"
  ) %>%
  mutate(Period = if_else(post == 1, "Post-2019", "Pre-2019")) %>%
  select(Region, Period, Mean, SD, N) %>%
  print()

cat("\n========== SUMMARY STATS (ulez_stats style) ==========\n")
cycling_stats <- df_did %>%
  mutate(ULEZ = if_else(post == 1, "Post-ULEZ", "Pre-ULEZ")) %>%
  group_by(Region, ULEZ) %>%
  summarise(
    Mean   = round(mean(CycleTraffic_bn,   na.rm = TRUE), 3),
    Median = round(median(CycleTraffic_bn, na.rm = TRUE), 3),
    SD     = round(sd(CycleTraffic_bn,     na.rm = TRUE), 3),
    IQR    = round(IQR(CycleTraffic_bn,    na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  rename(Mode = Region) %>%
  arrange(Mode, desc(ULEZ))

print(cycling_stats)

# Save all summary stats to file
sink("cycling_summary_stats.txt")
skim(df_did)
df_did %>%
  group_by(Region) %>%
  summarise(Mean = round(mean(CycleTraffic_bn, na.rm=TRUE), 3),
            SD   = round(sd(CycleTraffic_bn,   na.rm=TRUE), 3),
            Min  = round(min(CycleTraffic_bn,  na.rm=TRUE), 3),
            Max  = round(max(CycleTraffic_bn,  na.rm=TRUE), 3),
            N    = n()) %>% print()
df_did %>%
  group_by(Region, post) %>%
  summarise(Mean = round(mean(CycleTraffic_bn, na.rm=TRUE), 3),
            SD   = round(sd(CycleTraffic_bn,   na.rm=TRUE), 3),
            N    = n(), .groups="drop") %>%
  mutate(Period = if_else(post==1, "Post-2019", "Pre-2019")) %>%
  select(Region, Period, Mean, SD, N) %>% print()
print(cycling_stats)
sink()

# ============================================================
# 6. Parallel Trends Plot (Pre-treatment visual check)
# ============================================================

pre_trends <- df_did %>%
  filter(Year <= 2019) %>%
  ggplot(aes(x = Year, y = CycleTraffic_bn,
             colour = Region, group = Region)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  geom_vline(xintercept = 2019, linetype = "dashed", colour = "black") +
  annotate("text", x = 2019.3, y = 0.85,
           label = "ULEZ 2019", hjust = 0, size = 3.5) +
  scale_colour_manual(values = c("London"     = "#1f77b4",
                                 "North West" = "#ff7f0e",
                                 "Scotland"   = "#2ca02c")) +
  labs(title    = "Parallel Trends Check - Cycle Traffic (Pre-2019)",
       subtitle = "Billion vehicle miles by region, annual",
       x = "Year", y = "Cycle Traffic (Billion Vehicle Miles)",
       colour = "Region") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

print(pre_trends)
ggsave("cycling_parallel_trends.png", pre_trends, width = 9, height = 5.5, dpi = 150)

# ============================================================
# 7. Full Time Series Plot
# ============================================================

full_plot <- df_did %>%
  ggplot(aes(x = Year, y = CycleTraffic_bn,
             colour = Region, group = Region)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  geom_vline(xintercept = 2019, linetype = "dashed", colour = "black") +
  annotate("text", x = 2019.3,
           y = max(df_did$CycleTraffic_bn, na.rm = TRUE) * 0.95,
           label = "ULEZ 2019", hjust = 0, size = 3.5) +
  scale_colour_manual(values = c("London"     = "#1f77b4",
                                 "North West" = "#ff7f0e",
                                 "Scotland"   = "#2ca02c")) +
  labs(title    = "Cycle Traffic Miles by Region (1993-2024)",
       subtitle = "Billion vehicle miles, annual",
       x = "Year", y = "Cycle Traffic (Billion Vehicle Miles)",
       colour = "Region") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

print(full_plot)
ggsave("cycling_full_timeseries.png", full_plot, width = 9, height = 5.5, dpi = 150)

# ============================================================
# 8. DiD Regressions
# NOTE: Using HC3 (Jackknife) SE adjustment as recommended.
# HC3 penalizes high leverage points and is preferred over HC1
# for noisy data with small number of independent observations.
# CR2 clustered SEs at Region level included as robustness check.
# Caution: only 3 clusters - interpret clustered SEs carefully.
# ============================================================

# Model 1: Basic DiD (levels)
model1 <- lm(CycleTraffic_bn ~ treated + post + did, data = df_did)

# Model 2: DiD (log levels)
model2 <- lm(log_cycle ~ treated + post + did, data = df_did)

# Model 3: DiD levels excluding COVID years (2020-2021)
df_no_covid <- df_did %>% filter(!Year %in% c(2020, 2021))
model3 <- lm(CycleTraffic_bn ~ treated + post + did, data = df_no_covid)

# Model 4: DiD log excluding COVID
model4 <- lm(log_cycle ~ treated + post + did, data = df_no_covid)

# --- HC3 (Jackknife) Standard Errors ---
se1_hc3 <- sqrt(diag(vcovHC(model1, type = "HC3")))
se2_hc3 <- sqrt(diag(vcovHC(model2, type = "HC3")))
se3_hc3 <- sqrt(diag(vcovHC(model3, type = "HC3")))
se4_hc3 <- sqrt(diag(vcovHC(model4, type = "HC3")))

# --- CR2 Clustered SEs at Region level (robustness check) ---
cluster_se1 <- sqrt(diag(vcovCR(model1, cluster = df_did$Region,      type = "CR2")))
cluster_se2 <- sqrt(diag(vcovCR(model2, cluster = df_did$Region,      type = "CR2")))
cluster_se3 <- sqrt(diag(vcovCR(model3, cluster = df_no_covid$Region, type = "CR2")))
cluster_se4 <- sqrt(diag(vcovCR(model4, cluster = df_no_covid$Region, type = "CR2")))

# ============================================================
# 9. Print Results
# ============================================================

cat("\n========== MODEL 1: DiD Levels (HC3) ==========\n")
print(coeftest(model1, vcov = vcovHC(model1, type = "HC3")))

cat("\n========== MODEL 2: DiD Log Levels (HC3) ==========\n")
print(coeftest(model2, vcov = vcovHC(model2, type = "HC3")))

cat("\n========== MODEL 3: DiD Levels No COVID (HC3) ==========\n")
print(coeftest(model3, vcov = vcovHC(model3, type = "HC3")))

cat("\n========== MODEL 4: DiD Log Levels No COVID (HC3) ==========\n")
print(coeftest(model4, vcov = vcovHC(model4, type = "HC3")))

# Stargazer table - HC3 SEs
stargazer(model1, model2, model3, model4,
          se               = list(se1_hc3, se2_hc3, se3_hc3, se4_hc3),
          type             = "text",
          title            = "DiD Estimates - Cycling Demand (HC3 Jackknife SE)",
          dep.var.labels   = c("Cycle Traffic (bn miles)", "Log Cycle Traffic",
                               "Cycle Traffic (no COVID)", "Log Cycle Traffic (no COVID)"),
          covariate.labels = c("London (Treated)", "Post-2019",
                               "DiD (Treated x Post)", "Constant"),
          omit.stat        = c("f", "ser"),
          digits           = 4,
          out              = "cycling_did_results_HC3.txt")

# Stargazer table - CR2 Clustered SEs (robustness)
stargazer(model1, model2, model3, model4,
          se               = list(cluster_se1, cluster_se2, cluster_se3, cluster_se4),
          type             = "text",
          title            = "DiD Estimates - Cycling Demand (CR2 Clustered SE by Region)",
          dep.var.labels   = c("Cycle Traffic (bn miles)", "Log Cycle Traffic",
                               "Cycle Traffic (no COVID)", "Log Cycle Traffic (no COVID)"),
          covariate.labels = c("London (Treated)", "Post-2019",
                               "DiD (Treated x Post)", "Constant"),
          omit.stat        = c("f", "ser"),
          digits           = 4,
          out              = "cycling_did_results_CR2.txt")

# Comparison table: HC3 vs CR2
stargazer(model1, model1,
          se               = list(se1_hc3, cluster_se1),
          type             = "text",
          title            = "Comparison of SE Methods - DiD Levels Model",
          column.labels    = c("HC3 (Jackknife)", "CR2 (Clustered by Region)"),
          add.lines        = list(
            c("SE Method",     "HC3 Jackknife",  "CR2 Clustered"),
            c("Cluster Level", "None",           "Region (G=3)")
          ),
          covariate.labels = c("London (Treated)", "Post-2019",
                               "DiD (Treated x Post)", "Constant"),
          omit.stat        = c("f", "ser"),
          digits           = 4,
          out              = "cycling_se_comparison.txt")

# ============================================================
# 10. DiD Coefficient Plot (HC3)
# ============================================================

get_did_ci <- function(model, label) {
  coef_rob <- coeftest(model, vcov = vcovHC(model, type = "HC3"))
  did_row  <- coef_rob["did", ]
  tibble(
    Model    = label,
    Estimate = did_row["Estimate"],
    SE       = did_row["Std. Error"],
    CI_lo    = did_row["Estimate"] - 1.96 * did_row["Std. Error"],
    CI_hi    = did_row["Estimate"] + 1.96 * did_row["Std. Error"]
  )
}

did_coefs <- bind_rows(
  get_did_ci(model1, "Levels"),
  get_did_ci(model2, "Log Levels"),
  get_did_ci(model3, "Levels (No COVID)"),
  get_did_ci(model4, "Log Levels (No COVID)")
)

coef_plot <- ggplot(did_coefs,
                    aes(x = Model, y = Estimate,
                        ymin = CI_lo, ymax = CI_hi)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_pointrange(colour = "#1f77b4", size = 0.9, linewidth = 1) +
  labs(title    = "DiD Coefficient - Effect of ULEZ on London Cycling Demand",
       subtitle = "95% CI with HC3 Jackknife Standard Errors",
       x = "Model Specification",
       y = "DiD Estimate (Treated x Post)") +
  theme_minimal(base_size = 13)

print(coef_plot)
ggsave("cycling_did_coefplot.png", coef_plot, width = 8, height = 5, dpi = 150)

# ============================================================
# 11. Event Study (Year-by-year treatment effects, HC3)
# ============================================================

df_event <- df_did %>%
  mutate(year_fct = relevel(factor(Year), ref = "2018"))

model_event <- lm(CycleTraffic_bn ~ treated:year_fct + year_fct + treated,
                  data = df_event)

# HC3 Jackknife SEs for event study
event_vcov  <- vcovHC(model_event, type = "HC3")
event_coefs <- coeftest(model_event, vcov = event_vcov)

# Extract coefficients safely by position
event_df_raw <- data.frame(
  term     = rownames(event_coefs),
  Estimate = event_coefs[, 1],
  SE       = event_coefs[, 2],
  tval     = event_coefs[, 3],
  pval     = event_coefs[, 4],
  stringsAsFactors = FALSE
)

event_df <- event_df_raw %>%
  filter(grepl("treated:year_fct", term)) %>%
  mutate(
    Year  = as.integer(str_extract(term, "\\d{4}")),
    CI_lo = Estimate - 1.96 * SE,
    CI_hi = Estimate + 1.96 * SE
  ) %>%
  bind_rows(data.frame(term = "base", Year = 2018L,
                       Estimate = 0, SE = 0,
                       tval = NA_real_, pval = NA_real_,
                       CI_lo = 0, CI_hi = 0)) %>%
  arrange(Year)

event_plot <- ggplot(event_df,
                     aes(x = Year, y = Estimate,
                         ymin = CI_lo, ymax = CI_hi)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 2018.5, linetype = "dotted", colour = "red") +
  geom_ribbon(alpha = 0.15, fill = "#1f77b4") +
  geom_line(colour = "#1f77b4", linewidth = 1) +
  geom_point(colour = "#1f77b4", size = 2.5) +
  annotate("text", x = 2019.2,
           y = max(event_df$CI_hi, na.rm = TRUE) * 0.9,
           label = "ULEZ ->", colour = "red", size = 3.5) +
  labs(title    = "Event Study - London Cycling Demand vs Control",
       subtitle = "Year-by-year DiD coefficients (base year = 2018), HC3 Jackknife SE",
       x = "Year",
       y = "Estimated Treatment Effect (bn vehicle miles)") +
  theme_minimal(base_size = 13)

print(event_plot)
ggsave("cycling_event_study.png", event_plot, width = 10, height = 5.5, dpi = 150)

# ============================================================
# 12. Placebo Tests
# ============================================================

# --- Placebo 1: Fake policy year = 2015 (pre-2019 data only) ---
df_placebo_time <- df_did %>%
  filter(Year <= 2018) %>%
  mutate(
    post_placebo = if_else(Year >= 2015, 1L, 0L),
    did_placebo  = treated * post_placebo
  )

model_placebo_time <- lm(CycleTraffic_bn ~ treated + post_placebo + did_placebo,
                         data = df_placebo_time)

cat("\n========== PLACEBO 1: Fake Policy Year 2015 ==========\n")
print(coeftest(model_placebo_time, vcov = vcovHC(model_placebo_time, type = "HC3")))

# --- Placebo 2: Fake treatment = North West (exclude London) ---
df_placebo_space <- df_did %>%
  filter(Region != "London") %>%
  mutate(
    treated_placebo = if_else(Region == "North West", 1L, 0L),
    did_placebo     = treated_placebo * post
  )

model_placebo_space <- lm(CycleTraffic_bn ~ treated_placebo + post + did_placebo,
                          data = df_placebo_space)

cat("\n========== PLACEBO 2: Fake Treatment Region (North West) ==========\n")
print(coeftest(model_placebo_space, vcov = vcovHC(model_placebo_space, type = "HC3")))

cat("\nNote: Both placebo DiD coefficients should be insignificant.\n")
cat("If significant, this undermines the credibility of the main result.\n")

# ============================================================
# Done
# ============================================================

cat("\nDone! Files saved to:", getwd(), "\n")
cat("  - cycling_parallel_trends.png\n")
cat("  - cycling_full_timeseries.png\n")
cat("  - cycling_did_coefplot.png\n")
cat("  - cycling_event_study.png\n")
cat("  - cycling_did_results_HC3.txt\n")
cat("  - cycling_did_results_CR2.txt\n")
cat("  - cycling_se_comparison.txt\n")
cat("  - cycling_summary_stats.txt\n")
