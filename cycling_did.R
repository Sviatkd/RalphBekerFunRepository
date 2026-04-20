# ============================================================
# Cycling Demand - Difference-in-Differences (DiD) Analysis
# TRA0403: Pedal cycle traffic (billion vehicle miles)
# Treatment: London
# Controls: Greater Manchester FMC, West Midlands FMC,
#           Tyne and Wear FMC, West Yorkshire FMC,
#           South Yorkshire FMC, Scotland
# These match the air quality control cities exactly
# Policy cutoff: 2019 (ULEZ)
# SE: HC3 Jackknife (primary) + CR2 Region-clustered (robustness)
# ============================================================


# ============================================================
# 1. Libraries
# ============================================================

if (!require("readODS"))      install.packages("readODS")
if (!require("tidyverse"))    install.packages("tidyverse")
if (!require("lmtest"))       install.packages("lmtest")
if (!require("sandwich"))     install.packages("sandwich")
if (!require("clubSandwich")) install.packages("clubSandwich")
if (!require("stargazer"))    install.packages("stargazer")
if (!require("ggplot2"))      install.packages("ggplot2")
if (!require("skimr"))        install.packages("skimr")
if (!require("aod"))          install.packages("aod")
if (!require("car"))          install.packages("car")

library(readODS)
library(tidyverse)
library(lmtest)
library(sandwich)
library(clubSandwich)
library(stargazer)
library(ggplot2)
library(skimr)
library(aod)
library(car)


# ============================================================
# 2. Load Data
# ============================================================

setwd(".")  # Change to your folder path if needed

raw <- read_ods("tra0403-miles-pedal-cycle-traffic-by-region.ods",
                sheet = "TRA0403",
                skip  = 3)


# ============================================================
# 3. Clean & Reshape
# ============================================================

colnames(raw)[1] <- "Region"
colnames(raw)[2] <- "Subarea"
colnames(raw)[3] <- "Notes"
colnames(raw)[4] <- "Units"
colnames(raw)[5:ncol(raw)] <- as.character(1993:2024)

# Keep London (treated) + 5 FMC controls + Scotland
# FMCs match air quality control cities exactly
controls_keep <- c(
  "Greater Manchester FMC",
  "West Midlands FMC",
  "Tyne and Wear FMC",
  "West Yorkshire FMC",
  "South Yorkshire FMC",
  "Scotland All"
)

df_wide <- raw %>%
  filter(
    # London treated
    (Region == "London" & grepl("All", Subarea, ignore.case = TRUE)) |
    # FMC controls
    trimws(Subarea) %in% controls_keep
  ) %>%
  select(Region, Subarea, all_of(as.character(1993:2024))) %>%
  mutate(across(as.character(1993:2024), as.numeric)) %>%
  mutate(
    # Clean label for each unit
    Unit = case_when(
      Region == "London"                                   ~ "London",
      trimws(Subarea) == "Greater Manchester FMC"          ~ "Greater Manchester",
      trimws(Subarea) == "West Midlands FMC"               ~ "West Midlands",
      trimws(Subarea) == "Tyne and Wear FMC"               ~ "Tyne and Wear",
      trimws(Subarea) == "West Yorkshire FMC"              ~ "West Yorkshire",
      trimws(Subarea) == "South Yorkshire FMC"             ~ "South Yorkshire",
      trimws(Subarea) == "Scotland All"                    ~ "Scotland",
      TRUE ~ Subarea
    )
  ) %>%
  select(Unit, all_of(as.character(1993:2024)))

# Pivot to long format
df_long <- df_wide %>%
  pivot_longer(cols      = as.character(1993:2024),
               names_to  = "Year",
               values_to = "CycleTraffic_bn") %>%
  mutate(Year            = as.integer(Year),
         CycleTraffic_bn = as.numeric(CycleTraffic_bn))


# ============================================================
# 4. DiD Variables
# ============================================================

df_did <- df_long %>%
  mutate(
    treated   = if_else(Unit == "London", 1L, 0L),
    post      = if_else(Year >= 2019, 1L, 0L),
    did       = treated * post,
    log_cycle = log(CycleTraffic_bn)
  )


# ============================================================
# 5. Summary Statistics
# ============================================================

cat("\n========== BY UNIT ==========\n")
df_did %>%
  group_by(Unit) %>%
  summarise(
    Mean = round(mean(CycleTraffic_bn, na.rm = TRUE), 3),
    SD   = round(sd(CycleTraffic_bn,   na.rm = TRUE), 3),
    Min  = round(min(CycleTraffic_bn,  na.rm = TRUE), 3),
    Max  = round(max(CycleTraffic_bn,  na.rm = TRUE), 3),
    N    = n()
  ) %>% print()

cat("\n========== PRE vs POST 2019 ==========\n")
df_did %>%
  group_by(Unit, post) %>%
  summarise(
    Mean = round(mean(CycleTraffic_bn, na.rm = TRUE), 3),
    SD   = round(sd(CycleTraffic_bn,   na.rm = TRUE), 3),
    N    = n(), .groups = "drop"
  ) %>%
  mutate(Period = if_else(post == 1, "Post-2019", "Pre-2019")) %>%
  select(Unit, Period, Mean, SD, N) %>% print()

cat("\n========== ULEZ STYLE SUMMARY ==========\n")
cycling_stats <- df_did %>%
  mutate(ULEZ = if_else(post == 1, "Post-ULEZ", "Pre-ULEZ")) %>%
  group_by(Unit, ULEZ) %>%
  summarise(
    Mean   = round(mean(CycleTraffic_bn,   na.rm = TRUE), 3),
    Median = round(median(CycleTraffic_bn, na.rm = TRUE), 3),
    SD     = round(sd(CycleTraffic_bn,     na.rm = TRUE), 3),
    IQR    = round(IQR(CycleTraffic_bn,    na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  rename(Mode = Unit) %>%
  arrange(Mode, desc(ULEZ))
print(cycling_stats)


# ============================================================
# 6. Parallel Trends Plot (pre-2019)
# ============================================================

pre_trends <- df_did %>%
  filter(Year <= 2019) %>%
  ggplot(aes(x = Year, y = CycleTraffic_bn,
             colour = Unit, group = Unit)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2) +
  geom_vline(xintercept = 2019, linetype = "dashed", colour = "black") +
  annotate("text", x = 2019.3,
           y = max(df_did$CycleTraffic_bn, na.rm = TRUE) * 0.95,
           label = "ULEZ 2019", hjust = 0, size = 3.2) +
  labs(title    = "Parallel Trends Check - Cycle Traffic (Pre-2019)",
       subtitle = "London vs FMC control cities, billion vehicle miles",
       x = "Year", y = "Cycle Traffic (Billion Vehicle Miles)",
       colour = "Unit") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

print(pre_trends)
ggsave("cycling_parallel_trends_FMC.png", pre_trends,
       width = 10, height = 5.5, dpi = 150)


# ============================================================
# 7. Full Time Series Plot
# ============================================================

full_plot <- df_did %>%
  ggplot(aes(x = Year, y = CycleTraffic_bn,
             colour = Unit, group = Unit)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2) +
  geom_vline(xintercept = 2019, linetype = "dashed", colour = "black") +
  annotate("text", x = 2019.3,
           y = max(df_did$CycleTraffic_bn, na.rm = TRUE) * 0.95,
           label = "ULEZ 2019", hjust = 0, size = 3.2) +
  labs(title    = "Cycle Traffic Miles - London vs FMC Controls (1993-2024)",
       subtitle = "Billion vehicle miles, annual",
       x = "Year", y = "Cycle Traffic (Billion Vehicle Miles)",
       colour = "Unit") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

print(full_plot)
ggsave("cycling_full_timeseries_FMC.png", full_plot,
       width = 10, height = 5.5, dpi = 150)


# ============================================================
# 8. DiD Models
# HC3 Jackknife (primary) + CR2 clustered at Unit level (robustness)
# Note: 7 units now vs 3 before - clustering more reliable
# ============================================================

model1 <- lm(CycleTraffic_bn ~ treated + post + did, data = df_did)
model2 <- lm(log_cycle       ~ treated + post + did, data = df_did)

df_no_covid <- df_did %>% filter(!Year %in% c(2020, 2021))
model3 <- lm(CycleTraffic_bn ~ treated + post + did, data = df_no_covid)
model4 <- lm(log_cycle       ~ treated + post + did, data = df_no_covid)

se1_hc3 <- sqrt(diag(vcovHC(model1, type = "HC3")))
se2_hc3 <- sqrt(diag(vcovHC(model2, type = "HC3")))
se3_hc3 <- sqrt(diag(vcovHC(model3, type = "HC3")))
se4_hc3 <- sqrt(diag(vcovHC(model4, type = "HC3")))

cluster_se1 <- sqrt(diag(vcovCR(model1, cluster = df_did$Unit,      type = "CR2")))
cluster_se2 <- sqrt(diag(vcovCR(model2, cluster = df_did$Unit,      type = "CR2")))
cluster_se3 <- sqrt(diag(vcovCR(model3, cluster = df_no_covid$Unit, type = "CR2")))
cluster_se4 <- sqrt(diag(vcovCR(model4, cluster = df_no_covid$Unit, type = "CR2")))


# ============================================================
# 9. Print Results
# ============================================================

cat("\n========== MODEL 1: Levels (HC3) ==========\n")
print(coeftest(model1, vcov = vcovHC(model1, type = "HC3")))

cat("\n========== MODEL 2: Log Levels (HC3) ==========\n")
print(coeftest(model2, vcov = vcovHC(model2, type = "HC3")))

cat("\n========== MODEL 3: Levels No COVID (HC3) ==========\n")
print(coeftest(model3, vcov = vcovHC(model3, type = "HC3")))

cat("\n========== MODEL 4: Log Levels No COVID (HC3) ==========\n")
print(coeftest(model4, vcov = vcovHC(model4, type = "HC3")))

stargazer(model1, model2, model3, model4,
          se               = list(se1_hc3, se2_hc3, se3_hc3, se4_hc3),
          type             = "text",
          title            = "DiD - Cycling Demand with FMC Controls (HC3 SE)",
          dep.var.labels   = c("Cycle Traffic (bn miles)", "Log Cycle Traffic",
                               "Cycle Traffic (no COVID)", "Log Cycle Traffic (no COVID)"),
          covariate.labels = c("London (Treated)", "Post-2019",
                               "DiD (Treated x Post)", "Constant"),
          omit.stat        = c("f", "ser"),
          digits           = 4,
          out              = "cycling_did_FMC_HC3.txt")

stargazer(model1, model2, model3, model4,
          se               = list(cluster_se1, cluster_se2, cluster_se3, cluster_se4),
          type             = "text",
          title            = "DiD - Cycling Demand with FMC Controls (CR2 Clustered SE)",
          dep.var.labels   = c("Cycle Traffic (bn miles)", "Log Cycle Traffic",
                               "Cycle Traffic (no COVID)", "Log Cycle Traffic (no COVID)"),
          covariate.labels = c("London (Treated)", "Post-2019",
                               "DiD (Treated x Post)", "Constant"),
          omit.stat        = c("f", "ser"),
          digits           = 4,
          out              = "cycling_did_FMC_CR2.txt")


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
  labs(title    = "DiD Coefficient - ULEZ Effect on London Cycling (FMC Controls)",
       subtitle = "95% CI with HC3 Jackknife Standard Errors",
       x = "Model Specification",
       y = "DiD Estimate (Treated x Post)") +
  theme_minimal(base_size = 13)

print(coef_plot)
ggsave("cycling_did_coefplot_FMC.png", coef_plot, width = 8, height = 5, dpi = 150)


# ============================================================
# 11. Event Study
# ============================================================

df_event <- df_did %>%
  mutate(year_fct = relevel(factor(Year), ref = "2018"))

model_event <- lm(CycleTraffic_bn ~ treated:year_fct + year_fct + treated,
                  data = df_event)

event_vcov  <- vcovHC(model_event, type = "HC3")
event_coefs <- coeftest(model_event, vcov = event_vcov)

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
  labs(title    = "Event Study - London Cycling vs FMC Controls",
       subtitle = "Year-by-year DiD coefficients (base year = 2018), HC3 SE",
       x = "Year",
       y = "Estimated Treatment Effect (bn vehicle miles)") +
  theme_minimal(base_size = 13)

print(event_plot)
ggsave("cycling_event_study_FMC.png", event_plot, width = 10, height = 5.5, dpi = 150)


# ============================================================
# 12. Placebo Tests
# ============================================================

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

df_placebo_space <- df_did %>%
  filter(Unit != "London") %>%
  mutate(
    treated_placebo = if_else(Unit == "Greater Manchester", 1L, 0L),
    did_placebo     = treated_placebo * post
  )

model_placebo_space <- lm(CycleTraffic_bn ~ treated_placebo + post + did_placebo,
                          data = df_placebo_space)

cat("\n========== PLACEBO 2: Fake Treatment - Greater Manchester ==========\n")
print(coeftest(model_placebo_space, vcov = vcovHC(model_placebo_space, type = "HC3")))

p_placebo1 <- coeftest(model_placebo_time,
                       vcov = vcovHC(model_placebo_time, type = "HC3"))["did_placebo", 4]
p_placebo2 <- coeftest(model_placebo_space,
                       vcov = vcovHC(model_placebo_space, type = "HC3"))["did_placebo", 4]


# ============================================================
# 13. Wald Tests
# ============================================================

df_pre    <- df_did %>% filter(Year < 2019)
model_pre <- lm(CycleTraffic_bn ~ treated + Year + treated:Year,
                data = df_pre)

cat("\n========== WALD TEST 1: Pre-Trends ==========\n")
cat("H0: treated:Year = 0 (London trend same as FMC controls pre-2019)\n\n")
wald_pre <- linearHypothesis(model_pre, "treated:Year = 0",
                             vcov = vcovHC(model_pre, type = "HC3"))
print(wald_pre)

cat("\n========== WALD TEST 2: Post-Policy Effect ==========\n")
cat("H0: DiD = 0\n\n")
wald_post <- linearHypothesis(model1, "did = 0",
                              vcov = vcovHC(model1, type = "HC3"))
print(wald_post)

cat("\n--- Summary ---\n")
cat("Pre-trends Wald p-value:  ", round(wald_pre$`Pr(>F)`[2],  4), "\n")
cat("Post-policy Wald p-value: ", round(wald_post$`Pr(>F)`[2], 4), "\n")
cat("Placebo 1 (fake 2015) p:  ", round(p_placebo1, 4), "\n")
cat("Placebo 2 (fake GM) p:    ", round(p_placebo2, 4), "\n")


# ============================================================
# Done
# ============================================================

cat("\nDone! Files saved to:", getwd(), "\n")
cat("  - cycling_parallel_trends_FMC.png\n")
cat("  - cycling_full_timeseries_FMC.png\n")
cat("  - cycling_did_coefplot_FMC.png\n")
cat("  - cycling_event_study_FMC.png\n")
cat("  - cycling_did_FMC_HC3.txt\n")
cat("  - cycling_did_FMC_CR2.txt\n")
