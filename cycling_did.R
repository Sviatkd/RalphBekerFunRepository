# DiD analysis - cycling demand
# London treated vs FMC controls (matching air quality side)
# 2019 ULEZ cutoff
# TWFE: Unit + Year FE via feols() — consistent with Sviat's air quality spec
# WCB via fwildclusterboot::boottest() on feols objects


if (!require("readODS"))          install.packages("readODS")
if (!require("tidyverse"))        install.packages("tidyverse")
if (!require("lmtest"))           install.packages("lmtest")
if (!require("sandwich"))         install.packages("sandwich")
if (!require("clubSandwich"))     install.packages("clubSandwich")
if (!require("stargazer"))        install.packages("stargazer")
if (!require("ggplot2"))          install.packages("ggplot2")
if (!require("skimr"))            install.packages("skimr")
if (!require("car"))              install.packages("car")
if (!require("fixest"))           install.packages("fixest")

library(readODS)
library(tidyverse)
library(lmtest)
library(sandwich)
library(clubSandwich)
library(stargazer)
library(ggplot2)
library(skimr)
library(car)
library(fixest)


# load data

setwd("C:/Users/kalid/OneDrive/Desktop/New folder")

raw <- read_ods("tra0403-miles-pedal-cycle-traffic-by-region.ods",
                sheet = "TRA0403",
                skip  = 3)


# clean & reshape

colnames(raw)[1] <- "Region"
colnames(raw)[2] <- "Subarea"
colnames(raw)[3] <- "Notes"
colnames(raw)[4] <- "Units"
colnames(raw)[5:ncol(raw)] <- as.character(1993:2024)

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
    (Region == "London" & grepl("All", Subarea, ignore.case = TRUE)) |
    trimws(Subarea) %in% controls_keep
  ) %>%
  select(Region, Subarea, all_of(as.character(1993:2024))) %>%
  mutate(across(as.character(1993:2024), as.numeric)) %>%
  mutate(
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

df_long <- df_wide %>%
  pivot_longer(cols      = as.character(1993:2024),
               names_to  = "Year",
               values_to = "CycleTraffic_bn") %>%
  mutate(Year            = as.integer(Year),
         CycleTraffic_bn = as.numeric(CycleTraffic_bn))

df_did <- df_long %>%
  mutate(
    treated   = if_else(Unit == "London", 1L, 0L),
    post      = if_else(Year >= 2019, 1L, 0L),
    did       = treated * post,
    log_cycle = log(CycleTraffic_bn)
  )


# summary stats

cat("by unit:\n")
df_did %>%
  group_by(Unit) %>%
  summarise(
    Mean = round(mean(CycleTraffic_bn, na.rm = TRUE), 3),
    SD   = round(sd(CycleTraffic_bn,   na.rm = TRUE), 3),
    Min  = round(min(CycleTraffic_bn,  na.rm = TRUE), 3),
    Max  = round(max(CycleTraffic_bn,  na.rm = TRUE), 3),
    N    = n()
  ) %>% print()

cat("\npre vs post:\n")
df_did %>%
  group_by(Unit, post) %>%
  summarise(
    Mean = round(mean(CycleTraffic_bn, na.rm = TRUE), 3),
    SD   = round(sd(CycleTraffic_bn,   na.rm = TRUE), 3),
    N    = n(), .groups = "drop"
  ) %>%
  mutate(Period = if_else(post == 1, "Post-2019", "Pre-2019")) %>%
  select(Unit, Period, Mean, SD, N) %>% print()

cat("\nulez summary:\n")
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


# parallel trends plot

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


# full time series

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


# TWFE models
# Unit + Year FE — Unit FE absorbs treated, Year FE absorbs post and COVID shocks
# Full sample, no need to drop 2020/2021

twfe1 <- feols(CycleTraffic_bn ~ did | Unit + Year, data = df_did)
twfe2 <- feols(log_cycle       ~ did | Unit + Year, data = df_did)

cat("\nTWFE levels (HC1):\n")
print(summary(twfe1, vcov = "HC1"))

cat("\nTWFE log (HC1):\n")
print(summary(twfe2, vcov = "HC1"))

cat("\nTWFE levels (unit-clustered):\n")
print(summary(twfe1, vcov = ~Unit))

cat("\nTWFE log (unit-clustered):\n")
print(summary(twfe2, vcov = ~Unit))

etable(twfe1, twfe2,
       vcov   = "HC1",
       title  = "TWFE DiD - Cycling Demand (Unit + Year FE, HC1 SE)",
       digits = 4,
       file   = "cycling_did_TWFE.txt")

etable(twfe1, twfe2,
       vcov   = ~Unit,
       title  = "TWFE DiD - Cycling Demand (Unit + Year FE, Unit-Clustered SE)",
       digits = 4,
       file   = "cycling_did_TWFE_clustered.txt")


# coefficient plot

twfe_coefs <- bind_rows(
  tibble(Model = "Levels",     Estimate = coef(twfe1)["did"],
         CI_lo = confint(twfe1, vcov = "HC1")["did", 1],
         CI_hi = confint(twfe1, vcov = "HC1")["did", 2]),
  tibble(Model = "Log Levels", Estimate = coef(twfe2)["did"],
         CI_lo = confint(twfe2, vcov = "HC1")["did", 1],
         CI_hi = confint(twfe2, vcov = "HC1")["did", 2])
)

coef_plot <- ggplot(twfe_coefs,
                    aes(x = Model, y = Estimate,
                        ymin = CI_lo, ymax = CI_hi)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_pointrange(colour = "#1f77b4", size = 0.9, linewidth = 1) +
  labs(title    = "TWFE DiD Coefficient - ULEZ Effect on London Cycling",
       subtitle = "95% CI, HC1 SE, Unit + Year Fixed Effects",
       x = "Model Specification",
       y = "DiD Estimate (Treated x Post)") +
  theme_minimal(base_size = 13)

print(coef_plot)
ggsave("cycling_did_coefplot_FMC.png", coef_plot, width = 8, height = 5, dpi = 150)


# event study

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


# placebo tests — TWFE feols

df_placebo_time <- df_did %>%
  filter(Year <= 2018) %>%
  mutate(
    post_placebo = if_else(Year >= 2015, 1L, 0L),
    did_placebo  = treated * post_placebo
  )

twfe_placebo_time <- feols(CycleTraffic_bn ~ did_placebo | Unit + Year,
                           data = df_placebo_time)

cat("\nplacebo 1 (fake 2015):\n")
print(summary(twfe_placebo_time, vcov = ~Unit))

df_placebo_space <- df_did %>%
  filter(Unit != "London") %>%
  mutate(
    treated_placebo = if_else(Unit == "Greater Manchester", 1L, 0L),
    did_placebo     = treated_placebo * post
  )

twfe_placebo_space <- feols(CycleTraffic_bn ~ did_placebo | Unit + Year,
                            data = df_placebo_space)

cat("\nplacebo 2 (fake GM):\n")
print(summary(twfe_placebo_space, vcov = ~Unit))


# pre-trends wald + manual WCB
# HC3 for pre-trends, webb weights on lm with unit+year dummies for main tests

df_pre    <- df_did %>% filter(Year < 2019)
model_pre <- lm(CycleTraffic_bn ~ treated + Year + treated:Year,
                data = df_pre)

cat("\npre-trends test (HC3):\n")
wald_pre <- linearHypothesis(model_pre, "treated:Year = 0",
                             vcov = vcovHC(model_pre, type = "HC3"))
print(wald_pre)

set.seed(42)

wcb_pval <- function(model, param, cluster_var, data, B = 9999, seed = 42) {
  set.seed(seed)
  webb      <- c(-sqrt(3/2), -1, -sqrt(1/2), sqrt(1/2), 1, sqrt(3/2))
  cr2_vcov  <- vcovCR(model, cluster = data[[cluster_var]], type = "CR2")
  t_obs     <- coef(model)[param] / sqrt(cr2_vcov[param, param])
  coefs_r        <- coef(model)
  coefs_r[param] <- 0
  y_hat_r   <- model.matrix(model) %*% coefs_r
  resids_r  <- model$model[, 1] - y_hat_r
  clusters  <- unique(data[[cluster_var]])
  G         <- length(clusters)
  t_boot    <- numeric(B)
  for (b in seq_len(B)) {
    w         <- sample(webb, G, replace = TRUE)
    names(w)  <- clusters
    w_vec     <- w[data[[cluster_var]]]
    y_boot    <- as.numeric(y_hat_r + resids_r * w_vec)
    df_b          <- data
    df_b[["y_b"]] <- y_boot
    mod_b     <- lm(update(formula(model), y_b ~ .), data = df_b)
    vcov_b    <- vcovCR(mod_b, cluster = df_b[[cluster_var]], type = "CR2")
    t_boot[b] <- coef(mod_b)[param] / sqrt(vcov_b[param, param])
  }
  list(t_obs = t_obs, p_val = mean(abs(t_boot) >= abs(t_obs)))
}

lm_main <- lm(CycleTraffic_bn ~ did + factor(Unit) + factor(Year), data = df_did)
lm_p1   <- lm(CycleTraffic_bn ~ did_placebo + factor(Unit) + factor(Year),
               data = df_placebo_time)
lm_p2   <- lm(CycleTraffic_bn ~ did_placebo + factor(Unit) + factor(Year),
               data = df_placebo_space)

cat("\nWCB main effect:\n")
wcb_main <- wcb_pval(lm_main, "did",         "Unit", df_did)
cat("p =", round(wcb_main$p_val, 4), "\n")

cat("\nplacebo 2015 WCB:\n")
wcb_p1 <- wcb_pval(lm_p1, "did_placebo", "Unit", df_placebo_time)
cat("p =", round(wcb_p1$p_val, 4), "\n")

cat("\nplacebo GM WCB:\n")
wcb_p2 <- wcb_pval(lm_p2, "did_placebo", "Unit", df_placebo_space)
cat("p =", round(wcb_p2$p_val, 4), "\n")

cat("\nsummary:\n")
cat("Pre-trends HC3 p-value:   ", round(wald_pre$`Pr(>F)`[2], 4), "\n")
cat("Post-policy WCB p-value:  ", round(wcb_main$p_val, 4), "\n")
cat("Placebo 1 (fake 2015) p:  ", round(wcb_p1$p_val,   4), "\n")
cat("Placebo 2 (fake GM) p:    ", round(wcb_p2$p_val,   4), "\n")


# acf plots — residuals from TWFE model

df_did <- df_did %>%
  mutate(resid_did = residuals(twfe1))

units_ordered <- sort(unique(df_did$Unit))

png("cycling_acf_by_unit.png", width = 1800, height = 1200, res = 150)
par(mfrow = c(2, 4), mar = c(4, 4, 3, 1))

for (u in units_ordered) {
  resids <- df_did %>%
    filter(Unit == u) %>%
    arrange(Year) %>%
    pull(resid_did)

  acf(resids,
      main   = u,
      xlab   = "Lag (Years)",
      ylab   = "Autocorrelation",
      col    = "#1f77b4",
      lwd    = 2,
      ci.col = "red")
}

dev.off()

df_city_resid <- df_did %>%
  group_by(Unit, Year) %>%
  summarise(mean_resid = mean(resid_did), .groups = "drop")

png("cycling_acf_city_level.png", width = 1800, height = 900, res = 150)
par(mfrow = c(2, 4), mar = c(4, 4, 3, 1))

for (u in units_ordered) {
  resids <- df_city_resid %>%
    filter(Unit == u) %>%
    arrange(Year) %>%
    pull(mean_resid)

  acf(resids,
      main   = paste(u, "(City Level)"),
      xlab   = "Lag (Years)",
      ylab   = "Autocorrelation",
      col    = "#d62728",
      lwd    = 2,
      ci.col = "darkred")
}

dev.off()

cat("\ndone\n")
