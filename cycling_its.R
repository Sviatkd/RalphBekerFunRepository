# ============================================================
# Cycling Demand - Interrupted Time Series (ITS) Analysis
# London only - pre-policy period as its own counterfactual
# No control group, no parallel trends assumption needed
# Formula: Cycling = a + time + ULEZdummy + ULEZ:time
# time      = pre-existing trend in London
# ulez      = immediate level change at policy (2019)
# ulez_time = change in slope after policy
# SE: HC3 Jackknife
# Run with and without COVID years (Option 3 - present both)
# ============================================================


# ============================================================
# 1. Libraries
# ============================================================

if (!require("tidyverse"))    install.packages("tidyverse")
if (!require("readODS"))      install.packages("readODS")
if (!require("lmtest"))       install.packages("lmtest")
if (!require("sandwich"))     install.packages("sandwich")
if (!require("stargazer"))    install.packages("stargazer")
if (!require("ggplot2"))      install.packages("ggplot2")
if (!require("car"))          install.packages("car")

library(tidyverse)
library(readODS)
library(lmtest)
library(sandwich)
library(stargazer)
library(ggplot2)
library(car)


# ============================================================
# 2. Load & Prepare Data (London only)
# ============================================================

setwd(".")  # Change to your folder path if needed

raw <- read_ods("tra0403-miles-pedal-cycle-traffic-by-region.ods",
                sheet = "TRA0403",
                skip  = 3)

colnames(raw)[1] <- "Region"
colnames(raw)[2] <- "Subarea"
colnames(raw)[3] <- "Notes"
colnames(raw)[4] <- "Units"
colnames(raw)[5:ncol(raw)] <- as.character(1993:2024)

df_london <- raw %>%
  filter(Region == "London",
         grepl("All", Subarea, ignore.case = TRUE)) %>%
  select(Region, all_of(as.character(1993:2024))) %>%
  mutate(across(as.character(1993:2024), as.numeric)) %>%
  pivot_longer(cols      = as.character(1993:2024),
               names_to  = "Year",
               values_to = "CycleTraffic_bn") %>%
  mutate(
    Year            = as.integer(Year),
    CycleTraffic_bn = as.numeric(CycleTraffic_bn),
    time            = Year - 2019,
    ulez            = if_else(Year >= 2019, 1L, 0L),
    ulez_time       = if_else(Year >= 2019, Year - 2019, 0L),
    log_cycle       = log(CycleTraffic_bn)
  )

df_no_covid <- df_london %>% filter(!Year %in% c(2020, 2021))


# ============================================================
# 3. ITS Models
# ============================================================

# Full sample
its1 <- lm(CycleTraffic_bn ~ time + ulez + ulez_time, data = df_london)
its2 <- lm(log_cycle       ~ time + ulez + ulez_time, data = df_london)

# COVID excluded
its3 <- lm(CycleTraffic_bn ~ time + ulez + ulez_time, data = df_no_covid)
its4 <- lm(log_cycle       ~ time + ulez + ulez_time, data = df_no_covid)

# HC3 SEs
se1 <- sqrt(diag(vcovHC(its1, type = "HC3")))
se2 <- sqrt(diag(vcovHC(its2, type = "HC3")))
se3 <- sqrt(diag(vcovHC(its3, type = "HC3")))
se4 <- sqrt(diag(vcovHC(its4, type = "HC3")))


# ============================================================
# 4. Print Results
# ============================================================

cat("\n========== ITS MODEL 1: Levels - Full Sample (HC3) ==========\n")
print(coeftest(its1, vcov = vcovHC(its1, type = "HC3")))

cat("\n========== ITS MODEL 2: Log Levels - Full Sample (HC3) ==========\n")
print(coeftest(its2, vcov = vcovHC(its2, type = "HC3")))

cat("\n========== ITS MODEL 3: Levels - No COVID (HC3) ==========\n")
print(coeftest(its3, vcov = vcovHC(its3, type = "HC3")))

cat("\n========== ITS MODEL 4: Log Levels - No COVID (HC3) ==========\n")
print(coeftest(its4, vcov = vcovHC(its4, type = "HC3")))

# Full sample vs No COVID side by side (levels)
stargazer(its1, its3,
          se            = list(se1, se3),
          type          = "text",
          title         = "ITS - London Cycling: Full Sample vs No COVID (Levels, HC3)",
          column.labels = c("Full Sample", "COVID Excluded"),
          add.lines     = list(c("COVID years", "Included", "Excluded (2020-2021)")),
          covariate.labels = c("Time Trend", "ULEZ Level Change",
                               "ULEZ Slope Change", "Constant"),
          omit.stat     = c("f", "ser"),
          digits        = 4,
          out           = "cycling_its_levels_comparison.txt")

# Full sample vs No COVID side by side (log)
stargazer(its2, its4,
          se            = list(se2, se4),
          type          = "text",
          title         = "ITS - London Cycling: Full Sample vs No COVID (Log, HC3)",
          column.labels = c("Full Sample", "COVID Excluded"),
          add.lines     = list(c("COVID years", "Included", "Excluded (2020-2021)")),
          covariate.labels = c("Time Trend", "ULEZ Level Change",
                               "ULEZ Slope Change", "Constant"),
          omit.stat     = c("f", "ser"),
          digits        = 4,
          out           = "cycling_its_log_comparison.txt")


# ============================================================
# 5. Wald Tests (on full sample model)
# ============================================================

cat("\n========== WALD TEST 1: Joint Effect (level + slope) ==========\n")
cat("H0: ulez = 0 AND ulez_time = 0\n")
cat("Rejection = ULEZ had some effect on cycling\n\n")
wald_joint <- linearHypothesis(its1,
                               c("ulez = 0", "ulez_time = 0"),
                               vcov = vcovHC(its1, type = "HC3"))
print(wald_joint)

cat("\n========== WALD TEST 2: Slope Change ==========\n")
cat("H0: ulez_time = 0 (no change in trend after ULEZ)\n\n")
wald_slope <- linearHypothesis(its1, "ulez_time = 0",
                               vcov = vcovHC(its1, type = "HC3"))
print(wald_slope)

cat("\n========== WALD TEST 3: Level Change ==========\n")
cat("H0: ulez = 0 (no immediate jump at ULEZ)\n\n")
wald_level <- linearHypothesis(its1, "ulez = 0",
                               vcov = vcovHC(its1, type = "HC3"))
print(wald_level)

cat("\n--- Wald Test Summary ---\n")
cat("Joint (level + slope) p: ", round(wald_joint$`Pr(>F)`[2], 4), "\n")
cat("Slope change p:           ", round(wald_slope$`Pr(>F)`[2], 4), "\n")
cat("Level change p:           ", round(wald_level$`Pr(>F)`[2], 4), "\n")


# ============================================================
# 6. ITS Plot - Full Sample
# Solid line = fitted model
# Dashed line = counterfactual (pre-policy trend extrapolated)
# ============================================================

df_london <- df_london %>%
  mutate(
    fitted      = predict(its1),
    counterfact = coef(its1)["(Intercept)"] + coef(its1)["time"] * time
  )

its_plot <- ggplot(df_london, aes(x = Year)) +
  geom_point(aes(y = CycleTraffic_bn), colour = "#1f77b4", size = 2.5) +
  geom_line(aes(y = fitted), colour = "#1f77b4", linewidth = 1.1) +
  geom_line(aes(y = counterfact), colour = "grey50",
            linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = 2019, linetype = "dotted",
             colour = "red", linewidth = 0.9) +
  annotate("text", x = 2019.3,
           y = max(df_london$CycleTraffic_bn, na.rm = TRUE) * 0.95,
           label = "ULEZ 2019", colour = "red", hjust = 0, size = 3.5) +
  annotate("text", x = 1995,
           y = min(df_london$CycleTraffic_bn, na.rm = TRUE) * 1.08,
           label = "Dashed = counterfactual (pre-policy trend)",
           colour = "grey40", hjust = 0, size = 3.2) +
  labs(title    = "ITS Model - London Cycling Demand (Full Sample)",
       subtitle = "Solid = fitted ITS | Dashed = pre-policy trend extrapolated",
       x = "Year", y = "Cycle Traffic (Billion Vehicle Miles)") +
  theme_minimal(base_size = 13)

print(its_plot)
ggsave("cycling_its_plot.png", its_plot, width = 10, height = 5.5, dpi = 150)


# ============================================================
# 7. ITS Plot - COVID Excluded
# ============================================================

df_no_covid <- df_no_covid %>%
  mutate(
    fitted      = predict(its3),
    counterfact = coef(its3)["(Intercept)"] + coef(its3)["time"] * time
  )

its_plot_nc <- ggplot(df_no_covid, aes(x = Year)) +
  geom_point(aes(y = CycleTraffic_bn), colour = "#d62728", size = 2.5) +
  geom_line(aes(y = fitted), colour = "#d62728", linewidth = 1.1) +
  geom_line(aes(y = counterfact), colour = "grey50",
            linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = 2019, linetype = "dotted",
             colour = "red", linewidth = 0.9) +
  annotate("text", x = 2019.3,
           y = max(df_no_covid$CycleTraffic_bn, na.rm = TRUE) * 0.95,
           label = "ULEZ 2019", colour = "red", hjust = 0, size = 3.5) +
  annotate("text", x = 1995,
           y = min(df_no_covid$CycleTraffic_bn, na.rm = TRUE) * 1.08,
           label = "Dashed = counterfactual (pre-policy trend)",
           colour = "grey40", hjust = 0, size = 3.2) +
  labs(title    = "ITS Model - London Cycling Demand (COVID Years Excluded)",
       subtitle = "Solid = fitted ITS | Dashed = pre-policy trend extrapolated",
       x = "Year", y = "Cycle Traffic (Billion Vehicle Miles)") +
  theme_minimal(base_size = 13)

print(its_plot_nc)
ggsave("cycling_its_plot_nocovid.png", its_plot_nc, width = 10, height = 5.5, dpi = 150)


# ============================================================
# 8. Combined Comparison Plot - Full Sample vs COVID Excluded
# Both versions side by side in one figure using facet_wrap
# ============================================================

# Build combined dataset with a label for each version
df_combined <- bind_rows(
  df_london %>%
    mutate(Version = "Full Sample"),
  df_no_covid %>%
    mutate(Version = "COVID Excluded (2020-2021 removed)")
)

# Counterfactual for each version
cf_full    <- coef(its1)["(Intercept)"] + coef(its1)["time"] * df_london$time
cf_nocovid <- coef(its3)["(Intercept)"] + coef(its3)["time"] * df_no_covid$time

df_combined <- df_combined %>%
  mutate(counterfact = c(cf_full, cf_nocovid))

combined_plot <- ggplot(df_combined, aes(x = Year)) +
  geom_point(aes(y = CycleTraffic_bn), colour = "#1f77b4", size = 2) +
  geom_line(aes(y = fitted), colour = "#1f77b4", linewidth = 1.1) +
  geom_line(aes(y = counterfact), colour = "grey50",
            linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = 2019, linetype = "dotted",
             colour = "red", linewidth = 0.9) +
  annotate("text", x = 2019.3, y = 0.66,
           label = "ULEZ 2019", colour = "red", hjust = 0, size = 3.2) +
  facet_wrap(~Version) +
  labs(title    = "ITS Model - London Cycling: Full Sample vs COVID Excluded",
       subtitle = "Solid = fitted ITS | Dashed = pre-policy counterfactual trend",
       x = "Year", y = "Cycle Traffic (Billion Vehicle Miles)") +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold", size = 11))

print(combined_plot)
ggsave("cycling_its_comparison.png", combined_plot,
       width = 13, height = 5.5, dpi = 150)


# ============================================================
# Done
# ============================================================

cat("\nDone! Files saved to:", getwd(), "\n")
cat("  - cycling_its_levels_comparison.txt\n")
cat("  - cycling_its_log_comparison.txt\n")
cat("  - cycling_its_plot.png\n")
cat("  - cycling_its_plot_nocovid.png\n")
cat("  - cycling_its_comparison.png\n")

