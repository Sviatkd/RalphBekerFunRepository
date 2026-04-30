library(readxl)
library(readxl)
library(tidyverse)
library(stargazer)
library(lmtest)
library(sandwich)
library(coefplot)
library(clubSandwich)
library(aod)
library(modelsummary)
library(gt)
#Data Cleaning#
Rawdata <- read_excel("~/Rawdata.xlsx", sheet = "Sheet1")
View(Rawdata)
no2quarterly <- Rawdata %>%
  # Returns numeric year.quarter (e.g., 2024.1)
  mutate(yr_qtr = quarter(date, with_year = TRUE)) %>% 
  group_by(yr_qtr,Station) %>%
  summarise(mean_value = mean(pm10, na.rm = TRUE))
no2quarterly <- no2quarterly %>%
  mutate(
    pct_changepol = (mean_value - lag(mean_value)) / lag(mean_value) * 100
  )
# --- 1. SET PATH AND LIST FILES ---
# This assumes your CSV files are in your current working directory. 
# If they are in a subfolder, change "." to "C:/your/path"
file_list <- list.files(pattern = ".*-air-quality.*\\.csv$", full.names = TRUE)

# --- 2. IMPORT AND MERGE ALL CSVs ---
# We read each file, add a column for the Station name derived from the filename,
# and stack them all into one big data frame.
AllStationsData <- file_list %>%
  set_names() %>% 
  map_df(~read_csv(.x, show_col_types = FALSE), .id = "FileName") %>%
  mutate(
    # Clean the filename to get a nice Station name (e.g., "glasgow-high street")
    Station = str_remove(basename(FileName), "-air-quality\\.csv"),
    # Ensure date is in a proper format (adjust format string if your CSVs differ)
    date = as.POSIXct(date) 
  ) %>%
  select(-FileName)

# --- 3. COMBINE WITH ORIGINAL EXCEL DATA ---
FullDataset <- bind_rows(Rawdata, AllStationsData)
#Data wrangling#
# --- 4. CALCULATE QUARTERLY MEANS AND CHANGE ---
# CRITICAL: We group by Station BEFORE lagging to avoid mixing data between cities
no2quarterly <- FullDataset %>%
  mutate(yr_qtr = quarter(date, with_year = TRUE)) %>% 
  group_by(yr_qtr, Station) %>%
  summarise(mean_value = mean(pm10, na.rm = TRUE), .groups = "drop") %>%
  arrange(Station, yr_qtr) %>%
  group_by(Station) %>%
  mutate(
    pct_changepol = (mean_value - lag(mean_value)) / lag(mean_value) * 100
  ) %>%
  ungroup()
no2quarterly <- no2quarterly %>%
  mutate(City = case_when(
    str_detect(Station, regex("Birmingham", ignore_case = TRUE)) ~ "Birmingham",
    str_detect(Station, regex("Glasgow", ignore_case = TRUE))    ~ "Glasgow",
    str_detect(Station, regex("Leeds", ignore_case = TRUE))      ~ "Leeds",
    str_detect(Station, regex("London", ignore_case = TRUE))     ~ "London",
    str_detect(Station, regex("Manchester", ignore_case = TRUE)) ~ "Manchester",
    str_detect(Station, regex("Newcastle", ignore_case = TRUE))  ~ "Newcastle",
    str_detect(Station, regex("Sheffield", ignore_case = TRUE))  ~ "Sheffield",
    str_detect(Station, regex("Camden", ignore_case = TRUE))     ~ "London",
    str_detect(Station, regex("Salford", ignore_case = TRUE))    ~ "Salford",
    str_detect(Station, regex("Sunderland", ignore_case = TRUE)) ~ "Sunderland",
    str_detect(Station, regex("Westminster", ignore_case = TRUE)) ~ "London",
    TRUE ~ "Other" # Catch-all for anything missed
  ))
ggplot(no2quarterly,mapping = aes(x=yr_qtr,y=mean_value,colour = factor(Station)))+stat_summary(fun=mean,geom="line")+geom_vline(xintercept = 2019.2,linetype="dashed")+ ggtitle("Means with a policy breakpoint(No2)")
#Defining groups and Intervention Period#
intervention <- c(2019.2,2019.3,2019.4,2020.1,2020.2,2020.3,2020.4,2021.1,2021.2,2021.3,2021.4,2022.1,2022.2,2022.3,2022.4)
treatment <- c(
  "camden-- euston road", 
  "london-bloomsbury", 
  "LondonC", 
  "LondonMR"
)

control <- c(
  "glasgow-anderston",
  "glasgow-byres road",
  "glasgow-high street",
  "glasgow-townhead",
  "GlasgowK",
  "manchester-oxford road",
  "ManchesterPic",
  "NewcastleC",
  "salford-eccles",
  "salford-m60"
)
panel_dfind <- no2quarterly %>%
  mutate(
    intervention = ifelse(yr_qtr %in% intervention, 1,0),
    control = ifelse(Station %in% control, 1, 0),
    treatment= ifelse(Station %in% treatment,1,0),
    did= intervention*treatment)
treatment_control_df <- panel_dfind %>%
  filter(treatment + control == 1,yr_qtr < 2023.1 & yr_qtr > 2015.1)
#Creating The Station specific counterfactual trend plots"
# 1. Calculate the average "Control Trend" across all control stations
control_trend <- treatment_control_df %>%
  filter(treatment == 0) %>%
  group_by(yr_qtr) %>%
  summarise(control_mean = mean(mean_value, na.rm = TRUE)) %>%
  mutate(control_drift = control_mean - first(control_mean)) # Movement relative to start


# 2. Calculate the pre-intervention baseline for each Treated Station
treated_baselines <- treatment_control_df %>%
  filter(treatment == 1, yr_qtr < 2019.2) %>%
  group_by(Station) %>%
  summarise(baseline_val = mean(mean_value, na.rm = TRUE))
# 3. Create the Counterfactual Data
counterfactual_data <- treated_baselines %>%
  crossing(control_trend) %>% # Merge every station with every time point
  mutate(counterfactual_val = baseline_val + control_drift)

# 4. Plotting
ggplot() +
  # Actual Treated Data (Solid Lines)
  geom_line(data = filter(treatment_control_df, treatment == 1), 
            aes(x = yr_qtr, y = mean_value, color = Station), 
            size = 0.8, alpha = 0.6) +
  # Pseudo-Counterfactual Data (Dashed Lines)
  geom_line(data = counterfactual_data, 
            aes(x = yr_qtr, y = counterfactual_val, color = Station), 
            linetype = "dashed", size = 0.8) +
  # Policy Intervention Line
  geom_vline(xintercept = 2019.2, linetype = "dotted", color = "black", size = 1) +
  labs(
    title = "Treated Station Trajectories vs. Pseudo-Counterfactuals",
    subtitle = "Solid = Actual | Dashed = Control Trend applied to Station Baseline",
    x = "Year.Quarter",
    y = "NO2 Level",
    caption = "Counterfactual based on mean movement of all control stations"
  ) +
  theme_minimal()
ggplot() +
  # Actual Treated Data (Solid Lines)
  geom_line(data = filter(treatment_control_df, treatment == 1), 
            aes(x = yr_qtr, y = mean_value, color = Station), 
            size = 0.8, alpha = 0.6) +
  # Pseudo-Counterfactual Data (Dashed Lines)
  geom_line(data = counterfactual_data, 
            aes(x = yr_qtr, y = counterfactual_val, color = Station), 
            linetype = "dashed", size = 0.8) +
  # Policy Intervention Line
  geom_vline(xintercept = 2019.2, linetype = "dotted", color = "black", size = 1) +
  labs(
    title = "Treated Station Trajectories vs. Pseudo-Counterfactuals",
    subtitle = "Solid = Actual | Dashed = Control Trend applied to Station Baseline",
    x = "Year.Quarter",
    y = "NO2 Level",
    caption = "Counterfactual based on mean movement of all control stations"
  ) +
  theme_minimal()
#Creating a Smooth Version#
library(broom)
# 1. Fit a LOESS model to the Control Group to capture the "Natural" trend
control_loess_mod <- loess(mean_value ~ yr_qtr, 
                           data = filter(treatment_control_df, treatment == 0), 
                           span = 0.75)

# 2. Create a prediction grid for the Control Trend
qtr_grid <- seq(min(treatment_control_df$yr_qtr), max(treatment_control_df$yr_qtr), by = 0.1)
control_trend_shape <- data.frame(yr_qtr = qtr_grid) %>%
  mutate(pred_val = predict(control_loess_mod, newdata = .)) %>%
  mutate(drift = pred_val - first(pred_val)) # Calculate movement relative to T=0

# 3. Get the pre-intervention mean for each Treated Station to use as the "Shift"
treated_starts <- treatment_control_df %>%
  filter(treatment == 1, yr_qtr < 2019.2) %>%
  group_by(Station) %>%
  summarise(start_mean = mean(mean_value, na.rm = TRUE))

# 4. Create the Counterfactual Data by applying the Control Drift to each Treated Start
loess_counterfactuals <- treated_starts %>%
  crossing(control_trend_shape) %>%
  mutate(counterfactual_loess = start_mean + drift)

# 5. Plot
ggplot() +
  # Actual Treated Stations (Individual LOESS lines)
  geom_smooth(data = filter(treatment_control_df, treatment == 1), 
              aes(x = yr_qtr, y = mean_value, color = Station), 
              method = "loess", se = FALSE, size = 1) +
  # Pseudo-Counterfactuals (Shifted Control LOESS - Dashed)
  geom_line(data = loess_counterfactuals, 
            aes(x = yr_qtr, y = counterfactual_loess, color = Station), 
            linetype = "dashed", alpha = 0.7, size = 0.8) +
  # Visual Aids
  geom_vline(xintercept = 2019.2, linetype = "dotted", size = 1) +
  labs(
    title = "Treated Stations vs. Smooth Conditional mean  Counterfactuals",
    subtitle = "Solid = Actual LOESS trend | Dashed = Control Group LOESS shape shifted to station baseline",
    x = "Year.Quarter",
    y = "NO2 Level (Shifted)",
    caption = "The divergence after the dotted line represents the estimated policy impact."
  ) +
  theme_minimal()
ggplot(treatment_control_df, aes(x = mean_value, fill = factor(intervention), color = factor(intervention))) +
  geom_density(alpha = 0.3) +
  # Use facet_wrap to create the side-by-side view
  facet_wrap(~treatment, labeller = as_labeller(c("0" = "Control", "1" = "Treatment"))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_minimal() +
  labs(
    title = "Comparison of Densities(KDE)(NO2)",
    subtitle = "Side-by-side density plots of Treatment vs Control",
    x = "Mean Value",
    y = "Density",
    fill = "Intervention",
    color = "Intervention"
  ) +
  theme(legend.position = "bottom")
#Creating Density Plots#
ggplot(filter(treatment_control_df,treatment==1), aes(x = mean_value, fill = factor(intervention), color = factor(intervention))) +
  geom_density(alpha = 0.3) +
  # Use facet_wrap to create the side-by-side view
  facet_wrap(~Station) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_minimal() +
  labs(
    title = "Comparison of Densities(KDE)(NO2)",
    subtitle = "Side-by-side density plots of Treated London Stations",
    x = "Mean Value",
    y = "Density",
    fill = "Intervention",
    color = "Intervention"
  ) +
  theme(legend.position = "bottom")
ggplot(filter(treatment_control_df,treatment==0), aes(x = mean_value, fill = factor(intervention), color = factor(intervention))) +
  geom_density(alpha = 0.3) +
  # Use facet_wrap to create the side-by-side view
  facet_wrap(~Station) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_minimal() +
  labs(
    title = "Comparison of Densities(KDE)(NO2)",
    subtitle = "Side-by-side density plots of Control Stations",
    x = "Mean Value",
    y = "Density",
    fill = "Intervention",
    color = "Intervention"
  ) +
  theme(legend.position = "bottom")
#Estimating Models#
treatment_control_df$Station <- as.factor(treatment_control_df$Station)
treatment_control_df$yr_qtr <- as.factor(treatment_control_df$yr_qtr)
treatment_control_df$City <- as.factor(treatment_control_df$City)
treatment_control_df$did <-as.factor(treatment_control_df$did)
reducedmodel <- lm(mean_value~did+yr_qtr+Station,data=treatment_control_df)
hetmodel <- lm(mean_value~did:Station+yr_qtr+Station,data=treatment_control_df)
eventmodel <- lm(mean_value~treatment*yr_qtr+Station,data=treatment_control_df)
#Estimating clustered and Jacknife robust se matrices#
hcr_vcov <- vcovHC(eventmodel,type = "HC3")
hcr_vcovh <- vcovHC(hetmodel,type = "HC3")
hcr_vcovr <- vcovHC(reducedmodel,type = "HC3")
cluster_vcov <- vcovCR(eventmodel,type = "CR2",cluster = treatment_control_df$Station)
cluster_vcovh <- vcovCR(hetmodel,type = "CR2",cluster = treatment_control_df$Station)
cluster_vcovr <- vcovCR(reducedmodel,type = "CR2",cluster = treatment_control_df$Station)
#Event Study Plots#
all_coefs <- names(coef(eventmodel))
qtr_coefs <- all_coefs[grep("treatment:yr_qtr", all_coefs)]
qtr_coefs <- qtr_coefs[order(gsub(".*\\)", "", qtr_coefs))]
modelplot(eventmodel,vcov = hcr_vcov,coef_omit = "^(?!.*treatment:)",color="darkblue")+geom_vline(xintercept=0)+labs(title="Event(2019.2) No2 Study,HC3(Jacknife) Corrected SE")+theme_classic()+ geom_hline(yintercept ="treatment:yr_qtr2019.2",linetype = "dashed", colour = "red")+coord_flip()+theme(axis.text.x = element_text(angle = 90))
modelplot(hetmodel,vcov = hcr_vcovh,coef_omit = "^(?!did1:)",color="darkgreen")+geom_vline(xintercept=0)+labs(title="Heterogenous(Station) Average ULEZ Effects")
modelplot(eventmodel,vcov = cluster_vcov,coef_omit = "^(?!.*treatment:)",color="darkblue")+geom_vline(xintercept=0)+labs(title="Event(2019.2) No2 Study,Cluster(Station) Corrected SE")+theme_classic()+ geom_hline(yintercept ="treatment:yr_qtr2019.2",linetype = "dashed", colour = "red")+coord_flip()+theme(axis.text.x = element_text(angle = 90))
modelplot(hetmodel,vcov = cluster_vcovh,coef_omit = "^(?!did1:)",color="darkgreen")+geom_vline(xintercept=0)+labs(title="Heterogenous(Station) Average ULEZ Effects")
#Stargazer Tables with Jacknife and Clustered SE#
aic_row <- c("AIC", round(AIC(eventmodel), 2), round(AIC(hetmodel), 2), round(AIC(reducedmodel), 2))
bic_row <- c("BIC", round(BIC(eventmodel), 2), round(BIC(hetmodel), 2), round(BIC(reducedmodel), 2))
stargazer(hetmodel,reducedmodel,
          se=list(sqrt(diag(hcr_vcovh)),sqrt(diag(hcr_vcovr))),column.labels = c("Heterogeneous Model","Reduced Model"),add.lines = list(aic_row, bic_row),omit.stat = c("f", "ser"),
          header = FALSE,type = "text",title="Pollutant Levels,Station and Quarter FE,HC3(Jacknife)SE",omit = c("^Station", "yr_qtr"))
stargazer(hetmodel,reducedmodel,
          se=list(sqrt(diag(cluster_vcovh)),sqrt(diag(cluster_vcovr))),column.labels = c("Heterogeneous Model","Reduced Model"),add.lines = list(aic_row, bic_row),omit.stat = c("f", "ser"),
          header = FALSE,type = "text",title="Pollutant Levels,Station and Quarter FE,CR2 Cluster(Station) SE",omit = c("^Station", "yr_qtr"))
#Table comparing robust SE estimates#
stargazer(eventmodel, eventmodel,
          se = list(sqrt(diag(cluster_vcov)), sqrt(diag(hcr_vcov))), 
          column.labels = c("Small-Cluster Station (CR2)", "Leverage-Corr (HC3)"),
          add.lines = list(
            c("Clustering", "Station" ,"None (Jacknife)"),
            c("Robustness", "Small-G(15) Adj", "Leverage Adj"),
            # "n" omits the observation count
            omit.stat = c("f", "ser", "rsq", "adj.rsq", "n"),
            header = FALSE,
            type = "text",
            title = "Comparison of Robust Inference Methods (Event Model)",
            omit = c("^Station", "^yr_qtr")))
#Wald Tests#
# Get clean coefficients (no NAs)
betas_event <- na.omit(coef(eventmodel))
# 1. Identify Pre-Policy indices (2014.1 through 2019.1)
pre_indices <- which(grepl("treatment:yr_qtr201[4-8]", names(betas_event)) | 
                       grepl("treatment:yr_qtr2019.1", names(betas_event)))

# 1. Identify Pre-Policy indices (2014.1 through 2018.4)
pre_indicesant <- which(grepl("treatment:yr_qtr201[4-8]", names(betas_event)) | 
                          grepl("treatment:yr_qtr2018.4", names(betas_event)))
# 2. Identify Post-Policy indices (2019.2 through 2022.4)
post_indices <- which(grepl("treatment:yr_qtr2019.[2-4]", names(betas_event)) | 
                        grepl("treatment:yr_qtr202[0-2]", names(betas_event)))
library(aod)
# TEST 1: Pre-Trends (Null = No difference between groups before 2019.2)
wald_pre <- wald.test(Sigma = hcr_vcov, b = betas_event, Terms = pre_indices)
print(wald_pre)
# TEST 1.5: Pre-Trends including anticipation of 1 quarter (Null = No difference between groups before 2019.1)
wald_preanticipation <- wald.test(Sigma = hcr_vcov, b = betas_event, Terms = pre_indicesant)
print(wald_preanticipation)

# TEST 2: Post-Policy (Null = No impact after the start of 2019.2)
wald_post <- wald.test(Sigma = hcr_vcov, b = betas_event, Terms = post_indices)
print(wald_post)
betas <- coef(hetmodel)
betas_clean <- na.omit(betas)
target_indices <- which(grepl("did1:StationLondon", names(betas_clean)))
# 1. Get the names that exist in the vcov matrix
vcov_names <- colnames(hcr_vcov)

# 2. Get indices based on vcov names
post_indices_v <- grep("treatment:yr_qtr2019.[2-4]|treatment:yr_qtr202[0-2]", vcov_names)

# 3. Create matrix matching vcov dimensions
L_homog <- matrix(0, nrow = length(post_indices_v) - 1, ncol = ncol(hcr_vcov))

for(i in 1:(length(post_indices_v) - 1)) {
  L_homog[i, post_indices_v[i]] <- 1      
  L_homog[i, post_indices_v[i+1]] <- -1   
}

#Test 3 Time Homogeneity test(Null=Difference in coefficients=0)
wald_constanttimeeffects <- wald.test(Sigma = hcr_vcov, 
                                      b = coef(eventmodel)[vcov_names], 
                                      L = L_homog)

print(wald_constanttimeeffects)
# 1. Get the full vector of coefficients
betas <- coef(hetmodel)

# 2. Identify the numeric positions (indices) of the London terms
# Note: we use the coefficients that ARE NOT NA to match your hcr_vcovh dimensions
betas_clean <- na.omit(betas)
target_indices <- which(grepl("did1:", names(betas_clean)))
num_targets <- length(target_indices)

# Initialize a matrix of zeros
# Rows = number of restrictions (N-1), Cols = total coefficients
Lhomogens <- matrix(0, nrow = num_targets - 1, ncol = length(betas_clean))

# Fill the matrix to create pairwise differences
for (i in 1:(num_targets - 1)) {
  Lhomogens[i, target_indices[i]] <- 1      # Set current coef to 1
  Lhomogens[i, target_indices[i+1]] <- -1   # Set next coef to -1
}
# 3. Perform the Wald Test
# Sigma is your HC3 matrix, b is the vector of estimated coefficients
joint_test_hc3homogenousstationeffects <- wald.test(Sigma = hcr_vcovh, 
                                                    b = betas_clean, 
                                                    Terms = target_indices)
#Test 4 Individual Station Level Homogeneity test(Null=Difference in station coefficients=0)
joint_test_hc3homogenousstationhypothesis <- wald.test(Sigma = hcr_vcovh, 
                                                       b = betas_clean, 
                                                       L=Lhomogens)
print(joint_test_hc3homogenousstationeffects)
print(joint_test_hc3homogenousstationhypothesis)
print(wald_pre)
print(wald_preanticipation)
print(wald_post)
print(wald_constanttimeeffects)
# Helper function to extract stats from different wald test object types
extract_wald <- function(test_obj, label) {
  # For aod::wald.test objects
  if (inherits(test_obj, "wald.test")) {
    res <- test_obj$result$chi2
    return(tibble(
      Test = label,
      Statistic = res["chi2"],
      df = res["df"],
      p_val = res["P"]
    ))
  } 
  # For data.frame/tibble objects (like your previous waldboot results)
  else if (is.data.frame(test_obj)) {
    return(tibble(
      Test = label,
      Statistic = test_obj$statistic[1],
      df = NA, # boottest often doesn't report standard df
      p_val = test_obj$p.value[1]
    ))
  }
}

# Combine all tests into one table
results_df <- bind_rows(
  extract_wald(joint_test_hc3homogenousstationhypothesis, "Station Homogeneity (HC3)"),
  extract_wald(wald_pre, "Pre-Trend Joint Test"),
  extract_wald(wald_preanticipation, "Pre-Anticipation Effects"),
  extract_wald(wald_post, "Post-Treatment Effects"),
  extract_wald(wald_constanttimeeffects, "Constant Time Effects"))
gt(results_df)
pmmodelr = reducedmodel
pmhetmodel = hetmodel
pmvcovhchet = hcr_vcovh
pmvcovcr= cluster_vcovr
pmvcovhr= hcr_vcovr
stargazer(hetmodel,reducedmodel,pmhetmodel,pmmodelr,se=list(sqrt(diag(hcr_vcovh)),sqrt(diag(hcr_vcovr)),sqrt(diag(pmvcovhchet)),sqrt(diag(pmvcovcr))),column.labels = c("HeterogeneousNo2","ReducedNo2","HeterogeneousPM10","ReducedPM10"),omit.stat = c("f", "ser"),
          header = FALSE,type = "text",title="Pollutant Levels,Station and Quarter FE,HC3(Jacknife)SE",omit = c("^Station", "yr_qtr"))
