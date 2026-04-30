# Load required package for wild cluster bootstrap
# Core mandatory dependencies
install.packages(c("collapse", "dreamerr", "generics", "dqrng", 
                   "Matrix", "Formula", "Rcpp", "RcppArmadillo", 
                   "RcppEigen", "JuliaConnectoR"))

# Recommended modeling & utility dependencies
install.packages(c("fixest", "lfe", "ivreg", "data.table", "modelsummary"))
install.packages('fwildclusterboot', repos ='https://s3alfisc.r-universe.dev')
library(fwildclusterboot)
library(modelsummary)
library(broom)
library(gt)

# 1. Extract the correct coefficient names from the event model#
set.seed(124)
event_terms <- grep("treatment:yr_qtr", names(coef(eventmodel)), value = TRUE)

#Create a vector to store loop results#
wild_se_event <- setNames(rep(NA_real_, length(event_terms)), event_terms)

# Loop through each event term and run Webb‑weighted wild bootstrap
for (ev in event_terms) {
  # Run wild cluster bootstrap for a single coefficient#
  boot_out <- boottest(
    eventmodel,
    clustid = "City",
    param = ev,
    B = 9999,
    type = "webb"
  )
  
  # Store the bootstrap standard error
  if (!is.null(boot_out$se)) {
    wild_se_event[ev] <- boot_out$se
  } else {
    # Fallback: derive SE from the confidence interval
    ci <- boot_out$conf_int
    wild_se_event[ev] <- (ci[2] - ci[1]) / (2 * 1.96)
  }
}

# 4. Verify results – print only the event‑study coefficients
print(wild_se_event[event_terms])
modelplot(eventmodel,vcov = wild_se_event,coef_omit = "^(?!.*treatment:)",color="darkblue")+geom_vline(xintercept=0)+labs(title="NO2 Event(2019.2) Study,WildBootCluster(City) Corrected SE")+theme_classic() + geom_hline(yintercept ="treatment:yr_qtr2019.2",linetype = "dashed", colour = "red")+coord_flip()+theme(axis.text.x = element_text(angle = 90))
# Identify the post-treatment terms#
all_names <- names(coef(eventmodel))
post_treatment_names <- grep("treatment:yr_qtr2019\\.[2-4]|treatment:yr_qtr202[0-2]", 
                             all_names, value = TRUE)
block1_vars <- c("treatment:yr_qtr2019.2", "treatment:yr_qtr2019.3", "treatment:yr_qtr2019.4", 
                 "treatment:yr_qtr2020.1", "treatment:yr_qtr2020.2", "treatment:yr_qtr2020.3", 
                 "treatment:yr_qtr2020.4")

block2_vars <- c("treatment:yr_qtr2021.1", "treatment:yr_qtr2021.2", "treatment:yr_qtr2021.3", 
                 "treatment:yr_qtr2021.4", "treatment:yr_qtr2022.1", "treatment:yr_qtr2022.2", 
                 "treatment:yr_qtr2022.3", "treatment:yr_qtr2022.4")

# Combine into one parameter vector#
comparison_vars <- c(block1_vars, block2_vars)

# Create the R vector (Averages)
# We use 1/7 for the first block and -1/8 for the second block#
r_vector <- c(rep(1/7, 7), rep(-1/8, 8))

# Run the boottest where the null is that the difference in pre and post 2020.4 effects is zero#
boot_split <- boottest(
  object = eventmodel,
  clustid = "City",
  param = comparison_vars,
  R = r_vector,
  r = 0,
  B = 9999,
  type = "webb"
)
joint_test <- boottest(
  object = eventmodel,
  clustid = "City",
  param = post_treatment_names,
  R= r_vector,
  B = 9999,
  type = "webb"
)
#  Extract all pre-treatment interaction names#
pre_treatment_terms <- grep("treatment:yr_qtr2014|treatment:yr_qtr2015|treatment:yr_qtr2016|treatment:yr_qtr2017|treatment:yr_qtr2018|treatment:yr_qtr2019\\.1", 
                            names(coef(eventmodel)), value = TRUE)

#  Joint Wild Bootstrap Test for Pre-Trends
# We test the Null Hypothesis: All pre-treatment coefficients = 0
pre_trend_test <- boottest(
  object = eventmodel,
  clustid = "City",
  param = pre_treatment_terms, 
  B = 9999,
  type = "webb"
)
reduced_test <- boottest(
  object = reducedmodel,
  clustid = "City",
  param = "did1", 
  B = 9999,
  type = "webb"
)
# Get ALL coefficient names from the model#
all_coefs <- coef(hetmodel)

# 2. Filter for the 'did1:Station' terms that were NOT dropped (not NA)
estimated_station_effects <- names(all_coefs)[
  grepl("did1:", names(all_coefs)) & !is.na(all_coefs)
]

# 3. Run the joint test again with the cleaned list
station_het_test <- boottest(
  hetmodel,
  clustid = "City",
  param = estimated_station_effects, 
  B = 9999,
  type = "webb"
)
testpre <- tidy(pre_trend_test)
testreduced <- tidy(reduced_test)
testpre$term <- "Allprecoefs=0"
print(joint_test)
print(station_het_test)
testconsttime <- tidy(boot_split)
testconsttime$term <- "DiDpre2020-DiDpost2020=0"
#Making A Summary table with gt#
waldboot <- rbind(testpre, testconsttime, testreduced)
print(waldboot)
gt(waldboot) %>%
  tab_header(title = "Wild(Cluster=City) Bootstrap Inference Results") %>%
  cols_label(term = "Test Specification", 
             estimate = "Estimate", 
             p.value = "p-value") %>%
  # 3. Format numbers (3 decimal places)
  fmt_number(columns = 2:6, decimals = 3) %>%
  # 4. Merge CI into one column (optional but looks nice)
  cols_merge(columns = c(conf.low, conf.high), pattern = "[{1}, {2}]") %>%
  cols_label(conf.low = "95% Conf. Interval")
