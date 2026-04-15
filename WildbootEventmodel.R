# Load required package for wild cluster bootstrap#
# Core mandatory dependencies#
install.packages(c("collapse", "dreamerr", "generics", "dqrng", 
                   "Matrix", "Formula", "Rcpp", "RcppArmadillo", 
                   "RcppEigen", "JuliaConnectoR"))

# Recommended modeling & utility dependencies#
install.packages(c("fixest", "lfe", "ivreg", "data.table", "modelsummary"))
install.packages('fwildclusterboot', repos ='https://s3alfisc.r-universe.dev')
library(fwildclusterboot)

# 1. Extract the correct coefficient names from the event model#
set.seed(124)
event_terms <- grep("treatment:yr_qtr", names(coef(eventmodel)), value = TRUE)

# 2. Initialize a named vector to store the wild bootstrap standard errors
#    (start with NA to ensure we only fill in the bootstrapped values)
wild_se_event <- setNames(rep(NA_real_, length(event_terms)), event_terms)

# 3. Loop through each event term and run Webb‑weighted wild bootstrap
for (ev in event_terms) {
  # Run wild cluster bootstrap for a single coefficient
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
coefplot::coefplot(eventmodel, 
                   coefficients = qtr_coefs,innerCI=0, vcov=wild_se_event,horizontal = TRUE) + scale_y_discrete(labels = function(x) {
                     clean <- gsub("yr_qtr|treatment:yr_qtr", "", x)
                     ifelse(grepl("\\.1$", clean), clean, "")
                   }) + 
  
  labs(title="NO2 Event(2019.2) Study (Yearly Ticks),WildCluster(City) Corrected SE",
       x="Estimate", y="Year (Q1 Marked)") +
  theme_gray()
stargazer(eventmodel,eventmodel, eventmodel,eventmodel,
          se = list(wild_se_event,sqrt(diag(cluster_vcov)), sqrt(diag(hcr_vcov)),sqrt(diag(cluster_vcovc))), 
          column.labels = c("WildBoot City","Small-Cluster Station (CR2)", "Leverage-Corr (HC3)","City(CR2)"),
          add.lines = list(
            c("Clustering", "City","Station" ,"None (Jacknife)", "City"),
            c("Robustness", "Webbweights","Small-G(15) Adj", "Leverage Adj", "Small-G(5) Adj")),
          # "n" omits the observation count
          omit.stat = c("f", "ser", "rsq", "adj.rsq", "n"),
          header = FALSE,
          type = "text",
          title = "Comparison of Robust Inference Methods (Event Model)",
          omit = c("^Station", "^yr_qtr"))
# Extract all coefficient names in the model
all_coef_names <- names(coef(eventmodel))

# CR2 SE clustered by station
se_cr2_full <- sqrt(diag(cluster_vcov))   # Already a full vector

# HC3 (leverage‑corrected)
se_hc3_full <- sqrt(diag(hcr_vcov))

# CR2 SE clustered by city
se_city_full <- sqrt(diag(cluster_vcovc))

# Compute p‑values using a normal approximation (common with wild bootstrap)
coef_vals <- coef(eventmodel)
p_wild_full <- 2 * pnorm(abs(coef_vals / se_wild_full), lower.tail = FALSE)
stargazer(eventmodel, eventmodel, eventmodel, eventmodel,
          se = list(se_wild_full, se_cr2_full, se_hc3_full, se_city_full),
          p  = list(p_wild_full, 
                    2 * pnorm(abs(coef_vals / se_cr2_full), lower.tail = FALSE),
                    2 * pnorm(abs(coef_vals / se_hc3_full), lower.tail = FALSE),
                    2 * pnorm(abs(coef_vals / se_city_full), lower.tail = FALSE)),
          column.labels = c("WildBoot City", "Small-Cluster Station (CR2)", 
                            "Leverage-Corr (HC3)", "City(CR2)"),
          add.lines = list(
            c("Clustering", "City", "Station", "None (Jacknife)", "City"),
            c("Robustness", "Webb weights", "Small-G(15) Adj", "Leverage Adj", "Small-G(5) Adj")
          ),
          omit.stat = c("f", "ser", "rsq", "adj.rsq", "n"),
          header = FALSE,
          type = "text",
          title = "Comparison of Robust Inference Methods (Event Model)",
          omit = c("^Station", "^yr_qtr"))
# 1. Identify the post-treatment terms
all_names <- names(coef(eventmodel))
post_treatment_names <- grep("treatment:yr_qtr2019\\.[2-4]|treatment:yr_qtr202[0-2]", 
                             all_names, value = TRUE)

# 2. To test if they are ALL EQUAL, we test if the differences are ZERO.
# We will test each quarter against the first post-treatment quarter (2019.2)
baseline <- post_treatment_names[1]
others <- post_treatment_names[-1]

# 3. Use the 'wald_test' function logic or a loop to get the joint significance
# Since boottest() is struggling with the joint syntax, we can use a 'pooled' 
# approach to see if the post-period as a whole differs from a constant.

# Define the Null: All coefficients are equal to their mean
# We can do this by testing the coefficients against 0 and checking the 
# variance of the estimates. 

# RE-TRY with the most basic 'boottest' joint syntax 
# (Wrapping the whole vector in the 'param' argument)
joint_test <- boottest(
  object = eventmodel,
  clustid = "City",
  param = post_treatment_names, # Just the names
  B = 9999,
  type = "webb"
)

print(joint_test)
# 1. Extract all pre-treatment interaction names
pre_treatment_terms <- grep("treatment:yr_qtr2014|treatment:yr_qtr2015|treatment:yr_qtr2016|treatment:yr_qtr2017|treatment:yr_qtr2018|treatment:yr_qtr2019\\.1", 
                            names(coef(eventmodel)), value = TRUE)

# 2. Run the Joint Wild Bootstrap Test for Pre-Trends
# We test the Null Hypothesis: All pre-treatment coefficients = 0
pre_trend_test <- boottest(
  object = eventmodel,
  clustid = "City",
  param = pre_treatment_terms, 
  B = 9999,
  type = "webb"
)

print(pre_trend_test)
reduced_test <- boottest(
  object = reducedmodel,
  clustid = "City",
  param = "did1", 
  B = 9999,
  type = "webb"
)

print(reduced_test)
# 1. Get ALL coefficient names from the model
all_coefs <- coef(het_model_fixed)

# 2. Filter for the 'did1:Station' terms that were NOT dropped (not NA)
estimated_station_effects <- names(all_coefs)[
  grepl("did1:Station", names(all_coefs)) & !is.na(all_coefs)
]

# 3. Run the joint test again with the cleaned list
station_het_test <- boottest(
  het_model_fixed,
  clustid = "City",
  param = estimated_station_effects, 
  B = 9999,
  type = "webb"
)
print(pre_trend_test)
print(reduced_test)
print(joint_test)
print(station_het_test)
