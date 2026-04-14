# Load required package for wild cluster bootstrap
library(fwildclusterboot)

# 1. Extract the correct coefficient names from the event model
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

# Start with the model's default SEs as a fallback (or NA if you prefer)
se_wild_full <- summary(eventmodel)$coefficients[, "Std. Error"]
names(se_wild_full) <- all_coef_names

# Overwrite the event‑study coefficients with the bootstrapped SEs
# (wild_se_event is a named vector from your bootstrap loop)
se_wild_full[names(wild_se_event)] <- wild_se_event
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
