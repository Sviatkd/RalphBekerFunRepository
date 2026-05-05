set.seed(124)
library(fwildclusterboot)
library(gt)
#Converting to TWFE#
df_did$Unit <- as.factor(df_did$Unit)
df_did$Year <- as.factor(df_did$Year)
model1 <- lm(CycleTraffic_bn ~ Unit+Year+did, data = df_did)
simpleboot <- boottest(model1,clustid = "Unit",param = "did",type="webb",B=9999)
model_event <- lm(CycleTraffic_bn ~ treated*Year+Unit,data = df_did)
postnames <- c("treated:Year2019","treated:Year2020","treated:Year2021","treated:Year2022","treated:Year2023","treated:Year2024")
allevent <- grep("treated:Year",names(coef(model_pre)),value = TRUE)
prenames <- allevent[!(allevent %in% postnames)]
pretest <- boottest(model_event,clustid = "Unit",param = prenames,type="webb",B=9999)
postCovid <- c("treated:Year2022","treated:Year2023","treated:Year2024")
preandCovid <- c("treated:Year2019","treated:Year2020","treated:Year2021")
#Making a Summary Table#
consttime <- boottest(model_event,clustid = "Unit",param = c(preandCovid,postCovid),R=c(1,1,1,-1,-1,-1),type="webb",B=9999)
reducedmodeltest <- tidy(simpleboot)
pretrendtest <- tidy(pretest)
constanttimeeffecttest <- tidy(consttime)
pretrendtest$term <- "sum(Allprecoefs)=0"
constanttimeeffecttest$term <- "DiDpreandcovid-DiDpostcovid=0"
bootresultscycling <- rbind(pretrendtest,constanttimeeffecttest,reducedmodeltest)
gt(bootresultscycling) %>%
  tab_header(title = "Wild(Cluster=Unit/Region) Bootstrap Inference Results") %>%
  cols_label(term = "Test Specification", 
             estimate = "Estimate", 
             p.value = "p-value") %>%
  # 3. Format numbers (3 decimal places)
  fmt_number(columns = 2:6, decimals = 3) %>%
  # 4. Merge CI into one column (optional but looks nice)
  cols_merge(columns = c(conf.low, conf.high), pattern = "[{1}, {2}]") %>%
  cols_label(conf.low = "95% Conf. Interval")
