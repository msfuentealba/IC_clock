library(tidyverse)
library(methylclock)
library(survival)
library(survminer)
library(ggpmisc)
library(viridis)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggh4x)
library(ggsci)

raw <- read_rds("./input/fhs/epigenomics.rds") 
clocks <- DNAmAge(raw, clocks = c("Horvath", "Hannum", "Levine"), cell.count = FALSE, toBetas = FALSE)
clocks <- clocks[,c("id","Horvath","Hannum","Levine")] %>% set_names("sample","Horvath","Hannum","PhenoAge")
raw <- raw %>% t()
model <- read_rds("./input/best_model.rds")
data <- raw[,model$gene[2:nrow(model)]] 
data_adjusted <- data %>% t()

predictions <- apply(data_adjusted, 2, function(x) sum(x*model$coefficient[2:nrow(model)])+model$coefficient[1]) %>% 
  enframe %>%
  set_names("sample","ic") %>% 
  data.frame
predictions <- predictions %>% left_join(clocks)
ages <- read_tsv("./input/fhs/phs000007.v33.pht003099.v8.p14.c1.vr_dates_2019_a_1175s.HMB-IRB-MDS.txt", comment = "#")
ages <- ages %>% dplyr::select(shareid, age8, sex) %>% na.omit %>% set_names("sample","age","sex")
ages$sex <- ifelse(ages$sex==2,"female","male")
ages$sample <- as.character(ages$sample)
predictions <- predictions %>% left_join(ages) %>% na.omit

predictions$ic_acc <- resid(lm(1-ic~age+sex, data = predictions))
predictions$horvath_acc <- resid(lm(Horvath~age+sex, data = predictions))
predictions$hannum_acc <- resid(lm(Hannum~age+sex, data = predictions))
predictions$phenoage_acc <- resid(lm(PhenoAge~age+sex, data = predictions))
predictions <- predictions[,c("sample","age","ic","ic_acc","horvath_acc","hannum_acc","phenoage_acc")]

hazard_ratios <- tibble()
ages <- read_tsv("./input/fhs/phs000007.v33.pht003099.v8.p14.c1.vr_dates_2019_a_1175s.HMB-IRB-MDS.txt.gz", comment = "#")

#mortality
mortality <- read_tsv("./input/fhs/phs000007.v33.pht003317.v10.p14.c1.vr_survdth_2019_a_1337s.HMB-IRB-MDS.txt.gz", comment = "#")
mortality <- mortality %>% filter(IDTYPE==1)
mortality$status <- ifelse(!is.na(mortality$DATEDTH),1,0)
mortality$date <- ifelse(mortality$status==1,mortality$DATEDTH,mortality$LASTCON)
mortality$sample <- mortality$shareid
mortality <- mortality %>% left_join(ages %>% dplyr::select(shareid, age1, age8) %>% set_names("sample","age1","age8"))
mortality$age_censor <- mortality$age1+(mortality$date/365)
mortality <- mortality %>% filter((!is.na(age8))&(age_censor>age8))
mortality$sample <- as.character(mortality$sample)
mortality$time <- mortality$age_censor-mortality$age8
mortality <- predictions %>% left_join(mortality) 
mortality <- mortality %>% left_join(ages[,c("shareid","sex")] %>% set_names("sample","sex") %>% mutate(sample = as.character(sample), sex = ifelse(sex==2,"female","male")) %>% unique)
mortality$ic_group <- ifelse(mortality$ic_acc<quantile(mortality$ic_acc,0.2),"accelerated",ifelse(mortality$ic_acc>quantile(mortality$ic_acc,1-0.2),"delayed",NA))
hazard_ratios <- rbind(hazard_ratios, do.call("rbind",lapply(c("ic_acc","horvath_acc","hannum_acc","phenoage_acc"), function(x) {
  df <- mortality[,c("sex","time","status",x)] %>% set_names("sex","time","status","clock")
  cox_model <- coxph(Surv(time, status) ~ scale(clock), data = df)  
  tibble(clock = x, hr = exp(coef(cox_model))[1], ci_lower = exp(coef(cox_model) - 1.96 * sqrt(diag(vcov(cox_model))))[1], ci_upper = exp(coef(cox_model) + 1.96 * sqrt(diag(vcov(cox_model))))[1], p = summary(cox_model)$coefficients[1,5])
})) %>% mutate(outcome = "mortality"))

#atrial fibrillation
af <- read_tsv("./input/fhs/phs000007.v33.pht003315.v10.p14.c1.vr_survaf_2019_a_1338s.HMB-IRB-MDS.txt.gz", comment = "#")
af <- af %>% filter(IDTYPE==1)
af$status <- af$AFIX
af$date <- af$AFIDATE
af$sample <- af$shareid
af <- af %>% left_join(ages %>% dplyr::select(shareid, age1, age8) %>% set_names("sample","age1","age8"))
af$age_censor <- af$age1+(af$date/365)
af <- af %>% filter((!is.na(age8))&(age_censor>age8))
af$sample <- as.character(af$sample)
af$time <- af$age_censor-af$age8
af <- predictions %>% left_join(af) 
af <- af %>% left_join(ages[,c("shareid","sex")] %>% set_names("sample","sex") %>% mutate(sample = as.character(sample), sex = ifelse(sex==2,"female","male")) %>% unique)
af$ic_group <- ifelse(af$ic_acc<quantile(af$ic_acc,0.2),"accelerated",ifelse(af$ic_acc>quantile(af$ic_acc,1-0.2),"delayed",NA))
hazard_ratios <- rbind(hazard_ratios, do.call("rbind",lapply(c("ic_acc","horvath_acc","hannum_acc","phenoage_acc"), function(x) {
  df <- af[,c("sex","time","status",x)] %>% set_names("sex","time","status","clock")
  cox_model <- coxph(Surv(time, status) ~ scale(clock), data = df)  
  tibble(clock = x, hr = exp(coef(cox_model))[1], ci_lower = exp(coef(cox_model) - 1.96 * sqrt(diag(vcov(cox_model))))[1], ci_upper = exp(coef(cox_model) + 1.96 * sqrt(diag(vcov(cox_model))))[1], p = summary(cox_model)$coefficients[1,5])
})) %>% mutate(outcome = "af"))

#cvd
cvd <- read_tsv("./input/fhs/phs000007.v33.pht003316.v10.p14.c1.vr_survcvd_2019_a_1334s.HMB-IRB-MDS.txt.gz", comment = "#")
cvd <- cvd %>% filter(idtype==1)
cvd$status <- cvd$cvd
cvd$date <- cvd$cvddate
cvd$sample <- cvd$shareid
cvd <- cvd %>% left_join(ages %>% dplyr::select(shareid, age1, age8) %>% set_names("sample","age1","age8"))
cvd$age_censor <- cvd$age1+(cvd$date/365)
cvd <- cvd %>% filter((!is.na(age8))&(age_censor>age8))
cvd$sample <- as.character(cvd$sample)
cvd$time <- cvd$age_censor-cvd$age8
cvd <- predictions %>% left_join(cvd) 
cvd <- cvd %>% left_join(ages[,c("shareid","sex")] %>% set_names("sample","sex") %>% mutate(sample = as.character(sample), sex = ifelse(sex==2,"female","male")) %>% unique)
cvd$ic_group <- ifelse(cvd$ic_acc<quantile(cvd$ic_acc,0.2),"accelerated",ifelse(cvd$ic_acc>quantile(cvd$ic_acc,1-0.2),"delayed",NA))
hazard_ratios <- rbind(hazard_ratios, do.call("rbind",lapply(c("ic_acc","horvath_acc","hannum_acc","phenoage_acc"), function(x) {
  df <- cvd[,c("sex","time","status",x)] %>% set_names("sex","time","status","clock")
  cox_model <- coxph(Surv(time, status) ~ scale(clock), data = df)  
  tibble(clock = x, hr = exp(coef(cox_model))[1], ci_lower = exp(coef(cox_model) - 1.96 * sqrt(diag(vcov(cox_model))))[1], ci_upper = exp(coef(cox_model) + 1.96 * sqrt(diag(vcov(cox_model))))[1], p = summary(cox_model)$coefficients[1,5])
})) %>% mutate(outcome = "cvd"))

#chf
chf <- read_tsv("./input/fhs/phs000007.v33.pht003316.v10.p14.c1.vr_survcvd_2019_a_1334s.HMB-IRB-MDS.txt.gz", comment = "#")
chf <- chf %>% filter(idtype==1)
chf$status <- chf$chf
chf$date <- chf$chfdate
chf$sample <- chf$shareid
chf <- chf %>% left_join(ages %>% dplyr::select(shareid, age1, age8) %>% set_names("sample","age1","age8"))
chf$age_censor <- chf$age1+(chf$date/365)
chf <- chf %>% filter((!is.na(age8))&(age_censor>age8))
chf$sample <- as.character(chf$sample)
chf$time <- chf$age_censor-chf$age8
chf <- predictions %>% left_join(chf)
chf <- chf %>% left_join(ages[,c("shareid","sex")] %>% set_names("sample","sex") %>% mutate(sample = as.character(sample), sex = ifelse(sex==2,"female","male")) %>% unique)
chf$ic_group <- ifelse(chf$ic_acc<quantile(chf$ic_acc,0.2),"accelerated",ifelse(chf$ic_acc>quantile(chf$ic_acc,1-0.2),"delayed",NA))
hazard_ratios <- rbind(hazard_ratios, do.call("rbind",lapply(c("ic_acc","horvath_acc","hannum_acc","phenoage_acc"), function(x) {
  df <- chf[,c("sex","time","status",x)] %>% set_names("sex","time","status","clock")
  cox_model <- coxph(Surv(time, status) ~ scale(clock), data = df)  
  tibble(clock = x, hr = exp(coef(cox_model))[1], ci_lower = exp(coef(cox_model) - 1.96 * sqrt(diag(vcov(cox_model))))[1], ci_upper = exp(coef(cox_model) + 1.96 * sqrt(diag(vcov(cox_model))))[1], p = summary(cox_model)$coefficients[1,5])
})) %>% mutate(outcome = "chf"))

#chd
chd <- read_tsv("./input/fhs/phs000007.v33.pht003316.v10.p14.c1.vr_survcvd_2019_a_1334s.HMB-IRB-MDS.txt.gz", comment = "#")
chd <- chd %>% filter(idtype==1)
chd$status <- chd$chd
chd$date <- chd$chddate
chd$sample <- chd$shareid
chd <- chd %>% left_join(ages %>% dplyr::select(shareid, age1, age8) %>% set_names("sample","age1","age8"))
chd$age_censor <- chd$age1+(chd$date/365)
chd <- chd %>% filter((!is.na(age8))&(age_censor>age8))
chd$sample <- as.character(chd$sample)
chd$time <- chd$age_censor-chd$age8
chd <- predictions %>% left_join(chd) 
chd <- chd %>% left_join(ages[,c("shareid","sex")] %>% set_names("sample","sex") %>% mutate(sample = as.character(sample), sex = ifelse(sex==2,"female","male")) %>% unique)
chd$ic_group <- ifelse(chd$ic_acc<quantile(chd$ic_acc,0.2),"accelerated",ifelse(chd$ic_acc>quantile(chd$ic_acc,1-0.2),"delayed",NA))
hazard_ratios <- rbind(hazard_ratios, do.call("rbind",lapply(c("ic_acc","horvath_acc","hannum_acc","phenoage_acc"), function(x) {
  df <- chd[,c("sex","time","status",x)] %>% set_names("sex","time","status","clock")
  cox_model <- coxph(Surv(time, status) ~ scale(clock), data = df)  
  tibble(clock = x, hr = exp(coef(cox_model))[1], ci_lower = exp(coef(cox_model) - 1.96 * sqrt(diag(vcov(cox_model))))[1], ci_upper = exp(coef(cox_model) + 1.96 * sqrt(diag(vcov(cox_model))))[1], p = summary(cox_model)$coefficients[1,5])
})) %>% mutate(outcome = "chd"))

#stroke/tia
stroketia <- read_tsv("./input/fhs/phs000007.v33.pht006024.v5.p14.c1.vr_svstktia_2019_a_1336s.HMB-IRB-MDS.txt.gz", comment = "#")
stroketia <- stroketia %>% filter(idtype==1)
stroketia$status <- stroketia$stroke_tia
stroketia$date <- stroketia$stroketiadate
stroketia$sample <- stroketia$shareid
stroketia <- stroketia %>% left_join(ages %>% dplyr::select(shareid, age1, age8) %>% set_names("sample","age1","age8"))
stroketia$age_censor <- stroketia$age1+(stroketia$date/365)
stroketia <- stroketia %>% filter((!is.na(age8))&(age_censor>age8))
stroketia$sample <- as.character(stroketia$sample)
stroketia$time <- stroketia$age_censor-stroketia$age8
stroketia <- predictions %>% left_join(stroketia) 
stroketia <- stroketia %>% left_join(ages[,c("shareid","sex")] %>% set_names("sample","sex") %>% mutate(sample = as.character(sample), sex = ifelse(sex==2,"female","male")) %>% unique)
stroketia$ic_group <- ifelse(stroketia$ic_acc<quantile(stroketia$ic_acc,0.2),"accelerated",ifelse(stroketia$ic_acc>quantile(stroketia$ic_acc,1-0.2),"delayed",NA))
hazard_ratios <- rbind(hazard_ratios, do.call("rbind",lapply(c("ic_acc","horvath_acc","hannum_acc","phenoage_acc"), function(x) {
  df <- stroketia[,c("sex","time","status",x)] %>% set_names("sex","time","status","clock")
  cox_model <- coxph(Surv(time, status) ~ scale(clock), data = df)  
  tibble(clock = x, hr = exp(coef(cox_model))[1], ci_lower = exp(coef(cox_model) - 1.96 * sqrt(diag(vcov(cox_model))))[1], ci_upper = exp(coef(cox_model) + 1.96 * sqrt(diag(vcov(cox_model))))[1], p = summary(cox_model)$coefficients[1,5])
})) %>% mutate(outcome = "stroketia"))

hazard_ratios$outcome <- ifelse(hazard_ratios$outcome=="mortality","All-cause",
                                ifelse(hazard_ratios$outcome=="af","Atrial fibrillation",
                                       ifelse(hazard_ratios$outcome=="chf","Congestive heart failure",
                                              ifelse(hazard_ratios$outcome=="chd","Coronary heart disease",
                                                     ifelse(hazard_ratios$outcome=="stroketia","Stroke/TIA",
                                                            ifelse(hazard_ratios$outcome=="cvd","Cardiovascular disease",NA))))))
hazard_ratios$outcome <- factor(hazard_ratios$outcome, levels = c("All-cause","Atrial fibrillation","Cardiovascular disease","Congestive heart failure","Coronary heart disease","Stroke/TIA"))
hazard_ratios$clock <- ifelse(hazard_ratios$clock=="ic_acc","Intrinsic capacity",
                              ifelse(hazard_ratios$clock=="phenoage_acc","PhenoAge",
                                     ifelse(hazard_ratios$clock=="horvath_acc","Horvath",
                                            ifelse(hazard_ratios$clock=="hannum_acc","Hannum",NA))))

hazard_ratios$clock <- factor(hazard_ratios$clock, levels = c("Intrinsic capacity","PhenoAge","Horvath","Hannum") %>% rev)
subset <- hazard_ratios %>% filter(outcome%in%c("All-cause","Cardiovascular disease","Congestive heart failure","Stroke/TIA"))
a <- ggplot(data=subset, aes(x=hr, y=clock)) +
  geom_errorbarh(aes(xmin=ci_lower,xmax=ci_upper), height = 0.2) +
  geom_point()+
  xlim(0.8,2)+
  geom_vline(xintercept = 1, lty = 2)+
  theme_pubr(border = TRUE)+
  geom_text(aes(x=1.65, y=clock, label=paste(round(hr,2))), size = 3.5, hjust="inward") +
  geom_text(aes(x=2, y=clock, label=paste(format(p, digits = 3, scientific = TRUE), sep="")), size = 3.5, hjust="inward")+
  geom_text(aes(x=1.6, y = 4.4), label = "HR", stat = "unique")+
  geom_text(aes(x=1.9, y = 4.4), label = "p", stat = "unique")+
  labs(x = "Hazard ratio (95% CI)", y = "Epigenetic clock")+
  facet_wrap(.~outcome, nrow = 1)

#figure 4b
fit_mortality <- survfit(Surv(age_censor, status) ~ ic_group, data = mortality)
fit_af <- survfit(Surv(age_censor, status) ~ ic_group, data = af)
fit_cvd <- survfit(Surv(age_censor, status) ~ ic_group, data = cvd)
fit_chf <- survfit(Surv(age_censor, status) ~ ic_group, data = chf)
fit_chd <- survfit(Surv(age_censor, status) ~ ic_group, data = chd)
fit_stroketia <- survfit(Surv(age_censor, status) ~ ic_group, data = stroketia)

summary_fit <- summary(fit_mortality)
round(summary_fit$table[2,"median"]-summary_fit$table[1,"median"],2)

p1 <- ggsurvplot(fit_mortality,
                 censor = FALSE,
                 legend="left",
                 pval = TRUE,
                 pval.coord = c(55,0.1),
                 legend.labs = c("High IC", "Low IC"),
                 palette = c("#E7B800", "#2E9FDF"),
                 xlim = c(min(mortality$age_censor, na.rm = TRUE),max(mortality$age_censor, na.rm = TRUE)),
                 xlab = "Age")

p2 <- ggsurvplot(fit_af,
                 censor = FALSE,
                 legend="none",
                 ylab = "",
                 pval = TRUE,
                 pval.coord = c(55,0.1),
                 legend.labs = c("High IC", "Low IC"),
                 palette = c("#E7B800", "#2E9FDF"),
                 xlim = c(min(mortality$age_censor, na.rm = TRUE),max(mortality$age_censor, na.rm = TRUE)),
                 xlab = "Age")

p3 <- ggsurvplot(fit_cvd,
                 censor = FALSE,
                 legend="none",
                 ylab = "",
                 pval = TRUE,
                 pval.coord = c(55,0.1),
                 legend.labs = c("High IC", "Low IC"),
                 palette = c("#E7B800", "#2E9FDF"),
                 xlim = c(min(mortality$age_censor, na.rm = TRUE),max(mortality$age_censor, na.rm = TRUE)),
                 xlab = "Age")

p4 <- ggsurvplot(fit_chf,
                 censor = FALSE,
                 legend="none",
                 ylab = "",
                 pval = TRUE,
                 pval.coord = c(55,0.1),
                 legend.labs = c("High IC", "Low IC"),
                 palette = c("#E7B800", "#2E9FDF"),
                 xlim = c(min(mortality$age_censor, na.rm = TRUE),max(mortality$age_censor, na.rm = TRUE)),
                 xlab = "Age")

p5 <- ggsurvplot(fit_chd,
                 censor = FALSE,
                 legend="none",
                 ylab = "",
                 pval = TRUE,
                 pval.coord = c(55,0.1),
                 legend.labs = c("High IC", "Low IC"),
                 palette = c("#E7B800", "#2E9FDF"),
                 xlim = c(min(mortality$age_censor, na.rm = TRUE),max(mortality$age_censor, na.rm = TRUE)),
                 xlab = "Age")

p6 <- ggsurvplot(fit_stroketia,
                 censor = FALSE,
                 legend="none",
                 ylab = "",
                 pval = TRUE,
                 pval.coord = c(55,0.1),
                 palette = c("#E7B800", "#2E9FDF"),
                 xlim = c(min(mortality$age_censor, na.rm = TRUE),max(mortality$age_censor, na.rm = TRUE)),
                 xlab = "Age") 


b <- ggarrange(p1$plot + ggtitle("All-cause") + theme(plot.title = element_text(hjust = 0.5), plot.margin = unit(c(1, 0, 1, 0), "lines")),
                 p3$plot + ggtitle("Cardiovascular disease") + theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(1, 0, 1, 0), "lines")),
                 p4$plot + ggtitle("Congestive heart failure") + theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(1, 0, 1, 0), "lines")),
                 p6$plot + ggtitle("Stroke/TIA") + theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(1, 0.5, 1, 0), "lines")),
                 nrow = 1, widths = c(1.55,1,1,1))

#figure 4c
extract_data_dict <- function(file){
  library(xml2)
  library(rvest)
  library(stringr)
  results <- read_xml(file)
  recs <- xml_find_all(results, "//variable")
  recs <- recs[3:length(recs)]
  df <- data.frame(description = xml_text(xml_find_first(recs, "./description")), stringsAsFactors = FALSE)
  rownames(df) <- xml_text(xml_find_first(recs, "./name"))
  df$description <- str_to_sentence(df$description)
  df$description <- gsub(" for exam 8 callback","",df$description)
  return(df %>% rownames_to_column(var = "feature"))
}

dates <- read_tsv("./input/fhs/phs000007.v33.pht003099.v8.p14.c1.vr_dates_2019_a_1175s.HMB-IRB-MDS.txt", comment = "#")
dates <- dates %>% dplyr::select(shareid, date8) %>% na.omit %>% set_names("sample","exam_date")
within_range <- function(sample,date){
  df <- tibble(sample,date)
  df$id <- 1:nrow(df)
  df <- df %>% left_join(dates)
  df$dif <- abs(df$exam_date-df$date)
  df <- df %>% na.omit %>% group_by(sample) %>% summarise(subject = id[which.min(dif)])
  return(df$subject)
}

raw <- read_rds("./input/fhs/epigenomics.rds") 
raw <- raw %>% t()
model <- read_rds("./input/best_model.rds")
data <- raw[,model$gene[2:nrow(model)]] 
predictions <- apply(data, 1, function(x) sum(x*model$coefficient[2:nrow(model)])+model$coefficient[1]) %>% 
  enframe %>%
  set_names("sample","ic") %>% 
  data.frame
ages <- read_tsv("./input/fhs/phs000007.v33.pht003099.v8.p14.c1.vr_dates_2019_a_1175s.HMB-IRB-MDS.txt", comment = "#")
ages <- ages %>% dplyr::select(shareid, age8, sex, date8) %>% na.omit %>% set_names("sample","age","sex","date")
ages$sample <- as.character(ages$sample)
ages$sex <- ifelse(ages$sex==1,"male","female")
predictions <- predictions %>% left_join(ages) %>% na.omit

all_main <- tibble()

#clinic lab assay
lab <- read_tsv("./input/fhs/phs000007.v33.pht000742.v8.p14.c1.fhslab1_8s.HMB-IRB-MDS.txt.gz", comment = "#")
lab <- lab %>% filter(IDTYPE==1&FASTING==1)
lab <- lab %>% column_to_rownames(var = "shareid")
lab <- lab[,4:11]
directions <- tibble(feature = colnames(lab), dir = c(-1,1,-1,-1,-1,-1,-1,-1))
lab <- reshape2::melt(lab %>% as.matrix) %>% mutate(sample = as.character(Var1), feature = as.character(Var2)) %>% dplyr::select(sample,feature,value)
lab <- lab %>% mutate(dataset = "Clinical lab assay") %>% left_join(directions)
all_main <- rbind(all_main,lab)

#inflammatory markers
markers <- read_tsv("./input/fhs/phs000007.v33.pht012882.v1.p14.c1.l_inflamm_ex08_1b_1106s.HMB-IRB-MDS.txt.gz", comment = "#")
markers <- markers %>% filter(idtype==1)
markers <- markers %>% column_to_rownames(var = "shareid")
markers <- markers[,3:12]
directions <- tibble(feature = colnames(markers), dir = rep(-1,10))
markers <- reshape2::melt(markers %>% as.matrix) %>% mutate(sample = as.character(Var1), feature = as.character(Var2)) %>% dplyr::select(sample,feature,value)
markers <- markers %>% mutate(dataset = "Inflammatory markers") %>% left_join(directions)
all_main <- rbind(all_main,markers)

#iAge
iage <- read_rds("./input/fhs/iage.rds")
iage <- tibble(sample = iage$sample, feature = "iAge", value = iage$iage, dataset = "Inflammatory markers", dir = -1)
all_main <- rbind(all_main,iage)

#risk factors
risk <- read_tsv("./input/fhs/phs000007.v33.pht006027.v4.p14.c1.vr_wkthru_ex09_1_1001s.HMB-IRB-MDS.txt.gz", comment = "#")
risk <- risk %>% filter(IDTYPE==1)
risk <- risk %>% column_to_rownames(var = "shareid")
risk <- risk[,c("BG8","BMI8","CALC_LDL8","CPD8","CREAT8","CURRSMK8","DBP8","DLVH8","FASTING_BG8","HDL8","SBP8","TC8","TRIG8","VENT_RT8","DMRX8","HRX8","LIPRX8")]
directions <- tibble(feature = colnames(risk), dir = rep(-1,17))
risk <- reshape2::melt(risk %>% as.matrix) %>% mutate(sample = as.character(Var1), feature = as.character(Var2)) %>% dplyr::select(sample,feature,value)
risk <- risk %>% mutate(dataset = "Cardiovascular risk factors") %>% left_join(directions)
all_main <- rbind(all_main,risk)

#insulin
insulin <- read_tsv("./input/fhs/phs000007.v33.pht003901.v5.p14.c1.l_insulin_2008_m_0704s.HMB-IRB-MDS.txt.gz", comment = "#")
insulin <- insulin %>% filter(IDTYPE==1&INSULIN_I==2)
insulin <- insulin %>% column_to_rownames(var = "shareid")
insulin <- insulin[,3,drop=FALSE] %>% set_names("Insulin (plasma)")
directions <- tibble(feature = colnames(insulin), dir = rep(1,1))
insulin <- reshape2::melt(insulin %>% as.matrix) %>% mutate(sample = as.character(Var1), feature = as.character(Var2)) %>% dplyr::select(sample,feature,value)
insulin <- insulin %>% mutate(dataset = "Insulin") %>% left_join(directions)
all_main <- rbind(all_main,insulin)

#tau
tau <- read_tsv("./input/fhs/phs000007.v33.pht009757.v3.p14.c1.l_tau_2011_a_1022s.HMB-IRB-MDS.txt.gz", comment = "#")
tau <- tau %>% filter(IDTYPE==1)
tau <- tau %>% column_to_rownames(var = "shareid")
tau <- tau[,3,drop=FALSE] %>% set_names("Total Tau (plasma)")
directions <- tibble(feature = colnames(tau), dir = rep(1,1))
tau <- reshape2::melt(tau %>% as.matrix) %>% mutate(sample = as.character(Var1), feature = as.character(Var2)) %>% dplyr::select(sample,feature,value)
tau <- tau %>% mutate(dataset = "Tau") %>% left_join(directions)
all_main <- rbind(all_main,tau)

#Pulmonary Function Test 
pft <- read_tsv("./input/fhs/phs000007.v33.pht000832.v8.p14.c1.pft1_8s.HMB-IRB-MDS.txt.gz", comment = "#")
pft <- pft %>% filter(idtype==1)
pft <- pft %>% column_to_rownames(var = "shareid")
pft <- pft[,c("fv1_8_1","fv6_8_1","fv3_8_1","fvc_8_1","pf_8_1","mmf_8_1","ff1_8_1","ff2_8_1","ff3_8_1")] 
directions <- tibble(feature = colnames(pft), dir = c(1,1,1,1,1,1,1,1,1))
pft <- reshape2::melt(pft %>% as.matrix) %>% mutate(sample = as.character(Var1), feature = as.character(Var2)) %>% dplyr::select(sample,feature,value)
pft <- pft %>% mutate(dataset = "Pulmonary Function Test") %>% left_join(directions)
all_main <- rbind(all_main,pft)

#Bone mineral density
bmd <- read_tsv("./input/fhs/phs000007.v33.pht003096.v4.p14.c1.t_bmdhs_2008_1_0748s.HMB-IRB-MDS.txt.gz", comment = "#")
bmd <- bmd %>% filter(idtype==1)
bmd <- bmd %>% column_to_rownames(var = "shareid")
bmd <- bmd[,c("f8cbnbmd","f8cbtobmd","f8cbtrbmd","s8cbl2bd","s8cbl3bd","s8cbl4bd","s8cbl24bd")] 
directions <- tibble(feature = colnames(bmd), dir = rep(1,7))
bmd <- reshape2::melt(bmd %>% as.matrix) %>% mutate(sample = as.character(Var1), feature = as.character(Var2)) %>% dplyr::select(sample,feature,value)
bmd <- bmd %>% mutate(dataset = "Bone mineral density") %>% left_join(directions)
all_main <- rbind(all_main,bmd)

ex8 <- read_tsv("./input/fhs/phs000007.v33.pht000747.v8.p14.c1.ex1_8s.HMB-IRB-MDS.txt.gz", comment = "#")

#Main - Rosow-Breslau Questions
main <- ex8 %>% column_to_rownames(var = "shareid")
main <- main[,c("H468","H469","H470")] #%>% na.omit
directions <- tibble(feature = colnames(main), dir = c(1,1,1))
main <- reshape2::melt(main %>% as.matrix) %>% mutate(sample = as.character(Var1), feature = as.character(Var2)) %>% dplyr::select(sample,feature,value)
main <- main %>% mutate(dataset = "Rosow-Breslau Questions") %>% left_join(directions)
all_main <- rbind(all_main,main)

#Main - CES-D Scale
main <- ex8 %>% column_to_rownames(var = "shareid")
main <- main[,c("H471","H472")] #%>% na.omit
directions <- tibble(feature = colnames(main), dir = c(-1,-1))
main <- reshape2::melt(main %>% as.matrix) %>% mutate(sample = as.character(Var1), feature = as.character(Var2)) %>% dplyr::select(sample,feature,value)
main <- main %>% mutate(dataset = "CES-D Scale") %>% left_join(directions)
all_main <- rbind(all_main,main)

#Main - Katz Activities of Daily Living Scale
main <- ex8 %>% column_to_rownames(var = "shareid")
main <- main[,c("H474","H475","H476","H477","H478")] #%>% na.omit
directions <- tibble(feature = colnames(main), dir = c(-1,-1,-1,-1,-1))
main <- reshape2::melt(main %>% as.matrix) %>% mutate(sample = as.character(Var1), feature = as.character(Var2)) %>% dplyr::select(sample,feature,value)
main <- main %>% mutate(dataset = "Katz Activities of Daily Living Scale") %>% left_join(directions)
all_main <- rbind(all_main,main)

#Main - Physical Activity Questionnaire
main <- ex8 %>% column_to_rownames(var = "shareid")
main <- main[,c("H481","H482","H483","H484","H485","H486")] #%>% na.omit #removed "H480" u-shape
directions <- tibble(feature = colnames(main), dir = c(-1,1,1,1,1,1))
main <- reshape2::melt(main %>% as.matrix) %>% mutate(sample = as.character(Var1), feature = as.character(Var2)) %>% dplyr::select(sample,feature,value)
main <- main %>% mutate(dataset = "Physical Activity Questionnaire") %>% left_join(directions)
all_main <- rbind(all_main,main)

#Main - Activities past year (2 week)
main <- ex8 %>% column_to_rownames(var = "shareid")
main <- main[,c("H489","H494","H499","H504","H509","H514","H519","H524","H529","H534","H539","H544","H549","H554","H559","H564")] %>% na.omit
directions <- tibble(feature = colnames(main), dir = rep(1,16))
main <- reshape2::melt(main %>% as.matrix) %>% mutate(sample = as.character(Var1), feature = as.character(Var2)) %>% dplyr::select(sample,feature,value)
main <- main %>% mutate(dataset = "Physical Activity Questionnaire") %>% left_join(directions)
all_main <- rbind(all_main,main)

#Main - Activities past year
main <- ex8 %>% column_to_rownames(var = "shareid")
main <- main[,c("H488","H518","H523","H548","H538","H553","H493","H498","H503","H508","H513","H528","H533","H543","H558","H563")] #%>% na.omit
main[main==8] <- NA
directions <- tibble(feature = colnames(main), dir = rep(1,16))
main <- reshape2::melt(main %>% as.matrix) %>% mutate(sample = as.character(Var1), feature = as.character(Var2)) %>% dplyr::select(sample,feature,value)
main <- main %>% mutate(dataset = "Physical Activity Questionnaire") %>% left_join(directions)
all_main <- rbind(all_main,main)

#Main - Nagi questions
main <- ex8 %>% column_to_rownames(var = "shareid")
main <- main[,c("H569","H570","h571","h572","H573","H574","H575","H576","H577")] %>% na.omit
directions <- tibble(feature = colnames(main), dir = rep(-1,9))
main <- reshape2::melt(main %>% as.matrix) %>% mutate(sample = as.character(Var1), feature = as.character(Var2)) %>% dplyr::select(sample,feature,value)
main <- main %>% mutate(dataset = "Nagi questions") %>% left_join(directions)
all_main <- rbind(all_main,main)

#MMSE (higher better)
mmse <- read_tsv("./input/fhs/phs000007.v33.pht005174.v4.p14.c1.vr_mmse_ex09_1b_0943s.HMB-IRB-MDS.txt.gz", comment = "#")
mmse <- mmse %>% filter(exam==8&idtype==1)
mmse <- mmse %>% column_to_rownames(var = "shareid")
mmse <- mmse[,5:20] 
directions <- tibble(feature = colnames(mmse), dir = rep(1,16))
mmse <- reshape2::melt(mmse %>% as.matrix) %>% mutate(sample = as.character(Var1), feature = as.character(Var2)) %>% dplyr::select(sample,feature,value)
mmse <- mmse %>% mutate(dataset = "Mini-Mental State Examination") %>% left_join(directions)
all_main <- rbind(all_main,mmse)

#Factors affecting cognition
main <- ex8 %>% column_to_rownames(var = "shareid")
main <- main[,c("H600","H601","H602")] %>% na.omit
main <- main %>% mutate(across(everything(), ~ recode(., `1` = 2, `2` = 1)))
directions <- tibble(feature = colnames(main), dir = rep(-1,3))
main <- reshape2::melt(main %>% as.matrix) %>% mutate(sample = as.character(Var1), feature = as.character(Var2)) %>% dplyr::select(sample,feature,value)
main <- main %>% mutate(dataset = "Factors affecting cognition") %>% left_join(directions)
all_main <- rbind(all_main,main)

#Main - Hand grip strength
main <- ex8 %>% column_to_rownames(var = "shareid")
main <- main[,c("H605","H606","H607","H608","H609","H610")] %>% na.omit
directions <- tibble(feature = colnames(main), dir = rep(1,6))
main <- reshape2::melt(main %>% as.matrix) %>% mutate(sample = as.character(Var1), feature = as.character(Var2)) %>% dplyr::select(sample,feature,value)
main <- main %>% mutate(dataset = "Hand grip strength") %>% left_join(directions)
all_main <- rbind(all_main,main)

#Main - Walk time
main <- ex8 %>% column_to_rownames(var = "shareid")
main <- main[,c("H617","H620","H623")] %>% na.omit
directions <- tibble(feature = colnames(main), dir = rep(-1,3))
main <- reshape2::melt(main %>% as.matrix) %>% mutate(sample = as.character(Var1), feature = as.character(Var2)) %>% dplyr::select(sample,feature,value)
main <- main %>% mutate(dataset = "Walk time") %>% left_join(directions)
all_main <- rbind(all_main,main)

#Main - SF-12® Health Survey (Standard) Self-administered
main <- ex8 %>% column_to_rownames(var = "shareid")
main <- main[,c("H714","H715","H716","H717","H718","H719","H720","H721","H722","H723","H724","H725")] %>% na.omit
directions <- tibble(feature = colnames(main), dir = c(1,-1,-1,-1,-1,-1,-1,-1,1,1,-1,-1))
main <- reshape2::melt(main %>% as.matrix) %>% mutate(sample = as.character(Var1), feature = as.character(Var2)) %>% dplyr::select(sample,feature,value)
main <- main %>% mutate(dataset = "SF-12 Health Survey") %>% left_join(directions)
all_main <- rbind(all_main,main)

head(all_main)
dir <- all_main[,c("feature","dir")] %>% unique %>% set_names("feature","age_dir")
mat <- reshape2::acast(all_main, sample~feature, value.var = "value") %>% na.omit
mat <- mat[,apply(mat, 2, function(x) length(unique(x %>% na.omit)))>=3]
dim(mat)
dir <- dir %>% filter(feature%in%colnames(mat))
mat <- mat[,dir$feature]

scaled_data <- scale(mat)
adjusted_scaled_data <- sweep(scaled_data, 2, dir$age_dir, `*`)
pca_result <- prcomp(adjusted_scaled_data, center = TRUE, scale. = FALSE)
pca_data <- data.frame(PC1 = pca_result$x[, 1], PC2 = pca_result$x[, 2]) %>% rownames_to_column(var = "sample") %>% left_join(predictions[,c("sample","ic","age")])

variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
pc1_var <- round(variance_explained[1] * 100, 2)
pc2_var <- round(variance_explained[2] * 100, 2)

c <- ggplot(pca_data, aes(x = PC1, y = PC2, color = ic)) +
  geom_point() +
  scale_color_viridis() +
  xlab(paste0('PC1 (', pc1_var, '% variance)')) +
  ylab(paste0('PC2 (', pc2_var, '% variance)')) +
  theme_pubr()+
  labs(color = "DNAm IC")+
  theme(legend.position = c(0.2, 0.2))

#figure 4d
raw <- read_rds("./input/fhs/epigenomics.rds") 
raw <- raw %>% t()
model <- read_rds("./input/best_model.rds")
data <- raw[,model$gene[2:nrow(model)]] 
predictions <- apply(data, 1, function(x) sum(x*model$coefficient[2:nrow(model)])+model$coefficient[1]) %>% 
  enframe %>%
  set_names("sample","ic") %>% 
  data.frame
ages <- read_tsv("./input/fhs/phs000007.v33.pht003099.v8.p14.c1.vr_dates_2019_a_1175s.HMB-IRB-MDS.txt", comment = "#")
ages <- ages %>% dplyr::select(shareid, age8, sex, date8) %>% na.omit %>% set_names("sample","age","sex","date")
ages$sample <- as.character(ages$sample)
ages$sex <- ifelse(ages$sex==1,"male","female")
predictions <- predictions %>% left_join(ages) %>% na.omit

all_main$dir <- NULL
i <- 0.2
predictions$ic_acc <- resid(lm(ic~age, data = predictions))
predictions$ic_group <- ifelse(predictions$ic_acc<quantile(predictions$ic_acc,0+i),0,ifelse(predictions$ic_acc>quantile(predictions$ic_acc,1-i),1,NA))
combined <- rbind(all_main) %>% left_join(predictions) %>% filter(!is.na(feature)) 

model_list <- lapply(unique(combined$feature), function(x) {
  print(x)
  subset_data <- filter(combined, feature == x) %>% na.omit
  model <- glm(ic_group ~ value, data = subset_data, family = binomial)
  cor_age <- cor(subset_data$value,subset_data$age)
  summary(model)$coefficients %>% data.frame %>% rownames_to_column(var = "variable") %>% mutate(feature=x, dataset = unique(subset_data$dataset), dir = unique(subset_data$dir), cor_age) %>% filter(variable=="value")
})

model_list <- do.call("rbind",model_list) %>% mutate(score = -log10(`Pr...z..`), dir = sign(Estimate))

dict <- rbind(extract_data_dict("./input/fhs/phs000007.v33.pht000747.v8.ex1_8s.data_dict.xml"),
              extract_data_dict("./input/fhs/phs000007.v33.pht000832.v8.pft1_8s.data_dict.xml"),
              extract_data_dict("./input/fhs/phs000007.v33.pht003096.v4.t_bmdhs_2008_1_0748s.data_dict.xml"),
              extract_data_dict("./input/fhs/phs000007.v33.pht005174.v4.vr_mmse_ex09_1b_0943s.data_dict.xml"),
              extract_data_dict("./input/fhs/phs000007.v33.pht000742.v8.fhslab1_8s.data_dict.xml"),
              extract_data_dict("./input/fhs/phs000007.v33.pht012882.v1.l_inflamm_ex08_1b_1106s.data_dict.xml"),
              extract_data_dict("./input/fhs/phs000007.v33.pht006027.v4.vr_wkthru_ex09_1_1001s.data_dict.xml"))

model_list <- model_list %>% left_join(dict)
model_list$description <- gsub("\\s*\\([^\\)]+\\)", "", model_list$description)
model_list$description <- gsub("During the past year -","During the past year,",model_list$description)
model_list$description <- str_to_sentence(sapply(model_list$description, function(x) strsplit(x, split = " - ")[[1]][1]))
model_list$description[model_list$dataset=="Insulin"] <- "Insulin"
model_list$description[model_list$dataset=="Tau"] <- "Tau"
model_list$description[model_list$feature=="iAge"] <- "iAge"
model_list <- model_list %>% filter(!description%in%c("During the past year, other","During the past year, in a typical 2 week period of time how often do you other?"))
model_list$description <- gsub("\\, exam 8$","",model_list$description)

color_mapping <- data.frame(
  type_test = c(
    "Physical Health and Function", "Physical Health and Function", "Physical Health and Function", 
    "Physical Health and Function", "Physical Health and Function", 
    "Mental Health and Cognitive Function", "Mental Health and Cognitive Function", 
    "Daily Living and Functional Activities", "Daily Living and Functional Activities", 
    "Daily Living and Functional Activities", 
    "Health Surveys and Quality of Life", 
    "Biomarkers and Clinical Measurements", "Biomarkers and Clinical Measurements", 
    "Biomarkers and Clinical Measurements", "Biomarkers and Clinical Measurements"
  ),
  dataset = c(
    "Bone mineral density", "Hand grip strength", "Physical Activity Questionnaire", 
    "Pulmonary Function Test", "Walk time", 
    "CES-D Scale", "Mini-Mental State Examination", 
    "Katz Activities of Daily Living Scale", "Nagi questions", 
    "Rosow-Breslau Questions", 
    "SF-12 Health Survey", 
    "Tau", "Clinical lab assay", "Inflammatory markers", 
    "Cardiovascular risk factors"
  ),
  color = c(
    "#1f78b4", "#6baed6", "#9ecae1", "#c6dbef", "#deebf7", 
    "#31a354", "#74c476", 
    "#fd8d3c", "#fdae6b", "#fdd0a2", 
    "#756bb1", 
    "#a50f15", "#de2d26", "#fc9272", "#fee0d2"
  ),
  stringsAsFactors = FALSE
)

model_list <- model_list %>% left_join(color_mapping)
model_list <- model_list %>%
  filter(score > (-log10(0.05))) %>%
  mutate(description_num = as.numeric(factor(description)))
model_list$description <- factor(model_list$description, levels = model_list %>% arrange(score) %>% pull(description) %>% unique)

color_direct_mapping <- setNames(model_list$color, model_list$dataset)
model_list$dataset <- factor(model_list$dataset, levels = color_mapping$dataset)
model_list$dir <- ifelse(model_list$dir==1,"Positive association", "Negative association")
model_list$dir <- factor(model_list$dir, levels = c("Positive association", "Negative association"))

d <- ggplot(model_list, aes(x = description, y = score)) +
  geom_tile(aes(fill = cor_age), width = 0.9, height = Inf, alpha = 0.1) +
  geom_segment(aes(x = description, xend = description, y = 0, yend = score, color = dataset)) +
  geom_point(size = 3.5, color = "black") +
  geom_point(size = 3, aes(color = dataset)) +
  theme_pubr(border = TRUE) +
  scale_color_manual(values = color_direct_mapping) +
  facet_grid(dir~., scale="free", space="free")+
  coord_flip() +
  scale_fill_gradient2(low = "#2166ac",mid = "white",high = "#b2182b")+
  theme(legend.position = "bottom") +
  labs(x = "", y = "-log10 (P-value)", fill = "Age correlation (R)", color = "Dataset")

#figure 4e
extract_data_dict <- function(cols,file){
  library(xml2)
  library(rvest)
  library(stringr)
  results <- read_xml(file)
  recs <- xml_find_all(results, "//variable")
  recs <- recs[3:length(recs)]
  df <- data.frame(description = xml_text(xml_find_first(recs, "./description")), stringsAsFactors = FALSE)
  rownames(df) <- xml_text(xml_find_first(recs, "./name"))
  df$description <- gsub(" for exam 8 callback","",df$description)
  df$description <- gsub("sample type: ","",df$description)
  return(df[cols,])
}

raw <- read_rds("./input/fhs/epigenomics.rds") 
raw <- raw %>% t()
model <- read_rds("./input/best_model.rds")
data <- raw[,model$gene[2:nrow(model)]] 
data_adjusted <- data %>% t()
predictions <- apply(data_adjusted, 2, function(x) sum(x*model$coefficient[2:nrow(model)])+model$coefficient[1]) %>% 
  enframe %>%
  set_names("sample","ic") %>% 
  data.frame
ages <- read_tsv("./input/fhs/phs000007.v33.pht003099.v8.p14.c1.vr_dates_2019_a_1175s.HMB-IRB-MDS.txt", comment = "#")
ages <- ages %>% dplyr::select(shareid, age8, sex, date8) %>% na.omit %>% set_names("sample","age","sex","date")
ages$sample <- as.character(ages$sample)
ages$sex <- ifelse(ages$sex==1,"male","female")
predictions <- predictions %>% left_join(ages) %>% na.omit
predictions$ic_acc <- resid(lm(ic~age, data = predictions))

npg_colors <- pal_npg("nrc")(10)

#flavonoid
flavonoids <- read_tsv("./input/fhs/phs000007.v33.pht012913.v1.p14.c1.vr_flavind_ex08_1_1309s.HMB-IRB-MDS.txt.gz", comment = "#")
flavonoids <- flavonoids %>% filter(IDTYPE==1&EXAM==8)
flavonoids <- flavonoids %>% column_to_rownames(var = "shareid")
flavonoids <- flavonoids[,5:44]
colnames(flavonoids) <- extract_data_dict(colnames(flavonoids),"./input/fhs/phs000007.v33.pht012913.v1.vr_flavind_ex08_1_1309s.data_dict.xml")
colnames(flavonoids) <- gsub("\\, USDA|\\, USDA\\, 2007| USDA 2003|\\, USDA 2007|\\(revised March 2013 for no banana contribution\\)|\\, 2007|\\(corrected udlp 2013\\)","",colnames(flavonoids))
flavonoids <- reshape2::melt(flavonoids %>% as.matrix) %>% mutate(sample = as.character(Var1), feature = as.character(Var2)) %>% dplyr::select(sample,feature,value)
flavonoids <- flavonoids %>% left_join(predictions[,c("sample","ic_acc")]) %>% na.omit
flavonoids <- flavonoids %>% group_by(feature) %>% summarise(r = cor.test(value,ic_acc, method = "s")$estimate, p = cor.test(value,ic_acc, method = "s")$p.value)
flavonoids$fdr <- p.adjust(flavonoids$p)

p1 <- ggplot(flavonoids, aes(r, -log10(p)))+
  geom_point(fill = "grey", shape = 21, color = "grey")+
  geom_point(data = flavonoids %>% filter(p<0.05&r>0), fill = npg_colors[8], shape = 21)+
  geom_point(data = flavonoids %>% filter(p<0.05&r<0), fill = npg_colors[4], shape = 21)+
  geom_point(data = flavonoids %>% filter(fdr<0.05&r>0), fill = npg_colors[8], shape = 23, size = 3)+
  geom_point(data = flavonoids %>% filter(fdr<0.05&r<0), fill = npg_colors[4], shape = 23, size = 3)+
  geom_text_repel(data = flavonoids %>% filter(p<0.05&r>0) %>% arrange(p), aes(label = feature), max.overlaps = 100, segment.size = 0.1, size = 3, xlim = c(1, 0.07), direction = "y")+
  coord_cartesian(expand = T, clip = "off")+
  theme_pubr(border = TRUE)+
  ggtitle(label = "Flavonoids (Derived intake)")+
  theme(plot.margin = unit(c(1, 5, 1, 5), "lines"))+
  geom_hline(yintercept = 0, lty = 2)+
  geom_vline(xintercept = 0, lty = 2)+
  labs(x = "Pearson correlation coefficient (R)", y = "-log10 (P-value)")

#blood cells fatty acids
fatty_acids <- read_tsv("./input/fhs/phs000007.v33.pht009756.v3.p14.c1.l_rbcfae_ex08_1b_1072s.HMB-IRB-MDS.txt.gz", comment = "#")
fatty_acids <- fatty_acids %>% filter(idtype==1&exam==8)
fatty_acids <- fatty_acids %>% column_to_rownames(var = "shareid")
fatty_acids <- fatty_acids[,2:29,drop=FALSE]
colnames(fatty_acids) <- extract_data_dict(colnames(fatty_acids),"~/data/FHS/clinical_data/HMB-IRB-MDS/phs000007.v33.pht009756.v3.l_rbcfae_ex08_1b_1072s.data_dict.xml")
colnames(fatty_acids) <- gsub(" \\(Sample type: red blood cells\\)|  expanded \\(Sample type: red blood cells\\)|   expanded \\(Sample type: red blood cells\\)","",colnames(fatty_acids))
fatty_acids <- reshape2::melt(fatty_acids %>% as.matrix) %>% mutate(sample = as.character(Var1), feature = as.character(Var2)) %>% dplyr::select(sample,feature,value)
fatty_acids <- fatty_acids %>% left_join(predictions[,c("sample","ic_acc")]) %>% na.omit
fatty_acids <- fatty_acids %>% group_by(feature) %>% summarise(r = cor.test(value,ic_acc, method = "s")$estimate, p = cor.test(value,ic_acc, method = "s")$p.value)
fatty_acids$fdr <- p.adjust(fatty_acids$p)

#diet
food <- read_tsv("./input/fhs/phs000007.v33.pht002350.v7.p14.c1.vr_ffreq_ex08_1_0615s.HMB-IRB-MDS.txt.gz", comment = "#")
food <- food %>% filter(idtype==1)
food <- food %>% column_to_rownames(var = "shareid")
food <- food[,3:ncol(food),drop=FALSE]

results <- read_xml("./input/fhs/phs000007.v33.pht002350.v7.vr_ffreq_ex08_1_0615s.data_dict.xml")
recs <- xml_find_all(results, "//variable")
recs <- recs[4:length(recs)]
annotation <- tibble()
for (id in 1:length(recs)) {
  print(id)
  variable_name <- xml_text(xml_find_first(recs[id], ".//name"))
  values <- xml_find_all(recs[id], ".//value")
  if(length(values)>1){
    value_codes <- xml_attr(values, "code")
    value_descriptions <- xml_text(values)
    value_df <- data.frame(code = value_codes, description = value_descriptions, stringsAsFactors = FALSE)
    value_df$name <- variable_name
    annotation <- rbind(annotation,value_df)  
  }
}
table(annotation$description) %>% sort
to_empty <- c("Invalid data or other (missing data)","Blank","Don't know")
to_empty <- tibble(name = annotation$name[annotation$description%in%to_empty],
                   description = annotation$description[annotation$description%in%to_empty], 
                   code = annotation$code[annotation$description%in%to_empty]) %>% filter(code!=".")
sapply(1:nrow(to_empty), function(x) {
  name <- to_empty$name[[x]]
  rows <- which(food[,name]==to_empty$code[[x]])
  food[rows,name] <- NA
})
colnames(food) <- extract_data_dict(colnames(food),"./input/fhs/phs000007.v33.pht002350.v7.vr_ffreq_ex08_1_0615s.data_dict.xml")
colnames(food) <- gsub(" USDA 2003","",colnames(food))
add <- ifelse(grepl("FFQ: ",colnames(food)), "(Q)",ifelse(grepl("Derived field: Nutrient value: |Derived field: Nutrient value - |DERIVED FIELD:",colnames(food)), "(D)", ""))
colnames(food) <- gsub("FFQ: |Derived field: Nutrient value: |Derived field: Nutrient value - |DERIVED FIELD:","",colnames(food))
colnames(food) <- str_to_sentence(colnames(food))
colnames(food) <- paste(colnames(food),add)
food <- food[,apply(food,2,function(x)length(unique(x%>%na.omit))) %>% enframe %>% filter(value>1) %>% pull(name)]
food <- food[,colMeans(is.na(food))<0.5]
food_q <- food[,grepl("\\(Q\\)",colnames(food))]
food_q$`Multivitamin brand (Q)` <- NULL
food_q <- reshape2::melt(food_q %>% as.matrix) %>% mutate(sample = as.character(Var1), feature = as.character(Var2)) %>% dplyr::select(sample,feature,value)
food_q <- food_q %>% left_join(predictions[,c("sample","ic_acc")]) %>% na.omit
food_q <- food_q %>% group_by(feature) %>% summarise(r = cor.test(value,ic_acc, method = "s")$estimate, p = cor.test(value,ic_acc, method = "s")$p.value)
food_q <- food_q %>% mutate(fdr = p.adjust(p))
food_q$feature <- gsub(" \\(Q\\)","",food_q$feature)

food_d <- food[,grepl("\\(D\\)",colnames(food))]
food_d <- reshape2::melt(food_d %>% as.matrix) %>% mutate(sample = as.character(Var1), feature = as.character(Var2)) %>% dplyr::select(sample,feature,value)
food_d <- food_d %>% left_join(predictions[,c("sample","ic_acc")]) %>% na.omit
food_d <- food_d %>% group_by(feature) %>% summarise(r = cor.test(value,ic_acc, method = "s")$estimate, p = cor.test(value,ic_acc, method = "s")$p.value)
food_d <- food_d %>% mutate(fdr = p.adjust(p))
food_d$feature <- gsub(" \\(D\\)","",food_d$feature)

p4 <- ggplot(food_d, aes(r, -log10(p)))+
  geom_point(fill = "grey", shape = 21, color = "grey")+
  geom_point(data = food_d %>% filter(p<0.05&r>0), fill = npg_colors[8], shape = 21)+
  geom_point(data = food_d %>% filter(p<0.05&r<0), fill = npg_colors[4], shape = 21)+
  geom_point(data = food_d %>% filter(fdr<0.05&r>0), fill = npg_colors[8], shape = 23, size = 3)+
  geom_point(data = food_d %>% filter(fdr<0.05&r<0), fill = npg_colors[4], shape = 23, size = 3)+
  geom_text_repel(data = food_d %>% filter(p<0.05&r>0) %>% arrange(p) %>% slice_head(n = 10), aes(label = feature), max.overlaps = 10, segment.size = 0.1, size = 3, xlim = c(1, 0.1), direction = "y", force = 100)+
  geom_text_repel(data = food_d %>% filter(p<0.05&r<0) %>% arrange(p) %>% slice_head(n = 10), aes(label = feature), max.overlaps = 10, segment.size = 0.1, hjust = 1, size = 3, xlim = c(-1, -0.12), direction = "y")+
  coord_cartesian(expand = T, clip = "off")+
  theme_pubr(border = TRUE)+
  ggtitle(label = "Food frequency (Derived)")+
  theme(plot.margin = unit(c(1, 5, 1, 5), "lines"))+
  geom_hline(yintercept = 0, lty = 2)+
  geom_vline(xintercept = 0, lty = 2)+
  labs(x = "Pearson correlation coefficient (R)", y = "-log10 (P-value)")

#diet guidelines
diet <- read_tsv("./input/fhs/phs000007.v33.pht007774.v3.p14.c1.vr_dgai2010_ex08_1_1009s.HMB-IRB-MDS.txt.gz", comment = "#")
diet <- diet %>% filter(idtype==1)
diet <- diet %>% column_to_rownames(var = "shareid")
diet <- diet[,3:30,drop=FALSE]
colnames(diet) <- extract_data_dict(colnames(diet),"./input/fhs/phs000007.v33.pht007774.v3.vr_dgai2010_ex08_1_1009s.data_dict.xml")
colnames(diet) <- gsub(" COMPONENT SCORE","",colnames(diet))
colnames(diet) <- gsub("^ES |^HC ","",colnames(diet))
colnames(diet) <- str_to_sentence(colnames(diet))
diet <- diet[,colMeans(is.na(diet))<0.5]
diet <- reshape2::melt(diet %>% as.matrix) %>% mutate(sample = as.character(Var1), feature = as.character(Var2)) %>% dplyr::select(sample,feature,value)
diet <- diet %>% left_join(predictions[,c("sample","ic_acc")]) %>% na.omit
diet <- diet %>% group_by(feature) %>% summarise(r = cor.test(value,ic_acc, method = "s")$estimate, p = cor.test(value,ic_acc, method = "s")$p.value)
diet$fdr <- p.adjust(diet$p)

combined <- rbind(food_q %>% mutate(group="Food frequency (Questionnaire)"), flavonoids %>% mutate(group="Flavonoids (Derived intake)"))
combined$group <- factor(combined$group, levels = c("Food frequency (Questionnaire)","Flavonoids (Derived intake)"))
p1 <- ggplot(combined, aes(r, -log10(p)))+
  geom_point(fill = "grey", shape = 21, color = "grey")+
  geom_point(data = combined %>% filter(p<0.05&r>0), fill = npg_colors[8], shape = 21)+
  geom_point(data = combined %>% filter(p<0.05&r<0), fill = npg_colors[4], shape = 21)+
  geom_point(data = combined %>% filter(fdr<0.05&r>0), fill = npg_colors[8], shape = 23, size = 3)+
  geom_point(data = combined %>% filter(fdr<0.05&r<0), fill = npg_colors[4], shape = 23, size = 3)+
  geom_text_repel(data = combined %>% filter(group=="Food frequency (Questionnaire)"&p<0.05&r>0) %>% arrange(p) %>% slice_head(n = 10), aes(label = feature), max.overlaps = 100, segment.size = 0.1, size = 3, xlim = c(1, 0.17), direction = "y", force = 200)+
  geom_text_repel(data = combined %>% filter(group=="Food frequency (Questionnaire)"&p<0.05&r<0) %>% arrange(p) %>% slice_head(n = 10), aes(label = feature), max.overlaps = 100, segment.size = 0.1, hjust = 1, size = 3, xlim = c(-1, -0.17), direction = "y")+
  geom_text_repel(data = combined %>% filter(group=="Flavonoids (Derived intake)"&p<0.05&r>0) %>% arrange(p) %>% slice_head(n = 10), aes(label = feature), max.overlaps = 100, segment.size = 0.1, size = 3, xlim = c(1, 0.17), direction = "y", force = 200)+
  geom_text_repel(data = combined %>% filter(group=="Flavonoids (Derived intake)"&p<0.05&r<0) %>% arrange(p) %>% slice_head(n = 10), aes(label = feature), max.overlaps = 100, segment.size = 0.1, hjust = 1, size = 3, xlim = c(-1, -0.17), direction = "y")+
  coord_cartesian(expand = T, clip = "off")+
  theme_pubr(border = TRUE)+
  facet_wrap(.~group, ncol = 1)+
  xlim(-0.15,0.15)+
  ylim(0,15)+
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+
  geom_hline(yintercept = 0, lty = 2)+
  geom_vline(xintercept = 0, lty = 2)+
  coord_axes_inside(labels_inside = TRUE, clip = "off")+
  labs(x = "", y = "")+
  theme(plot.margin = unit(c(0, 2, -0.3, 2), "cm"))

combined2 <- rbind(fatty_acids %>% mutate(group="Fatty acids in red blood cells"), diet %>% mutate(group="Dietary guidelines adherence index"))
combined21 <- combined2 %>% filter(group=="Fatty acids in red blood cells")
p21 <- ggplot(combined21, aes(r, -log10(p)))+
  geom_point(fill = "grey", shape = 21, color = "grey")+
  geom_point(data = combined21 %>% filter(p<0.05&r>0), fill = npg_colors[8], shape = 21)+
  geom_point(data = combined21 %>% filter(p<0.05&r<0), fill = npg_colors[4], shape = 21)+
  geom_point(data = combined21 %>% filter(fdr<0.05&r>0), fill = npg_colors[8], shape = 23, size = 3)+
  geom_point(data = combined21 %>% filter(fdr<0.05&r<0), fill = npg_colors[4], shape = 23, size = 3)+
  geom_text_repel(data = combined21 %>% filter(group=="Fatty acids in red blood cells"&p<0.05&r>0) %>% arrange(p) %>% slice_head(n = 10), aes(label = feature), max.overlaps = 100, segment.size = 0.1, size = 3, xlim = c(1, 0.17), direction = "y", force = 200)+
  geom_text_repel(data = combined21 %>% filter(group=="Fatty acids in red blood cells"&p<0.05&r<0) %>% arrange(p) %>% slice_head(n = 10), aes(label = feature), max.overlaps = 100, segment.size = 0.1, hjust = 1, size = 3, xlim = c(-1, -0.17), direction = "y")+
  coord_cartesian(expand = T, clip = "off")+
  theme_pubr(border = TRUE)+
  facet_wrap(.~group, ncol = 1)+
  xlim(-0.15,0.15)+
  ylim(0,15)+
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+
  geom_hline(yintercept = 0, lty = 2)+
  geom_vline(xintercept = 0, lty = 2)+
  coord_axes_inside(labels_inside = TRUE, clip = "off")+
  labs(x = "", y = "")+
  theme(plot.margin = unit(c(0, 2, -0.3, 2), "cm"))

combined22 <- combined2 %>% filter(group=="Dietary guidelines adherence index")
p22 <- ggplot(combined22, aes(r, -log10(p)))+
  geom_point(fill = "grey", shape = 21, color = "grey")+
  geom_point(data = combined22 %>% filter(p<0.05&r>0), fill = npg_colors[8], shape = 21)+
  geom_point(data = combined22 %>% filter(p<0.05&r<0), fill = npg_colors[4], shape = 21)+
  geom_point(data = combined22 %>% filter(fdr<0.05&r>0), fill = npg_colors[8], shape = 23, size = 3)+
  geom_point(data = combined22 %>% filter(fdr<0.05&r<0), fill = npg_colors[4], shape = 23, size = 3)+
  geom_text_repel(data = combined22 %>% filter(group=="Dietary guidelines adherence index"&p<0.05&r>0) %>% arrange(p) %>% slice_head(n = 10), aes(label = feature), max.overlaps = 100, segment.size = 0.1, size = 3, xlim = c(1, 0.17), direction = "y", force = 200)+
  geom_text_repel(data = combined22 %>% filter(group=="Dietary guidelines adherence index"&p<0.05&r<0) %>% arrange(p) %>% slice_head(n = 10), aes(label = feature), max.overlaps = 100, segment.size = 0.1, hjust = 1, size = 3, xlim = c(-1, -0.17), direction = "y")+
  coord_cartesian(expand = T, clip = "off")+
  theme_pubr(border = TRUE)+
  facet_wrap(.~group, ncol = 1)+
  xlim(-0.15,0.15)+
  ylim(0,15)+
  geom_hline(yintercept = 0, lty = 2)+
  geom_vline(xintercept = 0, lty = 2)+
  coord_axes_inside(labels_inside = TRUE, clip = "off")+
  labs(x = "Spearman correlation coefficient (R)", y = "")+
  theme(plot.margin = unit(c(0, 2, -0.3, 2), "cm"))#, axis.title.x = element_text(margin = margin(t = 20)), axis.text.x = element_text(margin = margin(t = 10)))

e <- grid.arrange(p1,p21,p22, ncol = 1, heights = c(2,1,1.2))

pdf(file = "./output/figures/figure4.pdf",width=20, height=20)
plot_grid(plot_grid(
  plot_grid(a, plot_grid(NULL,b, nrow = 1, rel_widths = c(0.01,1)), nrow = 2, labels = c("a","b"), rel_heights = c(0.8,1)),
  plot_grid(c, ncol = 1, labels = c("c")),
  ncol = 2,
  rel_widths = c(2,1),
  label_size = 12,
  align = "v"
), plot_grid(d,NULL,e,NULL, rel_widths = c(2,0.05,1,0.05), ncol = 4, labels = c("d","e")), 
nrow = 2, rel_heights = c(0.65,1.35))
dev.off()
