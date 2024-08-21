

library(dplyr) # data management 
library(survival) # coxph
library(survminer) # coxph
library(gtools) # function/macro
library(Hmisc) 
library(mets) # for estimate
library(rms) # for splines
library(ggcorrplot) # plot 
library(corrplot) # plot
library(lmtest)
library(tidyverse)


# setwd("U:\\P0_UE_hairproduct\\Rscript\\Upload to FTP\\")
dat <- read_csv("hairUEca.csv")


# Table 2 ------------------
cov <- c('age','race2','educ', 'BMI', 'baseline_meno', 'parity', 'BCduration', 
         'HRT2', 'smokingcat2', 'alcohol', 'menarchecat', 'PA',"UECa_event", "eof_age", "PermDye_yn")

cox_mod1 <- 
  defmacro(event, EOFage, exposure_yn, data, expr = {
    n_yn <- data %>% 
      filter_at(c('age','race2','educ',
                   'BMI', 'baseline_meno', 'parity', 'BCduration', 'HRT2', 'smokingcat2', 'alcohol', 
                    'menarchecat', 'PA'), all_vars(!is.na(.)))%>%
      filter(!is.na(event) & !is.na(EOFage) & !is.na(exposure_yn))%>%
      group_by(exposure_yn)%>%
      summarise(py = as.integer(sum(py, na.rm=TRUE)), 
                n = n(), 
                case_n = length(event[event==1]))%>%
      mutate(n_per = round(n/sum(n)*100,1), 
             case_n_per = round(case_n/sum(case_n)*100,1))
  
    # Crude
    datX <- data
    datX <- datX %>%
      filter_at(c('age','race2','educ',
                  'BMI', 'baseline_meno', 'parity', 'BCduration', 'HRT2', 'smokingcat2', 'alcohol', 
                  'menarchecat', 'PA'), all_vars(!is.na(.)))%>%
      filter(!is.na(event) & !is.na(EOFage) & !is.na(exposure_yn))
                  
    mod_yn <- coxph(Surv(age, EOFage, event) ~ factor(exposure_yn) + cluster(HH_PSID), ties='breslow', data = datX)
    
    HR0 <- exp(coef(mod_yn)[[1]])
    lCI0 <- summary(mod_yn)$conf.int[1,3] 
    uCI0 <- summary(mod_yn)$conf.int[1,4]
    
    # Adjusted 
    ad_mod_yn <- coxph(Surv(age, EOFage, event) ~ factor(exposure_yn) + 
                         factor(race2) + factor(educ) + rcs(BMI, quantile(BMI, c(0, .5, .35, .65, .95, 1), na.rm = TRUE)) + factor(baseline_meno) + 
                         factor(parity) + factor(BCduration) + factor(HRT2) + 
                         factor(smokingcat2) + factor(alcohol) + 
                         menarchecat + PA + cluster(HH_PSID), ties='breslow', data = datX)
    
    # ties='breslow', to match default SAS options
    ad_HR0 <- exp(coef(ad_mod_yn)[[1]])
    ad_lCI0 <- summary(ad_mod_yn)$conf.int[1,3] 
    ad_uCI0 <- summary(ad_mod_yn)$conf.int[1,4]
    
    final <- data.frame(
      n_yn[1:2,2:6], 
      
      cHR = c("Ref", round(HR0,2)), 
      clCI = c("Ref", round(lCI0,2)),
      cuCI = c("Ref", round(uCI0,2)),
      cCI = c("Ref", paste0("(", round(lCI0,2), ", ",  round(uCI0,2),")")),
      
      adjHR = c("Ref", round(ad_HR0,2)), 
      adjlCI = c("Ref", round(ad_lCI0,2)),
      adjuCI = c("Ref", round(ad_uCI0,2)),
      adjCI = c("Ref", paste0("(", round(ad_lCI0,2), ", ",  round(ad_uCI0,2),")")))
    return(final)
  })


product_yn <- c("PermDye_yn", "PermDye_other_yn", "SemiPermDye_yn", "SemiPermDye_other_yn", 
                "HairColorRinse_yn", "Bleach_yn", "Highlight_yn", "Straightener_yn", "Straightener_other_yn", "HairPerm_yn")
m1 <- c()

for(i in 1:length(product_yn)){
  res <- cox_mod1(event = UECa_event, EOFage = eof_age, 
                  exposure_yn = get(product_yn[i]), data = dat)
  res <- cbind(res, rep(product_yn[i],2), c("Never", "Ever"))
  colnames(res)[14] <- "product"
  colnames(res)[15] <- "freq"
  m1 <- rbind(m1, res)
}

m1 <- apply(m1, 2, as.character)

table2 <- data.frame(m1)


### Table 3 ###
product = c("PermDye","SemiPermDye","Straightener","HairPerm","HairPerm2")

cox_mod2 <- 
  defmacro(event, EOFage, exposure, data, expr = {
    n <- data %>% 
      filter_at(c('age','race2','educ',
                  'BMI', 'baseline_meno', 'parity', 'BCduration', 'HRT2', 'smokingcat2', 'alcohol', 
                  'menarchecat', 'PA'), all_vars(!is.na(.)))%>%
      filter(!is.na(event) & !is.na(EOFage) & !is.na(exposure))%>%
      group_by(exposure)%>%
      summarise(py = as.integer(sum(py, na.rm=TRUE)), 
                n = n(), 
                case_n = length(event[event==1]))%>%
      mutate(n_per = round(n/sum(n)*100,1), 
             case_n_per = round(case_n/sum(case_n)*100,1))
    
    # Crude
    datX <- data
    datX <- datX %>%
      filter_at(c('age','race2','educ',
                  'BMI', 'baseline_meno', 'parity', 'BCduration', 'HRT2', 'smokingcat2', 'alcohol', 
                  'menarchecat', 'PA'), all_vars(!is.na(.)))%>%
      filter(!is.na(event) & !is.na(EOFage) & !is.na(exposure))
    
    mod <- coxph(Surv(age, EOFage, event) ~ factor(exposure) + cluster(HH_PSID), ties='breslow', data = datX)
    
    HR1 <- exp(coef(mod)[[1]])
    lCI1 <- summary(mod)$conf.int[1,3] 
    uCI1 <- summary(mod)$conf.int[1,4]
    
    HR2 <- exp(coef(mod)[[2]])
    lCI2 <- summary(mod)$conf.int[2,3]  
    uCI2 <- summary(mod)$conf.int[2,4]
    
    # Adjusted 
    ad_mod <- coxph(Surv(age, EOFage, event) ~ factor(exposure) + 
                      factor(race2) + factor(educ) + rcs(BMI, quantile(BMI, c(0, .5, .35, .65, .95, 1), na.rm = TRUE)) + factor(baseline_meno) + 
                      factor(parity) + factor(BCduration) + factor(HRT2) + 
                      factor(smokingcat2) + factor(alcohol) + 
                      menarchecat + PA + cluster(HH_PSID), ties='breslow', data = data)
    
    ad_HR1 <- exp(coef(ad_mod)[[1]])
    ad_lCI1 <- summary(ad_mod)$conf.int[1,3] 
    ad_uCI1 <- summary(ad_mod)$conf.int[1,4]
    
    ad_HR2 <- exp(coef(ad_mod)[[2]])
    ad_lCI2 <- summary(ad_mod)$conf.int[2,3]  
    ad_uCI2 <- summary(ad_mod)$conf.int[2,4]
    
    final <- data.frame(
      n[1:3,2:6], 
      
      cHR = c("Ref", round(HR1,2), round(HR2,2)), 
      clCI = c("Ref", round(lCI1,2), round(lCI2,2)),
      cuCI = c("Ref", round(uCI1,2), round(uCI2,2)),
      cCI = c("Ref", paste0("(", round(lCI1,2), ", ",  round(uCI1,2),")"), 
              paste0("(", round(lCI2,2), ", ",  round(uCI2,2),")")),
      
      adjHR = c("Ref", round(ad_HR1,2), round(ad_HR2,2)), 
      adjlCI = c("Ref", round(ad_lCI1,2), round(ad_lCI2,2)),
      adjuCI = c("Ref",  round(ad_uCI1,2), round(ad_uCI2,2)),
      adjCI = c("Ref", paste0("(", round(ad_lCI1,2), ", ",  round(ad_uCI1,2),")"), 
                paste0("(", round(ad_lCI2,2), ", ",  round(ad_uCI2,2),")")))
    return(final)
  })

# run loop for multiple exposures
m1 <- c()
for(i in 1:length(product)){
  res <- cox_mod2(event = UECa_event, EOFage = eof_age, 
                  exposure = get(product[i]), data = dat)
  res <- cbind(res, rep(product[i],3), c("Never", "Less frequent", "More frequent"))
  colnames(res)[14] <- "product"
  colnames(res)[15] <- "freq"
  m1 <- rbind(m1, res)
}
m1 <- apply(m1, 2, as.character)

table3 <- data.frame(m1)

# Table 2 ------------------
## Color -------------------
# Make it consistently with table 1, 
# assuming that the frequency question of none user (In the past 12 months, how frequently have you used XXX?) is correct.

dat2 = dat %>%
  mutate(PermDye_dark = replace(PermDye_dark, PermDye_yn == 0, 0))%>%
  mutate(PermDye_light = replace(PermDye_light, PermDye_yn == 0, 0))%>%
  mutate(SemiPermDye_dark = replace(SemiPermDye_dark, SemiPermDye_yn == 0, 0))%>%
  mutate(SemiPermDye_light = replace(SemiPermDye_light, SemiPermDye_yn == 0, 0))


cox_mod3 <- 
  defmacro(event, EOFage, exposure1, exposure2, data, expr = {
    n1 <- data %>% 
      filter_at(c('age','race2','educ',
                  'BMI', 'baseline_meno', 'parity', 'BCduration', 'HRT2', 'smokingcat2', 'alcohol', 
                  'menarchecat', 'PA'), all_vars(!is.na(.)))%>%
      filter(!is.na(event) & !is.na(EOFage) & !is.na(exposure1) & !is.na(exposure2))%>%
      summarise(py = as.integer(sum(py[exposure1==1], na.rm=TRUE)), 
                n = length(exposure1[exposure1==1]), 
                case_n = length(exposure1[event==1 & exposure1==1]))
    
    n2 <- data %>% 
      filter_at(c('age','race2','educ',
                  'BMI', 'baseline_meno', 'parity', 'BCduration', 'HRT2', 'smokingcat2', 'alcohol', 
                  'menarchecat', 'PA'), all_vars(!is.na(.)))%>%
      filter(!is.na(event) & !is.na(EOFage) & !is.na(exposure1) & !is.na(exposure2))%>%
      summarise(py = as.integer(sum(py[exposure2==1], na.rm=TRUE)), 
                n = length(exposure2[exposure2==1]), 
                case_n = length(exposure2[event==1 & exposure2==1]))

    # Crude
    datX <- data
    datX <- datX %>%
      filter_at(c('age','race2','educ',
                  'BMI', 'baseline_meno', 'parity', 'BCduration', 'HRT2', 'smokingcat2', 'alcohol', 
                  'menarchecat', 'PA'), all_vars(!is.na(.)))%>%
      filter(!is.na(event) & !is.na(EOFage) & !is.na(exposure1)&!is.na(exposure2))
    
    mod <- coxph(Surv(age, EOFage, event) ~ factor(exposure1) + factor(exposure2) + cluster(HH_PSID), ties='breslow', data = datX)
    
    HR1 <- exp(coef(mod)[[1]])
    lCI1 <- summary(mod)$conf.int[1,3] 
    uCI1 <- summary(mod)$conf.int[1,4]
    
    HR2 <- exp(coef(mod)[[2]])
    lCI2 <- summary(mod)$conf.int[2,3]  
    uCI2 <- summary(mod)$conf.int[2,4]
    
    # Adjusted 
    ad_mod <- coxph(Surv(age, EOFage, event) ~ factor(exposure1) + factor(exposure2) + 
                      factor(race2) + factor(educ) + rcs(BMI, quantile(BMI, c(0, .5, .35, .65, .95, 1), na.rm = TRUE)) + factor(baseline_meno) + 
                      factor(parity) + factor(BCduration) + factor(HRT2) + 
                      factor(smokingcat2) + factor(alcohol) + 
                      menarchecat + PA + cluster(HH_PSID), ties='breslow', data = data)
    
    ad_HR1 <- exp(coef(ad_mod)[[1]])
    ad_lCI1 <- summary(ad_mod)$conf.int[1,3] 
    ad_uCI1 <- summary(ad_mod)$conf.int[1,4]
    
    ad_HR2 <- exp(coef(ad_mod)[[2]])
    ad_lCI2 <- summary(ad_mod)$conf.int[2,3]  
    ad_uCI2 <- summary(ad_mod)$conf.int[2,4]  
    
    final <- data.frame(
      rbind(n1[1,1:3], n2[1,1:3]),
      
      cHR = c(round(HR1,2), round(HR2,2)), 
      clCI = c(round(lCI1,2), round(lCI2,2)),
      cuCI = c(round(uCI1,2), round(uCI2,2)),
      cCI = c(paste0("(", round(lCI1,2), ", ",  round(uCI1,2),")"), 
              paste0("(", round(lCI2,2), ", ",  round(uCI2,2),")")), 
      adHR = c(round(ad_HR1,2), round(ad_HR2,2)), 
      adlCI = c(round(ad_lCI1,2), round(ad_lCI2,2)),
      aduCI = c(round(ad_uCI1,2), round(ad_uCI2,2)),
      adCI = c(paste0("(", round(ad_lCI1,2), ", ",  round(ad_uCI1,2),")"),
               paste0("(", round(ad_lCI2,2), ", ",  round(ad_uCI2,2),")"))) 
    
    return(final)
  })

mod1 <- cox_mod3(event = UECa_event, EOFage = eof_age, 
                 exposure1 = PermDye_light, exposure2 = PermDye_dark, data = dat2)

mod2 <- cox_mod3(event = UECa_event, EOFage = eof_age, 
                 exposure1 = SemiPermDye_light, exposure2 = SemiPermDye_dark, data = dat2)


## Duration -------------------
product = c("PermDye_duration", # Permanent hair dye (how many years)
            "SemiPermDye_duration") # Semi - permanent hair dye (how many years)

cox_mod4 <- 
  defmacro(event, EOFage, exposure, data, expr = {
    n <- data %>% 
      filter_at(c('age','race2','educ',
                  'BMI', 'baseline_meno', 'parity', 'BCduration', 'HRT2', 'smokingcat2', 'alcohol', 
                  'menarchecat', 'PA'), all_vars(!is.na(.)))%>%
      filter(!is.na(event) & !is.na(EOFage) & !is.na(exposure))%>%
      group_by(exposure)%>%
      summarise(py = as.integer(sum(py, na.rm=TRUE)), 
                n = n(), 
                case_n = length(event[event==1]))%>%
      mutate(n_per = round(n/sum(n)*100,1), 
             case_n_per = round(case_n/sum(case_n)*100,1))
    
    # Crude
    datX <- data
    datX <- datX %>%
      filter_at(c('age','race2','educ',
                  'BMI', 'baseline_meno', 'parity', 'BCduration', 'HRT2', 'smokingcat2', 'alcohol', 
                  'menarchecat', 'PA'), all_vars(!is.na(.)))%>%
      filter(!is.na(event) & !is.na(EOFage) & !is.na(exposure))
    
    mod <- coxph(Surv(age, EOFage, event) ~ factor(exposure) + cluster(HH_PSID), ties='breslow', data = datX)
    
    HR1 <- exp(coef(mod)[[1]])
    lCI1 <- summary(mod)$conf.int[1,3] 
    uCI1 <- summary(mod)$conf.int[1,4]
    
    HR2 <- exp(coef(mod)[[2]])
    lCI2 <- summary(mod)$conf.int[2,3]  
    uCI2 <- summary(mod)$conf.int[2,4]
    
    HR3 <- exp(coef(mod)[[3]])
    lCI3 <- summary(mod)$conf.int[3,3]  
    uCI3 <- summary(mod)$conf.int[3,4]
    
    # Adjusted 
    ad_mod <- coxph(Surv(age, EOFage, event) ~ factor(exposure) + 
                      factor(race2) + factor(educ) + rcs(BMI, quantile(BMI, c(0, .5, .35, .65, .95, 1), na.rm = TRUE)) + factor(baseline_meno) + 
                      factor(parity) + factor(BCduration) + factor(HRT2) + 
                      factor(smokingcat2) + factor(alcohol) + 
                      menarchecat + PA + cluster(HH_PSID), ties='breslow', data = data)
    
    ad_HR1 <- exp(coef(ad_mod)[[1]])
    ad_lCI1 <- summary(ad_mod)$conf.int[1,3] 
    ad_uCI1 <- summary(ad_mod)$conf.int[1,4]
    
    ad_HR2 <- exp(coef(ad_mod)[[2]])
    ad_lCI2 <- summary(ad_mod)$conf.int[2,3]  
    ad_uCI2 <- summary(ad_mod)$conf.int[2,4]   
    
    ad_HR3 <- exp(coef(ad_mod)[[3]])
    ad_lCI3 <- summary(ad_mod)$conf.int[3,3]  
    ad_uCI3 <- summary(ad_mod)$conf.int[3,4]
    
    final <- data.frame(
      n[1:4,2:6], 
      
      cHR = c("Ref", round(HR1,2), round(HR2,2), round(HR3,2)), 
      clCI = c("Ref", round(lCI1,2), round(lCI2,2), round(lCI3,2)), 
      ulCI = c("Ref", round(uCI1,2), round(uCI2,2), round(uCI3,2)), 
      cCI = c("Ref", 
              paste0("(", round(lCI1,2), ", ",  round(uCI1,2),")"), 
              paste0("(", round(lCI2,2), ", ",  round(uCI2,2),")"),
              paste0("(", round(lCI3,2), ", ",  round(uCI3,2),")")), 
      adHR = c("Ref", round(ad_HR1,2), round(ad_HR2,2), round(ad_HR3, 2)), 
      adlCI = c("Ref", round(ad_lCI1,2), round(ad_lCI2,2), round(ad_lCI3,2)), 
      adlCI = c("Ref", round(ad_uCI1,2), round(ad_uCI2,2), round(ad_uCI3,2)), 
      adCI = c("Ref", 
               paste0("(", round(ad_lCI1,2), ", ",  round(ad_uCI1,2),")"),
               paste0("(", round(ad_lCI2,2), ", ",  round(ad_uCI2,2),")"),
               paste0("(", round(ad_lCI3,2), ", ",  round(ad_uCI3,2),")")))
    return(final)
  })


# run loop for multiple exposures
m2 <- c()

for(i in 1:length(product)){
  res <- cox_mod4(event = UECa_event, EOFage = eof_age, 
                  exposure = get(product[i]), data = dat)
  res <- cbind(res, rep(product[i],4), c("Never", "Less than 5 years",
                                         "5-9 years", 
                                         ">10 years"))
  colnames(res)[14] <- "product"
  colnames(res)[15] <- "freq"
  m2 <- rbind(m2, res)
}

m2 <- apply(m2, 2, as.character)
duration <- data.frame(m2)



## p for trend -------------------
product = c("PermDye","SemiPermDye","Straightener","HairPerm","HairPerm2","PermDye_duration", "SemiPermDye_duration")

p_trend <- 
  defmacro(event, EOFage, exposure, data, expr = {
    datX <- data
    datX <- datX %>%
      filter_at(c('age','race2','educ',
                  'BMI', 'baseline_meno', 'parity', 'BCduration', 'HRT2', 'smokingcat2', 'alcohol', 
                  'menarchecat', 'PA'), all_vars(!is.na(.)))%>%
      filter(!is.na(event) & !is.na(EOFage) & !is.na(exposure))
    mod <- coxph(Surv(age, EOFage, event) ~ exposure + cluster(HH_PSID), ties='breslow', data = datX)
    ad_mod <- coxph(Surv(age, EOFage, event) ~ exposure + 
                      factor(race2) + factor(educ) + rcs(BMI, quantile(BMI, c(0, .5, .35, .65, .95, 1), na.rm = TRUE)) + factor(baseline_meno) + 
                      factor(parity) + factor(BCduration) + factor(HRT2) + 
                      factor(smokingcat2) + factor(alcohol) + 
                      menarchecat + PA + cluster(HH_PSID), ties='breslow', data = data)
    p1 = round(summary(mod)$coef[1,6],3)
    p2 = round(summary(ad_mod)$coef[1,6],3)
    output = cbind(p1,p2)
    return(output)
  })

# run loop for multiple exposures 
pmat <- c()

for(i in 1:length(product)){
  res <- p_trend(event = UECa_event, EOFage = eof_age, 
                 exposure = get(product[i]), data = dat)
  res <- cbind(res, product[i])
  colnames(res)[1] <- "p-crude"
  colnames(res)[2] <- "p-adj"
  colnames(res)[3] <- "product"
  pmat <- rbind(pmat, res)
}

pmat <- apply(pmat, 2, as.character)
pmat

# Table 4 --------
product = c("PermDye","SemiPermDye","Straightener","HairPerm2")
product_yn <- c("PermDye_yn", "SemiPermDye_yn", "Straightener_yn", "HairPerm_yn")

cox_subtype <- 
  defmacro(event, EOFage, exposure, exposure_yn, data, expr = {
    n_yn <- data %>% 
      filter_at(c('age','race2','educ',
                  'BMI', 'baseline_meno', 'parity', 'BCduration', 'HRT2', 'smokingcat2', 'alcohol', 
                  'menarchecat', 'PA'), all_vars(!is.na(.)))%>%
      filter(!is.na(event) & !is.na(EOFage) & !is.na(exposure_yn))%>%
      group_by(exposure_yn)%>%
      summarise(py = as.integer(sum(py, na.rm=TRUE)), 
                n = n(), 
                case_n = length(event[event==1]))%>%
      mutate(n_per = round(n/sum(n)*100,1), 
             case_n_per = round(case_n/sum(case_n)*100,1))
    
    n <- data %>% 
      filter_at(c('age','race2','educ',
                  'BMI', 'baseline_meno', 'parity', 'BCduration', 'HRT2', 'smokingcat2', 'alcohol', 
                  'menarchecat', 'PA'), all_vars(!is.na(.)))%>%
      filter(!is.na(event) & !is.na(EOFage) & !is.na(exposure))%>%
      group_by(exposure)%>%
      summarise(py = as.integer(sum(py, na.rm=TRUE)), 
                n = n(), 
                case_n = length(event[event==1]))%>%
      mutate(n_per = round(n/sum(n)*100,1), 
             case_n_per = round(case_n/sum(case_n)*100,1))
    
    # Adjusted 
    ad_mod_yn <- coxph(Surv(age, EOFage, event) ~ factor(exposure_yn) + 
                         factor(race2) + factor(educ) + rcs(BMI, quantile(BMI, c(0, .5, .35, .65, .95, 1), na.rm = TRUE)) + factor(baseline_meno) + 
                         factor(parity) + factor(BCduration) + factor(HRT2) + 
                         factor(smokingcat2) + factor(alcohol) + 
                         menarchecat + PA + cluster(HH_PSID), ties='breslow', data = data)
    
    ad_HR0 <- exp(coef(ad_mod_yn)[[1]])
    ad_lCI0 <- summary(ad_mod_yn)$conf.int[1,3] 
    ad_uCI0 <- summary(ad_mod_yn)$conf.int[1,4]
    
    ad_mod <- coxph(Surv(age, EOFage, event) ~ factor(exposure) + 
                      factor(race2) + factor(educ) + rcs(BMI, quantile(BMI, c(0, .5, .35, .65, .95, 1), na.rm = TRUE)) + factor(baseline_meno) + 
                      factor(parity) + factor(BCduration) + factor(HRT2) + 
                      factor(smokingcat2) + factor(alcohol) + 
                      menarchecat + PA + cluster(HH_PSID), ties='breslow', data = data)
    
    ad_HR1 <- exp(coef(ad_mod)[[1]])
    ad_lCI1 <- summary(ad_mod)$conf.int[1,3] 
    ad_uCI1 <- summary(ad_mod)$conf.int[1,4]
    
    ad_HR2 <- exp(coef(ad_mod)[[2]])
    ad_lCI2 <- summary(ad_mod)$conf.int[2,3]  
    ad_uCI2 <- summary(ad_mod)$conf.int[2,4]
    
    final <- data.frame(
      rbind(n_yn[1:2,2:6],n[2:3,2:6]), 
      
      adHR = c("Ref", round(ad_HR0,2), round(ad_HR1,2), round(ad_HR2,2)), 
      adlCI = c("Ref", round(ad_lCI0,2), round(ad_lCI1,2), round(ad_lCI2,2)), 
      aduCI = c("Ref", round(ad_uCI0,2), round(ad_uCI1,2), round(ad_uCI2,2)), 
      adCI = c("Ref", paste0("(", round(ad_lCI0,2), ", ",  round(ad_uCI0,2),")"), 
               paste0("(", round(ad_lCI1,2), ", ",  round(ad_uCI1,2),")"), 
               paste0("(", round(ad_lCI2,2), ", ",  round(ad_uCI2,2),")")))
    
    return(final)
  })

## medically confirmed -------------
m1 <- c()
for(i in 1:length(product)){
  res <- cox_subtype(event = UECa_event_MR, EOFage = eof_age, 
                     exposure = get(product[i]), exposure_yn=get(product_yn[i]), data = dat)
  res <- cbind(res, rep(product[i],4), c("Never", "Ever","Less frequent", "More frequent"))
  colnames(res)[10] <- "product"
  colnames(res)[11] <- "freq"
  m1 <- rbind(m1, res)
}

m1 <- apply(m1, 2, as.character)
table4_MR <- data.frame(m1)

## endometrial ---------------
m1 <- c()
for(i in 1:length(product)){
  res <- cox_subtype(event = UECa_event_endometrial, EOFage = eof_age, 
                     exposure = get(product[i]), exposure_yn=get(product_yn[i]), data = dat)
  res <- cbind(res, rep(product[i],4), c("Never", "Ever","Less frequent", "More frequent"))
  colnames(res)[10] <- "product"
  colnames(res)[11] <- "freq"
  m1 <- rbind(m1, res)
}
m1 <- apply(m1, 2, as.character)
table4_endo <- data.frame(m1)

## type 1 -------------------
m1 <- c()
for(i in 1:length(product)){
  res <- cox_subtype(event = UECa_event_type1, EOFage = eof_age, 
                     exposure = get(product[i]), exposure_yn=get(product_yn[i]), data = dat)
  res <- cbind(res, rep(product[i],4), c("Never", "Ever","Less frequent", "More frequent"))
  colnames(res)[10] <- "product"
  colnames(res)[11] <- "freq"
  m1 <- rbind(m1, res)
}
m1 <- apply(m1, 2, as.character)
table4_type1 <- data.frame(m1)

## type 2 ----------------
m1 <- c()
for(i in 1:length(product)){
  res <- cox_subtype(event = UECa_event_type2, EOFage = eof_age, 
                     exposure = get(product[i]), exposure_yn=get(product_yn[i]), data = dat)
  res <- cbind(res, rep(product[i],4), c("Never", "Ever","Less frequent", "More frequent"))
  colnames(res)[10] <- "product"
  colnames(res)[11] <- "freq"
  m1 <- rbind(m1, res)
}
m1 <- apply(m1, 2, as.character)
table4_type2 <- data.frame(m1)

## menopausal status -------------
dat_meno = dat %>%
  mutate(postmeno_time = case_when(baseline_meno == 0 & eof_meno == 0 ~ 0, 
                                   baseline_meno == 0 & eof_meno == 1 ~ 1, 
                                   baseline_meno == 1 ~ 1))  
pre = dat_meno %>%
  mutate(tstart = case_when(baseline_meno == 0 & eof_meno == 0 ~ age, 
                            baseline_meno == 0 & eof_meno == 1 ~ age, 
                            baseline_meno == 1 ~ NA_real_)) %>% 
  mutate(tstop = case_when(baseline_meno == 0 & eof_meno == 0 ~ eof_age, 
                           baseline_meno == 0 & eof_meno == 1  ~ menoage,
                           baseline_meno == 1 ~ NA_real_)) %>% 
  mutate(event3 = ifelse(UECa_event == 1 & postmeno_time == 0, 1, 0)) %>%
  mutate(prepost = 0)

post = dat_meno %>%
  mutate(tstart = case_when(baseline_meno == 0 & eof_meno == 0 ~ NA_real_, 
                            baseline_meno == 0 & eof_meno == 1 ~ menoage, 
                            baseline_meno == 1 ~ age)) %>% 
  mutate(tstop = case_when(baseline_meno == 0 & eof_meno == 0 ~ NA_real_, 
                           baseline_meno == 0 & eof_meno == 1  ~ eof_age,
                           baseline_meno == 1 ~ eof_age)) %>% 
  mutate(event3 = ifelse(UECa_event == 1 & postmeno_time == 1, 1, 0)) %>%
  mutate(prepost = 1)

newdat = rbind(pre,post)
newdat = newdat[!is.na(newdat$tstart),]

cox_meno <- 
  defmacro(event, baselineage, EOFage, exposure, exposure_yn, data, expr = {
    n_yn_pre <- data %>% 
      filter_at(c('age','race2','educ',
                  'BMI', 'parity', 'BCduration', 'HRT2', 'smokingcat2', 'alcohol', 
                  'menarchecat', 'PA'), all_vars(!is.na(.)))%>%
      filter(!is.na(event) & !is.na(EOFage) & !is.na(exposure_yn) & !is.na(baselineage) & prepost==0)%>%
      group_by(exposure_yn)%>%
      summarise(py = as.integer(sum(py, na.rm=TRUE)), 
                n = n(), 
                case_n = length(event[event==1]))%>%
      mutate(per = case_n/sum(case_n)*100)
    
    n_yn_post <- data %>% 
      filter_at(c('age','race2','educ',
                  'BMI', 'parity', 'BCduration', 'HRT2', 'smokingcat2', 'alcohol', 
                  'menarchecat', 'PA'), all_vars(!is.na(.)))%>%
      filter(!is.na(event) & !is.na(EOFage) & !is.na(exposure_yn) & !is.na(baselineage) & prepost==1)%>%
      group_by(exposure_yn)%>%
      summarise(py = as.integer(sum(py, na.rm=TRUE)), 
                n = n(), 
                case_n = length(event[event==1]))%>%
      mutate(per = case_n/sum(case_n)*100)
    
    n_pre <- data %>% 
      filter_at(c('age','race2','educ',
                  'BMI', 'parity', 'BCduration', 'HRT2', 'smokingcat2', 'alcohol', 
                  'menarchecat', 'PA'), all_vars(!is.na(.)))%>%
      filter(!is.na(event) & !is.na(EOFage) & !is.na(exposure) & !is.na(baselineage) & prepost==0)%>%
      group_by(exposure)%>%
      summarise(py = as.integer(sum(py, na.rm=TRUE)), 
                n = n(), 
                case_n = length(event[event==1]))%>%
      mutate(per = case_n/sum(case_n)*100)
    
    n_post <- data %>% 
      filter_at(c('age','race2','educ',
                  'BMI', 'parity', 'BCduration', 'HRT2', 'smokingcat2', 'alcohol', 
                  'menarchecat', 'PA'), all_vars(!is.na(.)))%>%
      filter(!is.na(event) & !is.na(EOFage) & !is.na(exposure) & !is.na(baselineage) & prepost==1)%>%
      group_by(exposure)%>%
      summarise(py = as.integer(sum(py, na.rm=TRUE)), 
                n = n(), 
                case_n = length(event[event==1]))%>%
      mutate(per = case_n/sum(case_n)*100)
    
    # Pre_yn
    pre_mod <- coxph(Surv(baselineage, EOFage, event) ~ factor(exposure_yn) + 
                       factor(race2) + factor(educ) + rcs(BMI, quantile(BMI, c(0, .5, .35, .65, .95, 1), na.rm = TRUE)) + 
                       factor(parity) + factor(BCduration) + factor(HRT2) + 
                       factor(smokingcat2) + factor(alcohol) + 
                       menarchecat + PA + cluster(HH_PSID), ties='breslow', data = data, subset = {prepost==0})
    
    pre_HR <- exp(coef(pre_mod)[[1]])
    pre_lCI <- summary(pre_mod)$conf.int[1,3] 
    pre_uCI <- summary(pre_mod)$conf.int[1,4]
    
    # Pre_freq
    pre_mod_freq <- coxph(Surv(baselineage, EOFage, event) ~ factor(exposure) + 
                            factor(race2) + factor(educ) + rcs(BMI, quantile(BMI, c(0, .5, .35, .65, .95, 1), na.rm = TRUE)) +  
                            factor(parity) + factor(BCduration) + factor(HRT2) + 
                            factor(smokingcat2) + factor(alcohol) + 
                            menarchecat + PA + cluster(HH_PSID), ties='breslow', data = data, subset = {prepost==0})
    
    pre_HR_freq <- exp(coef(pre_mod_freq)[[2]])
    pre_lCI_freq <- summary(pre_mod_freq)$conf.int[2,3] 
    pre_uCI_freq <- summary(pre_mod_freq)$conf.int[2,4]
    
    # Post 
    post_mod <- coxph(Surv(baselineage, EOFage, event) ~ factor(exposure_yn) + 
                        factor(race2) + factor(educ) + rcs(BMI, quantile(BMI, c(0, .5, .35, .65, .95, 1), na.rm = TRUE)) + 
                        factor(parity) + factor(BCduration) + factor(HRT2) + 
                        factor(smokingcat2) + factor(alcohol) + 
                        menarchecat + PA + cluster(HH_PSID), ties='breslow', data = data, subset = {prepost==1})
    
    post_HR <- exp(coef(post_mod)[[1]])
    post_lCI <- summary(post_mod)$conf.int[1,3] 
    post_uCI <- summary(post_mod)$conf.int[1,4]
    
    # Post_freq
    post_mod_freq <- coxph(Surv(baselineage, EOFage, event) ~ factor(exposure) + 
                             factor(race2) + factor(educ) + rcs(BMI, quantile(BMI, c(0, .5, .35, .65, .95, 1), na.rm = TRUE)) + 
                             factor(parity) + factor(BCduration) + factor(HRT2) + 
                             factor(smokingcat2) + factor(alcohol) + 
                             menarchecat + PA + cluster(HH_PSID), ties='breslow', data = data, subset = {prepost==1})
    
    post_HR_freq <- exp(coef(post_mod_freq)[[2]])
    post_lCI_freq <- summary(post_mod_freq)$conf.int[2,3] 
    post_uCI_freq <- summary(post_mod_freq)$conf.int[2,4]
    
    final <- data.frame(
      rbind(n_yn_pre[2,2:5], n_pre[3,2:5], n_yn_post[2,2:5], n_post[3,2:5]), 
      
      HR = c(round(pre_HR,2), round(pre_HR_freq,2), round(post_HR,2), round(post_HR_freq,2)), 
      lCI = c(round(pre_lCI,2), round(pre_lCI_freq,2), round(post_lCI,2), round(post_lCI_freq,2)),
      uCI =  c(round(pre_uCI,2), round(pre_uCI_freq,2), round(post_uCI,2), round(post_uCI_freq,2)),
      CI = c(paste0("(", round(pre_lCI,2), ", ",  round(pre_uCI,2),")"), 
             paste0("(", round(pre_lCI_freq,2), ", ",  round(pre_uCI_freq,2),")"), 
             paste0("(", round(post_lCI,2), ", ",  round(post_uCI,2),")"),
             paste0("(", round(post_lCI_freq,2), ", ",  round(post_uCI_freq,2),")")),
      prepost = c("pre","pre_freq","post","post_freq"))
    
    return(final)
  })


# run loop for multiple exposures 
prepostmat <- c()
for(i in 1:length(product_yn)){
  res <- cox_meno(event = event3, baselineage = tstart, 
                  EOFage = tstop, exposure_yn = get(product_yn[i]), 
                  exposure = get(product[i]),data = newdat)
  res <- cbind(res, product[i])
  colnames(res)[10] <- "product"
  prepostmat <- rbind(prepostmat, res)
}

prepostmat <- apply(prepostmat, 2, as.character)
meno <- data.frame(prepostmat)


# #p for trend --------------------
pmat <- c()
for(i in 1:length(product)){
  res <- p_trend(event = UECa_event_MR, EOFage = eof_age, 
                 exposure = get(product[i]), data = dat)
  res <- cbind(res, product[i])
  colnames(res)[1] <- "p-crude"
  colnames(res)[2] <- "p-adj"
  colnames(res)[3] <- "product"
  pmat <- rbind(pmat, res)
}

pmat <- apply(pmat, 2, as.character)
pmat

pmat <- c()
for(i in 1:length(product)){
  res <- p_trend(event = UECa_event_endometrial, EOFage = eof_age, 
                 exposure = get(product[i]), data = dat)
  res <- cbind(res, product[i])
  colnames(res)[1] <- "p-crude"
  colnames(res)[2] <- "p-adj"
  colnames(res)[3] <- "product"
  pmat <- rbind(pmat, res)
}

pmat <- apply(pmat, 2, as.character)
pmat


pmat <- c()
for(i in 1:length(product)){
  res <- p_trend(event = UECa_event_type1, EOFage = eof_age, 
                 exposure = get(product[i]), data = dat)
  res <- cbind(res, product[i])
  colnames(res)[1] <- "p-crude"
  colnames(res)[2] <- "p-adj"
  colnames(res)[3] <- "product"
  pmat <- rbind(pmat, res)
}

pmat <- apply(pmat, 2, as.character)
pmat

pmat <- c()
for(i in 1:length(product)){
  res <- p_trend(event = UECa_event_type2, EOFage = eof_age, 
                 exposure = get(product[i]), data = dat)
  res <- cbind(res, product[i])
  colnames(res)[1] <- "p-crude"
  colnames(res)[2] <- "p-adj"
  colnames(res)[3] <- "product"
  pmat <- rbind(pmat, res)
}

pmat <- apply(pmat, 2, as.character)
pmat


p_trend2 <- 
  defmacro(baselineage, EOFage, event, exposure, data, expr = {
    post_mod_freq <- coxph(Surv(tstart, tstop, event3) ~ exposure + 
                             factor(race2) + factor(educ) + rcs(BMI, quantile(BMI, c(0, .5, .35, .65, .95, 1), na.rm = TRUE)) + 
                             factor(parity) + factor(BCduration) + factor(HRT2) + 
                             factor(smokingcat2) + factor(alcohol) + 
                             menarchecat + PA + cluster(HH_PSID), ties='breslow', data = newdat, subset = {prepost==1})
    pv = round(summary(post_mod_freq)$coef[,"Pr(>|z|)"][1],3)
    return(pv)
  })

postmat <- c()
for(i in 1:length(product)){
  res <- p_trend2(event = event3, baselineage = tstart, 
                  EOFage = tstop, exposure = get(product[i]), data = newdat)
  res <- cbind(res, product[i])
  colnames(res)[2] <- "product"
  postmat <- rbind(postmat, res)
}
postmat

p_trend3 <- 
  defmacro(baselineage, EOFage, event, exposure, data, expr = {
    post_mod_freq <- coxph(Surv(tstart, tstop, event3) ~ exposure + 
                             factor(race2) + factor(educ) + rcs(BMI, quantile(BMI, c(0, .5, .35, .65, .95, 1), na.rm = TRUE)) + 
                             factor(parity) + factor(BCduration) + factor(HRT2) + 
                             factor(smokingcat2) + factor(alcohol) + 
                             menarchecat + PA + cluster(HH_PSID), ties='breslow', data = newdat, subset = {prepost==0})
    pv = round(summary(post_mod_freq)$coef[,"Pr(>|z|)"][1],2)
    return(pv)
  })

premat <- c()
for(i in 1:length(product)){
  res <- p_trend3(event = event3, baselineage = tstart, 
                  EOFage = tstop, exposure = get(product[i]), data = newdat)
  res <- cbind(res, product[i])
  colnames(res)[2] <- "product"
  premat <- rbind(premat, res)
}
premat


