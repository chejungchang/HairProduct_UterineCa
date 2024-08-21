

library(dplyr) # data management 
library(survival) # coxph
library(survminer) # coxph
library(gtools) # function/macro
library(Hmisc) 
library(mets) # for estimate
library(rms) # for splines
library(ggcorrplot) # plot 
library(corrplot)# plot
library(lmtest)


# setwd("U:\\P0_UE_hairproduct\\Rscript\\Upload to FTP\\")
dat <- read_csv("hairUEca.csv")


# Table 5 ----------------
## Race ------------------
dat3 = dat[dat$race!=2,]

product = c("PermDye","SemiPermDye","Straightener","HairPerm2")
product_yn <- c("PermDye_yn", "SemiPermDye_yn", "Straightener_yn", "HairPerm_yn")

cox_race <- 
  defmacro(event, EOFage, exposure, exposure_yn, data, expr = {
    n_yn <- data %>% 
      filter_at(c('age','race2','educ',
                  'BMI', 'baseline_meno', 'parity', 'BCduration', 'HRT2', 'smokingcat2', 'alcohol', 
                  'menarchecat', 'PA'), all_vars(!is.na(.)))%>%
      filter(!is.na(event) & !is.na(EOFage) & !is.na(exposure_yn))%>%
      group_by(race, exposure_yn)%>%
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
      group_by(race, exposure)%>%
      summarise(py = as.integer(sum(py, na.rm=TRUE)), 
                n = n(), 
                case_n = length(event[event==1]))%>%
      mutate(n_per = round(n/sum(n)*100,1), 
             case_n_per = round(case_n/sum(case_n)*100,1))
    
    ad_mod_yn <- coxph(Surv(age, EOFage, event) ~ factor(exposure_yn)*factor(race) + 
                         factor(race) + factor(educ) + rcs(BMI, quantile(BMI, c(0, .5, .35, .65, .95, 1), na.rm = TRUE)) + factor(baseline_meno) + 
                         factor(parity) + factor(BCduration) + factor(HRT2) + 
                         factor(smokingcat2) + factor(alcohol) + 
                         menarchecat + PA + cluster(HH_PSID), ties='breslow', data = data)
    
    covnum <- dim(summary(ad_mod_yn)$coefficients)[1]
    
    e1 <- estimate(ad_mod_yn, 1)
    HR_e1 <- exp(e1$coefmat[1])
    lCI_e1 <- exp(e1$coefmat[3])
    uCI_e1 <- exp(e1$coefmat[4])
    p_e1 <- e1$coefmat[5]
    
    e2 <- estimate(ad_mod_yn, function(p) c(p[1]+p[covnum]))
    HR_e2 <- exp(e2$coefmat[1])
    lCI_e2 <- exp(e2$coefmat[3])
    uCI_e2 <- exp(e2$coefmat[4])
    p_e2 <- e2$coefmat[5]
    
    p_inter <- estimate(ad_mod_yn, covnum)$coefmat[5]
    
    ad_mod_freq <- coxph(Surv(age, EOFage, event) ~ factor(exposure) * factor(race) +  
                           factor(race) + factor(educ) + rcs(BMI, quantile(BMI, c(0, .5, .35, .65, .95, 1), na.rm = TRUE)) + factor(baseline_meno) + 
                           factor(parity) + factor(BCduration) + factor(HRT2) + 
                           factor(smokingcat2) + factor(alcohol) + 
                           menarchecat + PA + cluster(HH_PSID), ties='breslow', data = data)
    
    covnum_freq <- dim(summary(ad_mod_freq)$coefficients)[1]
    
    e1_freq <- estimate(ad_mod_freq, 2)
    HR_e1_freq <- exp(e1_freq$coefmat[1])
    lCI_e1_freq <- exp(e1_freq$coefmat[3])
    uCI_e1_freq <- exp(e1_freq$coefmat[4])
    p_e1_freq <- e1_freq$coefmat[5]
    
    e2_freq <- estimate(ad_mod_freq, function(p) c(p[2]+p[covnum_freq]))
    HR_e2_freq <- exp(e2_freq$coefmat[1])
    lCI_e2_freq <- exp(e2_freq$coefmat[3])
    uCI_e2_freq <- exp(e2_freq$coefmat[4])
    p_e2_freq <- e2_freq$coefmat[5]
    
    p_inter_freq <- estimate(ad_mod_freq, covnum_freq)$coefmat[5]
   
    final <- data.frame(
      modifier = c(0,1,0,1),
      freq = c("ever","ever", "frequent","frequent"),
      rbind(n_yn[2,3:7], n_yn[4,3:7], n[3,3:7], n[6,3:7]),
      HR = c(round(HR_e1,2), round(HR_e2,2), round(HR_e1_freq,2),round(HR_e2_freq,2)),
      lCI = c(round(lCI_e1,2), round(lCI_e2,2), round(lCI_e1_freq,2),round(lCI_e2_freq,2)),
      uCI = c(round(uCI_e1,2), round(uCI_e2,2), round(uCI_e1_freq,2),round(uCI_e2_freq,2)),
      CI = c(paste0("(", round(lCI_e1,2), ", ",  round(uCI_e1,2),")"), 
             paste0("(", round(lCI_e2,2), ", ",  round(uCI_e2,2),")"),
             paste0("(", round(lCI_e1_freq,2), ", ",  round(uCI_e1_freq,2),")"), 
             paste0("(", round(lCI_e2_freq,2), ", ",  round(uCI_e2_freq,2),")")),
      
      p_inter = c(round(p_inter, 4), "", round(p_inter_freq, 4), ""))
      
    return(final)
  })

racemat <- c()
for(i in 1:length(product_yn)){
  res <- cox_race(event = UECa_event, EOFage = eof_age, 
                    exposure_yn = get(product_yn[i]), exposure = get(product[i]), data = dat3)
  res <- cbind(res, product[i])
  colnames(res)[13] <- "product"
  racemat <- rbind(racemat, res)
  
}

racemat <- apply(racemat, 2, as.character)
race <- data.frame(racemat)
race

### p for trend ----------------
p_trend <- 
  defmacro(event, EOFage, exposure, data, expr = {
    ad_mod <- coxph(Surv(age, EOFage, event) ~ exposure*factor(race) + 
                      factor(race) + factor(educ) + rcs(BMI, quantile(BMI, c(0, .5, .35, .65, .95, 1), na.rm = TRUE)) + factor(baseline_meno) + 
                      factor(parity) + factor(BCduration) + factor(HRT2) + 
                      factor(smokingcat2) + factor(alcohol) + 
                      menarchecat + PA + cluster(HH_PSID), ties='breslow', data = data)
    
    covnum <- dim(summary(ad_mod)$coefficients)[1]
    
    e1 <- estimate(ad_mod, 1)
    p_e1 <- e1$coefmat[5]
    
    e2 <- estimate(ad_mod, function(p) c(p[1]+p[covnum]))
    p_e2 <- e2$coefmat[5]
    
    output = cbind(p_e1,p_e2)
    return(output)
  })

p_trend(exposure = Straightener, EOFage = eof_age, event = UECa_event, data=dat3)


## BMI -------------------
dat4 <- dat %>%
  mutate(BMI_cut1 = case_when(BMIcat %in% c(0,1) ~ 0, 
                              BMIcat == 2 ~ 1))

cox_BMI <- 
  defmacro(event, EOFage, exposure, exposure_yn, data, expr = {
    n_yn <- data %>% 
      filter_at(c('age','race2','educ',
                  'BMI_cut1', 'baseline_meno', 'parity', 'BCduration', 'HRT2', 'smokingcat2', 'alcohol', 
                  'menarchecat', 'PA'), all_vars(!is.na(.)))%>%
      filter(!is.na(event) & !is.na(EOFage) & !is.na(exposure_yn))%>%
      group_by(BMI_cut1, exposure_yn)%>%
      summarise(py = as.integer(sum(py, na.rm=TRUE)), 
                n = n(), 
                case_n = length(event[event==1]))%>%
      mutate(n_per = round(n/sum(n)*100,1), 
             case_n_per = round(case_n/sum(case_n)*100,1))
    
    n <- data %>% 
      filter_at(c('age','race2','educ',
                  'BMI_cut1', 'baseline_meno', 'parity', 'BCduration', 'HRT2', 'smokingcat2', 'alcohol', 
                  'menarchecat', 'PA'), all_vars(!is.na(.)))%>%
      filter(!is.na(event) & !is.na(EOFage) & !is.na(exposure))%>%
      group_by(BMI_cut1, exposure)%>%
      summarise(py = as.integer(sum(py, na.rm=TRUE)), 
                n = n(), 
                case_n = length(event[event==1]))%>%
      mutate(n_per = round(n/sum(n)*100,1), 
             case_n_per = round(case_n/sum(case_n)*100,1))
    
     
    ad_mod_yn <- coxph(Surv(age, EOFage, event) ~ factor(exposure_yn)*factor(BMI_cut1) + 
                         factor(race2) + factor(educ) +  factor(baseline_meno) + 
                         factor(parity) + factor(BCduration) + factor(HRT2) + 
                         factor(smokingcat2) + factor(alcohol) + 
                         menarchecat + PA + cluster(HH_PSID), ties='breslow', data = data)
    
    covnum <- dim(summary(ad_mod_yn)$coefficients)[1]
    
    e1 <- estimate(ad_mod_yn, 1)
    HR_e1 <- exp(e1$coefmat[1])
    lCI_e1 <- exp(e1$coefmat[3])
    uCI_e1 <- exp(e1$coefmat[4])
    p_e1 <- e1$coefmat[5]
    
    e2 <- estimate(ad_mod_yn, function(p) c(p[1]+p[covnum]))
    HR_e2 <- exp(e2$coefmat[1])
    lCI_e2 <- exp(e2$coefmat[3])
    uCI_e2 <- exp(e2$coefmat[4])
    p_e2 <- e2$coefmat[5]
    
    p_inter <- estimate(ad_mod_yn, covnum)$coefmat[5]
    
    ad_mod_freq <- coxph(Surv(age, EOFage, event) ~ factor(exposure) * factor(BMI_cut1) +  
                           factor(race2) + factor(educ) +  factor(baseline_meno) + 
                           factor(parity) + factor(BCduration) + factor(HRT2) + 
                           factor(smokingcat2) + factor(alcohol) + 
                           menarchecat + PA + cluster(HH_PSID), ties='breslow', data = data)
    
    covnum_freq <- dim(summary(ad_mod_freq)$coefficients)[1]
    
    e1_freq <- estimate(ad_mod_freq, 2)
    HR_e1_freq <- exp(e1_freq$coefmat[1])
    lCI_e1_freq <- exp(e1_freq$coefmat[3])
    uCI_e1_freq <- exp(e1_freq$coefmat[4])
    p_e1_freq <- e1_freq$coefmat[5]
    
    e2_freq <- estimate(ad_mod_freq, function(p) c(p[2]+p[covnum_freq]))
    HR_e2_freq <- exp(e2_freq$coefmat[1])
    lCI_e2_freq <- exp(e2_freq$coefmat[3])
    uCI_e2_freq <- exp(e2_freq$coefmat[4])
    p_e2_freq <- e2_freq$coefmat[5]
    
    p_inter_freq <- estimate(ad_mod_freq, covnum_freq)$coefmat[5]
    
    final <- data.frame(
      modifier = c(0,1,0,1),
      freq = c("ever","ever", "frequent","frequent"),
      rbind(n_yn[2,3:7], n_yn[4,3:7], n[3,3:7], n[6,3:7]),
      HR = c(round(HR_e1,2), round(HR_e2,2), round(HR_e1_freq,2),round(HR_e2_freq,2)),
      lCI = c(round(lCI_e1,2), round(lCI_e2,2), round(lCI_e1_freq,2),round(lCI_e2_freq,2)),
      uCI = c(round(uCI_e1,2), round(uCI_e2,2), round(uCI_e1_freq,2),round(uCI_e2_freq,2)),
      CI = c(paste0("(", round(lCI_e1,2), ", ",  round(uCI_e1,2),")"), 
             paste0("(", round(lCI_e2,2), ", ",  round(uCI_e2,2),")"),
             paste0("(", round(lCI_e1_freq,2), ", ",  round(uCI_e1_freq,2),")"), 
             paste0("(", round(lCI_e2_freq,2), ", ",  round(uCI_e2_freq,2),")")),
      
      p_inter = c(round(p_inter, 4), "", round(p_inter_freq, 4), ""))
    
    return(final)
  })

BMImat <- c()
for(i in 1:length(product_yn)){
  res <- cox_BMI(event = UECa_event, EOFage = eof_age, 
                  exposure_yn = get(product_yn[i]), exposure = get(product[i]), data = dat4)
  res <- cbind(res, product[i])
  colnames(res)[13] <- "product"
  BMImat <- rbind(BMImat, res)
  
}

BMImat <- apply(BMImat, 2, as.character)
BMI <- data.frame(BMImat)
BMI


### p for trend ----------------
p_trend <- 
  defmacro(event, EOFage, exposure, data, expr = {
    ad_mod <- coxph(Surv(age, EOFage, event) ~ exposure* factor(BMI_cut1) + 
                      factor(race2) + factor(educ) +  
                      factor(parity) + factor(BCduration) + factor(HRT2) + 
                      factor(smokingcat2) + factor(alcohol) + 
                      menarchecat + PA + cluster(HH_PSID), ties='breslow', data = data)
    
    covnum <- dim(summary(ad_mod)$coefficients)[1]
    
    e1 <- estimate(ad_mod, 1)
    p_e1 <- e1$coefmat[5]
    
    e2 <- estimate(ad_mod, function(p) c(p[1]+p[covnum]))
    p_e2 <- e2$coefmat[5]
    
    output = cbind(p_e1,p_e2)
    return(output)
  })

p_trend(exposure = Straightener, EOFage = eof_age, event = UECa_event, data=dat4)



## Physical activity ------------------

cox_PA <- 
  defmacro(event, EOFage, exposure, exposure_yn, data, expr = {
    n_yn <- data %>% 
      filter_at(c('age','race2','educ',
                  'BMI', 'baseline_meno', 'parity', 'BCduration', 'HRT2', 'smokingcat2', 'alcohol', 
                  'menarchecat', 'PA'), all_vars(!is.na(.)))%>%
      filter(!is.na(event) & !is.na(EOFage) & !is.na(exposure_yn))%>%
      group_by(PAcat, exposure_yn)%>%
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
      group_by(PAcat, exposure)%>%
      summarise(py = as.integer(sum(py, na.rm=TRUE)), 
                n = n(), 
                case_n = length(event[event==1]))%>%
      mutate(n_per = round(n/sum(n)*100,1), 
             case_n_per = round(case_n/sum(case_n)*100,1))
    
    # yn 
    ad_mod_yn <- coxph(Surv(age, EOFage, event) ~ factor(exposure_yn)*factor(PAcat) + 
                         factor(race2) + factor(educ) + rcs(BMI, quantile(BMI, c(0, .5, .35, .65, .95, 1), na.rm = TRUE)) + factor(baseline_meno) + 
                         factor(parity) + factor(BCduration) + factor(HRT2) + 
                         factor(smokingcat2) + factor(alcohol) + 
                         menarchecat + cluster(HH_PSID), ties='breslow', data = data)
    
    covnum <- dim(summary(ad_mod_yn)$coefficients)[1]
    
    e1 <- estimate(ad_mod_yn, 1)
    HR_e1 <- exp(e1$coefmat[1])
    lCI_e1 <- exp(e1$coefmat[3])
    uCI_e1 <- exp(e1$coefmat[4])
    p_e1 <- e1$coefmat[5]
    
    e2 <- estimate(ad_mod_yn, function(p) c(p[1]+p[covnum]))
    HR_e2 <- exp(e2$coefmat[1])
    lCI_e2 <- exp(e2$coefmat[3])
    uCI_e2 <- exp(e2$coefmat[4])
    p_e2 <- e2$coefmat[5]
    
    p_inter <- estimate(ad_mod_yn, covnum)$coefmat[5]
    
    ad_mod_freq <- coxph(Surv(age, EOFage, event) ~ factor(exposure) * factor(PAcat) +  
                           factor(race2) + factor(educ) + rcs(BMI, quantile(BMI, c(0, .5, .35, .65, .95, 1), na.rm = TRUE)) + factor(baseline_meno) + 
                           factor(parity) + factor(BCduration) + factor(HRT2) + 
                           factor(smokingcat2) + factor(alcohol) + 
                           menarchecat + cluster(HH_PSID), ties='breslow', data = data)
    
    covnum_freq <- dim(summary(ad_mod_freq)$coefficients)[1]
    
    e1_freq <- estimate(ad_mod_freq, 2)
    HR_e1_freq <- exp(e1_freq$coefmat[1])
    lCI_e1_freq <- exp(e1_freq$coefmat[3])
    uCI_e1_freq <- exp(e1_freq$coefmat[4])
    p_e1_freq <- e1_freq$coefmat[5]
    
    e2_freq <- estimate(ad_mod_freq, function(p) c(p[2]+p[covnum_freq]))
    HR_e2_freq <- exp(e2_freq$coefmat[1])
    lCI_e2_freq <- exp(e2_freq$coefmat[3])
    uCI_e2_freq <- exp(e2_freq$coefmat[4])
    p_e2_freq <- e2_freq$coefmat[5]
    
    p_inter_freq <- estimate(ad_mod_freq, covnum_freq)$coefmat[5]
    
    final <- data.frame(
      modifier = c(0,1,0,1),
      freq = c("ever","ever", "frequent","frequent"),
      rbind(n_yn[2,3:7], n_yn[4,3:7], n[3,3:7], n[6,3:7]),
      HR = c(round(HR_e1,2), round(HR_e2,2), round(HR_e1_freq,2),round(HR_e2_freq,2)),
      lCI = c(round(lCI_e1,2), round(lCI_e2,2), round(lCI_e1_freq,2),round(lCI_e2_freq,2)),
      uCI = c(round(uCI_e1,2), round(uCI_e2,2), round(uCI_e1_freq,2),round(uCI_e2_freq,2)),
      CI = c(paste0("(", round(lCI_e1,2), ", ",  round(uCI_e1,2),")"), 
             paste0("(", round(lCI_e2,2), ", ",  round(uCI_e2,2),")"),
             paste0("(", round(lCI_e1_freq,2), ", ",  round(uCI_e1_freq,2),")"), 
             paste0("(", round(lCI_e2_freq,2), ", ",  round(uCI_e2_freq,2),")")),
      
      p_inter = c(round(p_inter, 4), "", round(p_inter_freq, 4), ""))
    
    return(final)
  })

PAmat <- c()
for(i in 1:length(product_yn)){
  res <- cox_PA(event = UECa_event, EOFage = eof_age, 
                 exposure_yn = get(product_yn[i]), exposure = get(product[i]), data = dat)
  res <- cbind(res, product[i])
  colnames(res)[13] <- "product"
  PAmat <- rbind(PAmat, res)
  
}

PAmat <- apply(PAmat, 2, as.character)
PA <- data.frame(PAmat)

### p for trend ----------------
p_trend <- 
  defmacro(event, EOFage, exposure, data, expr = {
    ad_mod <- coxph(Surv(age, EOFage, event) ~ exposure*factor(PAcat) + 
                      factor(race2) + factor(educ) + rcs(BMI, quantile(BMI, c(0, .5, .35, .65, .95, 1), na.rm = TRUE)) + factor(baseline_meno) + 
                      factor(parity) + factor(BCduration) + factor(HRT2) + 
                      factor(smokingcat2) + factor(alcohol) + 
                      menarchecat + cluster(HH_PSID), ties='breslow', data = data)
    
    covnum <- dim(summary(ad_mod)$coefficients)[1]
    
    e1 <- estimate(ad_mod, 1)
    p_e1 <- e1$coefmat[5]
    
    e2 <- estimate(ad_mod, function(p) c(p[1]+p[covnum]))
    p_e2 <- e2$coefmat[5]
    
    output = cbind(p_e1,p_e2)
    return(output)
  })

res = p_trend(exposure = Straightener, EOFage = eof_age, event = UECa_event, data=dat)
round(res,3)
