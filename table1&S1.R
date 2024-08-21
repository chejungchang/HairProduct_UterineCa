

library(tidyverse) # data management 
library(table1) # table 1 
library(expss)

# setwd("U:\\P0_UE_hairproduct\\Rscript\\Upload to FTP\\")
dat <- read_csv("hairUEca.csv")

# Table 1 and S1 -----------------

dat.table1 <- dat %>%
   select(UECa_event, py, menarche, menarchecat, menoage, age, race2, educ, BMI, BMIcat, PA,  PAcat, baseline_meno, parity, smoking, alcohol, BCduration, HRT2, Straightener_yn, hairdresserjob)


label(dat.table1$age) <- "Age (year)"
label(dat.table1$py) <- "Follow-up time (year)"
label(dat.table1$menarche) <- "Age at menarche (year)"
label(dat.table1$menoage) <- "Age at menopause (year)"

dat.table1$race2 <- factor(dat.table1$race2, levels = c(0,1,2,3), 
                          labels = c("Non-Hispanic White", "Black", "Hispanic", "Other"))
label(dat.table1$race2) <- "Race/ethnicity"

dat.table1$educ <- factor(dat.table1$educ, level = c(0,1,2), 
                          labels = c("High school or less", "Some college", "College or above"))
label(dat.table1$educ) <- "Education"

dat.table1$BMIcat <- factor(dat.table1$BMIcat, levels = c(0,1,2), 
                                labels = c("Underweight and Normal", "Overweight", "Obesity"))
label(dat.table1$BMI) <- "Body mass index (BMI)"

label(dat.table1$PA) <- "Total average MET (metabolic equivalent)-hours/week"

dat.table1$PAcat <- factor(dat.table1$PAcat, levels = c(0,1), 
                              labels = c("<33rd", ">=33rd"))
label(dat.table1$PAcat) <- "Physical Activity"

dat.table1$baseline_meno <- factor(dat.table1$baseline_meno, levels = c(0,1), 
                                   labels = c("Premenopausal", "Postmenopausal"))
label(dat.table1$baseline_meno) <-"Menopausal status" 

dat.table1$parity <- factor(dat.table1$parity, levels = c(0,1,2,3), 
                            labels = c("Nulliparous", "1", "2", ">=3"))
label(dat.table1$parity) <-"Parity"

dat.table1$smoking <- factor(dat.table1$smoking, levels = c(0,1,2), 
                                     labels = c("Never", "Past", "Current"))
label(dat.table1$smoking) <-"Smoking status"

dat.table1$alcohol <- factor(dat.table1$alcohol, levels = c(0,1,2), 
                             labels = c("Never or past", "Current <1 drink/day", "Current >=1 drinks/day"))
label(dat.table1$alcohol) <-"Alcohol consumption" 

dat.table1$BCduration <- factor(dat.table1$BCduration, levels = c(0,1,2,3), 
                                labels = c("None", "<2 years", "2-<10 years", ">=10 years"))
label(dat.table1$BCduration) <-"Oral contraception use" 

dat.table1$HRT2 <- factor(dat.table1$HRT2, levels = c(0,1,2), 
                         labels = c("None", "Estrogen alone", "Estrogen plus Progestin"))
label(dat.table1$HRT2) <-"Hormone replacement therapy use" 

dat.table1$menarchecat <- factor(dat.table1$menarchecat, levels = c(1,0), 
                                   labels = c("<13", ">=13 years old"))
label(dat.table1$menarchecat) <- "Age at menarche" 


dat.table1$hairdresserjob <- factor(dat.table1$hairdresserjob, levels = c(1,0), 
                                   labels = c("Yes", "No"))
label(dat.table1$hairdresserjob) <- "Hairdressing job" 


tableone = table1(~.|factor(UECa_event), data = dat.table1)
tableone

dat.table1v2 = dat.table1[!is.na(dat.table1$Straightener_yn),]
tableSone = table1(~.|factor(Straightener_yn), data = dat.table1v2)
tableSone

