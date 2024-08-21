
# Library 
library(tidyverse) # data cleaning 
library(haven) # import sas data

# Data import
setwd("U:\\P0_UE_hairproduct\\raw_data\\")
data1 <- read_sas('dr00287_04_01.sas7bdat')

#####  Check missing data #####
# Duplicates 
sum(duplicated(data1$PSID)) # n = 0

# Exclusion: pull missing data information from SAS file
data2 = data1 %>%
  filter(!is.na(FU_UECa_Event)) # n = 50884 --> 50432 (delete 452)
# withdrew 3, pre-baseline diagnosis 380, Uncertain diagnosis 10, uncertain timing of diagnosis 59

# Impute missing EOF age and Hyst age 
data3 = data2 %>%
  mutate(FU_UECa_EOFAgeExact_imp = FU_UECa_EOFAgeExact)%>%
  mutate(FU_UECa_EOFAgeExact_imp = case_when(
    is.na(FU_UECa_EOFAgeExact_imp)&is.na(FU_UECa_EOFAge) ~ (FU_UECa_DxAgeExactMin+FU_UECa_DxAgeExactMax)/2,
    is.na(FU_UECa_EOFAgeExact_imp)&!is.na(FU_UECa_EOFAge) ~ FU_UECa_EOFAge,
    !is.na(FU_UECa_EOFAgeExact_imp) ~ FU_UECa_EOFAgeExact_imp))%>%
  mutate(HZ_HR_HystAgeExact_imp = ifelse(!is.na(HZ_HR_HystAgeExact), HZ_HR_HystAgeExact, HZ_HR_HystAge))


##### Remove those with hysterectomies before baseline ####
data4 = data3 %>%
  mutate(row_to_delet = case_when(HR34 == 1 ~ 1, 
                                  HZ_HR_Hyst == 1 & HZ_HR_HystAgeExact_imp<AgeExact_Baseline ~ 1, 
                                  TRUE ~ 0)) %>% 
  filter(row_to_delet == 0) 
# n = 50432 to 34847 --> 15585 hysterectomy
# for those have hyst but with unknown age, we keep them in the analysis. 

#### Fix EOF age and event for those who had hysterectomy during follow-up ####
data5 = data4 %>%
  mutate(FU_UECa_EOFAgeExact_imp_hyst = ifelse(HZ_HR_Hyst == 1 & 
                                                 HZ_HR_HystAgeExact_imp < FU_UECa_EOFAgeExact_imp & 
                                                 !is.na(HZ_HR_HystAgeExact_imp),  
                                               HZ_HR_HystAgeExact_imp, FU_UECa_EOFAgeExact_imp)) %>%
  mutate(FU_UECa_Event_hyst = ifelse(HZ_HR_Hyst == 1 & HZ_HR_HystAgeExact_imp < FU_UECa_EOFAgeExact_imp & 
                                       !is.na(HZ_HR_HystAgeExact_imp), 0, FU_UECa_Event))

##### Transform exposure variables and remove missing exposure variables #####
# Adjust Vanguard variables 
data5 = data5 %>%
  mutate(P129 = case_when(P128==0 ~ 1, P128==1 & is.na(P129) ~ 2, P128==1 & !is.na(P129) ~ P129, is.na(P128) ~ P129, TRUE ~ P129), 
         P134 = case_when(P133==0 ~ 1, P133==1 & is.na(P134) ~ 2, P133==1 & !is.na(P134) ~ P134, is.na(P133) ~ P134, TRUE ~ P134), 
         P121 = case_when(P120==0 ~ 1, P120==1 & is.na(P121) ~ 2, P120==1 & !is.na(P121) ~ P121, is.na(P120) ~ P121, TRUE ~ P121), 
         P126 = case_when(P125==0 ~ 1, P125==1 & is.na(P126) ~ 2, P125==1 & !is.na(P126) ~ P126, is.na(P125) ~ P126, TRUE ~ P126), 
         P118 = case_when(P117==0 ~ 1, P117==1 & is.na(P118) ~ 2, P117==1 & !is.na(P118) ~ P118, is.na(P117) ~ P118, TRUE ~ P118), 
         P137 = case_when(P136==0 ~ 1, P136==1 & is.na(P137) ~ 2, P136==1 & !is.na(P137) ~ P137, is.na(P136) ~ P137, TRUE ~ P137), 
         P140 = case_when(P139==0 ~ 1, P139==1 & is.na(P140) ~ 2, P139==1 & !is.na(P140) ~ P140, is.na(P139) ~ P140, TRUE ~ P140), 
         P143 = case_when(P142==0 ~ 1, P142==1 & is.na(P143) ~ 2, P142==1 & !is.na(P143) ~ P143, is.na(P142) ~ P143, TRUE ~ P143), 
         P146 = case_when(P145==0 ~ 1, P145==1 & is.na(P146) ~ 2, P145==1 & !is.na(P146) ~ P146, is.na(P145) ~ P146, TRUE ~ P146), 
         P149 = case_when(P148==0 ~ 1, P148==1 & is.na(P149) ~ 2, P148==1 & !is.na(P149) ~ P149, is.na(P148) ~ P149, TRUE ~ P149))


mutate.var <- function(x)(ifelse(x == 1, 0, # did not use 
                                 ifelse(x %in% c(2,3), 1, # < 1-2 times a year; every 3-4 months; # less than once a month and 1-3 times per month
                                        ifelse(x %in% c(4,5,6),2, NA)))) # every 5-8 weeks; once a month; more than once a month # 1-5 times per weeks and above

# Product use
product = c("PC85","PC89","PC91","PC95","PC97",
            "PC99","PC101","PC103","PC105","PC107", # Main study
            
            "P129","P134","P121","P126",
            "P118","P137","P140","P143",
            "P146","P149") # Vanguard


# Keep hair perm and body wave before transformation --> for different frequency category in the analysis
data5$PC107_perm = data5$PC107
data5$P149_perm = data5$P149
data5$PC103_relaxer = data5$PC103
data5$P143_relaxer = data5$P143

# Transform exposure variables to never, less, and more frequent user
data6 = as_tibble(data5) %>%
  mutate_at(product, mutate.var) %>% 
  
  # Combine Vanguard and Main study
  mutate(mix_PC85 = ifelse(is.na(PC85), P129, PC85),
         mix_PC89 = ifelse(is.na(PC89), P134, PC89),
         mix_PC91 = ifelse(is.na(PC91), P121, PC91),
         mix_PC95 = ifelse(is.na(PC95), P126, PC95),
         mix_PC97 = ifelse(is.na(PC97), P118, PC97),
         mix_PC99 = ifelse(is.na(PC99), P137, PC99),
         mix_PC101 = ifelse(is.na(PC101), P140, PC101),
         mix_PC103 = ifelse(is.na(PC103), P143, PC103), 
         mix_PC105 = ifelse(is.na(PC105), P146, PC105),
         mix_PC107 = ifelse(is.na(PC107), P149, PC107)) 



# Select exposure variables in this study 
exposureTS = c("mix_PC85", "mix_PC89", "mix_PC91", "mix_PC95", "mix_PC97", 
               "mix_PC99", "mix_PC101", "mix_PC103", "mix_PC105", "mix_PC107",
               "PC87", "PC93", 
               "PC86dark", "PC86light", "PC92dark", "PC92light")


# Create never/ever exposure variables 
data7 = tibble(data6) %>%
  # Remove those who did not answer any hair product use questions
  filter_at(exposureTS, any_vars(!is.na(.)))
# n = 34847 to 34111 --> 736 

##### Remove pre-baseline diagnosis from age (or not contribute any follow-up time) #####
data8 = data7 %>%
  filter(FU_UECa_EOFAgeExact_imp_hyst > AgeExact_Baseline) 
# n = 34111 to 33947 --> 164
33947-34111

##### Other variables ######
data9 = data8 %>%
  mutate(educ = case_when(SE18 %in% c(1,2,3,4,5) ~ 0,  # HS or less 
                          SE18 %in% c(6,7) ~ 1,    # some college
                          SE18 %in% c(8,9,10) ~ 2)) %>%  # college or higher 
  
  mutate(parity = case_when(PG_MedParity == 0 ~ 0, # 0
                            PG_MedParity == 1 ~ 1, # 1
                            PG_MedParity == 2 ~ 2, # 2 
                            PG_MedParity >= 3 ~ 3)) %>% # more
  
  mutate(BMIcat = case_when(EX_BMI_CDC_final %in% c(1,2) ~ 0, 
                             EX_BMI_CDC_final == 3 ~ 1, 
                             EX_BMI_CDC_final %in% c(4,5,6) ~ 2),
         BMI = EX_BMI_final)%>% # obese cat 3
  
  mutate(alcohol = case_when(AL_DrinkCat6 %in% c(0,1) ~ 0, # never and past 
                             AL_DrinkCat6 == 2 ~ 1, # current <1 drink  
                             AL_DrinkCat6 %in% c(3,4,5) ~ 2))%>% # current >=1 drinks 
  
  mutate(# for models 
    race = case_when(SE_RACE_ETH == 0 ~ 0, # non-Hispanic white
                     SE_RACE15 %in% c(3,4,11,12) ~ 1, # all Black, 
                     !is.na(SE_RACE15) ~ 2), # other
    # for table 1
    race2 = case_when (SE_RACE_ETH == 0 ~ 0, # non-Hispanic white
                       SE_RACE15 %in% c(3,4,11,12) ~ 1, # all Black, 
                       SE_RACE15 %in% c(2,6,8,10,14) ~ 2,
                       !is.na(SE_RACE15) ~ 3)) %>%
  mutate(PA = PH_CurrentTotMETHrsPerWeek)%>%
    
  mutate(PAcat = case_when(PH_CurrentTotMETHrsPerWeek<quantile(data8$PH_CurrentTotMETHrsPerWeek, probs = c(.33), na.rm=TRUE)~0,
                               PH_CurrentTotMETHrsPerWeek<=quantile(data8$PH_CurrentTotMETHrsPerWeek, probs = c(1), na.rm=TRUE)~1))%>%
  
  mutate(HRT_type = paste0(HR_HRTEstrOnly_Ever, HR_HRTProgOnly_Ever, HR_HRTEstrProg_Ever),
         HRT2 = case_when(HRT_type=="000"~0,
                          HRT_type %in% c("001","010","011")~2,# Estrogen plus progestin 
                          HRT_type %in% c("100","101","111","110")~1)) %>% # Estrogen only 
  
  mutate(BCduration = case_when(HR_BCpill_Years == 0 | HR_BCpill_Ever == 0 ~ 0, # None
                                0 < HR_BCpill_Years & HR_BCpill_Years < 2 ~ 1, # <2
                                10 > HR_BCpill_Years & HR_BCpill_Years >= 2 ~ 2, # >=2-10
                                HR_BCpill_Years >= 10 ~ 3))%>% # >=10 
  
  mutate(smoking = SM_SmokeStatusN, 
         smokingcat2 = ifelse(SM_SmokeStatusN %in% c(1,2), 1, SM_SmokeStatusN))%>% # never, past or current
  
  mutate(menarchecat = case_when(PG_MenarcheAge>0 & PG_MenarcheAge<13 ~ 1, 
                              PG_MenarcheAge>=13 ~ 0),
         menarche = PG_MenarcheAge) %>%
  
  mutate(py = FU_UECa_EOFAgeExact_imp_hyst - AgeExact_Baseline, 
         age = AgeExact_Baseline,
         eof_age = FU_UECa_EOFAgeExact_imp_hyst) %>%
  mutate(hairdresserjob = ifelse(!is.na(HM1)|OC48a==1, 1, 0),)


data9.1 = data9 %>%
    mutate(menoage = HZ_HR_MenopauseAgeExact)%>%
    mutate(menoage = case_when(
      is.na(menoage) & !is.na(HZ_HR_MenopauseAge) & HZ_HR_MenopauseAge != as.integer(AgeExact_Baseline) ~ HZ_HR_MenopauseAge, 
      is.na(menoage) & !is.na(HZ_HR_MenopauseAge) & HZ_HR_MenopauseAge == as.integer(AgeExact_Baseline) & HR_MenopauseStatus %in% c(6,7,11) ~ HZ_HR_MenopauseAge+0.9,
      is.na(menoage) & !is.na(HZ_HR_MenopauseAge) & HZ_HR_MenopauseAge == as.integer(AgeExact_Baseline) & (HR_MenopauseStatus %in% c(1,3,4,5,8,9,10)|is.na(HR_MenopauseStatus)) ~ HZ_HR_MenopauseAge, 
      TRUE ~ menoage)) %>%
  mutate(HZ_HR_MenopauseAgeExact_imp = ifelse(!is.na(HZ_HR_MenopauseAgeExact), HZ_HR_MenopauseAgeExact, HZ_HR_MenopauseAge))
      


##### Impute menopausal age and status ##### 

data10 = data9.1 %>% mutate(
  # menopausal status
    baseline_meno = case_when(
        menoage <= AgeExact_Baseline ~ 1, 
        menoage > AgeExact_Baseline ~ 0, 
        is.na(menoage) & HR_MenopauseStatus %in% c(6,7,11,9) ~ 0, 
        is.na(menoage) & HR_MenopauseStatus %in% c(1,2,3,4,5,10) ~ 1,
        is.na(menoage) & (HR_MenopauseStatus == 8 | is.na(HR_MenopauseStatus)) & HZ_HR_MenopauseStatus %in% c(6,7,11) ~ 0),
   # impute both menoage if HZ_menoage is missing; 
    menoage = case_when(
      is.na(menoage) & baseline_meno==0 & HZ_HR_MenopauseStatus %in% c(1,2,3,5,8) & AgeExact_Baseline<55 & FU_UECa_EOFAgeExact_imp_hyst>=55 ~ 55, 
      is.na(menoage) & baseline_meno==0 & HZ_HR_MenopauseStatus %in% c(1,2,3,5,8) & AgeExact_Baseline>=55 ~ AgeExact_Baseline+0.1, 
      is.na(menoage) & baseline_meno==0 & HZ_HR_MenopauseStatus %in% c(1,2,3,5) & FU_UECa_EOFAgeExact_imp_hyst<55 ~ FU_UECa_EOFAgeExact_imp_hyst-0.1, 
      is.na(menoage) & baseline_meno==0 & HZ_HR_MenopauseStatus==8 & FU_UECa_EOFAgeExact_imp_hyst<55 ~ NA_real_, 
      TRUE ~ menoage), 
   # menopausal status at EOF 
    eof_meno = case_when(
      baseline_meno==1 ~ 1, 
      baseline_meno==0 & menoage > FU_UECa_EOFAgeExact_imp_hyst ~ 0, 
      baseline_meno==0 & menoage <= FU_UECa_EOFAgeExact_imp_hyst ~ 1,
      baseline_meno==0 & HZ_HR_MenopauseStatus %in% c(6,7,11) ~ 0), 
    eof_meno = replace(eof_meno, HZ_HR_MenopauseStatus==9, 0),
    eof_meno = ifelse(baseline_meno==0 & HZ_HR_MenopauseStatus==8 & FU_UECa_EOFAgeExact_imp_hyst<55, 0, eof_meno))


##### Medical record #####
data11 = data10 %>%
  # medically confirmed UECa cases (DxType 6 self-report; 9 next of kin)
  # if not medically confirmed cases; coded them as non-cases, but censored at the original age (either diagnosis or hysterectomy)
  mutate(UECa_event = FU_UECa_Event_hyst, 
         UECa_event_MR = ifelse(FU_UECa_Event_hyst == 1 & FU_UECa_DxType_Source %in% c(6,9), 0, FU_UECa_Event_hyst), 
         # endometrial
         UECa_event_type1 = ifelse(UECa_event_MR == 1 & FU_UECa_DxHist01 %in% c(8140,8262,8380,8382,8480,8560,8570), 1, 0), 
         UECa_event_type2 = ifelse(UECa_event_MR == 1 & FU_UECa_DxHist01 %in% c(8310,8323,8441,8460,8950,8980), 1, 0),
         
         # Assuming the subjects with endometrial cancer related histology code are also endometrial 
         UECa_event_endometrial = ifelse(UECa_event_MR==1 & (FU_UECa_DxType == 11.2|
                                                                             FU_UECa_DxHist01 %in% c(8140,8262,8380,8382,8480,8560,8570,8310,8323,8441,8460,8950,8980)), 1, 0)) 



data12 <- data11 %>% rename("PermDye"="mix_PC85", 
                            "PermDye_other"="mix_PC89",
                            "SemiPermDye"="mix_PC91",
                            "SemiPermDye_other"="mix_PC95",
                            "HairColorRinse"="mix_PC97",
                            "Bleach"="mix_PC99",
                            "Highlight"="mix_PC101", 
                            "Straightener"="mix_PC103", 
                            "Straightener_other"="mix_PC105", 
                            "HairPerm"="mix_PC107",
                            "PermDye_duration"="PC87",
                            "SemiPermDye_duration"="PC93", 
                            "PermDye_dark"="PC86dark",
                            "PermDye_light"="PC86light",
                            "SemiPermDye_dark"="PC92dark",
                            "SemiPermDye_light"="PC92light")

exposure <- c("PermDye", "PermDye_other", "SemiPermDye", "SemiPermDye_other", 
              "HairColorRinse", "Bleach", "Highlight", "Straightener", "Straightener_other", "HairPerm")

mutate.yn <- function(x)(ifelse(x == 0, 0, # did not use
                                ifelse(x %in% c(1,2), 1, NA))) # ever

data13 <- data12 %>% mutate_at(exposure, funs(yn = mutate.yn(.))) %>% 
  mutate(HairPerm2 = ifelse(!is.na(PC107_perm), PC107_perm, P149_perm), 
         HairPerm2 = case_when(HairPerm2==1~0, 
                              HairPerm2%in%c(2)~1, 
                              HairPerm2%in%c(3,4,5,6)~2))


data14 <- data13 %>% select(UECa_event, UECa_event_MR, UECa_event_endometrial, UECa_event_type1, UECa_event_type2, eof_age, 
                            PermDye, PermDye_other, SemiPermDye, SemiPermDye_other, HairColorRinse, Bleach, Highlight, Straightener, Straightener_other, HairPerm, HairPerm2, 
                            PermDye_yn, PermDye_other_yn, SemiPermDye_yn, SemiPermDye_other_yn, HairColorRinse_yn, Bleach_yn, Highlight_yn, Straightener_yn, Straightener_other_yn, HairPerm_yn, 
                            PermDye_duration, SemiPermDye_duration, PermDye_dark, PermDye_light, SemiPermDye_dark, SemiPermDye_light, 
                            age, race, race2, educ, BMI, BMIcat, PA, PAcat, baseline_meno, menoage, eof_meno, 
                            parity, smoking, smokingcat2, alcohol, BCduration, HRT2, menarche, menarchecat, py, hairdresserjob, HH_PSID)
                 

setwd("U:\\P0_UE_hairproduct\\Rscript\\Upload to FTP\\")
write_csv(data14, "hairUEca.csv")

