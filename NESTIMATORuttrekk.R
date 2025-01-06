#This code doesn't produce a clean data set. Check missingness produced after each step.
library(haven)
library(dplyr)

pdir <- "C:/Users/kazem/Documents/Data/tdToTSD/td/"
tdds <- pdir %>% paste0("tdds.dta") %>% read_dta()
#Make STATA factor variables to R character
for (c in c("trtgrp", "wdrreas", "center")) {
  tdds[[c]] <- names(attributes(tdds[[c]])$labels)[
    match(tdds[[c]], attributes(tdds[[c]])$labels)]
}
#Save the full treatment name in a separate column
tdds$trtfull <- tdds$trtgrp
#Rename the treatments with the substance short name
myMeds <- c("Inf", "Eta", "Gol", "Cer", "Ada")
for(m in myMeds){
  tdds$trtgrp[grepl(m, tdds$trtgrp)] <- m
}
#Filter out treatment courses not in my meds
tdds <- tdds %>% filter((trtgrp %in% myMeds) & (bldt > as.Date("31012010", "%d%m%Y")) & (bldt < as.Date("01022024", "%d%m%Y")))
#Note: Missing trtgrps would disappear once we filter on RA later.
#Note: Missing bldts have no following visit.

my_tcs <- tdds$tcid

#Select useful columns
tdds <- tdds %>% select(tcid, patid, dob, center, trtgrp, trtfull, bldt, lvno, lvdt,
                        wdrdt, wdrreas, wdrreasc, trtdiscdt)

#Read more data
tddm <- pdir %>% paste0("tddm.dta") %>% read_dta() %>% filter(tcid %in% my_tcs)
#Make STATA factor variables to R character
for (c in c("diaggrp", "ethnicity", "smoker", "coffee", "edul")) {
  tddm[[c]] <- names(attributes(tddm[[c]])$labels)[
    match(tddm[[c]], attributes(tddm[[c]])$labels)]
}
#Select useful columns and filter out non RA
tddm <- tddm %>% select(tcid, diaggrp, diagdt, sympdebdt, sex, age,
                        smoker, coffee, edul, inv, nurse, rhmfact, 
                        anticcp) %>% filter(diaggrp == "RA")
#Note there are no missing diaggrp
#Merge tddm and tdds
MyTCData <- inner_join(tdds, tddm, by = "tcid")
my_tcs <- MyTCData$tcid

#Clean up
my_tcs <- MyTCData$tcid
remove(tdds,tddm)

#Fix missing prevbiol
tdtrt <- pdir %>% paste0("tdtrt.dta") %>% read_dta() %>% filter(tcid %in% my_tcs)
#Fix missing prevbiol
bla <- tdtrt %>% filter(is.na(prevbiol))
bla <- left_join(bla, MyTCData %>% filter(nchar(patid) > 2), by = "tcid")
temp <- MyTCData %>% filter(patid %in% bla$patid) %>% 
  left_join(tdtrt, by = "tcid") %>% arrange(bldt) %>% group_by(patid) %>% 
  mutate(regbio = 1:n()) %>% 
  select(tcid, patid, bldt, prevbiol, nprevbiol, regbio) %>% as.data.frame()
temp$nprevbiol[is.na(temp$nprevbiol)] <- temp$regbio[is.na(temp$nprevbiol)]-1
temp <- temp %>% group_by(patid) %>% 
  mutate(unobservedbio = (sum(nprevbiol) > (n()*(n()-1)/2))) %>% as.data.frame()
temp <- temp %>% mutate(prevbiolfix = 
                          (prevbiol == 1) |
                          (regbio > 1) | 
                          (nprevbiol > 0) |
                          unobservedbio) %>% select(tcid, prevbiolfix)
temp$prevbiolfix[is.na(temp$prevbiolfix)] <- FALSE
tdtrt <- left_join(tdtrt, temp, by = "tcid")
tdtrt$prevbiol[is.na(tdtrt$prevbiol)] <- 
  as.numeric(tdtrt$prevbiolfix[is.na(tdtrt$prevbiol)])

#Make STATA factor variables to R character
for (c in c("trtgrp", "cosdmard", "Biologicstype", "Treatmenttype", "Treatgrou",
            "groumain", paste0("M", c(1:34, 41,42,60,61), "Used"),
            paste0("M", c(1:34, 41,42,60,61), "Reas"))) {
  tdtrt[[c]] <- names(attributes(tdtrt[[c]])$labels)[
    match(tdtrt[[c]], attributes(tdtrt[[c]])$labels)]
}


#Select useful columns and filter out non bio naives
#Note prevmtx is in agreement with Methotrexate and methotrexate
tdtrt <- tdtrt %>% select(tcid, cosdmard, prevbiol, prevmtx) %>% 
  filter(prevbiol == 0)
MyTCData <- inner_join(MyTCData, tdtrt, by = "tcid")
my_tcs <- MyTCData$tcid

#Clean up
remove(tdtrt, bla, temp)

#Read more data
tdconmed <- pdir %>% paste0("tdconmed.dta") %>% read_dta() %>% filter((tcid %in% my_tcs) & (visit == 1))

for (c in c("medgroup", "medication")) {
  tdconmed[[c]] <- names(attributes(tdconmed[[c]])$labels)[
    match(tdconmed[[c]], attributes(tdconmed[[c]])$labels)]
}
tdconmed <- tdconmed %>%
  filter(medication %in% c("Sulfasalazine", "Methotrexate", "Prednisolone", "Prednisone")) %>% 
  select(tcid, medication) %>% 
  unique() %>% 
  mutate(taking = TRUE) %>% 
  tidyr::pivot_wider(id_cols = tcid, names_from = medication, 
                     values_from = taking,  names_prefix = "BL_") %>% 
  right_join(MyTCData %>% select(tcid), by = "tcid")
tdconmed[is.na(tdconmed)] <- FALSE
MyTCData <- left_join(MyTCData, tdconmed, by = "tcid")
remove(tdconmed)

#Read more data
tdsv <- pdir %>% paste0("tdsv.dta") %>% read_dta() %>% filter(tcid %in% my_tcs)

#Fix missing studyday
tdsv <- tdsv %>% 
  left_join(MyTCData %>% select(tcid, bldtfix = bldt), by = "tcid") %>% 
  mutate(bldt = ifelse(is.na(bldt), bldtfix, bldt),
         studyday = ifelse(is.na(studyday), visit_dt - bldt , studyday),
         visit_no = ifelse((studyday == 0) & 
                             (visit_no != 100) & 
                             (tcid != "205-00379"), 1, visit_no)) %>%
  select(-bldtfix)

tdsv <- tdsv %>% filter(studyday < 180)

#Create a visit ID
tdsv <- tdsv %>% mutate(visit_id = paste(tcid, visit_no))
my_visits <- unique(tdsv$visit_id)

#Read more data
tddisact <- pdir %>% paste0("tddisact.dta") %>% read_dta() %>% 
  mutate(visit_id = paste(tcid, visit_no)) %>% filter(visit_id %in% my_visits) %>% 
  select(-das28, -das28cat, -das28crp, -das28crpcat, -cdai, -cdaicat, -sdai,
         -sdaicat, -sdairem, -asdascat, -asdasrem) %>% 
  select(-all_of(paste0("basdai", 1:6)))

tdeq5d <- pdir %>% paste0("tdeq5d.dta") %>% read_dta() %>% 
  mutate(visit_id = paste(tcid, visit_no)) %>% 
  filter(visit_id %in% my_visits) %>% select(visit_id, EQ_index)

tdmhaq <- pdir %>% paste0("tdmhaq.dta") %>% read_dta() %>% 
  mutate(visit_id = paste(tcid, visit_no)) %>% 
  filter(visit_id %in% my_visits) %>% select(-all_of(paste0("mhaq", 1:8)))

tdraid <- pdir %>% paste0("tdraid.dta") %>% read_dta() %>% 
  mutate(visit_id = paste(tcid, visit_no)) %>% 
  filter(visit_id %in% my_visits) %>% select(visit_id, raid_score)
tdwpai <- pdir %>% paste0("tdwpai.dta") %>% read_dta() %>% 
  mutate(visit_id = paste(tcid, visit_no)) %>% filter(visit_id %in% my_visits) %>% 
  select(-wpai_2, -wpai_3, -wpai_4, -impairment, -absenteeism, -presenteeism)

my_joiner <- function(x, y, by = "visit_id"){
  commonCols <- colnames(y)[colnames(y) %in% colnames(x)]
  commonCols <- commonCols[commonCols != by]
  y <- y %>% select(-all_of(commonCols))
  y$MissingScore <- rowSums(is.na(y))
  y <- arrange(y, desc(MissingScore))
  out <- left_join(x, y, by = by, multiple = "last")
  return(out %>% select(-MissingScore))
}

MyVisitData <- tdsv %>% my_joiner(tddisact) %>% my_joiner(tdeq5d) %>% 
  my_joiner(tdmhaq) %>% my_joiner(tdraid) %>% my_joiner(tdwpai)

#Clean up
remove(tddisact, tdeq5d, tdmhaq, tdraid, tdwpai, tdsv)

#Make BL Data
MyVisitData <- MyVisitData %>% 
  select(-tcidn, -tcido, -timepoint, -visit_id, -bldt, -visit, -center, -vsource)
MyBLData <- MyVisitData %>% filter(visit_no == 1) %>% 
  select(-visit_no, -visit_dt, -studyday)
MyVisitData <- MyVisitData %>% filter(visit_no != 1) %>% 
  select(-visit_dt)

#Make target visit data
MyVisitData <- MyVisitData %>% filter(visit_no != 100)


#Remove duplicate rows
MyVisitData$MissingScore <- rowSums(is.na(MyVisitData))
MyVisitData <- MyVisitData %>% mutate(dupid = paste(tcid, studyday)) %>%
  group_by(dupid) %>% arrange(MissingScore) %>% filter(row_number() == 1) %>%
  as.data.frame() %>% select(-dupid, -MissingScore)


#Identify the visit closest to 3 months and within 180 days after baseline
MyVisitData <- MyVisitData %>% 
  group_by(tcid) %>% 
  mutate(studydaytargetdiff = min(abs(studyday - 90)),
         studydaytarget = ifelse((90+studydaytargetdiff) %in% studyday, 
                                 90+studydaytargetdiff, 
                                 90-studydaytargetdiff)) %>%
  filter(studyday == studydaytarget) %>% 
  mutate(visit_target = visit_no) %>% 
  select(tcid, visit_target) %>% 
  right_join(MyVisitData, by = "tcid") %>% 
  filter(visit_no == visit_target) %>%
  select(-visit_target)

#Read more data
comorbidity <- pdir %>% paste0("Comorbidity.dta") %>% read_dta() %>% 
  mutate(visit_id = paste0(2, pid, " ", visit_no)) %>% 
  filter(visit_id %in% my_visits)
my_comorbidities <- sapply(1:21, function(x) attributes(comorbidity[[paste0("comorb", x)]])$label)
comorbidity <- comorbidity %>% group_by(pid) %>% 
  summarise_if(is.numeric, sum, na.rm = TRUE) %>% select(-visit_no)
comorbidity[,2:22] <- (comorbidity[,2:22] > 0)
comorbidity <- comorbidity %>% mutate(tcid = paste0(2, pid)) %>% select(-pid)

MyTCData <- left_join(MyTCData, comorbidity, by = "tcid")
remove(comorbidity)

#Join all data
colnames(MyBLData)[2:ncol(MyBLData)] <- 
  paste("BL", colnames(MyBLData)[2:ncol(MyBLData)], sep = "_")

MyVisitData <- MyVisitData %>% select(-visit_no)
colnames(MyVisitData)[2:ncol(MyVisitData)] <- 
  paste("T", colnames(MyVisitData)[2:ncol(MyVisitData)], sep = "_")

MyData <- MyTCData %>% left_join(MyBLData, by = "tcid") %>% 
  left_join(MyVisitData, by = "tcid")

#Clean up
remove(MyBLData, MyTCData, MyVisitData)

#Fix missing target
MyData$Rem <- (MyData$T_das28crprem == 1)
MyData <- MyData %>% mutate(Rem = ifelse(!is.na(Rem), Rem, 
                                         T_das28rem == 1))
MyData <- MyData %>% mutate(Rem = ifelse(!is.na(Rem), Rem, 
                                         T_cdairem == 1))
MyData <- MyData %>% 
  mutate(Rem = ifelse(!is.na(Rem), Rem, 
                      ifelse((MyData$wdrreas %in% 
                               c("Adverse events", 
                                 "Lack of efficacy", 
                                 "Lack of efficacy and adverse events")) & 
                               ((lvdt - bldt)<365),
                             FALSE, NA)))

#Fix missing RF factor and antiCCP
library(readxl)
MyLabData <- read_excel(paste0(pdir,"ExtendedDataExtraction_2024-8-29.xlsx"), 
                        sheet = "DiagnosticTestsLabData") %>% 
  filter(DiagnosticTestLab %in% c("RF IgM", "CCP")) %>% 
  select(PatientID, DiagnosticTestLab, Pos_Neg_Inconclusive, 
         DiagnosticTestLabDate) %>%
  filter(!is.na(Pos_Neg_Inconclusive)) %>%
  group_by(PatientID, DiagnosticTestLab) %>% arrange(DiagnosticTestLabDate) %>%
  filter(row_number()==n()) %>% as.data.frame() %>% select(-DiagnosticTestLabDate) %>%
  tidyr::pivot_wider(id_cols = PatientID, names_from = DiagnosticTestLab, 
                     values_from = Pos_Neg_Inconclusive)
colnames(MyLabData) <- c("PatientID", "RF", "CCP")
MyData$PatientID <- sapply(strsplit(MyData$patid, "-"), function(x) as.numeric(x[2]))
MyData <- left_join(MyData, MyLabData, by = "PatientID")
MyData$rhmfact[is.na(MyData$rhmfact)] <- MyData$RF[is.na(MyData$rhmfact)]
MyData$anticcp[is.na(MyData$anticcp)] <- MyData$CCP[is.na(MyData$anticcp)]
MyData <- MyData %>% select(-PatientID, -RF, -CCP)
remove(MyLabData)


#Setting Z values
TLevels <- c("Inf", "Ada", "Cer", "Eta", "Gol")
Z <- list("2010" = c("Inf", "Gol", "Cer", "Eta", "Ada"),
          "2011" = c("Eta", "Inf", "Cer", "Gol", "Ada"),
          "2012" = c("Inf", "Eta", "Cer", "Gol", "Ada"),
          "2013" = c("Cer", "Inf", "Gol", "Eta", "Ada"),
          "2014" = c("Inf", "Cer", "Gol", "Ada", "Eta"),
          "2015" = c("Inf", "Cer", "Gol", "Eta", "Ada"),
          "2016" = c("Inf", "Cer", "Eta", "Gol", "Ada"),
          "2017" = c("Inf", "Eta", "Gol", "Cer", "Ada"),
          "2018" = c("Inf", "Eta", "Cer", "Ada", "Gol"),
          "2019" = c("Ada", "Inf", "Eta", "Cer", "Gol"),
          "2020" = c("Ada", "Inf", "Eta", "Cer", "Gol"),
          "2021" = c("Ada", "Inf", "Eta", "Cer", "Gol"),
          "2022" = c("Ada", "Inf", "Eta", "Cer", "Gol"),
          "2023" = c("Ada", "Inf", "Eta", "Cer", "Gol")
)
data <- MyData
remove(MyData)
data$Z_value <- NA
data$PeriodDay <- NA
Z_decider <- function(data, washout = 0) {
  for(y in 2010:2013){
    StartDate <- as.Date(paste0(y, "-01-31"), "%Y-%m-%d") + washout
    EndDate <- as.Date(paste0(y+1, "-02-01"), "%Y-%m-%d") + washout
    data$Z_value[(data$bldt > StartDate) & (data$bldt < EndDate)] <- y-2009
    data$PeriodDay[(data$bldt > StartDate) & (data$bldt < EndDate)] <- 
      data$bldt[(data$bldt > StartDate) & (data$bldt < EndDate)] - StartDate
  }
  y <- 2014
  StartDate <- as.Date(paste0(y, "-01-31"), "%Y-%m-%d") + washout
  EndDate <- as.Date(paste0(y+1, "-03-01"), "%Y-%m-%d") + washout
  data$Z_value[(data$bldt > StartDate) & (data$bldt < EndDate)] <- y-2009
  data$PeriodDay[(data$bldt > StartDate) & (data$bldt < EndDate)] <- 
    data$bldt[(data$bldt > StartDate) & (data$bldt < EndDate)] - StartDate
  y <- 2015
  StartDate <- as.Date(paste0(y, "-02-28"), "%Y-%m-%d") + washout
  EndDate <- as.Date(paste0(y+1, "-03-01"), "%Y-%m-%d") + washout
  data$Z_value[(data$bldt > StartDate) & (data$bldt < EndDate)] <- y-2009
  data$PeriodDay[(data$bldt > StartDate) & (data$bldt < EndDate)] <- 
    data$bldt[(data$bldt > StartDate) & (data$bldt < EndDate)] - StartDate
  y <- 2016
  StartDate <- as.Date(paste0(y, "-02-29"), "%Y-%m-%d") + washout
  EndDate <- as.Date(paste0(y+1, "-03-01"), "%Y-%m-%d") + washout
  data$Z_value[(data$bldt > StartDate) & (data$bldt < EndDate)] <- y-2009
  data$PeriodDay[(data$bldt > StartDate) & (data$bldt < EndDate)] <- 
    data$bldt[(data$bldt > StartDate) & (data$bldt < EndDate)] - StartDate
  y <- 2017
  StartDate <- as.Date(paste0(y, "-02-28"), "%Y-%m-%d") + washout
  EndDate <- as.Date(paste0(y+1, "-02-01"), "%Y-%m-%d") + washout
  data$Z_value[(data$bldt > StartDate) & (data$bldt < EndDate)] <- y-2009
  data$PeriodDay[(data$bldt > StartDate) & (data$bldt < EndDate)] <- 
    data$bldt[(data$bldt > StartDate) & (data$bldt < EndDate)] - StartDate
  for(y in 2018:2023){
    StartDate <- as.Date(paste0(y, "-01-31"), "%Y-%m-%d") + washout
    EndDate <- as.Date(paste0(y+1, "-02-01"), "%Y-%m-%d") + washout
    data$Z_value[(data$bldt > StartDate) & (data$bldt < EndDate)] <- y-2009
    data$PeriodDay[(data$bldt > StartDate) & (data$bldt < EndDate)] <- 
      data$bldt[(data$bldt > StartDate) & (data$bldt < EndDate)] - StartDate
  }
  return(data)
}

data <- Z_decider(data)
data <- data %>% filter(!is.na(Z_value))
data$leftover = NA
data$leftover[data$Z_value == 1] <- FALSE
data$leftover[data$Z_value > 1] <- 
  unname((data$trtgrp[data$Z_value > 1] == 
            sapply(Z[data$Z_value[data$Z_value > 1]-1], function(x) x[1])) & 
           (data$trtgrp[data$Z_value > 1] != 
              sapply(Z[data$Z_value[data$Z_value > 1]], function(x) x[1])))
hist((data %>% filter(leftover & (Z_value %in% c(3,5))))$PeriodDay, xlab = "Days since NDPC change",
     main = "Probability of picking previous winner", breaks = (0:40)*10, probability = T)

data <- data %>% select(-leftover, -PeriodDay)
#backup_data <- data
#data <- data %>% filter(Z_value < 12)
#leftover_prob <- c()
#for (washout in (0:10)*10) {
#  data$Z_value <- NA
#  data <- Z_decider(data, washout = washout)
#  leftover_data <- data %>% filter(!is.na(Z_value))
#  leftover_data$leftover = NA
#  leftover_data$leftover[leftover_data$Z_value == 1] <- FALSE
#  leftover_data$leftover[leftover_data$Z_value > 1] <- 
#    (leftover_data$trtgrp[leftover_data$Z_value > 1] == 
#       sapply(Z[leftover_data$Z_value[leftover_data$Z_value > 1]-1], 
#              function(x) x[1])) & 
#    (leftover_data$trtgrp[leftover_data$Z_value > 1] != 
#       sapply(Z[leftover_data$Z_value[leftover_data$Z_value > 1]], 
#              function(x) x[1]))
#  leftover_prob <- c(leftover_prob, mean(leftover_data$leftover))
#}
#plot((0:10)*10, leftover_prob)




#Handle the rest of the missingness
#Remove unnecessary variables and correct the type
data <- data %>% select(-tcid, -patid, -trtfull, -wdrreasc, -lvno, -lvdt,
                            -wdrdt, -diaggrp, -inv, -nurse, -prevbiol,
                            -BL_das28rem, -BL_das28crprem, -BL_cdairem, -BL_mcii,
                            -T_das28rem, -T_das28crprem, -T_cdairem, -T_mcii) %>% 
  mutate(center = as.factor(center),
         trtgrp = as.factor(trtgrp),
         wdrreas = as.factor(wdrreas),
         age = ifelse(!is.na(age), age, (bldt - dob)/365),
         T_treatmentdur = as.numeric(trtdiscdt-bldt),
         BL_diagdur = as.numeric(bldt-diagdt),
         BL_sympdur = as.numeric(bldt-sympdebdt),
         sex = as.factor(ifelse(sex %in% c("Male", "Female"), sex, NA)),
         smoker = as.factor(smoker),
         coffee = as.factor(coffee),
         edul = as.factor(edul),
         rhmfact = as.logical(as.numeric(rhmfact)),
         anticcp = as.logical(as.numeric(anticcp)),
         cosdmard = as.factor(cosdmard),
         prevmtx = as.logical(as.numeric(prevmtx)), 
         BL_ACREULARrem = as.logical(as.numeric(BL_ACREULARrem)),
         T_ACREULARrem = as.logical(as.numeric(T_ACREULARrem)), 
         BL_pass1 = as.logical(BL_pass1-1),
         T_pass1 = as.logical(T_pass1-1), 
         BL_wpai_1 = as.logical(as.numeric(BL_wpai_1)),
         T_wpai_1 = as.logical(as.numeric(T_wpai_1))) %>%
  select(-dob, -trtdiscdt, -diagdt, -sympdebdt, -bldt)
#Fix BL drugs
data$BL_Methotrexate[data$cosdmard == "MTX"] <- TRUE
data$BL_otherDmards <- data$cosdmard == "Other sDMARDs"
data$BL_Prednisolone <- data$BL_Prednisolone | data$BL_Prednisone
data <- data %>% select(-BL_Prednisone, -cosdmard)
data$BL_otherDmards[data$BL_Sulfasalazine | data$BL_Prednisolone] <- FALSE
#Remove variables with more than 70% missing
temp <- colMeans(is.na(data))
data <- data %>% select(-all_of(names(temp)[which(temp > 0.7)]))
remove(temp)

#data$bldt <- as.numeric(data$bldt)
data$Z_value <- as.factor(data$Z_value)

#Impute the rest of the missing data
MyMiceModel <- mice::mice(data = data %>% select(-trtgrp, -Z_value), m = 1, seed = 123, maxit = 1)
CompleteData <- mice::complete(MyMiceModel)
CompleteData$trtgrp <- data$trtgrp
CompleteData$Z_value <- data$Z_value
CompleteData <- CompleteData %>% select(Z_value, trtgrp, Rem, center, sex, age,
                                        smoker, rhmfact, anticcp, prevmtx, 
                                        all_of(colnames(CompleteData)
                                               [grepl("BL_", 
                                                      colnames(CompleteData))]),
                                        all_of(paste0("comorb", 1:21)))
#Remove ASDAS & BASDAI
CompleteData <- CompleteData %>% select(-BL_asdas, -BL_basdai)

colnames(CompleteData)[33:53] <- my_comorbidities


data <- data %>% select(Z_value, trtgrp, Rem, center, sex, age,
                                        smoker, rhmfact, anticcp, prevmtx, 
                                        all_of(colnames(data)
                                               [grepl("BL_", 
                                                      colnames(data))]),
                                        all_of(paste0("comorb", 1:21)))
#Remove ASDAS & BASDAI
data <- data %>% select(-BL_asdas, -BL_basdai)

colnames(data)[33:53] <- my_comorbidities

CompleteData <- CompleteData %>% mutate(BL_wpai_1 = as.logical(BL_wpai_1),
                                        prevmtx = as.logical(prevmtx),
                                        anticcp = as.logical(anticcp),
                                        rhmfact = as.logical(rhmfact),
                                        Rem = as.logical(Rem),
                                        BL_diagdur = BL_diagdur/365.25,
                                        BL_sympdur = BL_sympdur/365.25)
data <- data %>% mutate(BL_wpai_1 = as.logical(BL_wpai_1),
                                        prevmtx = as.logical(prevmtx),
                                        anticcp = as.logical(anticcp),
                                        rhmfact = as.logical(rhmfact),
                                        Rem = as.logical(Rem),
                                        BL_diagdur = BL_diagdur/365.25,
                                        BL_sympdur = BL_sympdur/365.25)


save(CompleteData, file = paste0(pdir, "CompleteData.Rdata"))
save(data, file = paste0(pdir, "data.Rdata"))


#The code after this is experimental









#Scale the numeric data
#Remove wdrreas
Main_data <- CompleteData %>% select(Z_value, trtgrp, Rem)
Cov_data <- CompleteData %>% select(-Z_value, -trtgrp, -Rem) %>% 
  fastDummies::dummy_cols(remove_most_frequent_dummy = TRUE,
                          remove_selected_columns = TRUE) %>% 
  mutate_all(as.numeric) %>% scale() %>% as.data.frame()

#Perform pca on covariate data and pick enough components to cover 90% variance
pca_model <- princomp(Cov_data)
summary(pca_model)
pca_data <- predict(pca_model, Cov_data)[,1:29] %>% as.data.frame()

#Tidy up
pca_data <- cbind(Main_data, pca_data)
pca_data$Rem <- as.logical(pca_data$Rem)
remove(Cov_data, Main_data, MyMiceModel, pca_model)

#GLM base model
my_base_model <- glm(Rem ~ ., data = pca_data %>% select(-Z_value), family = "binomial")
my_base_model %>% marginaleffects::avg_comparisons(variables = list(trtgrp = "pairwise"), type = "response", newdata = "marginalmeans")


#my_base_model <- glm(Rem ~ ., data = CompleteData %>% select(-Z_value), family = "binomial")
#my_base_model %>% marginaleffects::avg_comparisons(variables = list(trtgrp = "pairwise"), type = "response", newdata = "marginalmeans")

my_iptw <- ce_estimate(y = pca_data$Rem, 
                       x = pca_data %>% select(-Rem, -trtgrp, -Z_value), 
                       w = as.numeric(pca_data$trtgrp), 
                       method = "IPTW-SL", estimand = "ATE", 
                       sl_library =  c("SL.glm", "SL.glmnet", "SL.rpart"))
summary(my_iptw)

pca_data$w <- my_iptw$weight
w_q <- quantile(pca_data$w, probs = c(0.05,0.95))
pca_data <- pca_data %>% mutate(h_w = (w < w_q[1]) | (w > w_q[2]))

library(ggplot2)
ggplot(pca_data %>% filter(trtgrp == TLevels[5]), aes(x = Comp.1, y = Comp.2, colour = h_w)) + geom_point()

CompleteData$w <- my_iptw$weight
w_q <- quantile(CompleteData$w, probs = c(0.05,0.95))
CompleteData <- CompleteData %>% mutate(h_w = (w < w_q[1]) | (w > w_q[2]))

for (t in TLevels) {
  tree <- rpart::rpart(h_w ~ ., 
                       data = CompleteData %>% filter(trtgrp == t) %>% 
                         select(-w, -Z_value, -Rem, -trtgrp),
                       control=rpart::rpart.control(cp=.0001))
  best <- tree$cptable[which.min(tree$cptable[,"xerror"]),"CP"]
  pruned_tree <- rpart::prune(tree, cp=best)
  rpart.plot::prp(pruned_tree, main = t)
}

(CompleteData %>% group_by(Z_value) %>% summarise(Rem = mean(Rem)) %>% 
  mutate(Z_value = as.numeric(Z_value)) %>% ggplot(aes(x = Z_value, y = Rem))) +
  geom_point()




