#This code doesn't produce a clean data set. Check missingness produced after each step.
library(haven)
library(dplyr)

pdir <- "N:/Forskning/Nordmard/Nye NOR-DMARD/New Statistics/ndm/Datasets/td/"
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
#Filter out treatment courses before 2010
tdds <- tdds %>% filter((trtgrp %in% myMeds) & 
                          (tdds$bldt >= as.Date("01012010", "%d%m%Y")))

my_tcs <- tdds$tcid
#Read more data
tddm <- pdir %>% paste0("tddm.dta") %>% read_dta() %>% filter(tcid %in% my_tcs)
tdtrt <- pdir %>% paste0("tdtrt.dta") %>% read_dta() %>% filter(tcid %in% my_tcs)
#Select useful columns
tdds <- tdds %>% select(tcid, dob, center, trtgrp, trtfull, bldt, lvno, lvdt,
                        wdrdt, wdrreas, wdrreasc, trtdiscdt)
#Make STATA factor variables to R character
for (c in c("diaggrp", "ethnicity", "smoker", "coffee", "edul")) {
  tddm[[c]] <- names(attributes(tddm[[c]])$labels)[
    match(tddm[[c]], attributes(tddm[[c]])$labels)]
}
#Select useful columns and filter out non RA
tddm <- tddm %>% select(tcid, diaggrp, diagdt, sympdebdt, sex, age,
                        smoker, coffee, edul, inv, nurse, rhmfact, 
                        anticcp) %>% filter(diaggrp == "RA")
#Merge tddm and tdds
MyTCData <- tdds %>% merge(tddm, by = "tcid")

#Make STATA factor variables to R character
for (c in c("trtgrp", "cosdmard", "Biologicstype", "Treatmenttype", "Treatgrou",
            "groumain", paste0("M", c(1:34, 41,42,60,61), "Used"),
            paste0("M", c(1:34, 41,42,60,61), "Reas"))) {
  tdtrt[[c]] <- names(attributes(tdtrt[[c]])$labels)[
    match(tdtrt[[c]], attributes(tdtrt[[c]])$labels)]
}
#Select useful columns and filter out non bio naives
tdtrt <- tdtrt %>% select(tcid, trtgrp, cosdmard, prevbiol, prevmtx, 
                          methotrexate, Methotrexate,prednisolone, 
                          prednisolon_dose, sulfasalazine, Sulfasalazine, 
                          Biologicstype, Treatmenttype, Treatgrou, groumain) %>% 
  filter(prevbiol == 0)
MyTCData <- MyTCData %>% merge(tdtrt, by = "tcid")

#Clean up
remove(tdds, tddm, tdtrt)
my_tcs <- unique(MyTCData$tcid)
#Read more data
tdsv <- pdir %>% paste0("tdsv.dta") %>% read_dta() %>% filter(tcid %in% my_tcs)

#Identify the visit closest to 3 months and within 2 years after baseline
my_target_visits <- tdsv %>% filter((studyday != 0)) %>% 
  group_by(tcid) %>% 
  mutate(studydaytargetdiff = min(abs(studyday - 90)),
         studydaytarget = ifelse((90+studydaytargetdiff) %in% studyday, 
                                 90+studydaytargetdiff, 
                                 90-studydaytargetdiff)) %>%
  filter((studyday == studydaytarget) & (studyday < 730)) %>% 
  mutate(visit_target = visit_no) %>% 
  select(tcid, visit_target, studydaytargetdiff)

tdsv <- tdsv %>% filter(studyday < 730) %>% merge(my_target_visits, all.x = TRUE)
remove(my_target_visits)

#Create a visit ID
tdsv <- tdsv %>% mutate(visit_id = paste(tcid, visit_no))
my_visits <- unique(tdsv$visit_id)

#Read more data
tddisact <- pdir %>% paste0("tddisact.dta") %>% read_dta() %>% filter(paste(tcid, visit_no) %in% my_visits)
tdeq5d <- pdir %>% paste0("tdeq5d.dta") %>% read_dta() %>% filter(paste(tcid, visit_no) %in% my_visits)
tdmhaq <- pdir %>% paste0("tdmhaq.dta") %>% read_dta() %>% filter(paste(tcid, visit_no) %in% my_visits)
tdraid <- pdir %>% paste0("tdraid.dta") %>% read_dta() %>% filter(paste(tcid, visit_no) %in% my_visits)


MyVisitData <- merge(tdsv, tddisact, all.x = TRUE)
#remove(tddisact)

#Selecting and attaching eq5d, note some lack individual questions but not all
tdeq5d <- tdeq5d %>% select(-visit, -visit_dt)
MyVisitData <- MyVisitData %>% merge(tdeq5d, all.x = TRUE)
#remove(tdeq5d)

#Selecting and attaching mhaq
tdmhaq <- tdmhaq %>% select(tcid, visit_no, pass1, mcii, 
                            pain, jointpain, fatigue, mhaq)
MyVisitData <- MyVisitData %>% merge(tdmhaq, all.x = TRUE)
#remove(tdmhaq)

#selecting and attaching raid, note some lack individual questions but not all
tdraid <- tdraid %>% select(-visit, -visit_dt)
MyVisitData <- MyVisitData %>% merge(tdraid, all.x = TRUE)
#remove(tdraid)

#Read more data
tdconmed <- pdir %>% paste0("tdconmed.dta") %>% read_dta() %>% filter(paste(tcid, visit) %in% my_visits)






