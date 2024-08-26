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
tdds <- tdds %>% filter((trtgrp %in% myMeds) & (bldt >= as.Date("01012010", "%d%m%Y")) & (bldt < as.Date("01022024", "%d%m%Y")))
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
tddisact <- pdir %>% paste0("tddisact.dta") %>% read_dta() %>% mutate(visit_id = paste(tcid, visit_no)) %>% filter(visit_id %in% my_visits)
tdeq5d <- pdir %>% paste0("tdeq5d.dta") %>% read_dta() %>% mutate(visit_id = paste(tcid, visit_no)) %>% filter(visit_id %in% my_visits)
tdmhaq <- pdir %>% paste0("tdmhaq.dta") %>% read_dta() %>% mutate(visit_id = paste(tcid, visit_no)) %>% filter(visit_id %in% my_visits)
tdraid <- pdir %>% paste0("tdraid.dta") %>% read_dta() %>% mutate(visit_id = paste(tcid, visit_no)) %>% filter(visit_id %in% my_visits)
tdwpai <- pdir %>% paste0("tdwpai.dta") %>% read_dta() %>% mutate(visit_id = paste(tcid, visit_no)) %>% filter(visit_id %in% my_visits)

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
MyBLData <- MyVisitData %>% filter(visit_no == 1)
MyVisitData <- MyVisitData %>% filter(visit_no != 1)

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


temp <- colSums(is.na(MyTCData)) %>% as.data.frame()
temp2 <- MyTCData %>% filter(is.na(wdrdt) & (lvdt-bldt < 180) & (!(tcid %in% MyVisitData$tcid)))


