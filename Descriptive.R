#Paper2Experiments
load("C:/Users/kazem/Documents/Data/tdToTSD/td/CompleteData.Rdata")
load("C:/Users/kazem/Documents/Data/tdToTSD/td/data.Rdata")
setwd("C:/Users/kazem/Dropbox/NESTIMATOR PhD/R Scripts/nestimator2/Nestimator2")
source("Functions.R")

TLevels <- levels(data$trtgrp)

temp <- data$trtgrp
Summary_table_mat <- data %>% mutate_if(is_logical, as.integer) %>%
  select(-trtgrp) %>%
  fastDummies::dummy_cols(remove_most_frequent_dummy = TRUE,
                          remove_selected_columns = TRUE) %>% 
  mutate_if(is_integer, function(x) 100*x) %>% as.data.frame()
Summary_table_mat$trtgrp <- temp
Summary_table_missingness <- is.na(Summary_table_mat) %>% as.data.frame() %>% 
  mutate_all(function(x) 100*x)
Summary_table_missingness$trtgrp <- temp

temp <- Summary_table_mat %>% select(-trtgrp) %>%
    summarise_all(mean, na.rm = TRUE) %>% mutate(trtgrp = "All", .before = Rem)
Summary_table <- Summary_table_mat %>% group_by(trtgrp) %>%
  summarise_all(mean, na.rm = TRUE) %>% rbind(temp) %>% mutate_if(is.numeric, round)

temp <- Summary_table_mat %>% select(-trtgrp) %>%
  summarise_all(sd, na.rm = TRUE) %>% mutate(trtgrp = "All", .before = Rem)
Sd_table <- Summary_table_mat %>% group_by(trtgrp) %>%
  summarise_all(sd, na.rm = TRUE) %>% rbind(temp) %>% mutate_if(is.numeric, round)


temp <- Summary_table_missingness %>% select(-trtgrp) %>%
  summarise_all(mean, na.rm = TRUE) %>%
  mutate(trtgrp = "All", .before = Rem)
Missing_table <- Summary_table_missingness %>% group_by(trtgrp) %>%
  summarise_all(mean, na.rm = TRUE) %>%
  rbind(temp) %>% mutate_if(is.numeric, round)

main_cols <- c(1:3,68:71,4:9,28,10:15,20,21,19,22:24,26,27)
Summary_table_main <- Summary_table[,main_cols]
Sd_table_main <- Sd_table[, main_cols]
Missing_table_main <- Missing_table[, main_cols]
main_col_names <- 
  c("Treatment", "Remission at 3 months (%)", "Mean age (years)", "Male (%)", 
    "Smoker current (%)", "Smoker occationally (%)", "Smoker previously (%)", 
    "RF positive (%)", "Anti-CCP positive (%)", "Previous MTX user (%)", 
    "Methotrexate at BL (%)", "Prednisolone at BL (%)", "Sulfasalazine at BL (%)", 
    "Other DMARDs  at BL (%)", "Mean PGA at BL (0-100)", "Mean PhGA at BL (0-100)", 
    "Mean ESR at BL", "Mean CRP at BL", "Mean SJC28 at BL", "Mean TJC28 at BL",
    "Mean pain at BL (0-100)", "Mean joint pain at BL (0-100)", 
    "Mean EQ index at BL", "Mean MHAQ score at BL", "Mean RAID score at BL", 
    "Employment (%)", "Mean diagnosis duration (years)", 
    "Mean symptom duration (years)")
for (c in 2:length(main_cols)) {
  if(grepl("%", main_col_names[c])){
    Summary_table_main[[c]] <- 
      paste0(Summary_table_main[[c]], 
             " [", Missing_table_main[[c]], "]")
    
  }else{
    Summary_table_main[[c]] <- 
      paste0(Summary_table_main[[c]], " (", Sd_table_main[[c]], ")", 
             " [", Missing_table_main[[c]], "]")
  }
}

Summary_table_main$trtgrp <- 
  paste0(Summary_table_main$trtgrp, " (", 
         c(summary(data$trtgrp), nrow(data)), ")")


colnames(Summary_table_main) <- main_col_names
Summary_table_main <- t(Summary_table_main)
knitr::kable(Summary_table_main, "latex")

com_cols <- c(1, 29:49)
Summary_table_com <- Summary_table[,com_cols]
Missing_table_com <- Missing_table[, com_cols]
Missing_table_com_Missing <- rowMeans(Missing_table_com[, 2:22])

Summary_table_com$trtgrp <- 
  paste0(Summary_table_com$trtgrp, " (", 
         c(summary(data$trtgrp), nrow(data)), ")", " [", 
         Missing_table_com_Missing, "]")


Summary_table_com <- t(Summary_table_com)
knitr::kable(Summary_table_com, "latex")

Z_summary <- c("Year", paste0(1:5, ". place"), "N") %>%
  cbind(10:23 %>% rbind(data.frame(
          "2010" = c("Inf", "Gol", "Cer", "Eta", "Ada"),
          "2011" = c("Eta", "Inf", "Cer", "Gol", "Ada"),
          "2012" = c("Inf", "Eta", "Cer", "Gol", "Ada"),
          "2013" = c("Cer", "Inf", "Gol", "Eta", "Ada"),
          "2014" = c("Inf", "Cer", "Gol", "Ada", "Eta"),
          "2015" = c("Inf", "Cer", "Gol", "Eta", "Ada"),
          "2016" = c("Inf", "Cer", "Eta", "Gol", "Ada"),
          "2017" = c("Inf", "Eta", "Gol", "Cer", "Ada"),
          "2018" = c("Inf", "Eta", "Cer", "Ada", "Gol"),
          "2019" = c("Ada", "Inf", "Eta", "Cer", "Gol"),
          "2020" = c("Ada", "Inf", "Eta", "Cer", "Gol") ,
          "2021" = c("Ada", "Inf", "Eta", "Cer", "Gol") ,
          "2022" = c("Ada", "Inf", "Eta", "Cer", "Gol") ,
          "2023" = c("Ada", "Inf", "Eta", "Cer", "Gol")
), summary(data$Z_value)))

knitr::kable(Z_summary, "latex")

