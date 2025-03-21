#Paper2Experiments
load("C:/Users/kazem/Documents/Data/tdToTSD/td/CompleteData.Rdata")
load("C:/Users/kazem/Documents/Data/tdToTSD/td/data.Rdata")
setwd("C:/Users/kazem/Dropbox/NESTIMATOR PhD/R Scripts/nestimator2/Nestimator2")
source("Functions.R")
library(ggplot2)

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
Z <- list(
  "2010" = c("INX", "GOM", "CZP", "ETN", "ADM"),
  "2011" = c("ETN", "INX", "CZP", "GOM", "ADM"),
  "2012" = c("INX", "ETN", "CZP", "GOM", "ADM"),
  "2013" = c("CZP", "INX", "GOM", "ETN", "ADM"),
  "2014" = c("INX", "CZP", "GOM", "ADM", "ETN"),
  "2015" = c("INX", "CZP", "GOM", "ETN", "ADM"),
  "2016" = c("INX", "CZP", "ETN", "GOM", "ADM"),
  "2017" = c("INX", "ETN", "GOM", "CZP", "ADM"),
  "2018" = c("INX", "ETN", "CZP", "ADM", "GOM"),
  "2019" = c("ADM", "INX", "ETN", "CZP", "GOM"),
  "2020" = c("ADM", "INX", "ETN", "CZP", "GOM") ,
  "2021" = c("ADM", "INX", "ETN", "CZP", "GOM") ,
  "2022" = c("ADM", "INX", "ETN", "CZP", "GOM") ,
  "2023" = c("ADM", "INX", "ETN", "CZP", "GOM")
)

Z_summary <- data.frame(
  "2010" = c("INX", "GOM", "CZP", "ETN", "ADM"),
  "2011" = c("ETN", "INX", "CZP", "GOM", "ADM"),
  "2012" = c("INX", "ETN", "CZP", "GOM", "ADM"),
  "2013" = c("CZP", "INX", "GOM", "ETN", "ADM"),
  "2014" = c("INX", "CZP", "GOM", "ADM", "ETN"),
  "2015" = c("INX", "CZP", "GOM", "ETN", "ADM"),
  "2016" = c("INX", "CZP", "ETN", "GOM", "ADM"),
  "2017" = c("INX", "ETN", "GOM", "CZP", "ADM"),
  "2018" = c("INX", "ETN", "CZP", "ADM", "GOM"),
  "2019" = c("ADM", "INX", "ETN", "CZP", "GOM"),
  "2020" = c("ADM", "INX", "ETN", "CZP", "GOM") ,
  "2021" = c("ADM", "INX", "ETN", "CZP", "GOM") ,
  "2022" = c("ADM", "INX", "ETN", "CZP", "GOM") ,
  "2023" = c("ADM", "INX", "ETN", "CZP", "GOM")
)

Z_summary <- Z_summary %>% rbind(summary(CompleteData$Z_value))

Z_summary <- Z_summary %>% rbind(round((CompleteData %>% group_by(Z_value) %>% summarise(RR = mean(Rem)))$RR*100)/100)

Z_summary <- 10:23 %>% rbind(Z_summary)

Z_summary <- c("Year", paste0(1:5, ". place"), "N", "P(Y)") %>%
  cbind(Z_summary)

P_Z <- MakeP_Z(CompleteData, "Z_value", "trtgrp") %>% select(-Z_value) %>% t()
rownames(P_Z) <- c("ADM", "CZP", "ETN", "GOM", "INX")
counter <- 1
MyHeatMap <- data.frame("2010" = unname(P_Z[Z[[counter]],counter]))
for (counter in 2:14) {
  MyHeatMap <- MyHeatMap %>% cbind(unname(P_Z[Z[[counter]],counter]))
}
colnames(MyHeatMap) <- 2010:2023
dt2 <- MyHeatMap %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
dt3 <- Z_summary[2:8, 2:15] %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
dt3$rowname <- as.character(as.numeric(dt3$rowname) - 1)
dt3$colname <- substring(dt3$colname, 2)
dt3$rowname[dt3$rowname == "6"] <- "N"
dt3$rowname[dt3$rowname == "7"] <- "P(Y)"
ggplot(dt2) +
  geom_tile(aes(x = colname, y = rowname, fill = value)) + 
  scale_fill_gradient(low = "white", high = "red") +
  geom_text(data = dt3, aes(x = colname, y = rowname, label = value)) + 
  ylab("Position in LIS") + xlab("Year") + labs(fill = "Probability") + 
  geom_hline(yintercept = 2.5)+
  theme(panel.background = element_blank()) + 
  scale_y_discrete(limits = c("P(Y)", "N", as.character(5:1)))

knitr::kable(Z_summary, "latex")

data <- data %>% 
  mutate(BL_DAS28 = (0.56*sqrt(BL_tjc28)) + 
           (0.28*sqrt(BL_sjc28)) + 
           (0.014*BL_pga) + 
           (0.36*log(BL_crp+1)) + 0.96)
data %>% select(trtgrp, BL_DAS28) %>%group_by(trtgrp) %>% summarise_all(mean, na.rm = TRUE)
