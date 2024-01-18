#Paper2Experiments
load("adeff2.Rdata")
source("Functions2.R")
source("stata2.R")
TLevels <- levels(adeff$trtgrp)
Z <- list("2010" = c("Inf", "Gol", "Cer", "Eta", "Ada"),
          "2011" = c("Eta", "Inf", "Cer", "Gol", "Ada"),
          "2012" = c("Inf", "Eta", "Cer", "Gol", "Ada"),
          "2013" = c("Cer", "Inf", "Gol", "Eta", "Ada"),
          "2014" = c("Inf", "Cer", "Gol", "Ada", "Eta"),
          "2015" = c("Inf", "Cer", "Gol", "Eta", "Ada"),
          "2016" = c("Inf", "Cer", "Eta", "Gol", "Ada"),
          "2017" = c("Inf", "Eta", "Gol", "Cer", "Ada"),
          "2018" = c("Inf", "Eta", "Cer", "Ada", "Gol"),
          "2019" = c("Ada", "Inf", "Eta", "Cer", "Gol") 
)
adeff$Z_value <- factor(adeff$Z_value, levels = 1:10)
#Find the identifiables
nZ <- length(Z)
names(Z) <- 1:nZ
A_full <- GenerateA(TLevels)
R_full <- MakeR(A_full, Z, T_decider)
KB_full <- MakeKB(R_full, TLevels, 4)
b_full <- KbSolver(KB_full, 3)
Pis_full <- PiIdentifier(b_full)

A_limited <- list(
  A1 = c("Cer"),
  A2 = c("Inf"),
  A3 = c("Eta"),
#  A4 = c("Ada"),
  A5 = c("Gol", "Cer"),
#  A6 = c("Gol", "Ada"),
#  A7 = c("Ada", "Cer"),
  A8 = c("Inf", "Ada"),
#  A9 = c("Ada", "Eta"),
#  A10 = c("Ada", "Cer", "Gol"),
  A11 = c("Inf", "Ada", "Gol"),
  A12 = c("Ada", "Eta", "Inf"),
  A13 = c("Inf", "Ada", "Cer", "Gol"),
  A14 = c("Ada", "Cer", "Gol", "Eta"),
  A15 = c("Inf", "Eta", "Ada", "Cer", "Gol")
)
R_limited <- MakeR(A_limited, Z, T_decider)
KB_limited <- MakeKB(R_limited, TLevels, 4)
b_limited <- KbSolver(KB_limited, 3)
Pis_limited <- PiIdentifier(b_limited)

# Rename columns and create 'trtgrp' variable
P_Z_crude <- MakeP_Z(adeff, "Z_value", "trtgrp")
Q_Z_crude <- MakeQ_Z(adeff, "Z_value", "trtgrp", "das28crprem")
C_adj <- c("das28crp_BL", "pga_BL", "esr_BL", "crp_BL", "sjc28_BL",
           "tjc28_BL", "sjc32_BL", "tjc32_BL", "cdai_BL")
P_Z_adj <- MakeP_Z(adeff, "Z_value", "trtgrp", C_adj, TRUE)
Q_Z_adj <- MakeQ_Z(adeff, "Z_value", "trtgrp", "das28crprem", C_adj, TRUE, "binomial", P_Z_adj)

P_Z_crude$Z_value <- as.numeric(P_Z_crude$Z_value)
Q_Z_crude$Z_value <- as.numeric(Q_Z_crude$Z_value)
P_Z_adj$Z_value <- as.numeric(P_Z_adj$Z_value)
Q_Z_adj$Z_value <- as.numeric(Q_Z_adj$Z_value)

plot(1:10,rowSums(Q_Z_adj %>% select(-Z_value)))
plot(1:10,rowSums(Q_Z_crude %>% select(-Z_value)))
#Calculate adjusted P_Sigma and LATES
P_Sigma_adj_full <- P_SigmaIdentifier(P_Z_adj, KB_full, b_full)
LATEs_adj_full <- LATEIdentifier(Q_Z_adj, KB_full, b_full, P_Sigma_adj_full)


P_Sigma_crude_full <- P_SigmaIdentifier(P_Z_crude, KB_full, b_full)
LATEs_crude_full <- LATEIdentifier(Q_Z_crude, KB_full, b_full, P_Sigma_crude_full)


P_Sigma_adj_limited <- P_SigmaIdentifier(P_Z_adj, KB_limited, b_limited)
LATEs_adj_limited <- LATEIdentifier(Q_Z_adj, KB_limited, b_limited, P_Sigma_adj_limited)


P_Sigma_crude_limited <- P_SigmaIdentifier(P_Z_crude, KB_limited, b_limited)
LATEs_crude_limited <- LATEIdentifier(Q_Z_crude, KB_limited, b_limited, P_Sigma_crude_limited)