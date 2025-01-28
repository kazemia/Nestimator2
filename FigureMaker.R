library(dplyr)
library(ggplot2)
library(tidyr)
my_fig_data <- my_comparisons %>% 
  mutate(contrast = gsub("Inf", "INX", contrast),
         contrast = gsub("Cer", "CZP", contrast),
         contrast = gsub("Ada", "ADM", contrast),
         contrast = gsub("Eta", "ETN", contrast),
         contrast = gsub("Gol", "GOM", contrast)) %>%
  mutate(R_CI_low = sapply(strsplit(R_CI, ","), function(x) as.numeric(substring(x[1],2))),
         R_CI_high = sapply(strsplit(R_CI, ","), function(x) as.numeric(substring(x[2],1, nchar(x[2]) - 1))),
         IV_CI_low = sapply(strsplit(IV_CI, ","), function(x) as.numeric(substring(x[1],2))),
         IV_CI_high = sapply(strsplit(IV_CI, ","), function(x) as.numeric(substring(x[2],1, nchar(x[2]) - 1)))) %>%
  select(-R_CI, -IV_CI) %>%
  pivot_longer(
    cols = c(estimate, IV_estimate),
    names_to = "type",
    values_to = "value"
  ) %>%
  mutate(
    CI_low = ifelse(type == "estimate", R_CI_low, IV_CI_low),
    CI_high = ifelse(type == "estimate", R_CI_high, IV_CI_high),
    contrast = as.factor(contrast),
    y_position = ifelse(type == "estimate", as.numeric(contrast) + 0.1, as.numeric(contrast) - 0.1)
  ) %>% as.data.frame()

ggplot(my_fig_data, aes(x = value, y = y_position, color = type)) +
  geom_point(size = 3) +  # Points for estimates
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.1) +  # Horizontal error bars
  geom_vline(xintercept = 0, color = "black", linetype = "dashed", alpha = 0.5) +
  scale_y_continuous(
    breaks = seq_along(levels(my_fig_data$contrast)), 
    labels = levels(my_fig_data$contrast)
  ) +
  scale_x_continuous(
    n.breaks = 6  # Increase the number of x-axis breaks (ticks)
  ) +
  scale_color_manual(values = c("blue", "red"), labels = c("Regression", "IV")) +
  labs(
    x = "Difference in probability of remission at three months",
    y = "Contrast",
    color = "Analysis"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10)) + 
  annotate("text", x = 1, y = 5.5, label = "Favours left side", angle = 90) + 
  annotate("text", x = -0.5, y = 5.5, label = "Favours right side", angle = 90)
