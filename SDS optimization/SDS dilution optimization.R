### SDS dilution optimization
# Serial dilution was done for WT and mutant, but both gave similar results
# Only WT data is plotted for simplicity

#install.packages('tidyverse') #uncomment if not installed
library(tidyverse)

#read data
data <- read_csv('SDS optimization/all growth curves - SDS optimization.csv')

#clean data
data_clean <- data |>
  group_by(time) |>
  summarize(W5 = mean(W5.1:W5.3),
            W2.5 = mean(W2.5.1:W2.5.3),
            W1 = mean(W1.1:W1.3),
            W0.5 = mean(W0.5.1:W0.5.3),
            W0 = mean(W0.1:W0.3)) |>
  gather(concentration, od, -time)

legend <- c('W5' = '5%', 'W2.5' = '2.5%', 'W1' = '1%', 'W0.5' = '0.5%', 'W0' = '0%')
colours <- c('W5' = '#56B4E9', 'W2.5' = '#D55E00', 'W1' = '#009E73', 'W0.5' = '#0072B2', 'W0' = '#CC79A7')

#plot curves
growth_curve_sds <- ggplot(data_clean) +
  geom_line(aes(x = time, y = od, colour = concentration), linewidth = 0.9) +
  xlab('Time (minutes)') +
  ylab('Optical density at 600nm') +
  labs(colour = 'SDS concentration') +
  theme(text=element_text(size=14)) +
  scale_colour_discrete(type = colours,
                        labels = legend)
growth_curve_sds
