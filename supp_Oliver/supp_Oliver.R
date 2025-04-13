### supplemental growth curve - Oliver

#install and load necessary packages
#install.packages('tidyverse') #uncomment if not installed

library(tidyverse)

#set up plot aesthetics - keep consistent across plots
x_order <- c("PC", "E0", "E2", "E4")
x_names <- c("0% SDS + 0mM EDTA", "0.01% SDS + 0mM EDTA", "0.01% SDS + 0.2mM EDTA", "0.01% SDS + 0.4mM EDTA")
colours <-c("#0072B2", "#009E73", "#D55E00", "#CC79A7")

#read data
oliver <- read_csv('supp_Oliver/all growth curves - Oliver.csv')

#calculate technical replicate averages to get biological replicates
oliver_bioreps <- oliver |>
  group_by(time) |>
  summarize(PC1 = mean(PC1.1:PC1.2),
            PC2 = mean(PC2.1:PC2.2),
            PC3 = mean(PC3.1:PC3.2),
            E0.1 = mean(E0.1.1:E0.1.2),
            E0.2 = mean(E0.2.1:E0.2.2),
            E0.3 = mean(E0.3.1:E0.3.2),
            E2.1 = mean(E2.1.1:E2.1.2),
            E2.2 = mean(E2.2.1:E2.2.2),
            E2.3 = mean(E2.3.1:E2.3.2),
            E4.1 = mean(E4.1.1:E4.1.2),
            E4.2 = mean(E4.2.1:E4.2.2),
            E4.3 = mean(E4.3.1:E4.3.2))

#calculate biological replicate averages
oliver_avg <- oliver_bioreps |>
  group_by(time) |>
  summarize(PC = mean(PC1:PC3),
            E0 = mean(E0.1:E0.3),
            E2 = mean(E2.1:E2.3),
            E4 = mean(E4.1:E4.3)) |>
  gather(conc, od, -time) |>
  add_column(n = c(1:388)) #make column with unrepeated values for full_join to work

#calculate standard deviations
oliver_sd <- oliver_bioreps |>
  group_by(time) |>
  summarize(PC_sd = sd(c_across(PC1:PC3)),
            E0_sd = sd(c_across(E0.1:E0.3)),
            E2_sd = sd(c_across(E2.1:E2.3)),
            E4_sd = sd(c_across(E4.1:E4.3))) |>
  gather(conc_sd, sd, -time) |>
  add_column(n = c(1:388)) #make column with unrepeated values for full_join to work

#join od and sd and make clean dataframe
oliver_clean <- full_join(oliver_avg, oliver_sd) |>
  select(!(n:conc_sd))

#reorder data by EDTA concentration condition
oliver_clean$conc <- factor(oliver_clean$conc , levels = x_order)

#plot growth curves
growth_curve_oliver <- ggplot(oliver_clean) +
  geom_line(aes(x = time, y = od, colour = conc), linewidth = 0.9) +
  xlab('Time (minutes)') +
  ylab('Optical density at 600nm') +
  labs(colour = 'Experimental Condition') +
  theme(text=element_text(size=14)) +
  scale_colour_discrete(type = colours,
                        labels = x_names)
growth_curve_oliver