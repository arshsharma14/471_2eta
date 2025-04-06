### SDS-EDTA data analysis

library(tidyverse)
library(growthcurver)

data2 <- read_csv('all growth curves - 2mM.csv')
colnames(data2)[1] <- 'time'

#calculate technical replicate averages to get biological replicates
data2_bioreps <- data2 |>
  group_by(time) |>
  summarize(WT1 = mean(WT1.1:WT1.3), 
          WT2 = mean(WT2.1:WT2.2), 
          WT3 = mean(WT3.1:WT3.2), 
          WTBrkA1 = mean(WTBrkA1.1:WTBrkA1.3), 
          WTBrkA2 = mean(WTBrkA2.1:WTBrkA2.2), 
          WTBrkA3 = mean(WTBrkA3.1:WTBrkA3.2), 
          Oliver1 = mean(Oliver1.1:Oliver1.2), 
          Oliver2 = mean(Oliver2.1:Oliver2.2), 
          OliverBrkA1 = mean(OliverBrkA1.1:OliverBrkA1.2),
          OliverBrkA2 = mean(OliverBrkA2.1:OliverBrkA2.2), 
          Tropini1 = mean(Tropini1.1:Tropini1.3),
          Tropini2 = mean(Tropini2.1:Tropini2.2), 
          Tropini3 = mean(Tropini3.1:Tropini3.2),
          TropiniBrkA1 = mean(TropiniBrkA1.1:TropiniBrkA1.3), 
          TropiniBrkA2 = mean(TropiniBrkA2.1:TropiniBrkA2.2),
          TropiniBrkA3 = mean(TropiniBrkA3.1:TropiniBrkA3.2))

#calculate biological replicate averages
data2_avg <- data2_bioreps |>
  group_by(time) |>
  summarize(WT = mean(WT1:WT3),
            WTBrkA = mean(WTBrkA1:WTBrkA3),
            Oliver = mean(Oliver1:Oliver2),
            OliverBrkA = mean(OliverBrkA1:OliverBrkA2),
            Tropini = mean(Tropini1:Tropini3),
            TropiniBrkA = mean(TropiniBrkA1:TropiniBrkA3)) 

data2_avg_clean <- data2_avg |>
  gather(strain, od, -time) |>
  add_column(n = c(1:582)) #make column with unrepeated values for full_join to work

#calculate standard deviations
data2_sd <- data2_bioreps |>
  group_by(time) |>
  summarize(WT_sd = sd(c_across(WT1:WT3)),
         WTBrkA_sd = sd(c_across(WTBrkA1:WTBrkA3)),
         Oliver_sd = sd(c_across(Oliver1:Oliver2)),
         OliverBrkA_sd = sd(c_across(OliverBrkA1:OliverBrkA2)),
         Tropini_sd = sd(c_across(Tropini1:Tropini3)),
         TropiniBrkA_sd = sd(c_across(TropiniBrkA1:TropiniBrkA3))) |>
  gather(strain_sd, sd, -time) |>
  add_column(n = c(1:582)) #make column with unrepeated values for full_join to work

#join od and sd and make clean dataframe
data2_clean <- full_join(data2_avg_clean, data2_sd) |>
  select(!(n:strain_sd))

#playing around with Growthcurver
gc_bioreps <- SummarizeGrowthByPlate(data2_bioreps)

gc_rename <- mutate(gc_bioreps, sample = str_sub(sample, end = -2))

gc_avg <- gc_rename |>
  group_by(sample) |>
  summarize_at(vars(auc_l), list(mean = mean, sd = sd))

#plot auc
auc_plot <- ggplot(gc_avg) +
  geom_col(aes(x = sample, y = mean)) +
  geom_errorbar(aes(x = sample, ymin = mean-sd, ymax = mean+sd))
auc_plot

#plot growth curves
growth_curve_2mM <- ggplot(data2_clean) +
  geom_line(aes(x = time, y = od, colour = strain)) +
  xlab('Time (minutes)') +
  ylab('Optical density at 600nm') +
  labs(colour = 'Strain') +
  geom_errorbar(aes(x = time, ymin = od-sd, ymax = od+sd))
growth_curve_2mM
