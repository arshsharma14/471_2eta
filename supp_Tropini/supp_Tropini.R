### supplemental growth curve - Tropini



library(tidyverse)

#set up plot aesthetics - keep consistent across plots
x_order <- c("PC", "E0", "E2", "E4")
x_names <- c("0% SDS + 0mM EDTA", "0.01% SDS + 0mM EDTA", "0.01% SDS + 0.2mM EDTA", "0.01% SDS + 0.4mM EDTA")
colours <-c("#0072B2", "#009E73", "#D55E00", "#CC79A7")

#read data
tropini <- read_csv('all growth curves - Tropini.csv')

#calculate technical replicate averages to get biological replicates
tropini_bioreps <- tropini |>
  group_by(time) |>
  summarize(PC1 = mean(PC1.1:PC1.3),
            PC2 = mean(PC2.1:PC2.2),
            PC3 = mean(PC3.1:PC3.2),
            PC4 = mean(PC4.1:PC4.2),
            E0.1 = mean(E0.1.1:E0.1.3),
            E0.2 = mean(E0.2.1:E0.2.2),
            E0.3 = mean(E0.3.1:E0.3.2),
            E0.4 = mean(E0.4.1:E0.4.2),
            E2.1 = mean(E2.1.1:E2.1.3),
            E2.2 = mean(E2.2.1:E2.2.2),
            E2.3 = mean(E2.3.1:E2.3.2),
            E2.4 = mean(E2.4.1:E2.4.2),
            E4.1 = mean(E4.1.1:E4.1.3),
            E4.2 = mean(E4.2.1:E4.2.2),
            E4.3 = mean(E4.3.1:E4.3.2),
            E4.4 = mean(E4.4.1:E4.4.2))

#calculate biological replicate averages
tropini_avg <- tropini_bioreps |>
  group_by(time) |>
  summarize(PC = mean(PC1:PC4),
            E0 = mean(E0.1:E0.4),
            E2 = mean(E2.1:E2.4),
            E4 = mean(E4.1:E4.4)) |>
  gather(conc, od, -time) |>
  add_column(n = c(1:388)) #make column with unrepeated values for full_join to work

#calculate standard deviations
tropini_sd <- tropini_bioreps |>
  group_by(time) |>
  summarize(PC_sd = sd(c_across(PC1:PC4)),
            E0_sd = sd(c_across(E0.1:E0.4)),
            E2_sd = sd(c_across(E2.1:E2.4)),
            E4_sd = sd(c_across(E4.1:E4.4))) |>
  gather(conc_sd, sd, -time) |>
  add_column(n = c(1:388)) #make column with unrepeated values for full_join to work

#join od and sd and make clean dataframe
tropini_clean <- full_join(tropini_avg, tropini_sd) |>
  select(!(n:conc_sd))

#reorder data by EDTA concentration condition
tropini_clean$conc <- factor(tropini_clean$conc , levels = x_order)

#plot growth curves
growth_curve_tropini <- ggplot(tropini_clean) +
  
  geom_line(aes(x = time, y = od, colour = conc), linewidth = 1.1) +
  
  xlab('Time (minutes)') +
  ylab('Optical density at 600 nm') +
  labs(colour = 'Experimental Condition') +
  
  scale_colour_manual(values = colours, labels = x_names) +
  
  scale_y_continuous(
    breaks = seq(0, 1.25, by = 0.25),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.margin = margin(10, 10, 10, 10)
  )
growth_curve_tropini

