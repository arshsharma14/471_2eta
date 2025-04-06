### SDS-EDTA data analysis
install.packages(“growthcurver”)
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

library(tidyverse)

head(gc_rename)
kruskal.test(auc_l ~ sample, data = gc_rename)

install.packages("FSA")
library(FSA)
dunn_result <- dunnTest(auc_l ~ sample, 
                        data = gc_rename, 
                        method = "bh")  # Adjusts for multiple comparisons
print(dunn_result)
library(ggsignif)


library(ggplot2)
AUC_Strains <- ggplot(gc_rename, aes(x = sample, y = auc_l)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.6, color = "blue") +  # Shows replicates
  labs(title = "AUC Across Strains", x = "", y = "AUC")
AUC_Strains

significant_pairs <- list(c("TropiniBrkA", "WT"))  # The only pair with P.adj < 0.05


messing_around <- ggplot(gc_rename, aes(x = sample, y = auc_l)) +
  
  geom_boxplot(
    width = 0.6,                # Width of boxes
    outlier.shape = NA,         
    color = "black",            # Border color
    fill = "gray90",            # Fill color
    alpha = 0.7,                # Transparency
    lwd = 0.5                   #
  ) +
  
  geom_signif(
    comparisons = significant_pairs,
    annotations = "*",  # Or "p = 0.030" for exact value
    tip_length = 0.01,  # Length of the horizontal bars
    textsize = 5,       # Size of the asterisk
    vjust = -0.5,       # Vertical adjustment of the label
    y_position = max(gc_rename$auc_l) * 1.1  # Position above highest point
  ) +
  geom_jitter(
    width = 0.15,               
    height = 0,                 
    alpha = 0.7,                # Slightly more opaque
    size = 2.5,                 # Larger points
    color = "#3575b5"           # Nice blue color
  ) +
  
  # Axis and labels
  labs(
    y = "Area Under the Curve (AUC)", 
    x = NULL
  ) +
  
  
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),   
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5), # Add border
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11), # Rotated x-labels
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),
    plot.margin = margin(10, 10, 10, 10)     # Adjust plot margins
  ) +
  
  # Consistent y-axis expansion
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
messing_around

