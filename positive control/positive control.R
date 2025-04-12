### positive control growth curves

#install and load necessary packages
#install.packages('tidyverse') #uncomment if not installed
#install.packages('growthcurver') #uncomment if not installed
#install.packages('FSA') #uncomment if not installed
#install.packages('ggsignif') #uncomment if not installed
library(tidyverse)
library(growthcurver)
library(FSA)
library(ggsignif)

#set up plot aesthetics - keep consistent across plots
x_order <- c("WT", "WTBrkA", "Oliver", "OliverBrkA", "Tropini", "TropiniBrkA")
x_names <- c("WT" = "BW25113", "WTBrkA" = "BW25113+pPALMC1", 
             "Oliver" =  "JW0052-O", "OliverBrkA" = "JW0052-O+pPALMC1",
             "Tropini" =  "JW0052-T", "TropiniBrkA" = "JW0052-T+pPALMC1")
colours <-c("WT" = "navyblue", "WTBrkA" = "cyan",
            "Oliver" = "orange", "OliverBrkA" = "yellow",
            "Tropini" = "maroon", "TropiniBrkA" = "red")

#read data at 2mM EDTA
pc <- read_csv('positive control/all growth curves - PC.csv')

#calculate technical replicate averages to get biological replicates
pc_bioreps <- pc |>
  group_by(time) |>
  summarize(WT1 = mean(WT1.1:WT1.3), 
            WT2 = mean(WT2.1:WT2.2), 
            WT3 = mean(WT3.1:WT3.2), 
            WT4 = mean(WT4.1:WT4.2),
            WTBrkA1 = mean(WTBrkA1.1:WTBrkA1.3), 
            WTBrkA2 = mean(WTBrkA2.1:WTBrkA2.2), 
            WTBrkA3 = mean(WTBrkA3.1:WTBrkA3.2), 
            WTBrkA4 = mean(WTBrkA4.1:WTBrkA4.2),
            Oliver1 = mean(Oliver1.1:Oliver1.2), 
            Oliver2 = mean(Oliver2.1:Oliver2.2), 
            Oliver3 = mean(Oliver3.1:Oliver3.2),
            OliverBrkA1 = mean(OliverBrkA1.1:OliverBrkA1.2),
            OliverBrkA2 = mean(OliverBrkA2.1:OliverBrkA2.2), 
            OliverBrkA3 = mean(OliverBrkA3.1:OliverBrkA3.2),
            Tropini1 = mean(Tropini1.1:Tropini1.3),
            Tropini2 = mean(Tropini2.1:Tropini2.2), 
            Tropini3 = mean(Tropini3.1:Tropini3.2),
            Tropini4 = mean(Tropini4.1:Tropini4.2),
            TropiniBrkA1 = mean(TropiniBrkA1.1:TropiniBrkA1.3), 
            TropiniBrkA2 = mean(TropiniBrkA2.1:TropiniBrkA2.2),
            TropiniBrkA3 = mean(TropiniBrkA3.1:TropiniBrkA3.2))

#calculate biological replicate averages
pc_avg <- pc_bioreps |>
  group_by(time) |>
  summarize(WT = mean(WT1:WT4),
            WTBrkA = mean(WTBrkA1:WTBrkA4),
            Oliver = mean(Oliver1:Oliver3),
            OliverBrkA = mean(OliverBrkA1:OliverBrkA3),
            Tropini = mean(Tropini1:Tropini4),
            TropiniBrkA = mean(TropiniBrkA1:TropiniBrkA3))|>
  gather(strain, od, -time) |>
  add_column(n = c(1:582)) #make column with unrepeated values for full_join to work

#calculate standard deviations
pc_sd <- pc_bioreps |>
  group_by(time) |>
  summarize(WT_sd = sd(c_across(WT1:WT4)),
            WTBrkA_sd = sd(c_across(WTBrkA1:WTBrkA4)),
            Oliver_sd = sd(c_across(Oliver1:Oliver3)),
            OliverBrkA_sd = sd(c_across(OliverBrkA1:OliverBrkA3)),
            Tropini_sd = sd(c_across(Tropini1:Tropini4)),
            TropiniBrkA_sd = sd(c_across(TropiniBrkA1:TropiniBrkA3))) |>
  gather(strain_sd, sd, -time) |>
  add_column(n = c(1:582)) #make column with unrepeated values for full_join to work

#join od and sd and make clean dataframe
pc_clean <- full_join(pc_avg, pc_sd) |>
  select(!(n:strain_sd))

#reorder data by strain
pc_clean$strain <- factor(pc_clean$strain , levels = x_order)

#Growthcurver
gc_bioreps <- SummarizeGrowthByPlate(pc_bioreps)

gc_rename <- mutate(gc_bioreps, sample = str_sub(sample, end = -2))

#statistical analysis: Kruskal-Wallis and Dunn
kruskal_result <- kruskal.test(auc_l ~ sample, data = gc_rename)
print(kruskal_result)

dunn_result <- dunnTest(auc_l ~ sample, 
                        data = gc_rename, 
                        method = "bh")  # Adjusts for multiple comparisons
print(dunn_result)

significant_pairs1 <- list(c("TropiniBrkA", "WT"))
significant_pairs2 <- list(c("OliverBrkA", "WT"))

#reorder data by strain
gc_rename$sample <- factor(gc_rename$sample , levels = x_order)

auc_pc <- ggplot(gc_rename, aes(x = sample, y = auc_l, fill = sample)) +
  
  #boxplot  
  geom_boxplot(
    width = 0.6,                # Width of boxes
    outlier.shape = NA,         
    color = "black",            # Border color
    alpha = 0.7,                # Transparency
    lwd = 0.5                   
  ) +
  
  #significance bar (pair 1)  
  geom_signif(
    comparisons = significant_pairs1,
    annotations = "*",  
    tip_length = 0.05,  # Length of the horizontal bars
    textsize = 5,       # Size of the asterisk
    vjust = -0.5,       # Vertical adjustment of the label
    y_position = max(gc_rename$auc_l)*1.05  # Position above highest point
  ) +
  #significance bar (pair 2)  
  geom_signif(
    comparisons = significant_pairs2,
    annotations = "*",  
    tip_length = 0.05,  # Length of the horizontal bars
    textsize = 5,       # Size of the asterisk
    vjust = -0.5,       # Vertical adjustment of the label
    y_position = max(gc_rename$auc_l)*1.005  # Different position above highest point
  ) +
  geom_jitter(width = 0, height = 0) +
  
  # Axis and labels
  labs(
    y = "Area Under the Curve (AUC)", 
    x = "Strain",
    fill = "Strain"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_fill_manual(breaks = x_order,
                    values = colours,
                    labels = x_names) +
  
  #
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),   
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5), # Add border
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),
    plot.margin = margin(10, 10, 10, 10)     # Adjust plot margins
  ) 
auc_pc

#plot growth curves
growth_curve_pc <- ggplot(pc_clean) +
  geom_line(aes(x = time, y = od, colour = strain)) +
  xlab('Time (minutes)') +
  ylab('Optical density at 600nm') +
  labs(colour = 'Strain') +
  scale_colour_discrete(type = colours,
                        labels = x_names)
growth_curve_pc