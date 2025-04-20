library(dplyr)
library(ggplot2)
library(ggsignif)


# manually setting up small dataframe of uncleaved/cleaved ratios
df <- data.frame(
  replicate = rep(1:5, 3),
  strain = rep(c("BW25113", "JW0052-O", "JW0052-T"), each = 5),
  ratio = c(1.533743049, 1.086889208, 3.069795065, 1.064075768, 1.019896929,
            1.348154101, 1.178908676, 1.006219478, 1.713837556, 1.019843676,
            3.541619037, 4.542294511, 2.952085883, 2.785979036, 10.00331502)
)

#  indicating IPTG induction status
df <- df %>%
  mutate(IPTG = ifelse(replicate %in% c(4,5), "Yes", "No"))


print(df)

# Run two-way ANOVA (strain and IPTG as factors)
# strain = Do different strains have different mean ratios?
# IPTG = Does IPTG induction change the ratio, overall?	
# strain:IPTG = Does the effect of IPTG depend on the strain?	

anova_results <- aov(ratio ~ strain * IPTG, data = df)
summary(anova_results)

anova_strain <- aov(ratio ~ strain, data = df)
summary(anova_strain)

# Tukey's post-hoc test to identify specific differences:
TukeyHSD(anova_strain)

collapsed_df <- df %>%
  group_by(strain) %>%
  summarise(
    mean_ratio = mean(ratio),
    sd_ratio = sd(ratio),
    .groups = "drop"
  )

plot <- ggplot(collapsed_df, aes(x = strain, y = mean_ratio, fill = strain)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_errorbar(aes(ymin = mean_ratio - sd_ratio, ymax = mean_ratio + sd_ratio),
                width = 0.2) +
  labs(
    title = "",
    x = "Strain",
    y = "Uncleaved:Cleaved BrkA"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 11),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )

plot

plot_final <- plot +
  geom_signif(
    annotations = c("*", "*"),
    y_position = c(9.5, 11),  
    xmin = c(1, 2),
    xmax = c(3, 3),
    tip_length = 0.05,
    textsize = 5
  ) +
  coord_cartesian(ylim = c(0, 12), clip = "off") +
  scale_fill_manual(values = c(
    "#E69F00",  # orange
    "#56B4E9",  # sky blue
    "#009E73"   # bluish green
  )) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),   
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
    axis.text.x = element_text(color = "black"),
    axis.title.y = element_text(margin = margin(r = 10)),
    plot.margin = margin(10, 10, 10, 10),
    legend.position = "none"
  )


plot_final

df_2 <- df %>%
  mutate(group = paste(strain, IPTG, sep = "_"))

print(df_2)

anova_group <- aov(ratio ~ group, data = df_2)
summary(anova_group)

tukey_results <- TukeyHSD(anova_group)
print(tukey_results)

