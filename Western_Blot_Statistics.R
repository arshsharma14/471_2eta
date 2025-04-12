library(dplyr)

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
anova_results <- aov(ratio ~ strain * IPTG, data = df)
summary(anova_results)

anova_strain <- aov(ratio ~ strain, data = df)
summary(anova_strain)

# Tukey's post-hoc test to identify specific differences:
TukeyHSD(anova_strain)

