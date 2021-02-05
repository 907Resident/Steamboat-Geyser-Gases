library(tidyverse)
library(car)
library(gplots)
library(nortest)
library(EnvStats)
setwd("C:/Users/ayersj/OneDrive - Vanderbilt/Projects/VolcanicGases/Yellowstone/Data/EGM_CO2_Fluxes/All")
CO2Fluxes <- read_csv("EGM_CO2FluxesAllDataWithSoilTypePositiveFluxes.csv")
CO2Fluxes$Site <- factor(CO2Fluxes$Site)
hist(CO2Fluxes$fCO2, xlab = "CO2 Flux (g m-2 h-1)", main = "Histogram of CO2 Fluxes")
# Make Boxcox plots for fCO2 using EnvStats (see Millard 2013 pp. 75-77) to ID distribution type 
# CO2Fluxes$fCO2t <- replace(CO2Fluxes$fCO2, which(CO2Fluxes$fCO2 <= 0), NA)
boxcox.list <- boxcox(CO2Fluxes$fCO2)
plot(boxcox.list, main = "CO2 Flux (g m-2 h-1)")
# Peaks at a value of lambda = 0, indicating a log transform is appropriate.
hist(CO2Fluxes$log10fCO2, xlab = "log10 CO2 Flux (g m-2 h-1)", main = "Histogram of log10 CO2 Fluxes")
# Test for normality
lillie.test(CO2Fluxes$log10fCO2)
# P < 0.05, so reject the null hypothesis that log10fCO2 distribution is normal
CO2Flux.lnorm <- gofTest(CO2Fluxes$fCO2, dist = "lnorm")
plot(CO2Flux.lnorm, digits = 3)
# Since no other data transformation is better, use nonparametric tests on log transformed values.
myvars <- c("fCO2", "log10fCO2")
summary(CO2Fluxes[myvars])
CO2FluxesSumm <- CO2Fluxes %>%
  group_by(Site) %>%
  dplyr::summarize(
    median = round(median(log10fCO2, na.rm = TRUE), digits = 3),
    stdev = round(sd(log10fCO2, na.rm = TRUE), digits = 3),
    count = n()
  )
# Log10 CO2 Fluxes
CO2FluxesSumm
write.csv(CO2FluxesSumm, file = "SummaryLog10CO2Fluxes.csv")
stripChart(log10fCO2 ~ Site, data = CO2Fluxes, show.ci = FALSE, ylab = expression(paste(log[10], "[ fCO2 (g m-2 h-1)]")))
CO2Fluxes$Year <- factor(CO2Fluxes$Year)
p = ggplot(CO2Fluxes, aes(x = Site, y = log10fCO2, fill = Year, position = Year)) +
  geom_boxplot(na.rm = TRUE, alpha = 0.3) +
  theme_bw() +
  scale_fill_grey(start = 0, end = .9) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle=45, size=11, hjust = 1)) +
  labs(y = "log fCO2 (g m-2 h-1)") +
  stat_n_text()
print(p)
ggsave("log10fCO2Boxplot.png", plot = p)
# Use Levene's test to test if sites have the same variance for log10fCO2 (test for homoscedasticity). 
leveneTest(CO2Fluxes$log10fCO2, CO2Fluxes$Site)
# P < 0.05, so reject the null hypothesis that site log10fCO2 variances are equal
# ANOVA assumptions of normality and equal variance are invalid, so
# use nonparametric Kruskal-Wallis one-way analysis of variance
kruskal.test(log10fCO2 ~ Site, data = CO2Fluxes)
# P < 0.05, so reject the null hypothesis that site log10fCO2 medians are equal
plotmeans(CO2Fluxes$log10fCO2 ~ CO2Fluxes$Site, xlab = "Site", ylab = "log10 fCO2", main = "Mean plot\nwith 95% CI")
