#Individual Analysis
#Calibration Data
concentration <- c(0.5, 1, 2, 4, 6, 8, 10)
absorbance <- c(0.01, 0.011, 0.15, 0.201, 0.373, 0.599, 0.783)

#Linear Regression Model
cal_model <- lm(absorbance ~ concentration)
plot(concentration, absorbance,
           type = "b", 
           lty = 2, 
           pch = 19,
           col = "darkblue",
           main = "Aspartame Calibration Curve",
           xlab = "Concentration (mM)",
           ylab = "Absorbance")
abline(cal_model, col = "red",
               lwd = 2,
               lty = 1)
legend("topleft",
               legend = c("Concentration Data Curve", "Linear Regression Line"),
               col = c("darkblue", "red"),
               lty = c(2,1),
               pch = c(19, NA),
               lwd = c(1, 2))
summary(cal_model)
slope <- coef(cal_model)["concentration"]
se_slope <- summary(cal_model)$coefficients["concentration", "Std. Error"]
r2 <- summary(cal_model)$r.squared
molar_abs <- slope * 1000
molar_abs_se <- se_slope * 1000
molar_abs
molar_abs_se
r2

intercept <- round(coef(cal_model)[1], 4)
slope     <- round(coef(cal_model)[2], 4)
r2        <- round(summary(cal_model)$r.squared, 4)
eq_text   <- paste0("y = ", slope, "x + ", intercept, "\nR² = ", r2)

text(x = 10, y = 0.1, labels = eq_text, pos = 2, cex = 0.9)

#Equal unknown packets at 25 mg sweetner/mL
equal_abs <- c(0.630, 0.594, 0.727)
sweet_mg_per_mL <- 25

b <- coef(cal_model)[1] #intercept
m <- coef(cal_model)[2] #slope (abs per mM)

conc_mM <- (equal_abs - b) / m
asp_mg_per_mL <- conc_mM * 294.3 / 1000
asp_mg_per_g <- asp_mg_per_mL / sweet_mg_per_mL * 1000

asp_mg_per_g
mean(asp_mg_per_g)
sd(asp_mg_per_g)
t.test(asp_mg_per_g)$conf.int
t.test(asp_mg_per_g, mu = 37)

#Pooled Analysis (Harry Horrible)
concentration <- c(0.5, 1, 2, 4, 6, 8, 10)
group1 <- c(0.01, 0.011, 0.15, 0.201, 0.373, 0.599, 0.783)
group2 <- c(0, 0.201, 0.24, 0.457, 0.625, 0.68, 1.05)
group3 <- c(0.086, 0.164, 0.303, 0.524, 0.743, 0.909, 1.064)
group4 <- c(0.173, 0.176, 0.24, 0.543, 0.799, 1.019, 1.214)

concentration_all <- rep(concentration, 4)
absorbance_all <- c(group1, group2, group3, group4)

#Linear regression model
cal_model2 <- lm(absorbance_all ~ concentration_all)
summary(cal_model2)
confint(cal_model2)

slope2 <- coef(cal_model2)[2]
se_slope2 <- summary(cal_model2)$coefficients[2, 2]
r2_2 <- summary(cal_model2)$r.squared
molar_abs2 <- slope2 * 1000
molar_abs_se2 <- se_slope2 * 1000
molar_abs2
molar_abs_se2
r2_2

#Mean + SD
all_groups <- rbind(group1, group2, group3, group4)
mean_abs <- colMeans(all_groups)
se_abs <- apply(all_groups, 2, sd) / sqrt(4)

#Error Bars
plot(concentration, mean_abs,
     pch = 19,
     col = "darkblue",
     ylim = c(0, 1.2),
     main = "Aspartame Calibration Curve (Pooled)",
     xlab = "Concentration (mM)",
     ylab = "Absorbance")
arrows(concentration, mean_abs - se_abs,
       concentration, mean_abs + se_abs,
       angle = 90, code = 3, length = 0.05, col = "darkblue")
abline(cal_model2, col = "red", lwd = 2, lty = 1)
legend("topleft",
       legend = c("Mean Absorbance ± SE", "Linear Regression Line"),
       col = c("darkblue", "red"),
       pch = c(19, NA),
       lty = c(NA, 1),
       lwd = c(NA, 2))

#Regression equation
intercept2 <- round(coef(cal_model2)[1], 4)
slope2_r <- round(coef(cal_model2)[2], 4)
r2_r <- round(r2_2, 4)
eq_text2 <- paste0("y = ", slope2_r, "x + ", intercept2, "\nR² = ", r2_r)
text(x = 10, y = 0.1, labels = eq_text2, pos = 2, cex = 0.9)

plot(concentration_all, absorbance_all,
     pch = 19,
     col = "darkblue",
     main = "Pooled Calibration: Individual Points",
     xlab = "Concentration (mM)",
     ylab = "Absorbance")
abline(cal_model2, col = "red", lwd = 2, lty = 1)

legend("topleft",
       legend = c("Individual Data Points", "Linear Regression Line"),
       col = c("darkblue", "red"),
       pch = c(19, NA),
       lty = c(NA, 1),
       lwd = c(NA, 2))

text(x = 10, y = 0.1, labels = eq_text2, pos = 2, cex = 0.9)

#Grubbs test
library(outliers)
grubbs.test(c(group1[7], group2[7], group3[7], group4[7]))

#Recalculation
equal_abs <- c(0.630, 0.594, 0.727)
sweet_mg_per_mL <- 25

b2 <- coef(cal_model2)[1]   #intercept
m2 <- coef(cal_model2)[2]   #slope (abs per mM)

conc_mM2 <- (equal_abs - b2) / m2
asp_mg_per_mL2 <- conc_mM2 * 294.3 / 1000
asp_mg_per_g2 <- asp_mg_per_mL2 / sweet_mg_per_mL * 1000

asp_mg_per_g2
mean(asp_mg_per_g2)
sd(asp_mg_per_g2)
t.test(asp_mg_per_g2)$conf.int
t.test(asp_mg_per_g2, mu = 37)
