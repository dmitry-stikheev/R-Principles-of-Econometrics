library(PoEdata) #for PoE datasets
library(knitr) #for referenced tables with kable()
library(xtable) #makes data frame for kable
library(printr) #automatically prints output nicely
library(effects)
library(car) #for testing hypotheses
library(AER)
library(broom) #for tidy lm output and function glance()
library(stats)
library(dplyr)


## Testing Simultaneous Hypotheses
data("andy")

alpha <- 0.05
N <- nrow(andy)
K <- 4
J <- 2
fcr <- qf(1 - alpha, J, N - K)
mod1 <- lm(sales ~ price + advert + I(advert^2), data = andy)
SSEu <- anova(mod1)[4, 2]
mod2 <- lm(sales ~ price, data = andy)
SSEr <- anova(mod2)[2, 2]
f <- ((SSEr - SSEu) / J) / (SSEu / (N - K))

curve(df(x, J, N - K), 0, 9, xlab = "f", ylab = "")
abline(v = c(fcr, f), col = c("red", "blue"), lty = c(2, 3))
legend("topright", legend = c("fcr", "f"),
       lty = c(2, 3), col = c("red", "blue"))

Hnull <- c("advert = 0", "I(advert^2) = 0")
linearHypothesis(mod1, Hnull)

hyp <- c("advert + 3.8*I(advert^2) = 1",
         "(Intercept) + 6*price + 1.9*advert + 3.61*I(advert^2) = 80")
results <- tidy(linearHypothesis(mod1, hyp))
kable(results)


## Omitted Variable Bias
data("edu_inc")

mod1 <- lm(faminc ~ he + we, data = edu_inc)
mod2 <- lm(faminc ~ he, data = edu_inc)
coef(summary(mod1))
coef(summary(mod2))

b3 <- coef(mod1)[[3]]
cov_x2x3 <- cov(edu_inc$he, edu_inc$we)
var_x2 <- var(edu_inc$he)

bias_b2 <- b3*(cov_x2x3 / var_x2)
bias_b21 <- coef(mod2)[[2]] - coef(mod1)[[2]] # alternative approach


## Irrelevant variables
data("edu_inc")

mod3 <- lm(faminc ~ he + we + kl6, data = edu_inc)
mod4 <- lm(faminc ~ he + we + kl6 + xtra_x5 + xtra_x6,
           data = edu_inc)
correct_model <- coef(summary(mod3))
incorrect_model <- coef(summary(mod4))
correct_model
incorrect_model


## Model Selection Criteria
data("edu_inc")

mod1 <- lm(faminc~he, data=edu_inc)                       # Comparing 4 models
mod2 <- lm(faminc~he+we, data=edu_inc)
mod3 <- lm(faminc~he+we+kl6, data=edu_inc)
mod4 <- lm(faminc~he+we+kl6+xtra_x5+xtra_x6, data=edu_inc)

# Extracting R^2, adjusted R^2, AIC and SC (BIC)
r1 <- as.numeric(glance(mod1))
r2 <- as.numeric(glance(mod2))
r3 <- as.numeric(glance(mod3))
r4 <- as.numeric(glance(mod4))

selection_criteria <- data.frame(rbind(r1, r2, r3, r4))[, c(1, 2, 8, 9)]

rownames(selection_criteria) <- c("he","he, we","he, we, kl6",
   "he, we, kl6, xtra_x5, xtra_x6")
colnames(selection_criteria) <- c("Rsq", "Adj Rsq", "AIC", "BIC")
selection_criteria

# Alternative way of extraction criteria
library(stats)
smod1 <- summary(mod1)

Rsq <- smod1$r.squared
AdjRsq <- smod1$adj.r.squared
aic <- AIC(mod1)
bic <- BIC(mod1)
criteria <- c(Rsq, AdjRsq, aic, bic)


## RESET Test
mod3 <- lm(faminc ~ he + we + kl6, data = edu_inc)

resettest(mod3, power = 2, type = "fitted")
resettest(mod3, power = 2:3, type = "fitted")


## Collinearity
data("cars", package = "PoEdata")

mod1 <- lm(mpg ~ cyl, data = cars)
m1 <- coef(summary(mod1))

mod2 <- lm(mpg ~ cyl + eng + wgt, data = cars)
m2 <- coef(summary(mod2))
comparison <- rbind(m1, m2)
comparison


## Variance Inflation Factor (VIF)
mod2 <- lm(mpg ~ cyl + eng + wgt, data = cars) # soecifying the model

# Setting up artificial models
art1 <- lm(cyl ~ eng + wgt, data = cars)
art2 <- lm(eng ~ cyl + wgt, data = cars)
art3 <- lm(wgt ~ cyl + eng, data = cars)
sart1 <- summary(art1)
sart2 <- summary(art2)
sart3<- summary(art3)

# Extracting R squared values
Rsq1 <- sart1$r.squared
Rsq2 <- sart2$r.squared
Rsq3 <- sart3$r.squared

# Computing VIF
vif_cyl <- 1 / (1 - Rsq1)
vif_eng <- 1 / (1 - Rsq2)
vif_wgt <- 1 / (1 - Rsq3)

# Alternative way through the "car" library
vif_mod2 <- vif(mod2)


## Prediction
data("andy")
mod3 <- lm(sales ~ price + advert + I(advert^2), data = andy)
predpoint <- data.frame(price = 6, advert = 1.9)
predict(mod3, newdata = predpoint, interval = "prediction")


## ENJOY!!!
