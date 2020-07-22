library(AER) #for `ivreg()`
library(lmtest) #for `coeftest()` and `bptest()`.
library(PoEdata) #for PoE4 datasets
library(car) #for `hccm()` robust standard errors
library(sandwich)
library(stargazer)
library(dplyr)


## The instrumental variables method (IV)
data("mroz", package = "PoEdata")
mroz <- mroz %>%
  filter(lfp == 1)

# First-stage model
stage1 <- lm(educ ~ exper + I(exper^2) + mothereduc, data = mroz)
coeftest(stage1)

educ_hat <- stage1$fitted.values

wage.2sls <- lm(log(wage) ~ educ_hat + exper + I(exper^2), data = mroz)
coeftest(wage.2sls)

# The standard errors calculated by the lm() are incorrect as they use the explanatory variable estimated by the first-stage model
# The way to overcome this issue is to eploit the dedicated R library containing ivreg() function

wage.ols <- lm(log(wage) ~ educ + exper + I(exper^2), data = mroz)

wage.iv1 <- ivreg(log(wage) ~ educ + exper + I(exper^2) | exper + I(exper^2) + mothereduc, data = mroz)
wage.iv2 <- ivreg(log(wage) ~ educ + exper + I(exper^3) | exper + I(exper^2) + mothereduc + fathereduc, data = mroz)
coeftest(wage.ols)
coeftest(wage.2sls)
coeftest(wage.iv1)
coeftest(wage.iv2)

# Testing for weak instruments 
educ.ols <- lm(educ ~ exper + I(exper^2) + mothereduc + fathereduc, data = mroz)

b4 <- coeftest(educ.ols)[[4]]
b5 <- coeftest(educ.ols)[[5]]

restricted <- lm(educ ~ exper + I(exper^2), data = mroz)

SSE_un <- anova(educ.ols)[5, 2]
SSE_r <- anova(restricted)[3, 2]

df_n <- 2
df_d <- df.residual(educ.ols)

Fst <- ((SSE_r - SSE_un) / df_n) / (SSE_un / df_d) 

if(Fst > 10) {
  print("The instruments are strong")
} else {
  print("The instruments are weak")
}

## Hausman spesification test

# Initialize the first-stage regression
educ.ols <- lm(educ ~ exper + I(exper^2) + mothereduc + fathereduc, data = mroz)

# Compute fitted values and residuals from the first-stage model
# I denote residuals v_hat to avoid confusion
educ_hat <- educ.ols$fitted.values
v_hat <- educ.ols$residuals

# add v_hat to the original model
test_model <- lm(log(wage) ~ educ + exper + I(exper^2) + v_hat, data = mroz)
coeftest(test_model)

# If cov(educ_i, e_i) = 0, then the coefficient of v_hat should be equal to the slope of education
# as its effect is also exogenous. If educ is indeed endogenous, then v_hat will be biased and not
# converge to the slope of education in large samples (beta2 - slope of education, gamma - coef of v_hat)

# H0: beta2 - gamma = 0 (educ is exogenous), H1: beta2 - gamma != 0 (educ is endogenous)
b2 <- coeftest(test_model)[2, 1]
gamma <- coeftest(test_model)[5, 1]

varb2 <- vcov(test_model)[2, 2]
vargamma <- vcov(test_model)[5, 5]
cov_b2gamma <- vcov(test_model)[2, 5]
se <- sqrt(varb2 + vargamma + 2*cov_b2gamma)

alpha <- 0.05
df <- df.residual(test_model)
tc <- qt(1 - alpha / 2, df)

t <- (b2 - gamma) / se

if (abs(t) > tc) {
  print("Variable education is endogenous")
} else {
  print("Variable education is exogenous")
}
