library(lmtest) #for coeftest() and bptest().
library(broom) #for glance() and tidy()
library(PoEdata) #for PoE4 datasets
library(car) #for hccm() robust standard errors
library(dplyr)


## Spotting heteroscedasticity in scatter plots
data("food", package = "PoEdata")

mod1 <- lm(food_exp ~ income, data = food)
b1 <- coef(mod1)[[1]]
b2 <- coef(mod1)[[2]]

# Plotting the regression line
plot(food$income, food$food_exp, type = "p",
     xlab = "Income", ylab = "Food expenditure")
curve(b1 + b2*x, add = TRUE, col = "black")

# Plotting model residuals against explanatory variable
res <- mod1$residuals
plot(res ~ food$income, xlab = "Income", ylab = "residuals")
abline(h = 0)


## Heteroscedasticity tests
alpha <- 0.05

mod1 <- lm(food_exp ~ income, data = food)
ressq <- mod1$residuals^2

modres <- lm(ressq ~ income, data = food)
smodres <- summary(modres)

# F test
df_n <- glance(modres)$df - 1
df_d <- df.residual(modres)
Fc <- qf(1 - alpha, df_n, df_d)

SST <- anova(modres)[1, 2] + anova(modres)[2, 2]
SSE <- anova(modres)[2, 2]

Fst <- ((SST - SSE) / df_n) / (SSE / df_d)

if (Fst > Fc) {
  print("Heteroscedasticity is detected")
} else {
  print("Heteroscedasticity is undetected")
}

# Breusch Pagan test
N <- nobs(modres)
df <- glance(modres)$df - 1
Rsq <- smodres$r.squared

chisqcr <- qchisq(1 - alpha, df)
chisq <- N*Rsq

if (chisq > chisqcr) {
  print("Heteroscedasticity is detected")
} else {
  print("Heteroscedasticity is undetected")
}

# The White test adding the quadratic form into the model
modres <- lm(ressq ~ income + I(income^2), data = food)

N <- nobs(modres)
Rsq <- summary(modres)$r.squared

chisqcr <- qchisq(1 - alpha, 2)
chisq <- N*Rsq

if (chisq > chisqcr) {
  print("Heteroscedasticity is detected")
} else {
  print("Heteroscedasticity is undetected")
}

# The Goldfeld-Quandt test
data("cps2", package = "PoEdata")

alpha <- 0.05

m <- cps2 %>%
  filter(metro == 1)
r <- cps2 %>%
  filter(metro == 0)

wg1 <- lm(wage ~ educ + exper, data = m)
wg0 <- lm(wage ~ educ + exper, data = r)

df1 <- df.residual(wg1)
df0 <- df.residual(wg0)

sig1squared <- summary(wg1)$sigma^2
sig0squared <- summary(wg0)$sigma^2

Fst <- (sig1squared / df1) / (sig0squared / df0)

Flc <- qf(alpha / 2, df1, df0)
Fuc <- qf(1 - alpha / 2, df1, df0)

if (Fst > Fuc | Fst < Flc) {
  print("Heteroscedasticity is detected")
} else {
  print("Heteroscedasticity is undetected")
}

# Applying Goldfeld-Quandt test dividing the data into subsamples: < median(income), > median(income)
medianincome <- median(food$income)

hi <- food %>%
  filter(income >= medianincome)
li <- food %>%
  filter(income <= medianincome)

mhi <- lm(food_exp ~ income, data = hi)
mli <- lm(food_exp ~ income, data = li)

sigmasq_hi <- summary(mhi)$sigma^2
sigmasq_li <- summary(mli)$sigma^2

df_hi <- df.residual(mhi)
df_li <- df.residual(mli)

Fst <- (sigmasq_hi / df_hi) / (sigmasq_li / df_li)

Flc <- qf(alpha / 2 , df_hi, df_li)
Fuc <- qf(1 - alpha / 2, df_hi, df_li)

if (Fst > Fuc | Fst < Flc) {
  print("Heteroscedasticity is detected")
} else {
  print("Heteroscedasticity is undetected")
}

fl <- qf(0.05, 574, 420)
fu <- qf(0.95, 96, 96)
qchisq(0.95, 4)


## Heteroscedasticity consistent standard errors
foodeq <- lm(food_exp ~ income, data = food)
regular_se <- coef(summary(foodeq))

cov1 <- hccm(foodeq, type = "hc1")
food.HC1 <- coeftest(foodeq, vcov. = cov1) # coeftest is used instead of the summary function, as it allows to use a custom vcov matrix (heteroscedasticity robust in our case)
regular_se

# Testing hypothesis using heteroscedasticity robust standard errors
data("andy", package = "PoEdata")

andy.eq <- lm(sales ~ price + advert, data = andy)

bp <- bptest(andy.eq)

b2 <- coef(andy.eq)[[2]]
b3 <- coef(andy.eq)[[3]]

# H0: price + advert = 0, H1: != 0
cov1 <- hccm(andy.eq, type = "hc1")
se_b2b3 <- sqrt(cov1[2, 2] + cov1[3, 3] + 2*cov1[2, 3])

alpha <- 0.05
df <- df.residual(andy.eq)
tc <- qt(1 - alpha / 2, df)
c <- 0

t <- (b2 + b3 - c) / se_b2b3

if (abs(t) > tc) {
  print("Reject H0")
} else {
  print("Do not reject H0")
}


## GLS: known form of variance
w <- 1 / food$income
food.wls <- lm(food_exp ~ income, weights = w, data = food)

wls_estimates <- coeftest(food.wls)

## Grouped data
data("cps2", package = "PoEdata")

rural.eq <- lm(wage ~ educ + exper, data = cps2, subset = (metro == 0))
sigR <- summary(rural.eq)$sigma

metro.eq <- lm(wage ~ educ + exper, data = cps2, subset = (metro == 1))
sigM <- summary(metro.eq)$sigma

# Creating a vector of weights
cps2$wght <- rep(0, nrow(cps2))

for(i in 1:nrow(cps2)) {
  if (cps2$metro[i] == 0) {
    cps2$wght[i] <- 1 / sigR^2
  } else {
    cps2$wght[i] <- 1 / sigM^2
  }
}

wage.fgls <- lm(wage ~ educ + exper + metro, weights = wght, data = cps2)
coeftest(wage.fgls)


## GLS: unknown form of variance
data("food", package = "PoEdata")

food.ols <- lm(food_exp ~ income, data = food)
ehatsq <- food.ols$residuals^2

sighatsq.ols <- lm(log(ehatsq) ~ log(income), data = food)
vari <- exp(sighatsq.ols$fitted.values)

food.fgls <- lm(food_exp ~ income, weights = 1 / vari, data = food)
coeftest(food.fgls)


## Heteroscedasticity in the linear probability model
data("coke", package = "PoEdata")

coke.ols <- lm(coke ~ pratio + disp_coke + disp_pepsi, data = coke)
coke.hc1 <- coeftest(coke.ols, vcov. = hccm(coke.ols, type = "hc1"))
p <- coke.ols$fitted.values

# Truncate values < 0 & > 1
pt <- p
pt[pt < 0.001] <- 0.001
pt[pt > 0.999] <- 0.999

sigsq <- pt*(1 - pt)
wght <- 1 / sigsq

coke.gls.trunc <- lm(coke ~ pratio + disp_coke + disp_pepsi,
                     weights = 1 / sigsq, data = coke)
coeftest(coke.gls.trunc)

# Elimimate values < 0 & > 1
pe <- p
pe[pe < 0.001 | pe > 0.999] <- NA

sigsq <- pe*(1 - pe)
wght <- 1 / sigsq

coke.gls.omit <- lm(coke ~ pratio + disp_coke + disp_pepsi,
                    weights = 1 / sigsq, data = coke)
coeftest(coke.gls.omit)
