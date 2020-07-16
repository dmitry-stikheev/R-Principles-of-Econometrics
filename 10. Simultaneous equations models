library(systemfit)
library(broom) #for `glance(`) and `tidy()`
library(PoEdata) #for PoE4 dataset
library(knitr) #for kable()


## Model of the supply and demand of truffles, where
# Qd = alpha1 + alpha2*p + alpha3*ps + alpha4*di + ed
# Qs = beta1 + beta2*p + beta3*pf + es
data("truffles", package = "PoEdata")

D <- q ~ p + ps + di
S <- q ~ p + pf
sys <- list(D, S)
instr <- ~ ps + di + pf

truff.sys <- systemfit(sys, inst = instr, method = "2SLS", data = truffles)
summary(truff.sys)

# Modeling the quantity demanded and supplied through the reduced forms
Q.red <- lm(q ~ ps + di + pf, data = truffles)
P.red <- lm(p ~ ps + di + pf, data = truffles)
coeftest(Q.red)
coeftest(P.red)
# Testing for weak instruments, using the rule of thumb (F > 10)
SST <- var(truffles$q)*(nrow(truffles) - 1)
SSE <- anova(Q.red)[4, 2]
df_n <- 3
df_d <- df.residual(Q.red)

Fst <- ((SST - SSE) / df_n) / (SSE / df_d)

if (Fst > 10) {
  print("Intruments are strong")
} else {
  print("Intruments are weak")
}
head(truffles, 1)

## Another model of supply and demand
# log(quan) = alpha1 + alpha2*price + alpha3*mon + alpha4*tue + alpha5*wed + alpha6*thu + ed
# log(quan) = beta1 + beta2*log(price) + beta3*stormy + es
data("fultonfish", package = "PoEdata")

fishQ.ols <- lm(lquan ~ mon + tue + wed + thu + stormy, data = fultonfish)
coeftest(fishQ.ols)

fishP.ols <- lm(lprice ~ mon + tue + wed + thu + stormy, data = fultonfish)
coeftest(fishP.ols)

# Estimating structural equations
fish.D <- lquan ~ lprice + mon + tue + wed + thu
fish.S <- lquan ~ lprice + stormy
sys <- list(fish.D, fish.S)
ivs <- ~ mon + tue + wed + thu + stormy

fish.sys <- systemfit(sys, inst = ivs, method = "2SLS", data = fultonfish)
summary(fish.sys)

### The estimated supply equation is not reliable due to the insignificance of weekday coefficients
### which are in charge for shifting the supply curve. Thus, the model's IVs should be reconsidered.
