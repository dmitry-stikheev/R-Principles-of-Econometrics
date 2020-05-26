library(PoEdata)
library(ggplot2)
library(dplyr)


# Forecasting (predicting a particular value)
data("food")
alpha <- 0.05
x <- 20
xbar <- mean(food$income)
m1 <- lm(food_exp ~ income, data = food)
sm1 <- summary(m1)
b1 <- coef(sm1)[1, 1]
b2 <- coef(sm1)[2, 1]
yhat_x <- b1 + b2*x
df <- df.residual(m1)
tcr <- qt(1 - alpha / 2, df)
N <- nrow(food)
varb2 <- vcov(m1)[2, 2]
var_ehat <- sm1$sigma^2
varf <- var_ehat + (var_ehat / N) + ((x - xbar)^2)*varb2
sef <- sqrt(var_ehat)
lb <- yhat_x - tcr*sef
ub <- yhat_x + tcr*sef


sefv <- sqrt(var_ehat + var_ehat / N + (food$income - xbar)^2*varb2)
yhatv <- m1$fitted.values
lbv <- yhatv - tcr*sefv
ubv <- yhatv + tcr*sefv
xincome <- food$income
dplot <- data.frame(xincome, yhatv, lbv, ubv)
xmax <- max(dplot$xincome)
xmin <- min(dplot$xincome)
ymax <- max(dplot$ubv)
ymin <- min(dplot$lbv)
plot(dplot$xincome, dplot$yhatv, xlim = c(xmin, xmax), ylim = c(ymin, ymax),
     xlab = "income", ylab = "food expenditure", type = "l")
lines(dplot$ubv ~ dplot$xincome, lty = 2)
lines(dplot$lbv ~ dplot$xincome, lty = 2)
points(food$income, food$food_exp)

xmin <- min(food$income)
xmax <- max(food$income)
income <- seq(from = xmin, to = xmax)
ypredict <- predict(m1, newdata = data.frame(income), 
                    interval = "confidence")
yforecast <- predict(m1, newdata = data.frame(income),
                     interval = "prediction")
matplot(income, cbind(ypredict[, 1], ypredict[, 2], ypredict[, 3],
    yforecast[, 2], yforecast[, 3]), type = "l", lty = c(1, 2, 2, 3, 3),
    col = c("black", "red", "red", "blue", "blue"),
    xlab = "income", ylab = "food expenditure")
points(food$income, food$food_exp)
legend("topleft", 
       legend = c("E[y|x]", "lwr_pred", "upr_pred", "lwr_forcst", "upr_forcst"),
       lty = c(1, 2, 2, 3, 3), col = c("black", "red", "red", "blue", "blue"))


# Goodness of fit
rsq <- sm1$r.squared
rsq
anova(m1)


# Linear-log models
mod2 <- lm(food_exp ~ log(income), data = food)
smod2 <- summary(mod2)
coef(smod2)
b1 <- coef(smod2)[1, 1]
b2 <- coef(smod2)[2, 1]
confmod2 <- predict(mod2, 
        newdata = data.frame(income), interval = "confidence")
plot(food$income, food$food_exp, xlab = "income", ylab = "food expenditure")
lines(confmod2[, 1] ~ income, lty = 1, col = "black")
lines(confmod2[, 2] ~ income, lty = 2, col = "red")
lines(confmod2[, 3] ~ income, lty = 2, col = "red")

data("food")
m1 <- lm(food_exp ~ income, data = food)
sm1 <- summary(m1)
N <- nrow(food)
b1 <- coef(sm1)[1, 1]  
b2 <- coef(sm1)[2, 1]  
varb2 <- vcov(m1)[2, 2]  
vare <- sm1$sigma^2
alpha <- 0.05  
df <- df.residual(m1)  
tcr <- qt(1 - alpha / 2, df)  
xmin <- min(food$income)
xmax <- max(food$income)
xbar <- mean(food$income)
income <- seq(from = xmin, to = xmax)
yhat <- predict(m1, newdata = data.frame(income))
var_yhat <- vare / N + (income - xbar)^2*varb2   
se_yhat <- sqrt(var_yhat)  
lb_conf <- yhat - tcr*se_yhat
ub_conf <- yhat + tcr*se_yhat  
plot(income, yhat, xlim = c(xmin, xmax), ylim = c(min(lb_pr), max(ub_pr)), type = "l")
lines(lb_conf ~ income, col = "red", lty = 3)
lines(ub_conf ~ income, col = "red", lty = 3)  
points(food$income, food$food_exp)
varf <- vare + vare / N + (income - xbar)^2*varb2
sef <- sqrt(varf)
lb_pr <- yhat - tcr*sef
ub_pr <- yhat + tcr*sef
lines(lb_pr ~ income, col = "blue", lty = 2)
lines(ub_pr ~ income, col = "blue", lty = 2)


# Residuals and diagnostics
ehat <- mod2$residuals
plot(food$income, ehat, xlab = "income", ylab = "residuals", bg = "blue")


# Data generating process: Yi = 1 + Xi + Ei, where
# Xi is drawn from uniform distribution and Ei from standard normal
set.seed(12345)
x <- runif(300, 0, 10)
e <- rnorm(300, 0, 1)
y <- 1 + x + e
mod3 <- lm(y ~ x)
ehat <- mod3$residuals
plot(x, ehat, xlab = "x", ylab = "residuals")

set.seed(12345)
x <- runif(1000, -2.5, 2.5)
e <- rnorm(1000, 0, 4)
y <- 15 - 4*x^2 + e
mod3 <- lm(y ~ x)
ehat <- mod3$residuals
plot(x, ehat, xlab = "x", ylab = "residuals")

library(tseries)
mod1 <- lm(food_exp ~ income, data = food)
ehat <- residuals(mod1)
ebar <- mean(ehat)
sde <- sd(ehat)
hist(ehat, col = "grey", freq = FALSE, main = "",
     xlab = "ehat", ylab = "density")
curve(dnorm(x, ebar, sde), add = TRUE, col = "red")
jarque.bera.test(ehat)


# Polynomial models
data("wa_wheat")
mod1 <- lm(greenough ~ time, data = wa_wheat)
ehat <- mod1$residuals
plot(wa_wheat$time, ehat, xlab = "time", ylab = "residuals")

mod2 <- lm(greenough ~ I(time^3), data = wa_wheat)
ehat <- mod2$residuals
plot(wa_wheat$time, ehat, xlab = "time", ylab = "residuals")


# Log-linear models
mod4 <- lm(log(greenough) ~ time, data = wa_wheat)
smod4 <- summary(mod4)
coefficients(smod4)

data("cps4_small")
xeduc <- 12

mod5 <- lm(log(wage) ~ educ, data = cps4_small)
smod5 <- summary(mod5)
coef(smod5)
plot(cps4_small$educ, cps4_small$wage, col = "grey")
b1 <- coef(smod5)[1, 1]
b2 <- coef(smod5)[2, 1]
vare <- smod5$sigma^2
curve(exp(b1 + b2*x), add = TRUE, col = "blue", lwd = 2)
curve(exp(b1 + b2*x + vare / 2), add = TRUE, col = "red", lwd = 2)
legend("topleft", legend = c("predictor", "corrected predictor"),
      lty = c(1, 1) ,col = c("blue", "red"))

# Prediction interval for educ = 12
alpha <- 0.05
xeduc <- 12
xedbar <- mean(cps4_small$educ)
mod5 <- lm(log(wage) ~ educ, data = cps4_small)
smod5 <- summary(mod5)
b1 <- coef(smod5)[1, 1]
b2 <- coef(smod5)[2, 1]
df5 <- df.residual(mod5)
N <- nobs(mod5)
tcr <- qt(1 - alpha / 2, df5)
xmin <- min(cps4_small$educ)
xmax <- max(cps4_small$educ)
educ <- seq(from = xmin, to = xmax)
yhat <- exp(b1 + b2*educ)
varb2 <- vcov(mod5)[2, 2]
vare <- smod5$sigma^2
varf <- vare + vare / N + (educ - xedbar)^2*varb2
sef <- sqrt(varf)
lnyhat <- b1 + b2*educ
lb <- exp(lnyhat - tcr*sef)
ub <- exp(lnyhat + tcr*sef)
ymin <- min(lb)
ymax <- max(ub)
plot(educ, yhat, xlim = c(xmin, xmax), ylim = c(ymin, ymax), type = "l")
points(cps4_small$educ, cps4_small$wage, col = "grey")
lines(lb ~ educ, col = "blue", lty = 2)
lines(ub ~ educ, col = "blue", lty = 2)


data("newbroiler")
mod6 <- lm(log(q) ~ log(p), data = newbroiler)
smod6 <- summary(mod6)
b1 <- coef(smod6)[1, 1]
b2 <- coef(smod6)[2, 1]
vare <- smod6$sigma^2
coef(smod6)
plot(newbroiler$p, newbroiler$q, ylim = c(8, 60))
curve(exp(b1 + b2*log(x)), add = TRUE, col = "blue")
curve(exp(b1 + b2*log(x) + vare / 2), add = TRUE, col = "red")



# Prediction interval for log-log model
alpha <- 0.05
df <- df.residual(mod6)
tcr <- qt(1 - alpha / 2, df)
N <- nobs(mod6)
b1 <- coef(smod6)[1, 1]
b2 <- coef(smod6)[2, 1]
vare <- smod6$sigma^2
varb2 <- vcov(mod6)[2, 2]
xmin <- min(newbroiler$p)
xmax <- max(newbroiler$p)
p <- seq(from = xmin, to = xmax, by = 0.1)
lnp <- log(p)
ln_pbar <- mean(log(newbroiler$p))
varf <- vare + vare / N + (lnp - ln_pbar)^2*varb2
sef <- sqrt(varf)
lnyhat <- b1 + b2*lnp
yhat <- exp(lnyhat)
lb <- exp(lnyhat - tcr*sef)
ub <- exp(lnyhat + tcr*sef)
ymin <- min(lb)
ymax <- max(ub)
plot(p, yhat, type = "l", xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = "q")
lines(lb ~ p, lty = 2, col = "red")
lines(ub ~ p, lty = 2, col = "red")
points(newbroiler$p, newbroiler$q, col = "grey")
