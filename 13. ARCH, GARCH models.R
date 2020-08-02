library(FinTS) #for function `ArchTest()`
library(rugarch) #for GARCH models
library(tseries) # for `adf.test()`
library(dynlm) #for function `dynlm()`
library(vars) # for function `VAR()`
library(lmtest) #for `coeftest()` and `bptest()`.
library(PoEdata) #for PoE4 datasets
library(sandwich)
library(forecast) 


## The ARCH model
data("byd", package = "PoEdata")
r.ts <- ts(byd$r)
par(mfrow = c(2, 1))
ts.plot(r.ts)
hist(r.ts, main = "", breaks = 20, freq = FALSE, col = "grey")

# Testing for ARCH effects
byd.mean <- dynlm(r.ts ~ 1)
coeftest(byd.mean)

ehat.sq <- ts(byd.mean$residuals^2)
test.eq <- dynlm(ehat.sq ~ L(ehat.sq, 1))
acf(test.eq$residuals)

T <- nobs(byd.mean)
q <- length(coef(test.eq)) - 1
Rsq <- summary(test.eq)$r.squared
LM <- (T - q)*Rsq
alpha <- 0.05
chisqcr <- qchisq(1 - alpha, q)

if (LM > chisqcr) {
  print("Time series has ARCH effects")
} else {
  print("Time series has no ARCH effects")
}

# Estimating an ARCH model for returns volatility
byd.arch <- garch(r.ts, order = c(0, 1))
coeftest(byd.arch)
h.hat <- ts(byd.arch$fitted.values[-1, 1]^2)
ts.plot(h.hat)


## The GARCH model
garchSpec <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(0, 0)),
  distribution.model = "std")
garchFit <- ugarchfit(spec = garchSpec, data = r.ts)
coef(garchFit)

r.hat <- garchFit@fit$fitted.values
h.hat <- ts(garchFit@fit$sigma^2)

par(mfrow = c(1, 2))
ts.plot(r.hat)
ts.plot(h.hat)


## Threshhold GARCH model (TGARCH)
garchSpec <- ugarchspec(
  variance.model = list(model = "fGARCH",
                        garchOrder = c(1, 1),
                        submodel = "TGARCH"),
  mean.model = list(armaOrder = c(0, 0)),
  distribution.model = "std")
garchFit <- ugarchfit(garchSpec, data = r.ts)

r.hat <- garchFit@fit$fitted.values
h.hat <- ts(garchFit@fit$sigma^2)

ts.plot(r.hat)
ts.plot(h.hat)


## GARCH-in-mean
garchSpec <- ugarchspec(
  variance.model = list(model = "fGARCH",
                        garchOrder = c(1, 1),
                        submodel = "APARCH"),
  mean.model = list(armaOrder = c(0, 0),
                    include.mean = TRUE,
                    archm = TRUE,
                    archpow = 2),
  distribution.model = "std")
garchFit <- ugarchfit(garchSpec, data = r.ts)
coef(garchFit)

r.hat <- garchFit@fit$fitted.values
h.hat <- ts(garchFit@fit$sigma^2)

ts.plot(r.hat)
ts.plot(h.hat)
