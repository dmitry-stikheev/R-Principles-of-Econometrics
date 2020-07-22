library(dynlm) #for the `dynlm()` function
library(orcutt) # for the `cochrane.orcutt()` function
library(nlWaldTest) # for the `nlWaldtest()` function
library(pdfetch) # for retrieving data (just mentioned here)
library(lmtest) #for `coeftest()` and `bptest()`.
library(PoEdata) #for PoE4 datasets
library(car) #for `hccm()` robust standard errors
library(sandwich)
library(forecast) 


## Finite distributed lag model (FDL)
data("okun", package = "PoEdata")

check.ts <- is.ts(okun) # checking whether okun is a time series object or not
okun.ts <- ts(okun, start = c(1985, 2), end = c(2009, 3), frequency = 4)
?lag
okun.ts.tab <- cbind(okun.ts,
                 lag(okun.ts[, 2], k = -1),
                 diff(okun.ts[, 2], lag = 1),
                 lag(okun.ts[, 1], k = -1),
                 lag(okun.ts[, 1], k = -2),
                 lag(okun.ts[, 1], k = -3))
colnames(okun.ts.tab) <- c("g", "u", "uL1", "du", "gL1", "gL2", "gL3")

okunL3.dyn <- dynlm(d(u) ~ L(g, 0:3), data = okun.ts)
coeftest(okunL3.dyn)

okunL2.dyn <- dynlm(d(u) ~ L(g, 0:2), data = okun.ts)
coeftest(okunL2.dyn)


## Serial correlation
par(mfrow = c(1, 2))

# Detecting autocorrelation through scatterplots
plot(okun.ts[, "g"], ylab = "growth")
plot(okun.ts[, "u"], ylab = "unemployment")

ggL1 <- data.frame(cbind(okun.ts[, "g"], lag(okun.ts[, "g"], k = -1)))
colnames(ggL1) <- c("g", "gL1")

plot(ggL1)
abline(v = mean(ggL1$g, na.rm = TRUE), lty = 2)
abline(h = mean(ggL1$gL1, na.rm = TRUE), lty = 2)


ggL2 <- data.frame(cbind(okun.ts[, "g"], lag(okun.ts[, "g"], k = -2)))
colnames(ggL2) <- c("g", "gL2")

plot(ggL2)
abline(v = mean(ggL2$g, na.rm = TRUE), lty = 2)
abline(h = mean(ggL2$gL2, na.rm = TRUE), lty = 2)

# Testing for autocorrelation using correlogram
par(mfrow = c(1, 1))

growth_rate <- okun.ts[, "g"]
acf(growth_rate)

data("phillips_aus", package = "PoEdata")
phill.ts <- ts(phillips_aus,
            start = c(1987, 1),
            end = c(2009, 3),
            frequency = 4)
inflation <- phill.ts[, "inf"]
Du <- diff(phill.ts[, "u"], lag = 1)

par(mfrow = c(1, 2))
plot(inflation)
plot(Du)

phill.dyn <- dynlm(inf ~ diff(u), data = phill.ts)
coeftest(phill.dyn)

ehat <- phill.dyn$residuals

par(mfrow = c(1, 1))
plot(ehat)
abline(h = mean(ehat), lty = 2)

acf(ehat)

# Testing for autocorrelation through LM test
ehat <- phill.dyn$residuals
Lehat <- lag(phill.dyn$residuals, k = -1)

model <- dynlm(ehat ~ diff(phill.ts[, "u"]) + L(ehat))

Fst <- summary(model)$fstatistic[[1]]
df_n <- 2
df_d <- df.residual(model)
Fcr <- qf(0.95, df_n, df_d)

if (Fst > Fcr) {
  print("Autocorrelation of residuals is detected")
} else {
  print("Autocorrelation of residuals is undetected")
}


## Estimation with serially correlated errors
s0 <- coeftest(phill.dyn)
s1 <- coeftest(phill.dyn, vcov. = vcovHAC(phill.dyn))
s2 <- coeftest(phill.dyn, vcov. = NeweyWest(phill.dyn))
s3 <- coeftest(phill.dyn, vcov. = kernHAC(phill.dyn))

tbl <- cbind(s0[, 2], s1[, 2], s2[, 2], s3[, 2])
colnames(tbl) <- c("Incorrect", "vcovHAC", "NeweyWest", "kernHAC")


## Non-linear least squares estimation

# Initial linear model
phill.dyn <- dynlm(inf ~ diff(u), data = phill.ts)
head(phill.ts)
# Non-linear AR(1) model with Cochrane-Orcutt method
phill.ts.tab <- cbind(phill.ts[, "inf"],
                      phill.ts[, "u"],
                      lag(phill.ts[, "inf"], k = -1),
                      diff(phill.ts[, "u"], lag = 1),
                      lag(diff(phill.ts[, "u"], lag = 1), k = -1))
phill.dfr <- data.frame(phill.ts.tab)
colnames(phill.dfr) <- c("inf", "u", "Linf", "Du", "LDu")

phill.nls <- nls(inf ~ b1*(1 - rho) + b2*Du + rho*Linf - 
                      rho*b2*LDu,
                      data = phill.dfr,
                      start = c(rho = 0.5, b1 = 0.5, b2 = -0.5))
coeftest(phill.nls)


## Generalized least squares estimation
s.nls <- summary(phill.nls)

phill.gen <- dynlm(inf ~ L(inf) + diff(u) + L(diff(u)),
                   data = phill.ts)
s.gen <- summary(phill.gen)

nlW <- nlWaldtest(phill.gen, texts = "b[4] = -b[2]*b[3]")
nlW


## Autoregressive models
data("okun", package = "PoEdata")
okun.ts <- ts(okun)

okun.ar2 <- dynlm(g ~ L(g) + L(g, 2), data = okun.ts)
coeftest(okun.ar2)

res.ar2 <- okun.ar2$residuals
Acf(res.ar2, lag.max = 12)

# AR model specification
aics <- rep(0, 5)
bics <- rep(0, 5)
y <- okun.ts[, "g"]

for (i in 1:5) {
  ari <- dynlm(y ~ L(y, 1:i), start = i)
  aics[i] <- AIC(ari)
  bics[i] <- BIC(ari)
}

tbs <- data.frame(rbind(aics, bics))
colnames(tbs) <- c("AR(1)", "AR(2)", "AR(3)", "AR(4)", "AR(5)")
