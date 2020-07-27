library(tseries) # for `adf.test()`
library(dynlm) #for function `dynlm()`
library(nlWaldTest) # for the `nlWaldtest()` function
library(lmtest) #for `coeftest()` and `bptest()`.
library(PoEdata) #for PoE4 datasets
library(car) #for `hccm()` robust standard errors
library(sandwich)
library(forecast) 


## Estimating a VEC model
data("gdp", package = 'PoEdata')
gdp <- ts(gdp, start = c(1970, 1), end = c(2000, 4), frequency = 4)
ts.plot(gdp[, "usa"], gdp[, "aus"],
        lty = c(1, 2), col = c("black", "red"))
legend("topleft", legend = c("USA", "AUS"),
       lty = c(1, 2), col = c("black", "red"))

# Testing for stationarity and cointegration
gdp.ts.df <- ts.union(usa = gdp[, "usa"],
                      aus = gdp[, "aus"],
                      Dusa = diff(gdp[, "usa"]),
                      Daus = diff(gdp[, "aus"]),
                      dframe = TRUE)
test.eq <- dynlm(Dusa ~ time(usa) + L(usa, 1) + L(Dusa, 1:2), data = gdp.ts.df)
acf(test.eq$residuals)
tau <- coeftest(test.eq)[3, 3]

if (tau < -3.47) {
  print("USA is stationary")
} else {
  print("USA is non-stationary")
}

test.eq <- dynlm(Daus ~ time(aus) + L(aus, 1), data = gdp.ts.df)
acf(test.eq$residuals)
tau <- coeftest(test.eq)[3, 3]

if (tau < -3.47) {
  print("AUS is stationary")
} else {
  print("AUS is non-stationary")
}

plot(gdp.ts.df$Dusa)
test.eq <- dynlm(diff(Dusa) ~ L(Dusa, 1) + L(diff(Dusa), 1), data = gdp.ts.df)
acf(test.eq$residuals)
tau <- coeftest(test.eq)[2, 3]

if (tau < -2.86) {
  print("USA is I(1)")
} else {
  print("USA is not I(1)")
}

plot(gdp.ts.df$Daus)
test.eq <- dynlm(diff(Daus) ~ L(Daus, 1) + L(diff(Daus), 1:2), data = gdp.ts.df)
acf(test.eq$residuals)
tau <- coeftest(test.eq)[2, 3]

if (tau < -2.86) {
  print("AUS is I(1)")
} else {
  print("AUS is not I(1)")
}

aus.usa <- dynlm(aus ~ 0 + usa, data = gdp.ts.df)
coeftest(aus.usa)
e_hat <- aus.usa$residuals
plot(e_hat)

test.eq <- dynlm(diff(e_hat) ~ 0 + L(e_hat, 1))
acf(test.eq$residuals)
tau <- coeftest(test.eq)[1, 3]

if (tau < -2.76) {
  print("AUS and USA are cointegrated")
} else {
  print("Relationship between AUS and USA is spurious")
}

# Constructing the VEC model
vec.aus <- dynlm(diff(aus) ~ L(e_hat), data = gdp)
vec.usa <- dynlm(diff(usa) ~ L(e_hat), data = gdp)
coeftest(vec.aus)
coeftest(vec.usa)


## Estimating a VAR model
data("fred", package = "PoEdata")
fred <- ts(fred, start = c(1960, 1), end = c(2009, 4), frequency = 4)
fred.ts.df <- ts.union(c = fred[, "c"],
                       y = fred[, "y"],
                       Dc = diff(fred[, "c"]),
                       Dy = diff(fred[, "y"]),
                       dframe = TRUE)
ts.plot(fred[, "c"], fred[, "y"],
        lty = c(1, 2), col = c("black", "red"))
legend("topleft", legend = c("c", "y"),
       lty = c(1, 2), col = c("black", "red"))

# Assuming the series are I(1), we test for the cointegrating relationship
c.y <- dynlm(c ~ y, data = fred.ts.df)
e_hat <- c.y$residuals
plot(e_hat)

test.eq <- dynlm(diff(e_hat) ~ 0 + L(e_hat, 1) + L(diff(e_hat), 1))
acf(test.eq$residuals)
tau <- coeftest(test.eq)[1, 3]

if (tau < -3.37) {
  print("c and y are cointegrated")
} else {
  print("Relationship between c and y is spurious")
}

# Having confirmed that c and y are not cointegrated, we can proceed with a VAR model in first differences
library(vars) # for function "VAR()"
Dc <- na.omit(fred.ts.df$Dc)
Dy <- na.omit(fred.ts.df$Dy)

var.matrix <- as.matrix(cbind(Dc, Dy))
var.fit <- VAR(var.matrix)
summary(var.fit)


## Impulse responses and variance decompositions
impresp <- irf(var.fit)
plot(impresp)
plot(fevd(var.fit))
