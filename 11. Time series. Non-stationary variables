library(tseries) # for ADF unit root tests
library(dynlm)
library(lmtest) #for `coeftest()` and `bptest()`.
library(PoEdata) #for PoE4 datasets
library(sandwich)
library(forecast) 


## Binding time series data
data("usa", package = "PoEdata")

usa.ts <- ts(usa, start = c(1984, 1), end = c(2009, 4),
             frequency = 4)
Dgdp <- diff(usa.ts[, 1])
Dinf <- diff(usa.ts[, 2])
Df <- diff(usa.ts[, 3])
Db <- diff(usa.ts[, 4])
usa.ts.df <- ts.union(gdp = usa.ts[, 1],
                      inf = usa.ts[, 2],
                      f = usa.ts[, 3],
                      b = usa.ts[, 4],
                      Dgdp, Dinf, Df, Db,
                      dframe = TRUE)
par(mfrow = c(1, 2))

# GDP and its first difference
plot(usa.ts.df$gdp)
plot(usa.ts.df$Dgdp)

# Inf and its first difference
plot(usa.ts.df$inf)
plot(usa.ts.df$Dinf)

# Federal funds rate and its first difference
plot(usa.ts.df$f)
plot(usa.ts.df$Df)

# 3-year bond rate and its first difference
plot(usa.ts.df$b)
plot(usa.ts.df$Db)


## Simulating AR(1) process and random walks
N <- 1000
a <- 1
l <- 0.01
rho <- 0.7

set.seed(246810)
v <- ts(rnorm(N, 0, 1))

# AR(1) process
y <- ts(rep(0, N))

for (t in 2:N) {
  y[t] <- rho*y[t-1] + v[t]
}
plot(y, type = "l", ylab = "rho*y[t-1] + v[t]")
abline(h = 0)

par(mfrow = c(1, 1))

# AR(1) process with a constant term
y <- ts(rep(0, N))

for (t in 2:N) {
  y[t] <- a + rho*y[t-1] + v[t]
}
plot(y, type = "l", ylab = "a + rho*y[t-1] + v[t]")
abline(h = a/(1-rho))

# AR(1) process with constant and drift terms
y <- ts(rep(0, N))

for (t in 2:N) {
  y[t] <- a + l*time(y)[t] + rho*y[t-1] + v[t]
}
plot(y, type = "l", ylab = "a + l*time(y)[t] + rho*y[t-1] + v[t]")

# Simple random walk
y <- ts(rep(0, N))

for (t in 2:N) {
  y[t] <- y[t-1] + v[t]
}
plot(y, type = "l", ylab = "y[t-1] + v[t]")

# Random walk with a constant term
a <- 0.1
y <- ts(rep(0, N))

for (t in 2:N) {
  y[t] <- a + y[t-1] + v[t]
}
plot(y, type = "l", ylab = "a + y[t-1] + v[t]")

# Random walk with constand and drift terms
y <- ts(rep(0, N))

for (t in 2:N) {
  y[t] <- a + l*time(y)[t] + y[t-1] + v[t]
}
plot(y, type = "l", ylab = "a + l*time(y)[t] + y[t-1] + v[t]")


## Spurious regression

# Generating two random walks
T <- 1000
set.seed(1357)
y <- ts(rep(0, T))
vy <- ts(rnorm(T, 0, 1))

for (t in 2:T) {
  y[t] <- y[t-1] + vy[t]
}

set.seed(4365)
x <- ts(rep(0, T))
vx <- ts(rnorm(T, 0, 1))

for (t in 2:T) {
  x[t] <- x[t-1] + vx[t]
}
y <- ts(y[300:1000])
x <- ts(x[300:1000])
ts.plot(y, x, ylab = "y and x")

# Modeling the relationship
spurious.ols <- lm(y ~ x)
summary(spurious.ols)


## Unit root tests for stationarity
plot(usa.ts.df$f)

test.eq <- dynlm(diff(f) ~ L(f, 1) + L(diff(f), 1), data = usa.ts.df)
acf(test.eq$residuals)
tau <- coeftest(test.eq)[2, 3]

if (tau < 2.86) {
  print("The series is stationary")
} else {
  print("The series in nonstationary")
}

plot(usa.ts.df$b)

test.eq <- dynlm(diff(b) ~ L(b, 1) + L(diff(b), 1), data = usa.ts.df)
acf(test.eq$residuals)

tau <- coeftest(test.eq)[2, 3]
if (tau < 2.86) {
  print("The series is stationary")
} else {
  print("The series in nonstationary")
}


## Cointegration

# Testing whether federal funds rate and 3-year bond rate are I(1)
plot(usa.ts.df$Df)

test.eq <- dynlm(diff(Df) ~ 0 + L(Df, 1), data = usa.ts.df)
acf(test.eq$residuals)

tau <- coeftest(test.eq)[1, 3]
if (tau < 1.94) {
  print("The series is I(1)")
} else {
  print("The series is not I(1)")
}

plot(usa.ts.df$Db)

test.eq <- dynlm(diff(Db) ~ 0 + L(Db, 1), data = usa.ts.df)
acf(test.eq$residuals)

tau <- coeftest(test.eq)[1, 3]
if (tau < 1.94) {
  print("The series is I(1)")
} else {
  print("The series is not I(1)")
}

# Since both series are integrated of order 1, we can test whether they are cointegrated
# Test for cointegration is basically the same DF test, but applied to the regression errors
# However, since regression errors are unobservable, we estimate the model residuals, and
# therefore shall use special tau critical values for stationarity testing.
fb.dyn <- dynlm(b ~ f, data = usa.ts.df)
e_hat <- fb.dyn$residuals
plot(e_hat)

test.eq <- dynlm(diff(e_hat) ~ 0 + L(e_hat, 1) + diff(L(e_hat, 1)))
acf(test.eq$residuals)

tau <- coeftest(test.eq)[1, 3]
if (tau < -2.76) {
  print("Federal funds rate and 3-year bond rate are cointegrated")
} else {
  print("Regression with federal funds rate and 3-year bond rate is spurious")
}


## The error correction model for Federal funds and 3-year bond rates
b.ols <- dynlm(L(b) ~ L(f), data = usa.ts.df)
b1ini <- coef(b.ols)[[1]]
b2ini <- coef(b.ols)[[2]]

initial_ardl <- dynlm(b ~ L(b) + f + L(f), data = usa.ts.df)
aini <- 1 - coef(initial_ardl)[[2]]
d0ini <- coef(initial_ardl)[[3]]
d1ini <- coef(initial_ardl)[[4]]

b <- usa.ts.df$b
f <- usa.ts.df$f
Db <- diff(b)
Df <- diff(f)
Lb <- lag(b, k = -1)
Lf <- lag(f, k = -1)
LDf <- lag(Df, k = -1)
dataset <- data.frame(ts.union(cbind(b, f, Db, Df, Lb, Lf, LDf)))
formula <- Db ~ -a*(Lb - b1 - b2*Lf) + d0*Df + d1*LDf

ec_model <- nls(formula, na.action = na.omit, data = dataset,
                start = list(a = aini, b1 = b1ini, b2 = b2ini,
                             d0 = d0ini, d1 = d1ini))
coeftest(ec_model)
