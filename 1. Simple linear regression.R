library(PoEdata)
library(ggplot2)
library(dplyr)
data("cps_small")
head(cps_small)


# Plotting education agains wage
plot(cps_small$educ, cps_small$wage, 
    xlab = "Education", ylab = "wage")


# INCOME and FOOD_EXP data
data("food")

plot(food$income, food$food_exp, type = "p",
    xlab = "Weekly income in $100", ylab = "Weekly food expenditure in $",
      xlim = c(0, max(food$income)), ylim = c(0, max(food$food_exp)))

install.packages("ggplot2")
library(ggplot2)

ggplot(food, aes(x = food$income, y = food$food_exp)) +
  geom_point() + geom_smooth(method = lm, se = FALSE) +
  xlab("Income, in $hdr") + ylab("Food Expenditure")

# Estimating a linear regression
sd <- sqrt((food$food_exp - mean(food$food_exp))^2)
w <- 1 / sd
mod <- lm(food_exp ~ income, data = food, weights = w)
b1 <- coef(mod)[1]
b2 <- coef(mod)[2]
smod_1 <- summary(mod)
curve(b1 + b2*x, add = TRUE, col = "red", lwd = 2)
abline(b1, b2, col = "blue", lwd = 2)
names(mod)
mod_1$coefficients
smod_1$coefficients

food_exp_hat <- NULL
for (i in 1:length(food$food_exp)) {
  food_exp_hat[i] <- b1 + b2*food$income[i]
}
SST <- sum((food$food_exp - mean(food$food_exp))^2)
SSR <- sum((food_exp_hat - mean(food$food_exp))^2)
SSE <- sum((food$food_exp - food_exp_hat)^2)
r_sqr <- SSR/SST
smod_1
  

# Prediction with the linear regression model
new_x <- data.frame(income = c(20, 25, 27))
y_hat <- predict(mod_1, new_x)
names(y_hat) <- c("income = $2000", "$2500", "$2700")
(y_hat)


# Repeated samples to assess regression coefficients
N <- nrow(food)   # returns the number of observations in the dataset
C <- 50           # desired number of subsamples
S <- 38           # desired sample size
sum_b2 <- 0
for (i in 1:C) {  # a loop over the number of subsamples
  set.seed(3*i)   # a different seed for each subsample
  subsample <- food[sample(1:N, size = S, replace = TRUE), ]
  mod_2 <- lm(food_exp ~ income, data = subsample)
  sum_b2 <- sum_b2 + coef(mod_2)[[2]]
}
print(sum_b2/C, digits = 3)


# Estimated variances and covariances of the regression coefficients
var_b1_hat <- vcov(mod_1)[1, 1]
var_b2_hat <- vcov(mod_1)[2, 2]
cov_b1_b2_hat <- vcov(mod_1)[1, 2]


# Non-linear relationships
data("br")
mod_3 <- lm(price ~ I(sqft^2), data = br)
plot(br$sqft, br$price, 
     xlab = "Total square feet", ylab = "Sale price, $", col = "grey")
b1 <- coef(mod_3)[1]
b2 <- coef(mod_3)[2]
curve(b1 + b2*x^2, col = "red", add = TRUE)

par(mfrow = c(1, 1))
hist(br$price, col = "grey")
hist(log(br$price), col = "grey")
mod_4 <- lm(log(price) ~ sqft, data = br)
b1 <-coef(mod_4)[1]
b2 <- coef(mod_4)[2]
plot(br$sqft, br$price, 
    xlab = "Total square feet", ylab = "Sale price, $", col = "grey")
curve(exp(b1 + b2*x), col = "blue", add = TRUE)


# Using indicator variables in a regression
head(utown)
price0_bar <- mean(utown$price[which(utown$utown == 0)])
price1_bar <- mean(utown$price[which(utown$utown == 1)])
price_bar <- c(price0_bar, price1_bar)
names(price_bar) <- c("0", "1")
price_bar

mod_5 <- lm(price ~ utown, data = utown)
b1 <- coef(mod_5)[[1]]
b2 <- coef(mod_5)[[2]]


# Monte Carlo Simulation
N <- 40
x1 <- 10
x2 <- 20
b1 <- 100
b2 <- 10
mu <- 0
sig2_e <- 2500
sd_e <- sqrt(sig2_e)
y_hat1 <- b1 + b2*x1
y_hat2 <- b1 + b2*x2
curve(dnorm(x, mean = y_hat1, sd = sd_e), 0, 500, col = "blue")
curve(dnorm(x, mean = y_hat2, sd = sd_e), 0, 500, add = T, col = "red")
abline(v = y_hat1, col = "blue", lty = 2)
abline(v = y_hat2, col = "red", lty = 2)
legend("topright", legend = c("f(y|x=10)", "f(y|x=20"),
       lty = 1, col = c("blue", "red"))

x <- c(rep(x1, N/2), rep(x2, N/2))
x_bar <- mean(x)
sumx2 <- sum((x - x_bar)^2)
var_b2 <- sig2_e / sumx2
sd_b2 <- sqrt(var_b2)
left_lim <- b2 - 3*sd_b2
right_lim <- b2 + 3*sd_b2
curve(dnorm(x, mean = b2, sd = sd_b2, ), left_lim, right_lim)
abline(v = b2, lty = 2)

set.seed(12345)
y <- b1 + b2*x + rnorm(N, mean = 0, sd = sd_e)
mod_6 <- lm(y ~ x)
b1_hat <- coef(mod_6)[[1]]
b2_hat <- coef(mod_6)[[2]]
mod_6_summary <- summary(mod_6)
coef(mod_6_summary)
seb2_hat <- coef(mod_6_summary)[2, 2]

# enjoy!!!
