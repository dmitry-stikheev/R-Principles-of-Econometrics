library(PoEdata)
library(ggplot2)
library(dplyr)


# Confidence intervals in the food model
data("food")
alpha <- 0.05                # Specified significance level
mod_1 <- lm(food_exp ~ income, data = food)
b2 <- coef(mod_1)[[2]]
df <- df.residual(mod_1)     # degrees of freedom
smod_1 <- summary(mod_1)
se_b2 <- coef(smod_1)[2, 2]
tc <- qt(1 - alpha / 2, df)
low_b <- b2 - tc*se_b2
upp_b <- b2 + tc*se_b2

# Another approach through the built-in function
confint(mod_1)
lowb_b2 <- confint(mod_1)[2, 1]
uppb_b2 <- confint(mod_1)[2, 2]


# Confidence intervals in repeated samples
data("table2_2")
alpha <- 0.05
mod_1 <- lm(y1 ~ x, data = table2_2)    # just to determine df
df <- df.residual(mod_1)
tc <- qt(1 - alpha / 2, df)

# Initiating four vectors that will store the results
lowb_1 <- rep(0, 10)                     # repeat 0 ten times
uppb_1 <- rep(0, 10)                     # (alternatively, 'numeric(10)')
lowb_2 <- rep(0, 10)
uppb_2 <- rep(0, 10)

# One loop for each set of income
for (i in 2:ncol(table2_2)) {
  dat <- data.frame(cbind(table2_2[, 1], table2_2[, i]))
  names(dat) <- c("x", "y")
  mod_1 <- lm(y ~ x, data = dat)
  smod_1 <- summary(mod_1)
  b1 <- coef(mod_1)[[1]]
  b2 <- coef(mod_1)[[2]]
  se_b1 <- coef(smod_1)[1, 2]
  se_b2 <- coef(smod_1)[2, 2]
  lowb_1[i] <- b1 - tc*se_b1
  uppb_1[i] <- b1 + tc*se_b1
  lowb_2[i] <- b2 - tc*se_b2
  uppb_2[i] <- b2 + tc*se_b2
}
table <- data.frame(lowb_1, uppb_1, lowb_2, uppb_2)
kable(table, caption = "Confidence intervals for $b_{1}$ and $b_{2}$", 
      align = "c")


# Hypothesis tests
data("food")

alpha <- 0.05
mod_1 <- lm(food_exp ~ income, data = food)
smod_1 <- summary(mod_1)
coef(smod_1)
b2 <- coef(smod_1)[2, 1]
se_b2 <- sqrt(vcov(mod_1)[2, 2])
df <- df.residual(mod_1)
t <- b2 / se_b2
tc <- qt(1 - alpha / 2, df)
t > tc

# Plotting the density function and the values of t
curve(dt(x, df), -2.5*se_b2, 2.5*se_b2, xlab = "t", ylab = "")
abline(v = c(-tc, tc, t), col = c("red", "red", "blue"), lty = c(2, 2, 3))
legend("topleft", legend = c("-tc", "tc", "t"),
       lty = c(2, 2, 3), col = c("red", "red", "blue"))

# Testing the null hypothesis that beta2 is less than 5.5
c <- 5.5
aplha <- 0.05
t <- (b2 - c) / se_b2
tc <- qt(1 - aplha, df)
curve(dt(x, df), -2.5*se_b2, 2.5*se_b2, xlab = "t", ylab = "")
abline(v = c(tc, t), col = c("red", "blue"), lty = c(2, 3))
legend("topleft", legend = c("tc", "t"),
       lty = c(2, 3), col = c("red", "blue"))
t > tc

# Testing the null hypothesis that beta2 is greater than 15
c <- 15
alpha <- 0.05
t <- (b2 - c) / se_b2
tc <- qt(alpha, df)
curve(dt(x, df), -2.5*se_b2, 2.5*se_b2, xlab = "t", ylab = "")
abline(v = c(t, tc), col = c("blue", "red"), lty = c(3, 2))
legend("topleft", legend = c("t", "tc"),
       lty = c(3, 2), col = c("blue", "red"))


# The p-value

# Right-tail test,  H0: ß2 < c, A: ß2 > c
c <- 5.5
t <- (b2 - c) / se_b2

p <- 1 - pt(t, df)

# Left-tail test
c <- 15
t <- (b2 - c) / se_b2

p <- pt(t, df)

# Two-tail test
c <- 0
t <- (b2 - c) / se_b2

p <- 2*(1 - pt(abs(t), df))


# Testing linear combinations of parameters
alpha <- 0.05
x <- 20
mod_1 <- lm(food_exp ~ income, data = food)
tc <- qt(1 - alpha / 2, df)
df <- df.residual(mod_1)
b1 <- coef(mod_1)[[1]]
b2 <- coef(mod_1)[[2]]
var_b1 <- vcov(mod_1)[1, 1]
var_b2 <- vcov(mod_1)[2, 2]
cov_b1b2 <- vcov(mod_1)[1, 2]
L <- b1 + b2*x
var_L <- var_b1 + x^2*var_b2 + 2*x*cov_b1b2
se_L <- sqrt(var_L)
low_bL <- L - tc*se_L
upp_bL <- L + tc*se_L
low_bL
upp_bL

# enjoy!!!
