## Installing required packages
library(bookdown)
library(PoEdata)    #for PoE datasets
library(knitr)      #for referenced tables with kable()
library(xtable)     #makes data frame for kable
library(printr)     #automatically prints output nicely
library(effects)
library(car)
library(AER)
library(broom)      #for tidy lm output and function glance()
library(lmtest)     #for coeftest() and other test functions 
library(ggplot2)    #neat visualizing tool
library(dplyr)
## Indicator variables
data("utown")

# Transforming indicator variables into factors
utown$utown <- as.factor(utown$utown)
utown$pool <- as.factor(utown$pool)
utown$fplace <- as.factor(utown$fplace)
summary.data.frame(utown)

# Building the house prices model
mod4 <- lm(price ~ utown + sqft + utown*sqft + age + pool + fplace, data = utown)
smod4 <- summary(mod4)
coef(smod4)

bsqft <- 1000*coef(mod4)[["sqft"]]
bsqft1 <- bsqft + 1000*coef(mod4)[["utown1:sqft"]]

# Building the wage discrimination model
data("cps4_small")

mod5 <- lm(wage ~ educ + black + female + black*female, data = cps4_small)
smod5 <- summary(mod5)
coef(smod5)

delta1 <- coef(smod5)[3, 1]
delta2 <- coef(smod5)[4, 1]
gamma <- coef(smod5)[5, 1]

blfm <- delta1 + delta2 + gamma

# Testing the joint hypothesis whether at least race or sex indicator coefficient is significant
unrestricted <- lm(wage ~ educ + black + female + black*female, data = cps4_small)
sunrestricted <- summary(unrestricted)
SSEun <- sunrestricted$sigma^2*(nobs(unrestricted) - 5)

restricted <- lm(wage ~ educ, data = cps4_small)
srestricted <- summary(restricted)
SSEr <- srestricted$sigma^2*(nobs(restricted) - 2)

alpha <- 0.05
df_n <- 3
df_d <- df.residual(unrestricted)
Fc <- qf(1 - alpha, df_n, df_d)

Fst <- ((SSEr - SSEun) / df_n) / (SSEun / df_d)
if(Fst > Fc) {
  print("Reject Ho")
}  else {
  print("Do not reject Ho")
}


## Comparing two regressions: The Chow Test

# Dividing the data into two sumsamples based on a dummy variable
dnosouth <- cps4_small %>%
  filter(south == 0)
dsouth <- cps4_small %>%
  filter(south == 1)

# Constructing regression models
mod5ns <- lm(wage ~ educ + black + female + black*female, data = dnosouth)
mod5s <- lm(wage ~ educ + black + female + black*female, data = dsouth)
mod5 <- lm(wage ~ educ + black + female + black*female, data = cps4_small)

smod5ns <- summary(mod5ns)
smod5s <- summary(mod5s)
smod5 <- summary(mod5)

SSEns <- smod5ns$sigma^2*(nobs(mod5ns) - 5)
SSEs <- smod5s$sigma^2*(nobs(mod5s) - 5)
SSE <- smod5$sigma^2*(nobs(mod5) - 5)

Fc <- qf(0.955, 5, 990)
Fst <- ((SSE - SSEns - SSEs) / 5) / ((SSEns + SSEs) / (nobs(mod5) - 2*5))

if(Fst > Fc) {
  print("The data provides enough evidence towards constructing
        two regressions for south and non-south subsamples")
} else {
  print("There is not enough evidence to ensure that running
        two separate regressions will improve the fit significantly")
}


## Indicator variables in log-linear models
mod1 <- lm(log(wage) ~ educ + female, data = cps4_small)

# Computing approximate and exact percentage effects
approx <- 100*coef(mod1)[["female"]]
exact <- 100*(exp(coef(mod1)[["female"]]) - 1)


## The linear probability model

# Linear probability example
data("coke")

mod2 <- lm(coke ~ pratio + disp_coke + disp_pepsi, data = coke)
coef(summary(mod2))

# Graphing the linear probability model
b00 <- coef(mod2)[[1]]
b10 <- b00 + coef(mod2)[["disp_coke"]]
b01 <- b00 + coef(mod2)[["disp_pepsi"]]
b11 <- b01 + b10 - b00
b2 <- coef(mod2)[["pratio"]]

plot(coke$pratio, coke$coke, xlab = "Price ratio", ylab = "Pr[coke]")
abline(b00, b2, lty = 2, col = 2)
abline(b10, b2, lty = 3, col = 3)
abline(b11, b2, lty = 4, col = 4)
abline(b01, b2, lty = 5, col = 5)
legend("topright", c("00","10","11","01"),
       lty=c(2,3,4,5), col=c(2,3,4,5))


## Treatment effects

# Project STAR, an application of the simple difference estimator
data("star", package = "PoEdata")
star <- star[, c("totalscore", "small", "tchexper", "boy",
                         "freelunch", "white_asian", "tchwhite", "tchmasters",
                         "schurban", "schrural", "schid")]
starregular <- star %>%
  filter(small == 0)
starsmall <- star %>%
  filter(small == 1)

summary(starregular)
summary(starsmall)

mod3 <- lm(totalscore ~ small, data = star)
b2 <- coef(mod3)[["small"]]

# Including additional variables into regression to enhance precision
school <- as.factor(star$schid) # creating dummies for schools

mod4 <- lm(totalscore ~ small + tchexper, data = star)
mod5 <- lm(totalscore ~ small + tchexper + school, data = star)

b2_4 <- coef(mod4)[["small"]]
b2_5 <- coef(mod5)[["small"]]

anova(mod4, mod5)

# Regressing "small" on other regressor to identify whether the assignment for treatment and control groups was random
mod6 <- lm(small ~ boy + white_asian + tchexper + freelunch, data = star)
summary(mod6)


## The difference-in-difference estimator
data("njmin3", package = "PoEdata")

mod1 <- lm(fte ~ nj + d + nj*d, data = njmin3)
mod2 <- lm(fte ~ nj*d + kfc + roys + wendys + co_owned, data = njmin3)
mod3 <- lm(fte ~ nj*d + kfc + roys + wendys + co_owned + southj + centralj + pa1, 
           data = njmin3)
one <- coef(summary(mod1))
two <- coef(summary(mod2))
three <- coef(summary(mod3))
summary <- as.data.frame(rbind(one, two, three))

# Vizualizing the treatment effect
b1 <- coef(mod1)[[1]]
b2 <- coef(mod1)[["nj"]]
b3 <- coef(mod1)[["d"]]

delta <- coef(mod1)[["nj:d"]]

t_a <- b1 + b2 + b3 + delta  # treatment, after
c_a <- b1 + b3               # control, after
t_b <- b1 + b2               # treatment, before
c_b <- b1                    # control, before
D <- b1 + b2 + b3

plot(1, type = "n", xlab = "period", ylab = "fte", xaxt = "n",
     xlim = c(-0.01, 1.01), ylim = c(18, 24))

segments(x0 = 0, y0 = c_b, x1 = 1, y1 = c_a, lty = 1, col = "red")
segments(x0 = 0, y0 = t_b, x1 = 1, y1 = t_a, lty = 3, col = "green")
segments(x0 = 0, y0 = t_b, x1 = 1, y1 = D, lty = 4, col = "blue")

legend("topright", legend = c("control", "treated", "ccounterfactual"),
       lty = c(1, 3, 4), col = c(2, 3, 4))

axis(side = 1, at = c(0, 1), labels = NULL)
