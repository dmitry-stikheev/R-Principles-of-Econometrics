library(PoEdata)
library(knitr)
library(xtable)
library(printr)
library(effects)
library(car)
library(AER)
library(broom)

# Big Andy’s Hamburger Sales
data("andy")
s <- tidy(andy)[, c(1:5,8,9)]
s
mod1 <- lm(sales ~ price + advert, data = andy)
smod1 <- summary(mod1)
coef(smod1)
effprice <- effect("price", mod1)
plot(effprice)
alleffandy <- allEffects(mod1)
plot(alleffandy)


mod2 <- lm(sales ~ price + advert + I(advert^2), data = andy)
smod2 <- summary(mod2)
smod2
plot(effect("I(advert^2)", mod2))

mod1 <- lm(sales ~ price + advert, data = andy)
smod1 <- summary(mod1)
df <- df.residual(mod1)
N <- nobs(mod1)
b1 <- coef(smod1)[1, 1]
b2 <- coef(smod1)[2, 1]
b3 <- coef(smod1)[3, 1]
sighat2 <- smod1$sigma^2
anov <- anova(mod1)
anov
SSE <- anov[3, 2]
SST <- sum(anov[, 2])
SSR <- SST - SSE
vcov(mod1)
varb1 <- vcov(mod1)[1, 1]
varb2 <- vcov(mod1)[2, 2]
varb3 <- vcov(mod1)[3, 3]
covb1b2 <- vcov(mod1)[1, 2]
covb2b3 <- vcov(mod1)[2, 3]
covb1b3 <- vcov(mod1)[1, 3]
seb2 <- sqrt(varb2)
seb3 <- sqrt(varb3)

alpha <- 0.05
df <- df.residual(mod1)
N <- nobs(mod1)
tc <- qt(1 - alpha / 2, df)
lowb2 <- b2 - tc*seb2
upb2 <- b2 + tc*seb2
lowb3 <- b3 - tc*seb3
upb3 <- b3 + tc*seb3
cib2 <- cbind(lowb2, b2, upb2)
cib3 <- cbind(lowb3, b3, upb3)
cimatrix <- rbind(cib2, cib3)
rownames(cimatrix) <- c("b2", "b3")
colnames(cimatrix) <- c("lowb", "sample value", "upb")
cimatrix


# Testing hypotheses
c <- 0
alpha <- 0.05
df <- df.residual(mod1)
tc <- qt(1 - alpha / 2, df)
b2 <- coef(smod1)[2, 1]
seb2 <- sqrt(vcov(mod1)[2, 2])
t <- (b2 - c) / seb2 
t
curve(dt(x, df), -8, 5, xlab = "t", ylab = "")
abline(v = c(t, -tc, tc), col = c("blue", "red", "red"), lty = c(3, 2, 2))
legend("topleft", legend = c("t", "-tc", "tc"),
       lty = c(3, 2, 2), col = c("blue", "red", "red"))
pval <- 2*(1 - pt(abs(t), df))

c <- 0
alpha <- 0.05
df <- df.residual(mod1)
tc <- qt(1 - alpha / 2, df)
b3 <- coef(smod1)[3, 1]
seb3 <- sqrt(vcov(mod1)[3, 3])
t <- (b3 - c) / seb3
curve(dt(x, df), -5, 5, xlab = "t", ylab = "")
abline(v = c(-tc, tc, t), col = c("red", "red", "blue"), lty = c(2, 2, 3))
legend("topleft", legend = c("-tc", "tc", "t"),
       lty = c(2, 2, 3), col = c("red", "red", "blue"))
pval <- 2*(1 - pt(abs(t), df))
pval

alpha <- 0.05
c <- 0
df <- df.residual(mod1)
tc <- qt(1 - alpha, df)
b2 <- coef(mod1)[[2]]
seb2 <- sqrt(vcov(mod1)[2, 2])
t <- (b2 - c) / seb2
pval <- 1 - pt(abs(t), df)
pval

c <- 1
tc <- qt(1 - alpha, df)
t <- (b3 - c) / seb3
pval <- 1 - pt(abs(t), df)
pval

# H0: 0β1 − 0.2*β2 − 0.5*β3 ≤ 0 ,HA: 0β1 − 0.2*β2 − 0.5*β3 > 0
df <- df.residual(mod1)
A <- as.vector(c(0, -0.2, -0.5))
V <- vcov(mod1)
L <- as.numeric(t(A) %*% coef(mod1))
seL <- as.numeric(sqrt(t(A) %*% V %*% A))
t <- L / seL
pval <- 1 - pt(abs(t), df)
pval


# Polynomial regression models
mod2 <- lm(sales ~ price + advert + I(advert^2), data = andy)
smod2 <- summary(mod2)
smod2$r.squared

advlevels <- c(0.5, 2)
b3 <- coef(smod2)[3, 1]
b4 <- coef(smod2)[4, 1]
DsDa <- b3 + 2*b4*advlevels

alpha <- 0.05
df <- mod2$df.residual
tc <- qt(1 - alpha / 2, df)
g <- (1 - b3) / 2*b4
DgDb3 <- -1 / 2*b4
DgDb4 <- (b3 - 1) / 2*b4^2
varb3 <- vcov(mod2)[3, 3]
varb4 <- vcov(mod2)[4, 4]
covb3b4 <- vcov(mod2)[3, 4]
varg <- DgDb3^2*varb3 + DgDb4^2*varb4 + 2*DgDb3*DgDb4*covb3b4
seg <- sqrt(varg)
lowb <- g - tc*seg
upb <- g + tc*seg


# Interaction terms in linear regression
data("pizza4")
head(pizza4)
mod3 <- lm(pizza ~ age * income, data = pizza4)
smod3 <- summary(mod3)
inc <- c(25, 90)
b2 <- coef(smod3)[2, 1]
b4 <- coef(smod3)[4, 1]
DpDa <- b2 + b4*inc

data("cps4_small")

educ_bar <- mean(cps4_small$educ)
exper_bar <- mean(cps4_small$exper)
mod4 <- lm(log(wage) ~ educ * exper + I(exper^2), data = cps4_small)
smod4 <- summary(mod4)
b3 <- coef(mod4)[[3]]
b4 <- coef(mod4)[[4]]
b5 <- coef(mod4)[[5]]
pDwDex <- 100*(b3 + 2*b4*exper_bar + b5*educ_bar)


# Goodness of fit in multiple regression
mod1 <- lm(sales ~ price + advert, data = andy)
smod1 <- summary(mod1)
Rsq <- smod1$r.squared
anov <- anova(mod1)
tidy(anov)
