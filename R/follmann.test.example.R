library(survival)
library(Hmisc)

set.seed(1234)

dir.path = "~/research/projects/multitest"

source(file=file.path(dir.path,"glm.or.coxph.R"))
source(file=file.path(dir.path,"follmann.test.R"))
source(file=file.path(dir.path,"formula.interaction.R"))
source(file=file.path(dir.path,"multi.test.R"))
source(file=file.path(dir.path,"predict.R"))
source(file=file.path(dir.path,"plot.multi.R"))

#FOLLMAN LOGISTIC MODEL
#DATA GENERATION

N0 <- 1000             
pars <- c(-3, -2, 0.1, 0.4, -0.4, -0.7, 0.6)

Tx <- rbinom(N0, size=1, prob = 1/3)
age <- rlnorm(meanlog=log(64), N0, sdlog=0.1)
sex <- rbinom(N0, 1, prob=0.45)	# 1 FOR FEMALE
genhel <- c(t(rMultinom(rbind(c(0.1, 0.20, 0.35, 0.25, 0.10)), N0)))
lvef <- rbeta(N0, shape1=1.0, shape2=1.3)
Xmat <- as.matrix(cbind(age, sex, genhel, lvef) )

lpred <- (1 - Tx) * pars[1] + Tx * pars[2] +
    (1 + Tx * (pars[7] - 1)) * (age * pars[3] + sex * pars[4] + genhel * pars[5] + lvef * pars[6]) 
    
Y <- rbinom(N0, size=1, prob=plogis(lpred))

data <- data.frame(age, sex, genhel, lvef, Tx, Y)

#MODEL FIT
res1 <- multi.test(Y ~(age + sex + genhel + lvef)*Tx,data=data,family=binomial)
res1$fit
res1$T  #TEST STATISTIC
res1$a  #INTERACTION EFFECT MULTIPLIER beta.trt = a*beta.ctrl

#MEAN MODIFICATION WITH 1-UNIT CHANGE, MEAN SEVERITY UNDER NULL, WITH TRT, WITHOUT TRT
pre <- predict.multi(res1)
colMeans(pre)

#PLOT OF TREATMENT EFFECT VERSUS SEVERITY UNDER INTERACTION AND NULL HYPOTHESIS

plot.multi(res1)

#COX MODEL
#DATA GENERATION
pars <- c(-1.2, 0.5, -0.4, 0.8)
lam0 <- 0.5
n <- 500

trt <- rbinom(n, size=1, prob=0.4)
x1 <- rbinom(n, size=1, prob=0.3)
x2 <- rnorm(n, mean=1, sd=0.5)

linear.pred <- pars[1] * trt + (1 + trt*(pars[4]-1)) * (pars[2]*x1 + pars[3] * x2)

h <- lam0 * exp(linear.pred)
dt <- -log(runif(n))/h
Tcens <- quantile(dt, probs=0.7)

e <- ifelse(dt <= Tcens, 1, 0)
dt <- pmin(dt, Tcens)
data <- data.frame(time=dt, event=e, trt=trt, x1=x1, x2=x2)

#MODEL FIT
res1 <- multi.test(Surv(time,event)~(x1+x2)*trt,data=data,family="coxph")

res1$fit
res1$T
res1$a

plot.multi(res1,FUN=exp)

############################################################
# The exact approach for finding resrtricted parameters
# Simulation for Cox PH model
cox.fn <- function(a, data) {
f <- coxph(Surv(dt, e) ~ trt + I((1 + (a-1)*trt)*x1) + I((1 + (a-1)*trt)*x2), data=data)
-f$loglik[2]
}

a <- optimize(f=cox.fn, interval=c(-4,4), data=data)$min

f1 <- coxph(Surv(dt, e) ~ trt + I((1-trt)*x1 + a*trt*x1) + I((1-trt)*x2 + a*trt*x2), data=data)

f0 <- coxph(Surv(dt, e) ~ trt + x1 + x2, data=data)
T <- 2*(f1$loglik[2] - f0$loglik[2])

# a simulation study comparing approximate and exact approaches
pars <- c(-1.2, 0.5, -0.4, 0.5)
lam0 <- 0.5
n <- 500
nsim <- 100
T <- T2 <- a <- a2 <- rep(NA, nsim)

for (i in 1:nsim) {
print(i)
trt <- rbinom(n, size=1, prob=0.4)
x1 <- rbinom(n, size=1, prob=0.3)
x2 <- rnorm(n, mean=1, sd=0.5)

linear.pred <- pars[1] * trt + (1 + trt*(pars[4]-1)) * (pars[2]*x1 + pars[3] * x2)

h <- lam0 * exp(linear.pred)
dt <- -log(runif(n))/h
Tcens <- quantile(dt, probs=0.7)

e <- ifelse(dt <= Tcens, 1, 0)
dt <- pmin(dt, Tcens)
data <- data.frame(time=dt, event=e, trt=trt, x1=x1, x2=x2)

# approximate
res1 <- multi.test(Surv(time,event)~(x1+x2)*trt,data=data,family="coxph")

res1$fit
T[i] <- res1$T
a[i] <- 1/res1$a

# exact
a2[i] <- optimize(f=cox.fn, interval=c(-4,4), data=data)$min
f1 <- coxph(Surv(dt, e) ~ trt + I((1 + (a2[i]-1)*trt)*x1) + I((1 + (a2[i]-1)*trt)*x2), data=data)
f0 <- coxph(Surv(dt, e) ~ trt + x1 + x2, data=data)
T2[i] <- 2*(f1$loglik[2] - f0$loglik[2])
}

mean(a)  # biased
mean(a2) # good
mean(T)
mean(T2) # larger than T
##########################################################
# Simulation for logistic regression
logit.fn <- function(a, data) {
f <- glm(Y ~ Tx + I((1 + (a-1)*Tx)*age) + I((1 + (a-1)*Tx)*sex) + I((1 + (a-1)*Tx)*genhel) + I((1 + (a-1)*Tx)*lvef), data=data, family=binomial)
f$deviance
}

a <- optimize(f=logit.fn, interval=c(-4,4), data=data)$min

f1 <- glm(Y ~ Tx + I((1 + (a-1)*Tx)*age) + I((1 + (a-1)*Tx)*sex) + I((1 + (a-1)*Tx)*genhel) + I((1 + (a-1)*Tx)*lvef), data=data, family=binomial)

f0 <- glm(Y ~ Tx + age + sex + genhel + lvef, data=data, family=binomial)
lrt <- f0$deviance - f1$deviance
lrt

N0 <- 500             
pars <- c(-3, -2, 0.1, 0.4, -0.4, -0.7, 0.6)
nsim <- 100
T <- T2 <- a <- a2 <- rep(NA, nsim)

for (i in 1:nsim) {
print(i)
Tx <- rbinom(N0, size=1, prob = 1/3)
age <- rlnorm(meanlog=log(64), N0, sdlog=0.1)
sex <- rbinom(N0, 1, prob=0.45)	# 1 FOR FEMALE
genhel <- c(t(rMultinom(rbind(c(0.1, 0.20, 0.35, 0.25, 0.10)), N0)))
lvef <- rbeta(N0, shape1=1.0, shape2=1.3)
Xmat <- as.matrix(cbind(age, sex, genhel, lvef) )

lpred <- (1 - Tx) * pars[1] + Tx * pars[2] +
    (1 + Tx * (pars[7] - 1)) * (age * pars[3] + sex * pars[4] + genhel * pars[5] + lvef * pars[6]) 
    
Y <- rbinom(N0, size=1, prob=plogis(lpred))

data <- data.frame(age, sex, genhel, lvef, Tx, Y)

#MODEL FIT
res1 <- multi.test(Y ~(age + sex + genhel + lvef)*Tx,data=data,family=binomial)
res1$fit
T[i] <- res1$T  #TEST STATISTIC
a[i] <- 1/res1$a  #INTERACTION EFFECT MULTIPLIER beta.trt = a*beta.ctrl

a2[i] <- optimize(f=logit.fn, interval=c(-4,4), data=data)$min

f1 <- glm(Y ~ Tx + I((1 + (a2[i]-1)*Tx)*age) + I((1 + (a2[i]-1)*Tx)*sex) + I((1 + (a2[i]-1)*Tx)*genhel) + I((1 + (a2[i]-1)*Tx)*lvef), data=data, family=binomial)
f0 <- glm(Y ~ Tx + age + sex + genhel + lvef, data=data, family=binomial)
T2[i] <- f0$deviance - f1$deviance
}
mean(a)
mean(a2)
mean(T)
mean(T2)



