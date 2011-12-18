\name{mim}
\alias{mim}

\title{Perform a multivariate interaction test}

\description{

Fits and tests for a multivariate interaction with \code{coxph} or \code{glm} model.

}

\usage{
mim(formula,data,family=gaussian,exact=TRUE,select=FALSE,interval=c(-2,2),...)
}

\arguments{
  \item{formula}{multiple interaction formula for \code{glm} or \code{coxph}, see details}
  \item{data}{data.frame containing the variables of \code{formula}}
  \item{family}{family supplied to \code{family} for \code{glm}, or character "coxph" if \code{coxph} model is fit}
 \item{exact}{logical, If \code{TRUE} an exact estimate of the proportional interaction parameter is obtained}
  \item{select}{logical, If \code{TRUE}, a stepwise backward selection using AIC is applied to a model with interaction variables.}
  \item{interval}{supplied to \code{optimize} when \code{exact} estimation is used. Vector of min and max for the proportional interaction parameter.}
  \item{...}{additional arguments supplied to fit formula, \code{glm} or \code{coxph}}
}

\value{
Returns instance of \code{mim} class.
}

\details{
To test a multiple interaction between treatment variable indicator \code{z} (binary, 0 or 1) and variables \code{a}, \code{b}, with response \code{y} of a GLM model, \code{formula} must be \code{y~(a+b)*z}. If a Cox model with event time \code{time} and event indicator \code{event}, \code{formula} is \code{Surv(time,event)~(a+b)*z}.
}

\author{Ravi Varadhan, \email{rvaradhan@jhmi.edu} and S. Kovalchik \email{s.a.kovalchik@gmail.com}}

\seealso{\link{glm}, \link{coxph}}

\examples{

set.seed(11903)

#FOLLMAN LOGISTIC MODEL
#DATA GENERATION

N0 <- 1000             
pars <- c(-3, -1, 0.1, 0.4, -0.4, -0.7, 0.6)

Tx <- rbinom(N0, size=1, prob = 1/2)
age <- rlnorm(meanlog=log(64), N0, sdlog=0.1)
sex <- rbinom(N0, 1, prob=0.45)	# 1 FOR FEMALE
lvef <- rbeta(N0, shape1=1.0, shape2=1.3)
Xmat <- as.matrix(cbind(age, sex, lvef) )

lpred <- (1 - Tx) * pars[1] + Tx * pars[2] +
    (1 + Tx * (pars[7] - 1)) * (age * pars[3] + sex * pars[4] + lvef * pars[6]) 
    
Y <- rbinom(N0, size=1, prob=plogis(lpred))

data <- data.frame(age, sex, lvef, Tx, Y)

#FOLLMANN MODEL FIT
res1 <- mim(Y ~(age + sex + lvef)*Tx,data=data,family=binomial,exact=FALSE)
res1
coef(res1)
summary(res1)

#EXACT MODEL FIT
res1 <- mim(Y ~(age + sex + lvef)*Tx,data=data,family=binomial,exact=TRUE)
res1
coef(res1)
summary(res1)


#COX MODEL
#DATA GENERATION
pars <- c(-1.2, 0.5, -0.4, 0.8)
lam0 <- 0.5
n <- 1000

trt <- rbinom(n, size=1, prob=0.4)
x1 <- rbinom(n, size=1, prob=0.3)
x2 <- rnorm(n, mean=1, sd=0.5)
noise <- rep(c(1,0),length=n)

linear.pred <- pars[1] * trt + (1 + trt*(pars[4]-1)) * (pars[2]*x1 + pars[3] * x2)

h <- lam0 * exp(linear.pred)
dt <- -log(runif(n))/h
Tcens <- quantile(dt, probs=0.7)

e <- ifelse(dt <= Tcens, 1, 0)
dt <- pmin(dt, Tcens)
data <- data.frame(time=dt, event=e, trt=trt, x1=x1, x2=x2, noise=noise)

#FOLLMANN MODEL FIT
res1 <- mim(Surv(time,event)~(x1+x2+noise)*trt,data=data,family="coxph",exact=FALSE,select=TRUE)
res1
coef(res1)
summary(res1)

#EXACT MODEL FIT
res1 <- mim(Surv(time,event)~(x1+x2)*trt,data=data,family="coxph",exact=TRUE,select=TRUE)
res1
coef(res1)
summary(res1)

}


\references{
Follmann DA, Proschan MA. A multivariate test of interaction for use in clinical trials. \emph{Biometrics} 1999; 55(4):1151-1155
}
