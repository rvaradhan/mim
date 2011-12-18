\name{predict}
\docType{methods}
\alias{predict}
\alias{predict,mim-method}

\title{Get risk predictions for \code{mim} object.}

\description{

Computes linear predictors for the dataset used to fit the model.

}

\section{Methods}{
\describe{

 \item{predict}{\code{signature(object = "mim")}:
Linear predictors for dataset used to fit \code{mim} model. Separate fit with all observations supposed to be in the treatment group and all supposed to be in control group.
  }

 }
}


\author{S. Kovalchik \email{s.a.kovalchik@gmail.com}}

\keyword{methods}

\examples{

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

res1 <- mim(Y ~(age + sex + lvef)*Tx,data=data,family=binomial)

p <- predict(res1)

head(p)
}