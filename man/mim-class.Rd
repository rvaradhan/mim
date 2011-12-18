\name{mim-class}
\docType{class}
\alias{mim-class}

\title{Class "mim"}

\description{
Class for multiple interaction model (MIM).
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("mim", ...)}. 
}

\section{Slots}{
  \describe{
    \item{\code{coef}:}{Vector of coefficients for restricted interaction model.}
    \item{\code{vcov}:}{Matrix of variance-covaraince for restricted interaction model.}
    \item{\code{formula.object}:}{List of formula variables. See below.}
    \item{\code{fit.glm}:}{Restricted \code{glm} object for \code{glm} \code{mim} model.}
    \item{\code{fit.coxph}:}{Restricted \code{coxph} object for \code{coxph} \code{mim} model.}
    \item{\code{coxph}:}{Logical, whether \code{coxph} model.}
    \item{\code{exact}:}{Logical, whether exact method used for fitting.}
    \item{\code{data}:}{Data frame containing the outcomes and explanatory variables.}
    \item{\code{T}:}{Wald based test of a!=1 against no interaction model used when \code{exact} is \code{FALSE}}
    \item{\code{LRT}:}{Likelihood ratio test of a!=1 against no interaction model used when \code{exact} is \code{TRUE}}
    \item{\code{a}:}{Proportional interaction parameter}
}

The \code{formula.object} is a list with components:

\describe{
    \item{f}{formula of restricted fit}
    \item{beta0.index}{index of the control arm covariate parameters}
    \item{beta1.index}{index of the treatment arm covariate parameters}
    \item{a0.index}{index of the control baseline parameter}
    \item{a1.index}{index of the treatment baseline parameter}
    \item{trt}{name of the treatment variable}
    \item{selection.formula}{formula without interaction effect for running selection procedure}
    \item{cox}{logical indicator if model is a \code{coxph} model}
    \item{family}{character description of model family}
}

}

\section{Methods}{

  \describe{

\item{show}{\code{signature(object = "mim")}:
Display point estimates of \code{mim} object.
  }

\item{print}{\code{signature(x = "mim",...)}:
Display point estimates of \code{mim} object.
  }

\item{summary}{\code{signature(object = "mim",...)}:
List of estimates and convergence information.
  }

\item{coef}{\code{signature(object = "mim")}:
Extractor for MLEs. 
  }
  
\item{predict}{\code{signature(object = "mim")}:
Extractor for MLEs. 
  }

  \item{vcov}{\code{signature(object = "mim")}:
Extractor for variance-covariance of restricted \code{mim} object when \code{exact} method is used.
  }

 }
}

\keyword{classes}
