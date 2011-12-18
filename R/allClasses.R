setOldClass("coxph")

setClass("mim",
	representation(
		 coef="numeric",
		 vcov="matrix",
		 formula.object="list",
         fit.glm = "glm",
         fit.coxph = "coxph",
         coxph ="logical",
         exact = "logical",
		 data="data.frame",
		 T="numeric",
                 LRT = "numeric",
		 a="numeric"),
	prototype(
	     coef=numeric(),
		 vcov=matrix(),
		 formula.object=list(),
         fit.glm = glm(y~x,data=data.frame(y=runif(10),x=runif(10))),
         fit.coxph = coxph(Surv(t,e)~x,
                    data=data.frame(t=rgamma(10,1,1),
                    e=rep(1,10),x=runif(10))),
         coxph = logical(),
         exact = logical(),
		 data=data.frame(),
		 T=numeric(),
                 LRT = numeric(),
		 a=numeric()
))
