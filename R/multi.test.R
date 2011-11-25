mim <- function(formula,data,family=gaussian,select=FALSE,...){

	formula.interaction <- function(f=y~(a+b)*trt,cox=FALSE){
		
		pattern = "\\(.*\\+.*\\) \\*.*"
		
		if(length(as.character(f))!=3||length(grep(pattern,as.character(f)[3]))!=1)
		stop("Error: Multivariate interaction formula misspecified.")
		
		pattern = "(.*)(\\(.*\\))( \\* )([a-zA-Z]+)(.*)"
		formula = sub(pattern,"\\1\\2:\\4+\\2:I(1-\\4)+\\4\\5",as.character(f)[3])
		formula = formula(paste(as.character(f)[2],as.character(f)[1],formula,collapse=""))
		trt = sub(pattern,"\\4",as.character(f)[3])
		
		selection.formula = sub(pattern,"\\1\\2",as.character(f)[3])
		selection.formula = formula(paste(as.character(f)[2],as.character(f)[1],selection.formula,collapse=""))	
		
		ctrl = sub(pattern,"I(1 - \\4)",as.character(f)[3])
		ctrl.adj = sub(pattern,"~.-1+I(1-\\4)",as.character(f)[3])
		
		if(!cox) formula = update(formula,ctrl.adj)

		labels = attr(terms(formula),"term.labels")
		beta1.index = grep("(:I\\()|(\\):)",labels) #PROBLEM: OTHER INTERACTIONS WITH I()?
		beta0.index = grep(":",labels)
		beta0.index = beta0.index[-match(beta1.index,beta0.index)]
		
		a1.index = match(trt,labels)
		a0.index = match(ctrl,labels)
		
		list(f=formula,
			 beta0.index=beta0.index,
			 beta1.index=beta1.index,
			 a0.index=a0.index,
			 a1.index=a1.index,
			 trt=trt,
			 selection.formula=selection.formula) 
	}
	
	selection.update <- function(formula.old,formula.new){
		
		pattern = "(.*)\\(?(.*)\\)?( \\* )([a-zA-Z]+)(.*)"
		pattern.replace = "(.*)(\\(.*)\\))( \\* )([a-zA-Z]+)(.*)"
		f.replace = as.character(formula.new)[3]
		if(length(grep("\\(",f.replace))==0) f.replace = paste("(",f.replace,")",sep="",collapse="")
		f.replaced = sub("\\(.*\\)",f.replace,as.character(formula.old)[3])
		f.replaced = formula(paste(as.character(formula.old)[2],as.character(formula.old)[1],f.replaced,collapse=""))
		
		f.replaced
	}

	glm.or.coxph <-
	
	function(formula,data,family=gaussian,select=FALSE,...){
		
		cox <- FALSE
		
		if(is.character(family)) {
			if(family=="coxph") {
				cox <- TRUE
			}
			else{
				family <- get(family, mode = "function", envir = parent.frame())
			}
		}
		
		if (is.function(family)) 
        family <- family()
        
		if (is.null(family)|(is.list(family)&&is.null(family$family))) {
			print(family)
			stop("'family' not recognized")
		}
		
		formula.object <- formula.interaction(formula,cox)
		
		if(select){
	
			if(cox){				
				selection.fit <- coxph(formula.object$selection.formula,data=data,...)				
			}
			else{				
				selection.fit <- glm(formula.object$selection.formula,data=data,family=family,...)				
			}
			
			fit.interaction <- stepAIC(selection.fit,direction="backward")
			
			if(is.null(coef(fit.interaction))) stop("No terms of interaction selected.")
			
		    f <- selection.update(formula,fit.interaction$formula)
			
			formula.object <- formula.interaction(f,cox)
			
		
		}
										
		
		if(cox){
			
			fit <- coxph(formula.object$f,data=data,...)
			
		}
		else{
			
			fit <- glm(formula.object$f,data=data,family=family,...)
			
		}
		
		list(fit=fit,formula.object=formula.object)
	}
	
	
	follmann.test <- function(fit){
		
		V <- vcov(fit$fit)
		
		vmat1 <- V[fit$formula.object$beta0.index,fit$formula.object$beta0.index]
		vmat2 <- V[fit$formula.object$beta1.index,fit$formula.object$beta1.index]
		vmat <- (vmat1 + vmat2)/2
		vc <- solve(t(chol(vmat)))
		
		U0 <- c(vc %*% fit$fit$coef[fit$formula.object$beta0.index])
		U1 <- c(vc %*% fit$fit$coef[fit$formula.object$beta1.index])
		R <- sqrt(sum(U1*U1)/sum(U0*U0))
		cost <- sum(U0 * U1) / sqrt(sum(U1*U1) * sum(U0 * U0))
		a <- ((R - 1/R) + sqrt((R - 1/R)^2 + 4 * cost^2)) / (2 * cost)
		u0 <- (a * U1 + U0) / (1 + a^2)
		u0.norm <- sqrt(sum(u0*u0))/length(u0)
		u1 <- a * u0
		T <- sum((U1 - U0)^2)/2 - sum((U1 - u1)^2) - sum((U0 - u0)^2)
		
		list(T=T,a=a)
	}
	
       exact.follmann.test <- function(fit){

         #WRITING FUNCTION TO INPUT SCALING CONSTANT AND OPTIMIZE
         #TEST AGAINST MODEL WITH NO INTERACTION EFFECT
         
cox.fn <- function(a, data) {
f <- coxph(Surv(dt, e) ~ trt + I((1 + (a-1)*trt)*x1) + I((1 + (a-1)*trt)*x2), data=data)
-f$loglik[2]
}

a <- optimize(f=cox.fn, interval=c(-4,4), data=data)$min

f1 <- coxph(Surv(dt, e) ~ trt + I((1-trt)*x1 + a*trt*x1) + I((1-trt)*x2 + a*trt*x2), data=data)

f0 <- coxph(Surv(dt, e) ~ trt + x1 + x2, data=data)
T <- 2*(f1$loglik[2] - f0$loglik[2])

logit.fn <- function(a, data) {
f <- glm(Y ~ Tx + I((1 + (a-1)*Tx)*age) + I((1 + (a-1)*Tx)*sex) + I((1 + (a-1)*Tx)*genhel) + I((1 + (a-1)*Tx)*lvef), data=data, family=binomial)
f$deviance
}

a <- optimize(f=logit.fn, interval=c(-4,4), data=data)$min

f1 <- glm(Y ~ Tx + I((1 + (a-1)*Tx)*age) + I((1 + (a-1)*Tx)*sex) + I((1 + (a-1)*Tx)*genhel) + I((1 + (a-1)*Tx)*lvef), data=data, family=binomial)

f0 <- glm(Y ~ Tx + age + sex + genhel + lvef, data=data, family=binomial)
lrt <- f0$deviance - f1$deviance
lrt

       }

        
  fit <- glm.or.coxph(formula,data,family,select=select,...)
  test <- follmann.test(fit)
  coxph <- class(fit$fit)[1]=="coxph"
  
  
  if(coxph){
      new("mim",
		 coef=coef(fit$fit),
		 vcov=vcov(fit$fit),
		 formula.object=fit$formula.object,
                 fit.coxph = fit$fit,
                 coxph = coxph,
		 data=data,
		 T=test$T,
		 a=1/test$a
		 )
  }
  else{
      new("mim",
		 coef=coef(fit$fit),
		 vcov=vcov(fit$fit),
		 formula.object=fit$formula.object,
                 fit.glm = fit$fit,
                 coxph = coxph,
		 data=data,
		 T=test$T,
		 a=1/test$a
		 )
  }

}



