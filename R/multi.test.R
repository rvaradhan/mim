mim <- function(formula,data,family=gaussian,exact=TRUE,select=FALSE,interval=c(-2,2),...){

formula.interaction <- function(f=y~(a+b)*trt,cox=FALSE,family){
		
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
			 selection.formula=selection.formula,
                         cox=cox,
                         family=family) 
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


get.formula <- function(formula,data,family=gaussian,select=FALSE,...){
		
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
			stop("'family' not recognized")
		}

                       
		formula.object <- formula.interaction(formula,cox,family)
		
	if(select){
           if(formula.object$cox){				
		 selection.fit <- coxph(formula.object$selection.formula,data=data,...)		
			}
	    else{				
		 selection.fit <- glm(formula.object$selection.formula,data=data,family=formula.object$family,...)
               }
			
          fit.interaction <- stepAIC(selection.fit,direction="backward")			
          if(is.null(coef(fit.interaction))) stop("No terms of interaction selected.")	

           f <- selection.update(formula,fit.interaction$formula)			
	   formula.object <- formula.interaction(f,cox,family)			
		
		}
 formula.object
}

		
follmann.test <- function(fit,formula.object){
		
		V <- vcov(fit)
		coef <- coef(fit)
                
		vmat1 <- V[formula.object$beta0.index,formula.object$beta0.index]
		vmat2 <- V[formula.object$beta1.index,formula.object$beta1.index]
		vmat <- (vmat1 + vmat2)/2
		vc <- solve(t(chol(vmat)))
		
		U0 <- c(vc %*% coef[formula.object$beta0.index])
		U1 <- c(vc %*% coef[formula.object$beta1.index])
		R <- sqrt(sum(U1*U1)/sum(U0*U0))
		cost <- sum(U0 * U1) / sqrt(sum(U1*U1) * sum(U0 * U0))
		a <- ((R - 1/R) + sqrt((R - 1/R)^2 + 4 * cost^2)) / (2 * cost)
		u0 <- (a * U1 + U0) / (1 + a^2)
		u0.norm <- sqrt(sum(u0*u0))/length(u0)
		u1 <- a * u0
		T <- sum((U1 - U0)^2)/2 - sum((U1 - u1)^2) - sum((U0 - u0)^2)
		
		list(T=T,a=1/a,LRT=numeric(0))
	}

follmann.exact <- function(f.object,data,interval=c(-2,2),...){

for.optimizer <- function(prp.factor,data,f.object){
                    #FUNCTION OF PROPORTIONAL FACTOR USED IN EXACT OPTIMIZATION PROCEDURE
                    #FACTOR IS APPLIED TO ALL TREATMENT COVARIATES BEFORE FITTING
   prpdf <- model.matrix(f.object$f,data)
   if(f.object$cox&all(prpdf[,1]==1)) prpdf <- prpdf[,-1]
   prpdf[,f.object$beta0.index] <- prpdf[,f.object$beta0.index]*prp.factor
   prpdf[,f.object$beta0.index] <- prpdf[,f.object$beta0.index]+prpdf[,f.object$beta1.index]
   prpdf <- prpdf[,-f.object$beta1.index]

   var.names <- all.vars(f.object$f)
   response.index <- 1

   if(f.object$cox) response.index <- 1:2

   response <- var.names[response.index]
   covar.names <- var.names[-response.index]
   covar.names <- covar.names[covar.names!=f.object$trt]
  
   prpdf <- data.frame(prpdf)
   names(prpdf)[f.object$beta0.index] <- covar.names

   exact.formula <- paste("~-1+",paste(names(prpdf),collapse="+"),collapse="")
   exact.formula <- update(f.object$f,exact.formula)
     
   for(var in response){
    prpdf[[var]] <- data[,var]
  }
   
  if(f.object$cox){
    fit <- coxph(exact.formula,data=prpdf)
    -fit$loglik[2]
  }
  else{
   fit <- glm(exact.formula,data=prpdf,family=f.object$family)
   -as.numeric(logLik(fit)) 
 }
}

for.fit <- function(prp.factor,data,f.object){ #RETURN FINAL FIT WITH OPTIMAL PROP FACTOR

   prpdf <- model.matrix(f.object$f,data) 
   if(f.object$cox&all(prpdf[,1]==1)) prpdf <- prpdf[,-1]
   prpdf[,f.object$beta0.index] <- prpdf[,f.object$beta0.index]*prp.factor
   prpdf[,f.object$beta0.index] <- prpdf[,f.object$beta0.index]+prpdf[,f.object$beta1.index]
   prpdf <- prpdf[,-f.object$beta1.index]

   var.names <- all.vars(f.object$f)
   response.index <- 1
   if(f.object$cox) response.index <- 1:2
   response <- var.names[response.index]
   covar.names <- var.names[-response.index]
   covar.names <- covar.names[covar.names!=f.object$trt]
  
   prpdf <- data.frame(prpdf)
   names(prpdf)[f.object$beta0.index] <- covar.names
   
   exact.formula <- paste("~-1+",paste(names(prpdf),collapse="+"),collapse="")
   exact.formula <- update(f.object$f,exact.formula)
     
   for(var in response){
    prpdf[[var]] <- data[,var]
  }
   
  if(f.object$cox){
    fit <- coxph(exact.formula,data=prpdf)
    }
  else{
   fit <- glm(exact.formula,data=prpdf,family=f.object$family)
  }
   
fit
}

   result <- optimize(f=for.optimizer,interval=interval,data=data,f.object=f.object)
   loglik <- result$obj
   fit <- for.fit(result$min,data,f.object=f.object)

## LIKELIHOOD RATIO TEST FOR INTERACTION AGAINST NO INTERACTION MODEL

   response.index <- 1
   if(f.object$cox) response.index <- 1:2

   var.names <- all.vars(f.object$f)
   null.formula <- paste("~",paste(var.names[-response.index],collapse="+"),collapse="")
   null.formula <- update(f.object$f,null.formula)

   null.df <- data[var.names]
 
   if(f.object$cox){
    fit.null <- coxph(null.formula,data=null.df)
    loglik.null <- -fit.null$loglik[2]
   }
  else{
    fit.null <- glm(null.formula,data=null.df,family=f.object$family)
    loglik.null <- -as.numeric(logLik(fit.null))
  }

   list(
        fit=fit,
        test=list(
          a = result$min,
          T = numeric(0),
          LRT = 2*(loglik.null-loglik)
          )
        )
}

follmann.fit <- function(f.object,data,...){
  
    if(f.object$cox){
      
      fit <- coxph(f.object$f,data=data,...)
      test <- follmann.test(fit,f.object)
      
      }
    else{      
      fit <- glm(f.object$f,data=data,family=f.object$family,...)
      test <- follmann.test(fit,f.object)
     }
     
     list(fit=fit,test=test)
}

  formula.object <- get.formula(formula,data,family,select)
  
  if(!exact){
    fit <- follmann.fit(formula.object,data,...)
  }
  else{
    fit <- follmann.exact(formula.object,data,interval,...)
  }
   
  coxph <- class(fit$fit)[1] == "coxph"

 if (coxph) {
        new("mim", coef = coef(fit$fit), vcov = vcov(fit$fit), 
            formula.object = formula.object, fit.coxph = fit$fit,coxph = coxph,
            exact = exact, data = data, T = fit$test$T, LRT = fit$test$LRT, a = fit$test$a)
    }
    else {
        new("mim", coef = coef(fit$fit), vcov = vcov(fit$fit), 
            formula.object = formula.object, fit.glm = fit$fit, coxph = coxph,
            exact = exact, data = data, T = fit$test$T, LRT = fit$test$LRT, a =fit$test$a)
    }
}



