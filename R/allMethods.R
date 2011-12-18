setMethod("coef","mim",
function(object,...){
    object@coef
})


setMethod("vcov","mim",
function(object,...){
    object@vcov
})

setMethod("print","mim",
         function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\nCall:\n", paste(deparse(x@formula.object$f), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    
    if (length(coef(x))) {
        cat("Coefficients:\n")
        print.default(format(coef(x), digits = digits), print.gap = 2, 
            quote = FALSE)
    }
    else cat("No coefficients\n")
    cat("\n")

    cat("Multiple interaction parameter:\n")
    print.default(format(x@a, digits = digits), print.gap = 2, 
            quote = FALSE)
    cat("\n")
    invisible(x)
})


setMethod("show","mim",
          function(object){
            print(object)
          })


setMethod("summary","mim",
          function(object,digits = max(3, getOption("digits") - 3),...){
          
            if(object@coxph){
                fit <- object@fit.coxph
                }
            else{
                fit <- object@fit.glm
                }
                
            x <- summary(fit)
            x$test <- ifelse(length(object@T)==0,object@LRT,object@T)
            test.name <-  ifelse(length(object@T)==0,"LRT","t-test")
            x$a <- object@a
            
            print(summary(fit))
            cat("\n")
            
            cat("Multiple interaction parameter:\n\n")
            
            interaction.summary <- matrix(c(
                                    object@a,
                                    x$test,
                                    pchisq(x$test,df=1,lower=FALSE)
                                    ),ncol=1)
            
            rownames(interaction.summary) <- c(
                                        "Interaction parameter",
                                        test.name,
                                        "p-value"
                                        )
                                        
            colnames(interaction.summary) <- "a"
                                                                            
            print.default(format(interaction.summary, digits = digits),      
                                                    print.gap = 2, 
                                                    quote = FALSE)
                                                    
        invisible(x)
          })


setMethod("predict","mim",
    function(object,...){

    data.trt <- object@data
    data.ctrl <- object@data

    data.trt[[object@formula.object$trt]] <- 1
    data.ctrl[[object@formula.object$trt]] <- 0

    a <- ifelse(object@exact,object@a,1)

    if(object@coxph){
      f <- object@fit.coxph$formula
      }
    else{
      f <-object@fit.glm$formula
    }

    ctrl.var <- match(all.vars(f),names(data.trt)) #CONTROL VAR CORRECTION FOR GLM/EXACT
    ctrl <- any(is.na(ctrl.var))

    if(ctrl){
      data.trt[[all.vars(f)[which(is.na(ctrl.var))]]] <- 0
      data.ctrl[[all.vars(f)[which(is.na(ctrl.var))]]] <- 1
    }
    
    X.trt <- model.matrix(f,data.trt)
    X.ctrl <- model.matrix(f,data.ctrl)

    p.trt <- X.trt%*%(object@coef*a)
    p.ctrl <- X.ctrl%*%object@coef

    data.frame(pred.trt=p.trt,pred.ctrl=p.ctrl)
 })
    

setMethod("plot","mim",
    function(x,y,...){

    p <- predict(x)
    
    p0 <- p$pred.ctrl
    p1 <- p$pred.trt

    plot(y=p0,x=p0,xlab="prediction without modification",
         ylab="prediction with modification",...)
    points(y=p1,x=p0,col="blue")
    legend("topleft",pch=1,col=c("black","blue"),legend=c("Control","Treated"))

    })
    
    
    
