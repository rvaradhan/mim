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
            x$T <- object@T
            x$a <- object@a
            
            print(summary(fit))
            cat("\n")
            
            cat("Multiple interaction parameter:\n\n")
            
            interaction.summary <- matrix(c(
                                    object@a,
                                    object@T,
                                    pchisq(object@T,df=1,lower=FALSE)
                                    ),ncol=1)
            
            rownames(interaction.summary) <- c(
                                        "Interaction parameter",
                                        "t-value",
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
    
    if(object@coxph){
        predict(object@fit.coxph,...)
    }
    else{
        predict(object@fit.glm,...)
    }
    
    })
    

setMethod("plot","mim",
    function(x,y,...){
    
    control.data <- x@data
    control.data[,x@formula.object$trt] <- 0
    
    trt.data <- x@data
    trt.data[,x@formula.object$trt] <- 1
   
    beta <- rep(1,length(coef(x)))
    X <- model.matrix(control.data,x@formula.object$formula)
    S <- X%*%beta 
    
    p0 <- predict(y,newdata=control.data)
    p1 <- predict(y,newdata=trt.data)

    plot(y=p0,x=S,...)
    point(y=p1,x=S)

    })
    
    
    