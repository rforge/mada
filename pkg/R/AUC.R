AUC <- function(x, ...) UseMethod("AUC")

AUC.default <- function(x, fpr = 1:99/100, ...){
  sroc <- x
  stopifnot(is.vector(fpr), all(fpr <= 1), all(fpr >= 0), is.function(sroc))
  if(! length(formals(sroc)) == 1){stop("expected a function with exactly one argument")}
  n <- length(fpr)
  if(n < 10){stop("specify at least 10 FPR values!")}
  s <- numeric(n)
  for(i in 1: n){
    temp <- try(sroc(fpr[i]))
    if(class(temp) == "try-error"){stop(paste("calculation of sroc failed for value of FPR ", fpr[i]))}
    if(temp < 0 | temp > 1){stop("expected values of sroc to be >= 0 and <= 1, but this is not the case for FPR value ", fpr[i])}      
    s[i] <- temp
  } # end of loop 
  AUC <-  (s[1]/2 + sum(s[2:(n-1)]) + s[n]/2)/n
  ret <- list(AUC = AUC)
  class(ret) <- "AUC"
  return(ret)
}

AUC.phm <- function(x, level = 0.95, ...){
  theta <- coef(x)[1]
  ci <- confint(x, level = level)["theta",]
  AUC <- 1/(theta+1)
  ci <- 1/(ci+1)
  ret <- list(AUC = AUC, ci = ci)
  class(ret) <- "AUC"
  return(ret)
}

AUC.reitsma <- function(x, fpr = 1:99/100, ...){
  estimate <- x$coefficients
  alpha.sens <- 1
  alpha.fpr <- 1
  if(length(estimate) == 7){alpha.sens <- estimate[6]; alpha.fpr <- estimate[7]}
  mu1 <- estimate[1]
  mu2 <- estimate[2]
  sigma2 <- estimate[4]
  sigma <- estimate[5]  
  rsroc <- function(x){mada:::calc.sroc(x, alpha.sens, alpha.fpr, mu1, mu2, sigma2, sigma)}
  AUC.default(rsroc, fpr = fpr, ...)
}

