# function, output of class "madad"
madad <- function(x=NULL, TP, FN, FP, TN, level = 0.95, correction = 0.5, 
                         correction.control = "all", method = "wilson", yates = TRUE, ...){
  names <- x$names
  alpha<-1-level
	kappa<-qnorm((1-alpha/2))
  DNAME <- deparse(substitute(x))
  if(!is.null(x)){
    X <- as.data.frame(x)
    TP <- X$TP
    FN <- X$FN
    FP <- X$FP
    TN <- X$TN
    origdata <- data.frame(TP = TP, FN = FN, FP = FP, TN = TN)
  }
  if(is.null(x)){origdata <- data.frame(TP = TP, FN = FN, FP = FP, TN = TN)}

  checkdata(origdata)

  ## apply continuity correction to _all_ studies if one contains zero
  if(correction.control == "all"){if(any(c(TP,FN,FP,TN) == 0)){TP <- TP+correction;
    					FN <- FN + correction;
							FP <- FP + correction;
							TN <- TN + correction}}
  if(correction.control == "single"){
	  correction = ((((TP == 0)|(FN == 0))|(FP == 0))| (TN == 0))*correction
    TP <- correction + TP
	  FN <- correction + FN
	  FP <- correction + FP
	  TN <- correction + TN
	}
  
  
  number.of.pos<-TP+FN
  number.of.neg<-FP+TN
  sens.ci<-binomCIvector(TP,number.of.pos, conf.level  = level, method = method)
  colnames(sens.ci) <- c(paste(100*alpha/2, "%", sep="", collapse = ""), paste(100*(1-alpha/2), "%", sep="", collapse = ""))
  spec.ci<-binomCIvector(TN,number.of.neg, conf.level  = level, method = method)
  colnames(spec.ci) <- c(paste(100*alpha/2, "%", sep="", collapse = ""), paste(100*(1-alpha/2), "%", sep="", collapse = ""))
  fpr.ci<-binomCIvector(FP,number.of.neg, conf.level  = level, method = method)
  colnames(fpr.ci) <- c(paste(100*alpha/2, "%", sep="", collapse = ""), paste(100*(1-alpha/2), "%", sep="", collapse = ""))


  sens <- TP/number.of.pos
  spec <- TN/number.of.neg
  fpr <- FP/number.of.neg

  cor_da <- cor(as.matrix(data.frame(sens = sens, spec = spec, fpr = fpr)))
  
	posLR<-sens/fpr
	negLR<-(1-sens)/(1-fpr)
	
	se.lnposLR <- sqrt(1/TP + 1/FP - 1/number.of.pos - 1/number.of.neg)
	se.lnnegLR <- sqrt(1/TN + 1/FN - 1/number.of.pos - 1/number.of.neg)
	posLR.ci<-cbind(exp(-kappa*se.lnposLR),exp(kappa*se.lnposLR))
	posLR.ci<-posLR*posLR.ci
  colnames(posLR.ci) <- c(paste(100*alpha/2, "%", sep="", collapse = ""), paste(100*(1-alpha/2), "%", sep="", collapse = ""))

	negLR.ci<-cbind(exp(-kappa*se.lnnegLR),exp(kappa*se.lnnegLR))
	negLR.ci<-negLR*negLR.ci
	colnames(negLR.ci) <- c(paste(100*alpha/2, "%", sep="", collapse = ""), paste(100*(1-alpha/2), "%", sep="", collapse = ""))

  DOR<-posLR/negLR
  se.lnDOR<-sqrt(1/TP + 1/TN + 1/FN + 1/FP)
	DOR.ci<-cbind(exp(-kappa*se.lnDOR),exp(kappa*se.lnDOR))
	DOR.ci<-DOR*DOR.ci
  colnames(DOR.ci) <- c(paste(100*alpha/2, "%", sep="", collapse = ""), paste(100*(1-alpha/2), "%", sep="", collapse = ""))
  
  output <- list(sens = list(sens = sens, sens.ci = sens.ci),
                 spec = list(spec = spec, spec.ci = spec.ci),
                 fpr = list(fpr = fpr, fpr.ci = fpr.ci),
                 sens.htest = prop.test(TP, number.of.pos, correct = yates),
                 spec.htest = prop.test(TN, number.of.neg, correct = yates),
                 posLR = list(posLR = posLR, posLR.ci = posLR.ci, se.lnposLR = se.lnposLR),
                 negLR = list(negLR = negLR, negLR.ci = negLR.ci, se.lnnegLR = se.lnnegLR),
                 DOR = list(DOR = DOR, DOR.ci = DOR.ci, se.lnDOR = se.lnDOR),
                 cor_sens_fpr = cor(sens,fpr),
                 level = level, method = method, names = names,
                 nobs = nrow(origdata), data = origdata, data.name = DNAME,
                 correction = correction, correction.control = correction.control)
  class(output) <- "madad"
  output
}

print.madad <- function(x, digits = 3, ...){
  cat("Descriptive summary of", x$data.name, "with", x$nobs, "primary studies.\n")
  cat("Confidence level for all calculations set to", 100*x$level, "%\n")
  cat("Using a continuity correction of", x$correction, "if applicable \n")
  cat("\n")
  
  
  cat("Diagnostic accuracies \n")
  output1a <- round(cbind(x$sens$sens, x$sens$sens.ci, x$spec$spec, x$spec$spec.ci), digits)
  rownames(output1a) <- x$names
  colnames(output1a)[c(1,4)] <- c("sens", "spec")
  print(output1a)
  
  format.test.result<-function(x)
  {
	out <- character()
    if (!is.null(x$statistic)) 
        out <- c(out, paste(names(x$statistic), "=", format(round(x$statistic, 
            4))), ", ")
    if (!is.null(x$parameter)) 
        out <- c(out, paste(names(x$parameter), "=", format(round(x$parameter, 
            3))), ", ")
    if (!is.null(x$p.value)) {
        fp <- format.pval(x$p.value, digits = digits)
        out <- c(out, paste("p-value =", fp))}
    return(out)
	} # End of function format.test.result

  sens.result<-format.test.result(x$sens.htest)
	spec.result<-format.test.result(x$spec.htest)
  
	output1b<-paste(c("\nTest for equality of sensitivities: \n",sens.result, "\n",
                   "Test for equality of specificities: \n",spec.result, "\n\n"), collapse="")
	cat(output1b)
  
  cat("\n")
  cat("Diagnostic OR and likelihood ratios \n")
  output2 <- round(cbind(x$DOR$DOR, x$DOR$DOR.ci, x$posLR$posLR, x$posLR$posLR.ci,
                         x$negLR$negLR, x$negLR$negLR.ci), digits)
  rownames(output2) <- x$names
  colnames(output2)[c(1,4,7)] <- c("DOR", "posLR", "negLR")
  print(output2)
  
  cat("\n")
  cat("Correlation of sensitivities and false positive rates: \n")
  print(round(CIrho(x$cor_sens_fpr, x$nobs), digits)[1,])
  return(invisible(NULL))
  }

CIrho <- function(rho, N, level = 0.95){
  stopifnot(rho < 1, rho > -1, N > 3, round(N) == N)
  z <- atanh(rho)
  kappa <- qnorm(1-(1-level)/2)
  output <- matrix(cbind(rho, tanh(z - kappa*sqrt(1/(N-3))),tanh(z + kappa*sqrt(1/(N-3)))), ncol = 3)
  colnames(output) <- c("rho", paste(100*(1-level)/2, "%", collapse =""), 
                     paste(100*(1- (1-level)/2), "%", collapse =""))
  output
}

#calculate cochranes Q, expects (s.e. of) log transformed DOR, posLR, negLR 
cochran.Q<-function(x, weights)
{
	if(length(weights) != length(x)){stop("Length of weights does not match length of x")}				
	bar.x<-sum(weights*x)/sum(weights)
	Q<-sum(weights*((x-bar.x)^2))
	k<-length(x)
	if(k <= 1)stop("x needs to have length at least 2")
	p.value<-pchisq(Q, k-1, lower.tail = FALSE)
	#fp<-format(round(p.value, 4))
	#fQ<-format(round(Q,4))
	output<-c(Q,p.value,as.integer(k-1))
	names(output)<-c("Q", "p-value", "df")
	return(output)
	}



madauni <- function(x, type = "DOR", method = "DSL", suppress = TRUE){

if(suppress){x <- suppressWarnings(madad(x))
             }else{
               x <- madad(x)
             }  
  
TP <- x$data$TP
FP <- x$data$FP
FN <- x$data$FN
TN <- x$data$TN

number.of.pos<-TP+FN
number.of.neg<-FP+TN
nop <- number.of.pos
non <- number.of.neg
total <- nop + non

# from Cochran.Q and inverse variance weights calculate between study variance
naive.tausquared<-function(Q,weights)
{
k<-length(weights)
if(Q<(k-1)){return(0)}
else
return((Q-k+1)/(sum(weights)-(sum(weights^2)/sum(weights))))
}

if(! method %in% c("MH","DSL"))stop("method must be either \"MH\" or \"DSL\"")else
nobs <- x$nobs
theta <-switch(type, "DOR" = x$DOR$DOR, "posLR" = x$posLR$posLR, 
                 "negLR" = x$negLR$negLR)
  
if(method == "MH")
  {
  weights<-switch(type, "DOR" = FP*FN/total, "posLR" = FP*nop/total, 
                  "negLR" = TN*nop/total)
  coef <- log(sum(weights*theta)/sum(weights))
  CQ<-cochran.Q(theta, weights = weights)
  tau.squared <- NULL
  
  P <- sum((nop*non*(TP+FP) - TP*FP*total)/total^2)
  U <- sum(TP*non/total)
  V = sum(FN*nop/total)
  Uprime = sum(FP*non/total)
  Vprime = sum(TN*nop/total)
  R = sum(TP*TN/total)
  S = sum(FP*FN/total)
  E = sum((TP+TN)*TP*TN/(total^2))
  FF = sum((TP+TN)*FN*FP/(total^2))
  G = sum((FP+FN)*TP*TN/(total^2))
  H = sum((FP+FN)*FP*FN/(total^2))
  
  vcov <- switch(type, "DOR" = 0.5*(E/(R^2) + (F+G)/(R*S) + H/(S^2)),
                        "posLR" = P/(U*V), "negLR" = P/(Uprime*Vprime)) 
  }#end of method = "MH"

if(method == "DSL")
  {
  se.lntheta <- switch(type, "DOR" = x$DOR$se.lnDOR, "posLR" = x$posLR$se.lnposLR, 
                       "negLR" = x$negLR$se.lnnegLR)
  lntheta <- log(theta)
  CQ<-cochran.Q(lntheta, 1/se.lntheta^2)
  tau.squared<-naive.tausquared(CQ[1],1/(se.lntheta^2))
  weights<-1/(se.lntheta^2+tau.squared)
  #recalculate CQ based on new weights
  CQ<-cochran.Q(lntheta, weights)  
  coef <- sum(weights*lntheta)/sum(weights)
  vcov <- 1/sum(weights)
  }#end of method == "DSL"

names(coef) <- paste("ln", type, collapse ="", sep ="")

vcov <- matrix(vcov, nrow = 1, ncol = 1)
colnames(vcov) <- paste("ln", type, collapse ="", sep ="")
rownames(vcov) <- paste("ln", type, collapse ="", sep ="")


output <- list(coefficients = coef, vcov = vcov,  tau_sq = tau.squared, weights = weights,
               type = type, method = method, data = x$data, theta = theta, CQ = CQ, nobs = length(theta),
               call = match.call())
class(output) <- "madauni"
output
}# end of function madauni
  
print.madauni <- function(x, digits = 3, ...){
  cat("Call:\n")
  print(x$call)
  cat("\n")
  if(is.null(x$tau_sq)){
  ans <- exp(x$coefficients)
  names(ans) <- x$type
  print(ans)
      }else{
  ans <- c(exp(x$coefficients),x$tau_sq)
  names(ans) <- c(x$type, "tau^2")
  print(round(ans, digits))
  }  
}

vcov.madauni <- function(object){object$vcov}

summary.madauni <- function(object, level = .95, ...){
x <- object
# Calculate Higgins I^2, only if method is not "MH"  
Higgins.Isq<-function(T,df){return(max(0,100*(T-df)/T))}

if(object$method == "DSL"){
  Isq <- Higgins.Isq(x$CQ[1], x$CQ[3])}else{
  Isq <- NULL
}
  
CIcoef <- rbind(exp(cbind(coef(x), confint(x, level = level))),
                cbind(coef(x), confint(x, level = level)))
rownames(CIcoef) <- c(x$type,paste("ln",x$type, sep ="", collapse = ""))
colnames(CIcoef)[1] <- paste(x$method, "estimate", collapse ="")

Q <- function(tau2){sum(((log(x$theta) - coef(x)[1])^2)/(1/x$weights+tau2))}
CQ <- ifelse(is.null(x$tau_sq), Q(0), Q(x$tau_sq))

## Q-Profile Confidence interval for tau_sq like in Viechtbauer (2007)
if(!is.null(x$tau_sq)){
# browser()
  kappa_up <- qchisq(1-(1-level)/2, x$nobs - 1)
  kappa_low <- qchisq((1-level)/2, x$nobs - 1)
  if(Q(0) < kappa_up){lower <- 0}else{
      lower <-  uniroot(function(x){Q(x)-kappa_up}, lower = 0, upper = 10^4)$root}
  if(Q(0) < kappa_low){upper <- 0}else{
  upper <-  uniroot(function(x){Q(x)-kappa_low}, lower = 0, upper = 10^4)$root}
  CIcoef <- rbind(CIcoef, c(x$tau_sq, lower, upper), sqrt(c(x$tau_sq, lower, upper)))
  rownames(CIcoef)[3:4] <- c("tau^2","tau")
}

output <- list(x=object, Isq = Isq, CIcoef = CIcoef)
class(output) <- "summary.madauni"
output
  
}

print.summary.madauni <- function(x, digits = 3,...){

cat("Call:\n")
print(x$x$call)
cat("\nEstimates:\n")
print(round(x$CIcoef,digits))
cat("\nCochran's Q: ",round(x$x$CQ[1],digits),
    " (",round(x$x$CQ[3])," df, p = ", round(x$x$CQ[2], digits),")", sep = "")
if(!is.null(x$Isq)){cat("\nHiggins' I^2: ",round(x$Isq, digits),"%", sep ="")}
}




checkdata <- function(X, nrowwarn = 5){
  X <- as.data.frame(X)
  if(!all(c("TP","FN","FP","TN") %in% names(X))){
    stop("Data frame or matrix must have columns labelled TP, FN, FP and TN.")}
  if(!identical(round(X),X)){stop("Data must consist of counts.")}
  if(nrow(X) < nrowwarn){warning("There are very few primary studies!")}
  return(invisible(NULL))
}

sens <- function(x){x$TP/(x$TP + x$FN)}
fpr <- function(x){x$FP/(x$FP + x$TN)}
spec <- function(x){x$TN/(x$FP + x$TN)}



crosshair <- function(M, correction = 0.5, conf.level = 0.95, method = "wilson", pch = 1, add = FALSE, ...){
  x <- calc.SensSpec(M, FPR = TRUE, correction = correction, conf.level = conf.level, method = method)
  if(!add){plot(x[,4], x[,1], xlim = c(0,1), ylim =c(0,1), pch = pch, ...)}
  if(add){points(x[,4], x[,1], pch = pch, ...)}
  arrows(x[,4], x[,2], x[,4], x[,3], angle = 90, code = 3, length = 0.1, ...)
  arrows(x[,5], x[,1], x[,6], x[,1], angle = 90, code = 3, length = 0.1, ...)
  return(invisible(NULL))
  }

### ellipse from package ellipse

ROC.ellipse <- function(M, correction = 0.5, conf.level = 0.95, method = "wilson", pch = 1, add = FALSE, ...)
{
  x <- calc.SensSpec(M, FPR = TRUE, correction = correction, conf.level = 0.5, method = method)
  if(!add){plot(x[,4], x[,1], xlim = c(0,1), ylim =c(0,1), pch = pch, ...)}
  if(add){points(x[,4], x[,1], pch = pch, ...)}
  logit.x <- logit(x)
  half.conf.level <- 1-(1-conf.level)/2
  for(i in 1:dim(x)[1]){
    lines(inv.logit(ellipse(0,
          centre = c(logit.x[i,4],logit.x[i,1]),
          scale = c((logit.x[i,1]-logit.x[i,2])*qnorm(half.conf.level), (logit.x[i,4]-logit.x[i,5])*qnorm(half.conf.level)))),
          col ="grey")
  points(x[,4], x[,1], xlim = c(0,1), ylim =c(0,1), pch = pch, ... )
    }
  return(invisible(NULL))
}



## FROM package MKmisc, author: Matthias Kohl
## Confidence Intervals for Binomial Proportions
binomCI <- function(x, n, conf.level = 0.95, method = "wilson", rand = 123){
    if (!is.na(pmatch(method, "wilson")))
        method <- "wilson"

    METHODS <- c("wald", "wilson", "agresti-coull", "jeffreys", "modified wilson", 
                 "modified jeffreys", "clopper-pearson", "arcsine", "logit", "witting")
    method <- pmatch(method, METHODS)

    if (is.na(method))
        stop("invalid method")

    if (method == -1)
        stop("ambiguous method")

    if(length(x) != 1)
        stop("'x' has to be of length 1 (number of successes)")
    if(length(n) != 1)
        stop("'n' has to be of length 1 (number of trials)")
    if(length(conf.level) != 1)
        stop("'conf.level' has to be of length 1 (confidence level)")
    if(conf.level < 0.5 | conf.level > 1)
        stop("'conf.level' has to be in [0.5, 1]")

    alpha <- 1 - conf.level
    kappa <- qnorm(1-alpha/2)
    p.hat <- x/n
    q.hat <- 1 - p.hat

    if(method == 1){ # wald
        est <- p.hat
        term2 <- kappa*sqrt(p.hat*q.hat)/sqrt(n)
        CI.lower <- max(0, p.hat - term2)
        CI.upper <- min(1, p.hat + term2)
    }
    if(method == 2){ # wilson
        est <- p.hat
        term1 <- (x + kappa^2/2)/(n + kappa^2)
        term2 <- kappa*sqrt(n)/(n + kappa^2)*sqrt(p.hat*q.hat + kappa^2/(4*n))
        CI.lower <-  max(0, term1 - term2)
        CI.upper <- min(1, term1 + term2)
    }
    if(method == 3){ # agresti-coull
        x.tilde <- x + kappa^2/2
        n.tilde <- n + kappa^2
        p.tilde <- x.tilde/n.tilde
        q.tilde <- 1 - p.tilde
        est <- p.tilde
        term2 <- kappa*sqrt(p.tilde*q.tilde)/sqrt(n.tilde)
        CI.lower <- max(0, p.tilde - term2)
        CI.upper <- min(1, p.tilde + term2)
    }
    if(method == 4){ # jeffreys
        est <- p.hat
        if(x == 0)
            CI.lower <- 0
        else
            CI.lower <- qbeta(alpha/2, x+0.5, n-x+0.5)
        if(x == n)
            CI.upper <- 1
        else
            CI.upper <- qbeta(1-alpha/2, x+0.5, n-x+0.5)
    }
    if(method == 5){ # modified wilson
        est <- p.hat
        term1 <- (x + kappa^2/2)/(n + kappa^2)
        term2 <- kappa*sqrt(n)/(n + kappa^2)*sqrt(p.hat*q.hat + kappa^2/(4*n))
        if((n <= 50 & x %in% c(1, 2)) | (n >= 51 & n <= 100 & x %in% c(1:3)))
            CI.lower <- 0.5*qchisq(alpha, 2*x)/n
        else
            CI.lower <-  max(0, term1 - term2)

        if((n <= 50 & x %in% c(n-1, n-2)) | (n >= 51 & n <= 100 & x %in% c(n-(1:3))))
            CI.upper <- 1 - 0.5*qchisq(alpha, 2*(n-x))/n
        else
            CI.upper <- min(1, term1 + term2)
    }
    if(method == 6){ # modified jeffreys
        est <- p.hat
        if(x == 0)
            CI.lower <- 1 - (alpha/2)^(1/n)
        else{
            if(x == 1)
                CI.lower <- 0
            else
                CI.lower <- qbeta(alpha/2, x+0.5, n-x+0.5)
        }
        if(x == n)
            CI.upper <- (alpha/2)^(1/n)
        else{
            if(x == n-1)
                CI.upper <- 1
            else
                CI.upper <- qbeta(1-alpha/2, x+0.5, n-x+0.5)
        }
    }
    if(method == 7){ # clopper-pearson
        est <- p.hat
        CI.lower <- qbeta(alpha/2, x, n-x+1)
        CI.upper <- qbeta(1-alpha/2, x+1, n-x)
    }
    if(method == 8){ # arcsine
        p.tilde <- (x + 0.375)/(n + 0.75)
        est <- p.tilde
        CI.lower <- sin(asin(sqrt(p.tilde)) - 0.5*kappa/sqrt(n))^2
        CI.upper <- sin(asin(sqrt(p.tilde)) + 0.5*kappa/sqrt(n))^2
    }
    if(method == 9){ # logit
        est <- p.hat
        lambda.hat <- log(x/(n-x))
        V.hat <- n/(x*(n-x))
        lambda.lower <- lambda.hat - kappa*sqrt(V.hat)
        lambda.upper <- lambda.hat + kappa*sqrt(V.hat)
        CI.lower <- exp(lambda.lower)/(1 + exp(lambda.lower))
        CI.upper <- exp(lambda.upper)/(1 + exp(lambda.upper))
    }
    if(method == 10){ # witting
        set.seed(rand)
        x.tilde <- x + runif(1, min = 0, max = 1)
        pbinom.abscont <- function(q, size, prob){
            v <- trunc(q)
            term1 <- pbinom(v-1, size = size, prob = prob) 
            term2 <- (q - v)*dbinom(v, size = size, prob = prob)
            return(term1 + term2)
        }
        qbinom.abscont <- function(p, size, x){
            fun <- function(prob, size, x, p){
                pbinom.abscont(x, size, prob) - p
            }
            uniroot(fun, interval = c(0, 1), size = size, x = x, p = p)$root
        }
        est <- p.hat
        CI.lower <- qbinom.abscont(1-alpha, size = n, x = x.tilde)
        CI.upper <- qbinom.abscont(alpha, size = n, x = x.tilde)
    }

    CI <- c(CI.lower, CI.upper)
    attr(CI, "confidence level") <- conf.level
    return(list("estimate" = est, "CI" = CI))
}

# expects vector x of successes and vector n of trials, output
binomCIvector<-function(x, n, conf.level = 0.95, method="wilson"){
	N<-length(x)
	CI<-matrix(NA,nrow=N,ncol=2)
	colnames(CI)<-c("LL","UL")
	for(i in 1:N){
		CI[i,]<-binomCI(x[i],n[i],conf.level=conf.level,method=method)$CI
		}
		return(CI)
	}


