# function of the class "descript"
madaDescript <- function(x=NULL, TP, FN, FP, TN, level = 0.95, correction = 0.5, 
                         correction.control = "all", method = "wilson", ...){
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
                 posLR = list(posLR = posLR, posLR.ci = posLR.ci),
                 negLR = list(negLR = negLR, negLR.ci = negLR.ci),
                 DOR = list(DOR = DOR, DOR.ci = DOR.ci),
                 cor_sens_fpr = cor(sens,fpr),
                 level = level, method = method, names = names,
                 nobs = nrow(origdata), data = origdata, data.name = DNAME,
                 correction = correction, correction.control = correction.control)
  class(output) <- "madaDescript"
  output
}

print.madaDescript <- function(x, digits = 3){
  cat("Descriptive summary of", x$data.name, "with", x$nobs, "primary studies.\n")
  cat("Confidence level for all calculations set to", 100*x$level, "%\n")
  cat("Using a continuity correction of", x$correction, "if applicable \n")
  cat("\n")
  
  
  cat("Diagnostic accuracies \n")
  output1 <- round(cbind(x$sens$sens, x$sens$sens.ci, x$spec$spec, x$spec$spec.ci), digits)
  rownames(output1) <- x$names
  colnames(output1)[c(1,4)] <- c("sens", "spec")
  print(output1)
  
  cat("\n")
  cat("Diagnostic OR and likelihood Ratios \n")
  output2 <- round(cbind(x$DOR$DOR, x$DOR$DOR.ci, x$posLR$posLR, x$posLR$posLR.ci,
                         x$negLR$negLR, x$negLR$negLR.ci), digits)
  rownames(output2) <- x$names
  colnames(output2)[c(1,4,7)] <- c("DOR", "posLR", "negLR")
  print(output2)
  
  cat("\n")
  cat("Correlation of sensitivities and false positive rates: \n")
  print(round(CIrho(x$cor_sens_fpr, x$nobs), digits))
  return(invisible(NULL))
  }

CIrho <- function(rho, N, level = 0.95){
  stopifnot(rho < 1, rho > -1, N > 3, round(N) == N)
  z <- atanh(rho)
  kappa <- qnorm(1-(1-level)/2)
  output <- c(rho, tanh(z - kappa*sqrt(1/(N-3))),tanh(z + kappa*sqrt(1/(N-3))))
  names(output) <- c("rho", paste(100*(1-level)/2, "%", collapse =""), 
                     paste(100*(1- (1-level)/2), "%", collapse =""))
  output
}


madaNaive <- function(x, level = 0.95){
  
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

############################
# Naive Homogeneity checks #
############################

# Naively test for equality of proportions using R's prop.test
naive.prop.test<-function(M, conf.level = 0.95, correct = TRUE){
	
	M<-as.matrix(M)
	true.pos<-M[,1]
	false.neg<-M[,2]
	false.pos<-M[,3]
	true.neg<-M[,4]
	number.of.pos<-true.pos+false.neg
	number.of.neg<-false.pos+true.neg
	# First calculate the htests for sens and spec
	sens.htest<-prop.test(true.pos,number.of.pos,conf.level = conf.level, correct = correct)
	spec.htest<-prop.test(true.neg,number.of.neg,conf.level = conf.level, correct = correct)
	
	format.test.result<-function(x, digits = 4)
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
	sens.result<-format.test.result(sens.htest)
	spec.result<-format.test.result(spec.htest)
	output<-paste(c("Test for equality of sensitivities: \n",sens.result, "\n","Test for equality of specificities: \n",spec.result, "\n", "Both p-values for confidence level ",conf.level),collapse="")
	cat(output)
	}

#calculate cochranes Q, expects (s.e. of) log transformed DOR, posLR, negLR 
cochrane.Q<-function(x, sd.x, weights = NULL)
{
	if(is.null(weights)){
		if(length(weights != length(x))){stop("Length of weights does not match length of x")}				
		weights <- 1/(sd.x^2)}
	if(length(x) != length(sd.x)){stop("Length of x does not match length of sd.x")}
	bar.x<-sum(weights*x)/sum(weights)
	Q<-sum(weights*((x-bar.x)^2))
	k<-length(x)
	if(k <= 1)stop("x needs to have length at least 2")
	p.value<-pchisq(Q, k-1, lower.tail = FALSE)
	#fp<-format(round(p.value, 4))
	#fQ<-format(round(Q,4))
	output<-c(Q,p.value,as.integer(k-1))
	names(output)<-c("Cochran's Q", "p-value", "df")
	return(output)
	}
	
	
# Calculate Higgins I^2	
Higgins.Isq<-function(T,df){return(round(max(0,100*(T-df)/T),2))}

# naively pool proportions like sensitivity or specificity (in absence of any heterogeneity)
naive.prop.pool<-function(success, total){
	return(sum(success)/sum(total))	
	}

# from Cochran.Q and inverse variance weights calculate between study variance
naive.tausquared<-function(Q,weights)
{
k<-length(weights)
if(Q<(k-1)){return(0)}
else
return((Q-k+1)/(sum(weights)-(sum(weights^2)/sum(weights))))
}

	
# naively pool theta (i.e. LR or (D)OR) using either Mantel Haenszel or DerSimonian Laird method (method = c("MH","DSL"), type = c("posLR", "negLR", "DOR")
# TODO: method = "DSL" gives weird values! (very small!) look into that!

naive.pool<-function(M, type = "DOR", method="MH")
{
M<-as.matrix(M)
if(! method %in% c("MH","DSL"))stop("method must be either \"MH\" or \"DSL\"")else
total<-rowSums(M)
true.pos<-M[,1]
false.neg<-M[,2]
false.pos<-M[,3]
true.neg<-M[,4]
number.of.pos<-true.pos+false.neg
number.of.neg<-false.pos+true.neg
if(method == "MH")
{
weights<-switch(type, "DOR" = false.pos*false.neg/total, "posLR" = false.pos*number.of.pos/total, "negLR" = true.neg*number.of.pos/total)
output<-switch(type, "DOR" = sum(weights*calc.DOR(M)[,1])/sum(weights), "posLR" = sum(weights*calc.LR(M)[,1])/sum(weights), "negLR" = sum(weights*calc.LR(M)[,4])/sum(weights))
output<-cbind(output)
colnames(output)<-type
return(output)
}
if(method == "DSL")
{
ln.se<-sqrt(switch(type, "DOR"=(1/true.pos+1/true.neg+1/false.neg+1/false.pos),
		"posLR" = (1/true.pos+1/false.pos-1/number.of.pos-1/number.of.neg),
		"negLR" = (1/true.neg+1/false.neg-1/number.of.pos-1/number.of.neg)))

ln.theta<-log(switch(type, "DOR" = calc.DOR(M)[,1], "posLR" = calc.LR(M)[,1], "negLR" = calc.LR(M)[,4]))

cQ<-cochrane.Q(ln.theta,ln.se)[1]

tau.squared<-naive.tausquared(cQ,1/(ln.se^2))

weights.DL<-1/(ln.se^2+tau.squared)

output<-exp(sum(weights.DL*ln.theta)/sum(weights.DL))
output<-cbind(output)
colnames(output)<-type
return(output)
}
}