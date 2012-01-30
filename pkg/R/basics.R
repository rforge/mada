sens <- function(x){x$TP/(x$TP + x$FN)}
fpr <- function(x){x$FP/(x$FP + x$TN)}



# Given a Nx4 Matrix compute the Sen, Spe of each study, columns of the matrix are: TP, FN, FP, TN. Alternatively M can be a data frame consisting of 4 variables ordered correctly
# Calculates FPR instead of Specificity if FPR=TRUE
# Continuity Correction: Add 0.5 to all cells if one 0 cell is present
# TODO: maybe implement a function that checks the sanity of the data and the sanity of the continuity correction
#TODO: maybe implement other continuity corrections or even a function for corrections

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

calc.SensSpec<-function(M, FPR = FALSE, correction = 0.5, conf.level = 0.95, method = "wilson")
{
# convert M to matrix  if necessary
M<-as.matrix(M)
if(any(M == 0)){M<-M+correction}
true.pos<-M[,1]
false.neg<-M[,2]
false.pos<-M[,3]
true.neg<-M[,4]
number.of.pos<-true.pos+false.neg
number.of.neg<-false.pos+true.neg
sens<-true.pos/number.of.pos
sens.ci<-binomCIvector(true.pos,number.of.pos, conf.level  = conf.level, method = method)

if(FPR==FALSE)
{
spec<-true.neg/number.of.neg
spec.ci<-binomCIvector(true.neg,number.of.neg, conf.level  = conf.level, method = method)
output<-cbind(sens,sens.ci,spec,spec.ci)
colnames(output)<-c("sensitiviy","Se LL", "Se UL","specifity", "Sp LL", "Sp UL")
return(output)
}else
FPR<-false.pos/number.of.neg
FPR.ci<-binomCIvector(false.pos,number.of.neg, conf.level  = conf.level, method = method)
output<-cbind(sens,sens.ci,FPR,FPR.ci)
colnames(output)<-c("sensitiviy","Se LL", "Se UL","false positive rate", "FPR LL", "FPR UL")
return(output)
} # End of function calc.SensSpec


# Calculate positive and negative Likelihood Ratios
# Input as in calc.SensSpec
# conf.level gives the level of the confidence intervalls
calc.LR<-function(M, correction = 0.5, conf.level = 0.95)
{
	
	SenFPR<-calc.SensSpec(M, FPR = TRUE, correction = correction)
	posLR<-SenFPR[,1]/SenFPR[,4]
	negLR<-(1-SenFPR[,1])/(1-SenFPR[,4])
	
	# now for the confidence intervals: use normal approx for log-transformed LR and well-known formula
	true.pos<-M[,1]
	false.neg<-M[,2]
	false.pos<-M[,3]
	true.neg<-M[,4]
	number.of.pos<-true.pos+false.neg
	number.of.neg<-false.pos+true.neg
	se.lnposLR<-sqrt(1/true.pos+1/false.pos-1/number.of.pos-1/number.of.neg)
	se.lnnegLR<-sqrt(1/true.neg+1/false.neg-1/number.of.pos-1/number.of.neg)
	alpha<-1-conf.level
	kappa<-qnorm((1-alpha/2))
	posLR.ci<-cbind(exp(-kappa*se.lnposLR),exp(kappa*se.lnposLR))
	posLR.ci<-posLR*posLR.ci
	negLR.ci<-cbind(exp(-kappa*se.lnnegLR),exp(kappa*se.lnnegLR))
	negLR.ci<-negLR*negLR.ci
	
	
	output<-cbind(posLR, posLR.ci, negLR, negLR.ci)
	colnames(output)<-c("pos LR", "pos LR LL", "pos LR UL", "neg LR", "neg LR LL", "neg LR UL")
	return(output)
} # End of function calc.LR

# Calculate DOR
# Input as in calc.SensSpec
calc.DOR<-function(M, correction = 0.5, conf.level = 0.95)
{
	LR<-calc.LR(M, correction = correction)
	DOR<-LR[,1]/LR[,4]

	true.pos<-M[,1]
	false.neg<-M[,2]
	false.pos<-M[,3]
	true.neg<-M[,4]

	se.lnDOR<-sqrt(1/true.pos+1/true.neg+1/false.neg+1/false.pos)
	alpha<-1-conf.level
	kappa<-qnorm((1-alpha/2))
	DOR.ci<-cbind(exp(-kappa*se.lnDOR),exp(kappa*se.lnDOR))
	DOR.ci<-DOR*DOR.ci
	
	output<-cbind(DOR,DOR.ci)
	colnames(output)<-c("DOR", "DOR LL", "DOR UL")
	return(output)
} # End of function calc.DOR



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

# Naively check for cut-off value problem by Spearman's rho
# For this calculate the correlation between observed Sens and Spe or Sens and FPR
# method is passed on to cor()

naive.cor<-function(M, method="spearman", correction = 0.5, FPR = FALSE)
{
	SensSpec<-calc.SensSpec(M, FPR = FPR, correction = correction, conf.level = 0.95, method = "wilson")
	Sens<-SensSpec[,1]
	SpecFPR<-SensSpec[,4]
	output<-cor(Sens, SpecFPR, method = method)
	names(output)<-paste(c("Correlation of sensitivity and "),ifelse(FPR, "false positive rate", "specificity"),collapse="")
	return(output)
}

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