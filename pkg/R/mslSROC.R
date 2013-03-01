## MosesShapiroLittenberg
mslSROC <- function(data = NULL, subset=NULL,
  TP="TP", FN="FN", FP="FP", TN="TN", 
  fpr = 1:99/100,
  alphasens = 1, alphafpr = 1,
  correction = 0.5, correction.control = "all",
  add = FALSE, lty = 1, lwd = 1, col = 1, ...){
  
  stopifnot(is.numeric(correction), 0 <= correction,  
            correction.control %in% c("all", "single", "none"),
            0 <= alphasens, alphasens <= 2, 0 <= alphafpr, alphafpr <= 2,
            is.numeric(TP) | (is.character(TP) & length(TP) == 1),
            is.numeric(FP) | (is.character(FP) & length(FP) == 1),
            is.numeric(TN) | (is.character(TN) & length(TN) == 1),
            is.numeric(FN) | (is.character(FN) & length(FN) == 1))
  
  if(!is.null(data) & is.character(c(TP,FP,TN,FN))){
    X <- as.data.frame(data)
    origdata <- data
    TP <- getElement(X,TP)
    FN <- getElement(X,FN)
    FP <- getElement(X,FP)
    TN <- getElement(X,TN)
  }
  
  if(is.null(data) & !is.character(c(TP,FP,TN,FN))){
    origdata <- data.frame(TP = TP, FN = FN, FP = FP, TN = TN)
  }
  
  freqdata <- cbind(TP,FN,FP,TN)
  checkdata(freqdata)
  
  N <- length(TP)  
  
  ## apply continuity correction to _all_ studies if one contains zero
  if(correction.control == "all"){if(any(c(TP,FN,FP,TN) == 0))
  {TP <- TP + correction;
   FN <- FN + correction;
   FP <- FP + correction;
   TN <- TN + correction}}
  if(correction.control == "single"){
    correction = ((((TP == 0)|(FN == 0))|(FP == 0))| (TN == 0)) * 
      correction
    TP <- correction + TP
    FN <- correction + FN
    FP <- correction + FP
    TN <- correction + TN
  }
  
  number.of.pos <- TP + FN
  number.of.neg <- FP + TN
  sens<-TP/number.of.pos
  FPR <- FP/number.of.neg
  
  senstrafo <- function(x){return(mada:::talpha(alphasens)$linkfun(x))}
  fprtrafo <- function(x){return(mada:::talpha(alphafpr)$linkfun(x))}
  sensinvtrafo <- function(x){return(mada:::talpha(alphasens)$linkinv(x))}
  
  z <- fprtrafo(FPR)
  y <- senstrafo(sens)  
  D <- y - z
  S <- y + z
  fit <- lm(D~S) 
  A1 <- fit$coefficients[1]
  B1 <- fit$coefficients[2]

  if(add){lines(fpr, sensinvtrafo(-(-A1/(1-B1))*((1-fpr)/fpr)^((1+B1)/(1-B1))),
        col=col,lty=lty,lwd=lwd, ...)}
  if(!add){plot(fpr, sensinvtrafo(-(-A1/(1-B1))*((1-fpr)/fpr)^((1+B1)/(1-B1))),
               col=col,lty=lty,lwd=lwd, type = "l",...)}
  return(invisible(NULL))
}


