
library(survival)
library(mvtnorm)

# Function for computing importance of predictor variables in a survival model based on Breiman's permutation method.
# Note that only right censoring is allowed.

breiman <- function (model, time, event=rep(1,length(time)), covdata, nrep,
                     plot=TRUE){
  
  #calculate the concordance of model (will be used as benchmark)
  if (attr(model, "class")=="survreg") {
    cmod <- 1-survConcordance(Surv(cbind(time=time, event=event))~ predict(model, newdata=covdata))$concordance
    }  else if ("glmnet" %in% attr(model, "class")) {
      cmod <- survConcordance(Surv(cbind(time=time, event=event)) ~ predict(model, newdata=covdata))$concordance
    } else {
      cmod <- survConcordance(Surv(cbind(time=time, event=event))~ as.vector(predict(model, newdata=covdata)))$concordance
    }
  
  #compute variable importance
  ncov <- ncol(covdata)
  nobs <- nrow(covdata)
  
  vi <- sapply(1:ncov, function (i){
    simMat <- covdata[rep(seq(nobs), nrep), ]
    permutedCol <- as.vector(replicate(nrep, sample(covdata[,i], nobs, replace=FALSE)))
    indicatorNrep <- rep(1:nrep, each=nobs)
    simMat[,i] <- permutedCol
    if (attr(model, "class")=="survreg") {
      permuted.pred <- predict(model, newdata=simMat, type="response")
      scmod <- tapply(permuted.pred, indicatorNrep, function (j)
        1-survConcordance(Surv(cbind(time=time, event=event))~j)$concordance)
      } else if ("glmnet" %in% attr(model, "class")) {
        
      } else if ("CoxBoost" %in% attr(model, "class")) {
        permuted.pred <- as.vector(predict(model, newdata=simMat))
        scmod <- tapply(permuted.pred, indicatorNrep, function (j)
          survConcordance(Surv(cbind(time=time, event=event))~j)$concordance)
      } else {
        permuted.pred <- predict(model, newdata=simMat, type="risk")
      scmod <- tapply(permuted.pred, indicatorNrep, function (j)
        survConcordance(Surv(cbind(time=time, event=event))~j)$concordance)
      }
    
    mean((cmod-scmod)/cmod)})
  
  names(vi) <- colnames(covdata)
  
  #dotchart displaying variable importance
  if (plot==TRUE) {dotchart(vi, xlab="Relative importance")}
  
  return(vi)
}

