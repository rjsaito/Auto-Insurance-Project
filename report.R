setwd("C:/Users/Riki/Dropbox/UMN Courses/STAT 8051/Travelers/")
kangaroo <- read.csv("Kangaroo.csv")
kangtrain <- subset(kangaroo, split == "T")

#"NormalizedGini" is the other half of the metric. This function does most of the work, though
SumModelGini <- function(solution, submission) {
  df = data.frame(solution = solution, submission = submission)
  df <- df[order(df$submission, decreasing = TRUE),]
  df
  df$random = (1:nrow(df))/nrow(df)
  df
  totalPos <- sum(df$solution)
  df$cumPosFound <- cumsum(df$solution) # this will store the cumulative number of positive examples found (used for computing "Model Lorentz")
  df$Lorentz <- df$cumPosFound / totalPos # this will store the cumulative proportion of positive examples found ("Model Lorentz")
  df$Gini <- df$Lorentz - df$random # will store Lorentz minus random
  return(sum(df$Gini))
}

NormalizedGini <- function(solution, submission) {
  SumModelGini(solution, submission) / SumModelGini(solution, solution)
}

hist(kangtrain$claimcst0, main="Histogram of claimcst0", xlab="claimcst0")
max(kangtrain$claimcst0)

library(Hmisc)
plot(Ecdf(kangtrain$claimcst0), type="b", cex=.5, xlab="claimcst0", ylab="Cumulative Probability"); abline(h=1, lty=2)

table(kangtrain$clm)
table(kangtrain$numclaim)


#1 Full OLS
ols <-  lm(claimcst0 ~ veh_value+exposure+veh_body+veh_age+gender+
             area+agecat, data=kangtrain)
summary(ols)
library(car)
Anova(ols)

#2: GLM Tweedie distribution link
library(statmod)
tweed <- glm(claimcst0 ~ veh_value+exposure+veh_body+veh_age+gender+
  area+agecat, data=kangtrain,  family=tweedie(var.power=1.3,link.power=0))
summary(tweed)
Anova(tweed)



#cross validation with gini coefficients
cv <- function(fit, fit2 = NULL, data, data2 = NULL, K){
  cost = function(y, yhat) mean((y - yhat)^2)
  n = nrow(data)
  if(K > 1) s = sample(rep(1:K, ceiling(nrow(data)/K)),nrow(data)) else 
    if(K == 1) s = rep(1, nrow(data))
  glm.y <- fit$y
  cost.0 <- cost(glm.y, fitted(fit))
  ms <- max(s)
  call <- Call <- fit$call
  if(!is.null(fit2)) call2 <- Call2 <- fit2$call
  CV <- CV.coef <- NULL
  
  pb <- winProgressBar(title = "progress bar", min = 0, max = K, width = 300)
  Sys.time() -> start
  
  for (i in seq_len(ms)) {
    j.out <- seq_len(n)[(s == i)]
    if(K > 1) j.in <- seq_len(n)[(s != i)] else if (K==1) j.in = j.out
    Call$data <- data[j.in, , drop = FALSE]; 
    d.glm <- eval.parent(Call)
    pred.glm <- predict(d.glm, newdata=data[j.out,], type="response")
    if(!is.null(fit2) & !is.null(data2)){
      j2.out.data <- merge(data2, data[j.out,])
      if(K > 1) j2.in.data <- merge(data2, data[j.in,]) else if (K==1) j2.in.data = j2.out.data
      Call2$data <- j2.in.data
      d.glm2 <- eval.parent(Call2)
      pred.glm2 <- predict(d.glm2, newdata=data[j.out,], type="response")
    }
    if(!is.null(fit2)) CV$Fitted = rbind(CV$Fitted, cbind(j.out, pred.glm*pred.glm2)) else 
      CV$Fitted = rbind(CV$Fitted, cbind(j.out, pred.glm))
    CV.coef$coef <- rbind(CV.coef$coef, coef(d.glm))
    CV.coef$se <- rbind(CV.coef$se, coef(summary(d.glm))[,2])
    Sys.sleep(0.1); setWinProgressBar(pb, i, title=paste( round(i/K*100, 0),"% done"))
  }
  close(pb); Sys.time() -> end
  cat("Cross-Validation Time Elapsed: ", round(difftime(end, start, units="secs"),3) ,"seconds \n")
  Fitted <- CV$Fitted[order(CV$Fitted[,1]),2]
  Fitted
}

# bootstrap
library(boot)
bs <- function(formula, data, family, indices) {
  d <- data[indices,] # allows boot to select sample 
  fit <- glm(formula, family, data=d)
  return(coef(fit)) 
} 

#two part model (Count * Severity)
#model 1: count (offset(log(exposure))
pm.sub <- glm(numclaims ~ offset(log(exposure))+factor(agecat)+area+veh_value+veh_age+
                veh_value:veh_age+area:veh_value, family = poisson, data=subset(kangtrain))
summary(pm.sub)
Count = predict(pm.sub, newdata = kangaroo, type="response")

pm.bs <- boot(data=subset(kangtrain), statistic=bs, R=10, formula=formula(pm.sub), family = poisson, parallel="multicore")
cbind(coef(pm.sub),colMeans(pm.bs$t))

pm.sub.bs <- pm.sub
pm.sub.bs$coefficients <- colMeans(pm.bs$t)
Count.bs = predict(pm.sub.bs, newdata = kangaroo, type="response")

#model 2: severity inverse gaussian model
ivg.sub <- glm((claimcst0/numclaims) ~ gender + veh_age + agecat,
               family=inverse.gaussian(link="log"),data=subset(kangtrain, clm > 0))
Severity = predict(ivg.sub, newdata=kangaroo, type="response")

ivg.bs <- boot(data=subset(kangtrain, clm > 0 & veh_value>0), statistic=bs, R=10, formula=formula(ivg.sub), family=inverse.gaussian(link="log"), parallel="multicore")
cbind(coef(ivg.sub),colMeans(ivg.bs$t))

ivg.sub.bs <- ivg.sub
ivg.sub.bs$coefficients <- colMeans(ivg.bs$t)
Severity.bs = predict(ivg.sub.bs, newdata = kangaroo, type="response")






