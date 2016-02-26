setwd("C:/Users/Riki/Dropbox/UMN Courses/STAT 8051/Travelers/")
kangaroo <- read.csv("Kangaroo.csv")
library(dplyr)
kangtrain <- subset(kangaroo, split == "T")
kangtest <- subset(kangaroo, split != "T")

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


#data
boxplot(claimcst0~veh_body, data=kangtrain)
thsd = TukeyHSD(aov(claimcst0~veh_body, data=subset(kangtrain, clm>0)))
boxplot(claimcst0~gender,data=subset(kangtrain, clm>0))
boxplot(claimcst0~area, data=subset(kangtrain, clm>0))
boxplot(claimcst0~factor(agecat), data=subset(kangtrain, clm>0))
pairs(kangtrain[,c("numclaims","claimcst0","veh_value","exposure","veh_age","area","agecat")])
pairs(subset(kangtrain,clm>0)[,c("numclaims","claimcst0","veh_value","exposure","veh_age","area","agecat")])

jpeg("hist_claimcst0.jpeg")
hist(kangtrain$claimcst0, main="Histogram of claimcst0", xlab="claimcst0")
dev.off()

table(kangtrain$clm)
table(kangtrain$numclaim)

# models
#0 Null Model
NormalizedGini(solution = kangtrain$claimcst0, submission = mean(kangtrain$claimcst0))

#1 Full OLS
ev = c("veh_value","exposure","veh_body","veh_age","gender","area","agecat")
ols <-  lm(claimcst0 ~ veh_value+exposure+veh_body+veh_age+gender+
	area+agecat, data=kangtrain)
summary(ols)
library(car)
Anova(ols)
(ols.gini = NormalizedGini(solution = kangtrain$claimcst0, submission = predict(ols)))


#2: GLM Tweedie distribution link
library(statmod)
tweed <- glm((claimcst0/exposure) ~ (veh_value+area+veh_age+gender+factor(agecat))^2, 
	data=kangtrain,  family=tweedie(var.power=1.3,link.power=0))
summary(tweed)
(tweed.gini = NormalizedGini(solution = kangtrain$claimcst0, submission = predict(tweed,type="response")*kangtrain$exposure))



#3: two part model (Frequency * Severity)
# We will divide the data into two levels of severity
# offset as a weight (frequency) or offset (count)?

#model 1: count (offset(log(exposure))
pm.count = glm(numclaims ~ offset(log(exposure))+(veh_value+veh_age+gender+area+factor(agecat))^2, 
	family = poisson, data=kangtrain)

pm.count.sub = step(pm.count, test="Chi")
(dp<-sum(residuals(pm.count.sub,type="pearson")^2)/pm.count.sub$df.res)
summary(pm.count.sub)

pm.small = glm(numclaims ~ offset(log(exposure))+veh_value+area+factor(agecat)+
	veh_value:veh_age+veh_value:area+gender:area, family=poisson, data=kangtrain)
summary(pm.small)

# frequency ( exposure as weight ) 
pm.freq = glm(numclaims ~ (veh_value+veh_age+gender+area+agecat)^2, 
	weight = exposure, family = poisson, data=kangtrain)

pm.freq.sub = step(pm.freq, test="Chi")
(dp<-sum(residuals(pm.freq.sub,type="pearson")^2)/pm.freq.sub$df.res)
summary(pm.freq.sub)


# freq and count
Freq =  predict(pm.freq.sub, newdata = kangaroo, type="response")
Count = predict(pm.small, newdata = kangaroo, type="response")
plot(Freq, Count); abline(0,1)
kang.updated = mutate(kangaroo, freq = as.numeric(Freq), count = as.numeric(Count))

outlier = which(Freq == max(Freq))
kangaroo[outlier,] #this guy is a trouble maker



#model 2: severity
# gamma model (terrible!)
gam.freq <- glm((claimcst0/numclaims) ~ 0+I(freq^2)+(freq+veh_value+gender+area+factor(agecat))^2,
	family=Gamma(link="log"),data=subset(kang.updated,clm > 0 & split=="T"))
gam.freq.sub = step(gam.freq)
summary(gam.freq.sub)

gam.count <- glm((claimcst0/numclaims) ~ 0+I(count^2)+(count+veh_value+gender+area+factor(agecat))^2,
	family=Gamma(link="log"),data=subset(kang.updated,clm > 0 & split=="T"))
gam.count.sub = step(gam.count)
summary(gam.count.sub)

#inverse gaussian model
ivg <- glm((claimcst0/numclaims) ~ veh_value+gender+area+factor(agecat),
	family=inverse.gaussian(link="log"),data=subset(kang.updated, clm > 0 & split=="T"))
ivg.sub = step(ivg)
summary(ivg.sub)


# predictions
preds = data.frame(
	claimcst0 = kang.updated$claimcst0,
	numclaims = kang.updated$numclaims,
      tweedie = predict(tweed, newdata=kang.updated, type="response")*kang.updated$exposure,
	gam.freq = Freq*predict(gam.freq.sub, newdata=kang.updated, type="response"),
	gam.count = Count*predict(gam.count.sub, newdata=kang.updated, type="response"),
	ivg.freq = Freq*predict(ivg.sub, newdata=kang.updated, type="response"),
	ivg.count = Count*predict(ivg.sub, newdata=kang.updated, type="response"),
	split = kang.updated$split
)


# gini coefficient on train data
apply(subset(preds, split == "T")[,c(3,6:7)], 2, function(x) NormalizedGini(kangtrain$claimcst0, x))

# means of test data predictions
colMeans(subset(preds, split == "T")[,c(3,6:7)])
colMeans(subset(preds, split != "T")[,c(3,6:7)])

library(lattice)
library(latticeExtra)
ecdfplot(~ tweedie + ivg.freq + ivg.count, data=preds, xlim=c(0,20000), auto.key=list(space='right'))


#cross validation with gini coefficients
cv <- function(fit, fit2 = NULL, data, data2 = NULL, K, R){
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
  for (i in seq_len(ms)) {
    j.out <- seq_len(n)[(s == i)]
    if(K > 1) j.in <- seq_len(n)[(s != i)] else if (K==1) j.in = j.out
    Call$data <- data[j.in, , drop = FALSE]
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
  }
  CV$Fitted <- CV$Fitted[order(CV$Fitted[,1]),2]
  CV
}


cv.tweed <- cv(fit=tweed,data= kangtrain, K=10)
NormalizedGini(kangtrain$claimcst0, cv.tweed$Fitted*kangtrain$exposure)

cv.ivg.count <- cv(fit=pm.count.sub, fit2=ivg.sub, data = kangtrain, data2=subset(kangtrain, clm>0), K=10)
NormalizedGini(kangtrain$claimcst0, cv.ivg.count$Fitted)


