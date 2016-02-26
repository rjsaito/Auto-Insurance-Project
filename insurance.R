#Insurance program
#Author: Marcos Marin-Galiano, University of Dortmund
#Mail: marcos.marin-galiano@uni-dortmund.de
#Version 1.0
#This can be downloaded from:
#http://www.statistik.uni-dortmund.de/sfb475/berichte/insurance.zip
#Only for scientific use.
#NO WARRANTY.
insurance<-function(tvdat,varia,y,interval,seed=0,bigmean=FALSE,
bigprob=FALSE,prob.method="multinom",pred.method="svm",
multiteration=300,tratio=0.5,noiter=500,tolerance=10^(-6),
nooptim=FALSE,tablet=TRUE,criterion="MSE",pareto=FALSE,onlydata="no",
opti.method="Nelder-Mead",opti.mode="all-in-one",NA.hold=FALSE)
{
#Reading all required R-packages. Contribution of an output-list.
if (opti.mode!="all-in-one" & opti.mode!="single") stop("False value
chosen for opti.mode. Value must be all-in-one or single.")
library(MASS) ; library(e1071) ; library(nnet)
output<-NULL
#Erase all data with Missing Values in one of the independent
#variables (code: -9).

if (NA.hold==TRUE) ldata<-tvdat else
{if (length(varia)==1) ldata<-tvdat[which(tvdat[,varia]!=-9),]
else ldata<-tvdat[which(apply(tvdat[,varia]!=-9,1,prod)==1),]}
rm(tvdat) ; gc(TRUE)
#Computation of the number of classes (without "0"-class)
classes<-length(interval)+1
#Determine, to which class of claim size each observation belongs
#(according to interval).
classif<-rep(classes,dim(ldata)[1])
for (i in 1:(classes-1)) classif[which(ldata[,y]<
interval[classes-i])]<-classes-i
classif[which(ldata[,y]==0)]<-0
ldata<-cbind(ldata,classif)
if (length(unique(sort(classif)))!=classes+1) stop("Empty classes
present.")
rm(classif) ; gc(TRUE)
#Building a data frame, where all useless data will be omitted.
ldata<-ldata[,c(varia,y,dim(ldata)[2])]
#Re-determination of all column numbers.
varia<-1:length(varia) ; y<-length(varia)+1
classvar<-length(varia)+2
#Seperation of training and validation data.
set.seed(seed)
seperate<-sample(1:dim(ldata)[1],round(dim(ldata)[1]*tratio))
traindat<-ldata[seperate,] ; validat<-ldata[-seperate,]
rm(seperate) ; rm(ldata) ; gc(TRUE)
#Requesting, whether progamm was only used to generate a data set
#rather than analysing anything
if (onlydata=="training data") return(traindat)
if (onlydata=="validation data") return(validat)
if (onlydata!="no") stop()
#Estimating class probabilities with multinom.
if (bigprob==FALSE) classdat<-traindat[which(traindat[,classvar]<
classes),] else classdat<-traindat
if (bigprob==FALSE) pbig<-length(which(traindat[,y]>
interval[classes-1]))/dim(traindat)[1] else pbig<-NULL
classdat[,classvar]<-as.factor(classdat[,classvar])
for (i in 1:length(varia)) classdat[,i]<-factor(classdat[,i])
if (tablet==TRUE)
{
#Running multinom by using a frequency table.

freque<-as.data.frame(table(classdat[,-y]))
multiclass<-length(unique(sort(freque[,length(varia)+1])))
freque2<-matrix(freque[,length(varia)+2],ncol=multiclass)
classdat<-freque[1:(dim(freque)[1]/multiclass),varia]
if (mode(classdat)=="numeric")
{classdat<-as.data.frame(classdat)
names(classdat)<-names(freque)[varia]}
delete<-which(apply(freque2,1,sum)==0)
if(length(delete)!=0)
{freque2<-freque2[-delete,] ; classdat<-classdat[-delete,]}
if (mode(classdat)=="numeric")
{classdat<-as.data.frame(classdat)
names(classdat)<-names(freque)[varia]}
if (prob.method=="multinom") pmodel<-multinom(as.formula(paste
("freque2~",paste(names(classdat),collapse="+"))),classdat,
maxit=multiteration) else pmodel<-NULL
}
else
{
#Running multinom directly.
if (prob.method=="multinom") pmodel<-multinom(as.formula(paste
("classif~",paste(names(classdat)[varia],collapse="+"))),classdat,
maxit=multiteration) else pmodel<-NULL
multiclass<-NULL ; freque<-NULL ; freque2<-NULL
}
rm(classdat) ; rm(multiclass) ; rm(freque) ; rm(freque2) ; gc(TRUE)
#Applying pmodel on validation data.
classdat<-validat[,varia]
if (mode(classdat)=="numeric")
{classdat<-as.data.frame(classdat)
names(classdat)<-names(validat)[varia]}
for (i in 1:dim(classdat)[2]) classdat[,i]<-factor(classdat[,i])
if (prob.method=="multinom") pmat<-predict(pmodel,classdat,type=
"probs") else stop("Illegal probability model")
rm(classdat) ; gc(TRUE)
#Scaling data to [0,1] and generating the formula for the desired
#model.
scalemat<-matrix(0,nrow=length(varia),ncol=2)
for (i in 1:length(varia))
{scalemat[i,1]<-max(traindat[,i])
scalemat[i,2]<-min(traindat[,i])
traindat[,i]<-(traindat[,i]-scalemat[i,2])/
(scalemat[i,1]-scalemat[i,2])
validat[,i]<-(validat[,i]-scalemat[i,2])/
(scalemat[i,1]-scalemat[i,2])}
model.form<-as.formula(paste(c(names(traindat)[y],paste(names
(traindat)[varia],collapse="+")),collapse="~"))
dimnames(scalemat)[[2]]<-c("max","min")
dimnames(scalemat)[[1]]<-names(traindat)[varia]
#Seperation of traindat into partial data for each claim class.
traindat<-traindat[traindat[,classvar]>0,]
traindat<-split(traindat[,1:y],traindat[,classvar])
if (bigmean==FALSE) eff.class<-classes else eff.class<-classes-1
#Estimating the mean value for the highest claim by computating the
#mean or by using the expectation of a generalized Pareto distribution.
if (pareto==FALSE)
{
if (bigmean==TRUE) bigmeanvalue<-mean(traindat[[classes]][,y]) else
bigmeanvalue<-NULL
}
else
{
if (bigmean==TRUE)
{
#Two functions for estimating an appropriate generalized Pareto distribution.
ac.logdensityParetoAS<-function(x,alpha,sigma)
{
if (alpha!=0) {logdensity<- -sigma-(1+alpha)*log(1+x/(alpha*
exp(sigma)))}
if (alpha==0) {logdensity<- -sigma-x/exp(sigma)}
return(logdensity)
}
ac.MLE.GPD<-function(x,xi,beta,level=0.05)
{
my.data.x<-x ; my.p<-2
assign("my.data.x",my.data.x,env=.GlobalEnv)
# f1 is the function you want to minimize
f1<-function(theta)
{
logL<- -sum(ac.logdensityParetoAS(my.data.x,theta[1],theta[2]))
return(logL)
}
# Starting vector
theta0<-matrix(c(1/xi,log(beta)),ncol=1)
# Minimization.
out<-optim(theta0,f1,method="Nelder-Mead",control=list
(maxit=20000,reltol=1.0E-12))
theta<-out$par
xi<-1/theta[1] ; beta<-exp(theta[2])
mean.observed<-mean(my.data.x) ; expectation<-NA
if (xi<1) expectation<-beta/(1-xi)
out.cov<-matrix(c(0,0,0,0),ncol=2)
out.cov[1,1]<-(1+xi)^2 ; out.cov[1,2]<-beta*(1+xi)
out.cov[2,1]<-out.cov[1,2] ; out.cov[2,2] <- 2*beta*beta*(1+xi)
out.cov<-out.cov/length(x)
out.CIxi<-c(xi-qnorm(1-level/2)*sqrt(out.cov[1,1]),
xi+qnorm(1-level/2)*sqrt(out.cov[1,1]))
out.CIbeta<-c(beta-qnorm(1-level/2)*sqrt(out.cov[2,2]),
beta+qnorm(1-level/2)*sqrt(out.cov[2,2]))
return(out=list(xi=xi,beta=beta,CI.xi=out.CIxi,
CI.beta=out.CIbeta,cov=out.cov,CI.level=level,
expectation=expectation,mean.observed=mean.observed))
}
modell<-ac.MLE.GPD(traindat[[classes]][,y]-50000,10,2)
bigmeanvalue<-modell$expectation+50000}
else bigmeanvalue<-NULL
}
#Computing optimal parameters (c.f. Cherkassky/Ma (2004))
beta0<-rep(0,eff.class*3)
for (i in 1:eff.class)
{
var1<-mean(traindat[[i]][,y])+3*sqrt(var(traindat[[i]][,y]))
noise<-sqrt(mean((glm(model.form,family=Gamma(link="log"),
traindat[[i]],maxit=200)$fitted-traindat[[i]][,y])^2))
noise<-max(10^(-3),noise)
n<-dim(traindat[[i]])[1] ; var2<-3*noise*sqrt(log(n)/n)
var3<-1/(2*((0.35^(1/length(varia)))^2))
beta0[((i-1)*3+1):(i*3)]<-c(var1,var2,var3)
}
if (nooptim==FALSE)
{
#Optimizing the parameters of Cherkassky/Ma by using the
#Nelder-Mead-Algorithm.
#Function for estimating MSE of validat, when traindat and beta
#are used to fit the model.
MSE.func<-function(beta,traindat,validat,model.form,pred.method,
eff.class,y,pmat,bigprob,bigmean,pbig,bigmeanvalue,criterion,
opti.mode)
{
#Fitting svm model.
beta<-exp(beta)
costmat<-matrix(0,nrow=dim(validat)[1],ncol=eff.class)
if(opti.mode=="all-in-one")
{
if (pred.method=="svm")
for (i in 1:eff.class)
{model<-svm(model.form,traindat[[i]],type="eps-regression",
cost=beta[i*3-2],epsilon=beta[i*3-1],gamma=beta[i*3],scale=FALSE)
costmat[,i]<-predict(model,validat)}
else stop("Illegal prediction model")
#Computation of the estimated premia for different settings of
#bigprob and bigmean.
if (bigprob==TRUE & bigmean==TRUE)
{pmat2<-pmat[,-1] ; costmat2<-cbind(costmat,bigmeanvalue)
premia<-apply(pmat2*costmat2,1,sum)}
if (bigprob==TRUE & bigmean==FALSE)
{pmat2<-pmat[,-1] ; costmat2<-costmat
premia<-apply(pmat2*costmat2,1,sum)}
if (bigprob==FALSE & bigmean==TRUE)
{pmat2<-cbind(pmat[,-1]*(1-pbig),pbig)
costmat2<-cbind(costmat,bigmeanvalue)
premia<-apply(pmat2*costmat2,1,sum)}
if (bigprob==FALSE & bigmean==FALSE)
{pmat2<-cbind(pmat[,-1]*(1-pbig),pbig) ; costmat2<-costmat
premia<-apply(pmat2*costmat2,1,sum)}
}
else
{
if (pred.method=="svm")
{model<-svm(model.form,traindat,type="eps-regression",
cost=beta[1],epsilon=beta[2],gamma=beta[3],scale=FALSE)
premia<-predict(model,validat)}
else stop("Illegal prediction model")
}
#Computation of desired MSE according to the given criterion.
if (criterion=="trim.MSE" & opti.mode!="all-in-one") criterion<-"MSE"
if (criterion=="MSE") MSE<-mean((premia-validat[,y])^2)
if (criterion=="trim.MSE")
{trim.Punkte<-which(validat[,y]<50000)
MSE<-mean((premia[trim.Punkte]-validat[trim.Punkte,y])^2)}
if (criterion=="meanresid") MSE<-mean(abs(premia-validat[,y]))
if (criterion=="hosmer")
{
sortmat<-cbind(validat[,y],premia)
sortmat<-sortmat[order(sortmat[,2]),]
n<-dim(sortmat)[1] ; blocks<-max(round(n/1000),10)
observed<-rep(0,blocks) ; expected<-rep(0,blocks)
border<-round(n/blocks*(0:blocks))
for (i in 1:blocks) observed[i]<-
mean(sortmat[(border[i]+1):border[i+1],1])
for (j in 1:blocks) expected[j]<-
mean(sortmat[(border[j]+1):border[j+1],2])
MSE<-sum((observed-expected)^2/expected)}
return(MSE)
}
#Optimization step.
if (opti.mode=="all-in-one")
beta<-optim(log(beta0),MSE.func,method=opti.method,control=
list(maxit=noiter,reltol=tolerance),traindat=traindat,
validat=validat,model.form=model.form,pred.method=pred.method,
eff.class=eff.class,y=y,pmat=pmat,bigprob=bigprob,bigmean=bigmean,
pbig=pbig,bigmeanvalue=bigmeanvalue,criterion=criterion,
opti.mode=opti.mode)$par
else
{
beta<-rep(0,eff.class*3)
for (k in 1:eff.class) beta[((k-1)*3+1):(k*3)]<-
optim(log(beta0[((k-1)*3+1):(k*3)]),MSE.func,method="Nelder-Mead",
control=list(maxit=noiter,reltol=tolerance),traindat=traindat[[k]],
validat=validat[validat[,classvar]==k,],model.form=model.form,
pred.method=pred.method,eff.class=eff.class,y=y,pmat=pmat,
bigprob=bigprob,bigmean=bigmean,pbig=pbig,bigmeanvalue=bigmeanvalue,
criterion=criterion,opti.mode=opti.mode)$par
}
}
if (nooptim==FALSE) modelmat<-matrix(exp(beta),byrow=T,ncol=3)
else modelmat<-matrix(beta0,byrow=T,ncol=3)
dimnames(modelmat)[[2]]<-c("C","epsilon","gamma")
output<-list()
if (pred.method=="svm")
#Fitting the svm models for the optimized parameters. Writing of all
#results.
{for (i in 1:eff.class)
{output[[i]]<-svm(model.form,traindat[[i]],type="eps-regression",
cost=modelmat[i,1],epsilon=modelmat[i,2],gamma=modelmat[i,3],
scale=FALSE)
}}
beta0<-matrix(beta0,ncol=3,byrow=TRUE)
dimnames(beta0)[[2]]<-c("C","epsilon","gamma")
output<-list(pred.method=pred.method,pred.model=output,
bigprob=bigprob,pbig=pbig,prob.method=prob.method,pmodel=pmodel,
modelmat=modelmat,interval=interval,bigmean=bigmean,
bigmeanvalue=bigmeanvalue,scalemat=scalemat,beta0=beta0,
criterion=criterion,opti.mode=opti.mode, NA.hold=NA.hold)
class(output)<-c("insurance","predict")
output$call<-match.call()
return(output)
}